#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <fcntl.h>
#include "align_seq.h"

#define MAX_SEQ_LEN 6000 // max length 6000 bp for now
#define NCBI_API_KEY "d9d8a429deb9981e28976a848dd00a3fd708"
#define MAX_JSON 4096


/* Combine the NCBI and patient sequence into one string:
    * Given the NCBI sequence (wt_seq) and the patient sequence,
    * populate the empty string, combined. Since wt_seq and patient
    * strings are null-terminated (from implementation in fetch_seq.c),
    * and strncat() always adds a null terminator, we know that combined
    * is null-terminated.
 */
void combine_seq(char wt_seq[], char patient[], char combined[]) {
    combined[0] = '\0'; // strncat requires dest be null-terminated
    strncat(combined, wt_seq, strlen(wt_seq) + 1);
    combined[strlen(wt_seq) + 1] = '\n'; // ensure sequences separated by line
    strncat(combined, patient, strlen(patient) + 1);
}


/* Create the output file name:
    * Given the patient's file name, an empty string to populate,
    * and the required suffix, get the output file name. File
    * will have the same name as patient_filename with
    * "_aligned.fasta" suffix.
    * Alignment works best with a .fna or .txt patient file.
    * Returns the built string.
 */
char *get_outfile_name(char patient_filename[], char outfile[], const char *suffix) {
    // then when needed, can call "cat filename.fasta" to print to stdout
    char *dot_ptr = strchr(patient_filename, '.');
    // get the prefix of the patient's file name (without .fna for example)
    if (dot_ptr == NULL) {
        fprintf(stderr, "Error: patient file format is incorrect\n");
        exit(1);
    }
    int i = 0;
    while (patient_filename[i] != '.') {
        outfile[i] = patient_filename[i]; // add the prefix to outfile
        i++;
    }
    strncat(outfile, suffix, strlen(suffix) + 1); // add the suffix; null-terminated
    return outfile;
}


/* Align the nucleotide sequences via Clustal Omega:
    * Given the combined string of sequences, call clustalo to align
    * the sequences and produce a file with the resulting alignment.
    *
    * This function creates a pipe and forks a worker process to
    * run ./clustalo and all its arguments. Child process redirects
    * input to stdin as required by clustalo arguments. Parent will
    * write combined string into pipe so child can read. Then, the
    * worker/child process executes the command.
    *
    * Outputs a 0 on success, and 1 if any error occurred.
 */
int align(char combined[], char outfile[]) {
    // want parent to write combined into pipe, child reads
    int fd[2];
    if (pipe(fd) == -1) {
        perror("pipe");
        exit(1);
    }
    int pid = fork();
    if (pid < 0) {
        perror("fork");
        exit(1);
    }

    if (pid == 0) { // child
        /* example command:
            ./clustalo -i - -o sickle_aligned.fasta --outfmt=clustal --force */
        // -i - means let stdin be the input sequence file
        char *args[] = {"./clustalo", "-i", "-", "-o", outfile,
            "--outfmt=clustal", "--force", NULL};

        if (close(fd[1]) == -1) { // close child write
            perror("close write");
            exit(1);
        }
        if (dup2(fd[0], STDIN_FILENO) == -1) {
            perror("dup2");
            exit(1);
        }

        int devnull = open("/dev/null", O_WRONLY); // get fd for /dev/null
        if (dup2(devnull, STDERR_FILENO) == -1) {
            // ensure warning message from clustalo not printed to screen
            perror("dup2");
            exit(1);
        }
        if (close(devnull) == -1) {
            perror("close");
            exit(1);
        }

        if (close(fd[0]) == -1) { // read() blocks if pipe empty
            // read until EOF, when parent closes fd[1]
            perror("close read");
            exit(1);
        }
        execvp(args[0], args);
            // child process replaced; clustalo has internal read() call
        perror("execvp");
        exit(1);
    }
    // parent
    if (close(fd[0]) == -1) { // parent not reading
        perror("close read");
        exit(1);
    }
    // write combined into pipe
    if (write(fd[1], combined, strlen(combined) + 1) < 0) {
        perror("write");
        exit(1);
    }
    // close write to signal end
    if (close(fd[1]) == -1) {
        perror("close write");
        exit(1);
    }
    int status;
    if (wait(&status) == -1) {
        perror("wait");
        return 1;
    }
    if (!WIFEXITED(status)) { // if failed
        fprintf(stderr, "Error: Clustal Omega failed\n");
        return 1;
    }
    return 0;
}


/* Report mutations/changes between NCBI and patient sequence:
    * Given the name of the aligned output file. This function prints
    * to stdout a report of any base changes in the successfully
    * aligned files.
    *
    * This function traverses each line, storing in either wt_seq
    * or patient_seq. We CANNOT simply use the original NCBI sequence
    * or patient sequence, since the alignment optimizes the fit, potentially
    * altering positions. For example, a '-' dash might be introduced in
    * the sequence, representing either an addition (if in wild-type) or
    * deletion (if in patient).
    *
    * Outputs a 0 on success, and 1 if any error occurred.
 */
int report_mutations(char *alignment) {
    FILE *alignment_file = fopen(alignment, "r");
    if (alignment_file == NULL) {
        fprintf(stderr, "Error: could not open file\n");
        return 0;
    }
    char report[MAX_JSON] = "\0"; // will store all the mutations to print at once
    char line[MAX_JSON];
    char wt_seq[MAX_SEQ_LEN]; // will store just the wild-type line
    char patient_seq[MAX_SEQ_LEN]; // stores the patient line
    int wt_len = 0;
    int patient_len = 0;
    int which_line = 0;
        // tracks which line we are looking at
        // 0 will be for wild-type, 1 for patient, 2 for alignment matching

    // skip the header
    if (fgets(line, MAX_JSON, alignment_file) == NULL) {
        fprintf(stderr, "Error: could not read from file\n");
        fclose(alignment_file);
        return 1;
    }

    while (fgets(line, MAX_JSON, alignment_file) != NULL) {
        if (line[0] == '\n' || line[0] == '\r') {
            continue; // skip blank lines
        }
        if (which_line == 0) {
            // want pointer to last space char in gap before alignment
            char *seq_start = strrchr(line, ' ');
            if (seq_start == NULL) {
                // if character does not appear in string
                fprintf(stderr, "Error: incorrect alignment formatting\n");
                fclose(alignment_file);
                return 1;
            }
            seq_start++; // move to first base

            int i = 0;
            while (seq_start[i] != '\n' && seq_start[i] != '\r') {
                wt_seq[wt_len] = seq_start[i];
                wt_len++;
                i++;
            }
            which_line = 1;

        } else if (which_line == 1) {
            // on patient line
            char *seq_start = strrchr(line, ' ');
            if (seq_start == NULL) {
                // if character does not appear in string
                fprintf(stderr, "Error: incorrect alignment formatting\n");
                fclose(alignment_file);
                return 1;
            }
            seq_start++; // move to first base

            int j = 0;
            while (seq_start[j] != '\n' && seq_start[j] != '\r') {
                patient_seq[patient_len] = seq_start[j];
                patient_len++;
                j++;
            }
            which_line = 2;
        } else {
            which_line = 0; // alignment line
        }
    }
    wt_seq[wt_len] = '\0';
    patient_seq[patient_len] = '\0';
    fclose(alignment_file); // close file

    int mutations = 0; // track number of mutations
    int len;
    // want to use the shorter line
    if (wt_len < patient_len) {
        len = wt_len;
    } else {
        len = patient_len;
    }
    printf("\n==================================================\n");
    printf("Mutation report for %s:\n", alignment);
    for (int i = 0; i < len; i++) {
        if (wt_seq[i] != patient_seq[i]) {
            char change[5 * sizeof(char) + sizeof(int)] = "\0";
            // reps a single change, given some extra space
            if (snprintf(change, 3 * sizeof(char) + sizeof(int), "%c%d%c ",
                wt_seq[i], i + 1, patient_seq[i]) == -1) {
                perror("snprintf");
                exit(1);
            }
            // format the string to go into report
            strncat(report, change, strlen(change)); // auto null-terminates
            // format the string into report
            mutations++;
        }
    }
    if (mutations == 0) {
        printf("No mutations found\n");
        printf("==================================================\n");
    } else {
        printf("%s\n", report);
        printf("Total mutations: %d\n", mutations);
        printf("==================================================\n");
    }
    return 0;

}



