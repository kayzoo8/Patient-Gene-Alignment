#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h> // << what is this?
#include <sys/wait.h> // what are all these libraries

#define MAX_SEQ_LEN 6000 // max length 6000 bp for now
#define NCBI_API_KEY "d9d8a429deb9981e28976a848dd00a3fd708"
#define MAX_JSON 4096


/* Get the official gene id number from gene name:
    * Given the gene_name passed in by the user, populate
    * gene_id string with the corresponding gene ID.
    * For example, 3043 for HBB (Hemoglobin subunit beta).
    * Queries the NCBI Datasets command line tools to obtain id.
    *
    * Uses a pipe, where the child queries that database and
    * writes to the pipe, the parent reads from the pipe and
    * populates a summary_report string. The parent then extracts
    * the gene id number.
    *
    * Outputs a 0 on success, and 1 if any error occurred.
 */
int get_gene_id(const char *gene_name, char *gene_id) {
    // full command: ./datasets summary gene symbol HBB --report ids_only
    char *args[] = {"./datasets", "summary", "gene", "symbol", (char *)gene_name,
        "--report", "ids_only", NULL};

    char summary_report[MAX_JSON]; // 'buf'

    int fd[2];
    if (pipe(fd) == -1) {
        perror("pipe");
        return 1;
    }
    int pid = fork(); // fork a worker process
    if (pid < 0) {
        perror("fork");
        return 1;
    }

    // child process will run command
    if (pid == 0) {
        if (close(fd[0]) == -1) {
            perror("close child read");
            exit(1);
        }
        // want to redirect stdout to pipe write end since command prints to console
        if (dup2(fd[1], STDOUT_FILENO) == -1) {
            perror("dup2");
            exit(1);
        }

        if (close(fd[1]) == -1) { // close write end to prevent resource leaks
            perror("close child write");
            exit(1);
        }
        execvp(args[0], args);
        perror("execvp"); // shouldn't execute if no error
        exit(1);
    }

    // parent process now
    if (close(fd[1]) == -1) {
        perror("close parent write");
        return 1;
    }
    int total = 0;
    int num_bytes;
    // summary_report also rep pointer to first element of array
    while ((num_bytes = (int) read(fd[0], summary_report + total, MAX_JSON - total - 1)) > 0) {
        // - 1 to allow null-terminator
        total += num_bytes;
    } // read from pipe into summary_report

    summary_report[total] = '\0';
    if (close(fd[0]) == -1) {
        perror("close parent read");
        return 1;
    }

    int status;
    if (wait(&status) == -1) { // suspends parent process until child terminates
        perror("wait");
        return 1;
    }
    if (!WIFEXITED(status)) { // if failed
        // WIFEXITED(status) returns true if child exited normally
        fprintf(stderr, "Error: dataset summary retrieval failed\n");
        return 1;
    }
    // succeeded; have all of gene summary in summary_report
    // want to find the gene_id
    // formatted like: {"reports": [{"gene":{"gene_id":"3043"} ... only want first instance
    char *key = strstr(summary_report, "\"gene_id\":");
    if (key == NULL) { // returns null if str2 cannot be found in summary_report
        // error message prints in align.c
        return 1;
    }
    // else, gene_id was found
    key += strlen("\"gene_id\":\""); // pointing to first digit
    gene_id[0] = key[0];
    char *gene_id_end = key + 1; // offset, should eventually equal ending quote mark "
    int i = 1;
    while (*gene_id_end != '"') {
        gene_id[i] = *gene_id_end;
        i++;
        gene_id_end++;
    }
    gene_id[i] = '\0';
    return 0; // success
}


/* Get the nucleotide sequence from the gene id:
    * Given the gene_id number, this function will query
    * the NCBI database for the nucleotide sequence of
    * the specific gene. The length of the sequence is dictated
    * by the MAX_SEQ_LEN macro.
    *
    * Forked child queries the NCBI database, to download
    * the dataset related to that gene. A pipe is created
    * to unzip the resulting file and send the cDNA
    * sequence to the parent.
    *
    * Outputs a 0 on success, and 1 if any error occurred.
 */
int get_wt_ntseq(char *gene_id, char cdna_seq[]) {
    // example command: ./datasets download gene gene-id 3043 --include gene,protein,rna,cds
    int pid = fork();
    if (pid < 0) {
        perror("fork get_wt_ntseq 1");
        return -1;
    }
    if (pid == 0) { // child
        char *args[] = {"./datasets", "download", "gene", "gene-id", gene_id,
        "--include", "rna", "--no-progressbar", NULL};
        // only want rna file for now
        // after --include rna, can specify --filename gene_id.zip
        execvp(args[0], args);
        perror("execvp get_wt_ntseq 1");
        return -1;
    }
    // parent now
    int status;
    if (wait(&status) == -1) {
        perror("wait");
        return -1;
    }
    if (!WIFEXITED(status)) { // if failed
        fprintf(stderr, "Error: dataset download failed\n");
        return -1;
    }

    // now create a pipe and another process to unzip resulting file
    int fd[2];
    if (pipe(fd) == -1) {
        perror("pipe get_wt_ntseq");
        return -1;
    }

    int pid2 = fork();
    if (pid2 < 0) {
        perror("fork get_wt_ntseq 2");
        return -1;
    }

    if (pid2 == 0) { // child will write, close read
        char *args2[] = {"unzip", "-p", "-o", "ncbi_dataset.zip", "ncbi_dataset/data/rna.fna", NULL};
            // -o to overwrite existing files
        // command would print rna.fna file contents (sequence) to stdout

        if (close(fd[0]) == -1) {
            perror("close child read 2");
            exit(1);
        }
        if (dup2(fd[1], STDOUT_FILENO) == -1) { // redirect stdout to write end
            perror("dup2 get_wt_ntseq");
            exit(1);
        }

        if (close(fd[1]) == -1) {
            perror("close child write get_wt_ntseq");
            exit(1);
        }
        execvp(args2[0], args2);
        perror("execvp get_wt_ntseq");
        exit(1);
    }
    // now parent
    if (close(fd[1]) == -1) {
        perror("close parent write get_wt_ntseq");
        return 1;
    }

    int total = 0;
    int num_bytes;
    // read the data from pipe into cdna_seq
    while ((num_bytes = (int) read(fd[0], cdna_seq + total, MAX_SEQ_LEN - total - 1)) > 0) {
        total += num_bytes;
    }
    cdna_seq[total] = '\0';

    if (close(fd[0]) == -1) {
        perror("close parent read");
        return 1;
    }
    int status2;
    if (wait(&status2) == -1) {
        perror("wait");
        return 1;
    }
    if (!WIFEXITED(status2)) { // if failed
        fprintf(stderr, "Error: unzip failed\n");
        return 1;
    }

    return 0;
}


/* Get the nucleotide sequence of the patient from user input:
    * Given a file pointer to the patient file, and an empty string
    * patient_seq to read in the data. The function gets the length
    * of the file first, so it knows the approximate number of bytes
    * to read for fread(). The patient_seq string is then populated
    * and null-terminated.
    * All open files are closed within this function. This function
    * will keep any newline characters in the original patient file,
    * as it does not affect the alignment.
    *
    * Outputs a 0 on success, and 1 if any error occurred.
 */
int get_patient_seq(FILE *fasta, char patient_seq[]) {
    // get length of file for one fread call
    if (fseek(fasta, 0, SEEK_END) != 0) {
        fprintf(stderr, "Error: fseek\n");
        fclose(fasta); // close file after each error
        return 1;
    }
    // ftell gets current position of file pointer, returning number bytes from beginning of file
    int num_bytes = (int) ftell(fasta);
    if (num_bytes < 0) {
        fprintf(stderr, "Error: ftell\n");
        fclose(fasta);
        return 1;
    }
    if (fseek(fasta, 0, SEEK_SET) != 0) { // reset file pointer to beginning
        fprintf(stderr, "Error: fseek\n");
        fclose(fasta);
        return 1;
    }
    // read in all contents at once
    int bytes_read = (int) fread(patient_seq, sizeof(char), num_bytes, fasta);
    if (bytes_read == 0) {
        // print error message in align.c
        fclose(fasta);
        return 1;
    }
    patient_seq[bytes_read] = '\0';
    fclose(fasta);

    return 0;
}





