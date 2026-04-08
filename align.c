#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
#include "fetch_seq.h"
#include "align_seq.h"

#define MAX_SEQ_LEN 6000 // max length 6000 bp for now
#define NCBI_API_KEY "d9d8a429deb9981e28976a848dd00a3fd708"
#define MAX_JSON 4096

/* Main program the user calls:
    * Coordinates the functions in the other files to
    * perform the alignment.
 */
int main(int argc, char **argv) {
    if (argc < 3) { // expect at least 2 arguments + program name
        fprintf(stderr, "Fewer arguments than expected\n");
        return 1;
    }

    char gene_id[32];
    if (get_gene_id(argv[1], gene_id) != 0) { // null terminated in function
        fprintf(stderr, "Error: gene_id not found\n");
        return 1;
    } // else, then gene_id holds the string version of the gene id we want

    char cdna_seq[MAX_SEQ_LEN];
    if (get_wt_ntseq(gene_id, cdna_seq) != 0) {
        fprintf(stderr, "Error: cDNA sequence not found\n");
        return 1;
    }

    const int num_patients = argc - 2; // subtract program name and gene name
    char patient_seqs[num_patients][MAX_SEQ_LEN]; // to store each patient's sequence
    char combined_seqs[num_patients][2 * MAX_SEQ_LEN + 1]; // to store combined sequences

    // compute all output filenames
    char outfiles[num_patients][256];
    const char *suffix = "_aligned.fasta";
    for (int i = 0; i < num_patients; i++) {
        get_outfile_name(argv[i + 2], outfiles[i], suffix);
    }

    // store pids to match results to the right patient file
    pid_t pids[num_patients];

    for (int i = 0; i < num_patients; i++) {
        pids[i] = fork(); // fork a worker for each patient file
        if (pids[i] < 0) {
            perror("fork");
            exit(1);
        }

        if (pids[i] == 0) { // a child process
        // read patient file
            FILE *patient_fasta = fopen(argv[i + 2], "r");
            if (patient_fasta == NULL) {
                fprintf(stderr, "Error: patient file: %s\n", argv[i + 2]);
                exit(-1);
            }

            if (get_patient_seq(patient_fasta, patient_seqs[i]) != 0) { // retrieve patient sequence
                fprintf(stderr, "Error: patient sequence not found\n");
                exit(-1); // fclose already called in function
            }
            // now combine
            combine_seq(cdna_seq, patient_seqs[i], combined_seqs[i]);

            // align
            if (align(combined_seqs[i], outfiles[i]) != 0) {
                fprintf(stderr, "Error: alignment failed\n");
                exit(-1);
            } else {
                // show user which files were created
                printf("\n%s file created.\n", outfiles[i]);
                exit(0);
            }
        }
    }

    int success_count = 0;
    int j = 0;
    char success_files[num_patients][256]; // will store file names of successful alignments
    // parent should now wait for children
    for (int k = 0; k < num_patients; k++) {
        int status;
        if (waitpid(pids[k], &status, 0) == -1) {
            perror("waitpid");
            exit(1);
        }
        if (WIFEXITED(status)) { // if terminated normally
            if (WEXITSTATUS(status) == 0) {
                success_count++;
                strncpy(success_files[j], outfiles[k], 256);
                j++;
            }
        }
    }
    if (success_count == num_patients) { // all successful
        printf("\nAll alignments complete!\n");
    }
    if (success_count == 0) {
        printf("\nNo successful alignments were produced.\n");
        return 0;
    }

    // list available alignments
    printf("\nAvailable alignments:\n");
    for (int i = 0; i < success_count; i++) {
        printf("  [%d] %s\n", i + 1, success_files[i]);
    }

    // let user pick an alignment to print
    while (1) {
        printf("\nEnter a number to print that alignment (or -1 to skip): ");
        char answer[10];
        fgets(answer, sizeof(answer), stdin);
        int choice = atoi(answer);

        if (choice == -1) {
            printf("Skipped!\n");
            break;
        }
        if (choice >= 1 && choice <= success_count) {
            // fork a worker to cat the file to stdout
            int pid = fork();
            if (pid == 0) {
                char *args[] = {"cat", success_files[choice - 1], NULL};
                execvp("cat", args);
                perror("execvp cat");
                exit(1);
            }
            // parent
            int status;
            if (wait(&status) == -1) {
                perror("wait");
                exit(1);
            }
            if (!(WIFEXITED(status))) {
                fprintf(stderr, "Error: could not print alignment\n");
                // should hopefully still carry on if one worker runs into trouble
                exit(1);
            }
        } else {
            printf("Invalid alignment. Skipped!\n");
            break;
        }
    }

    // prompt user for mutations report
    char answer2[10]; // maybe avoid numbers
    while (1) {
        printf("\n Would you like a report of the mutations for each patient? [y/n]: ");
        fgets(answer2, sizeof(answer2), stdin);
        if (answer2[0] == 'n' || answer2[0] == 'N') {
            break; // exit early
        }
        if (answer2[0] == 'y' || answer2[0] == 'Y') {
            for (int i = 0; i < success_count; i++) {
                if (report_mutations(success_files[i]) != 0) { // get mutation report
                    // error message in function
                    exit(1);
                }
                // formatted in align_seq.c for an organized print
            }
            break;
        }
        fprintf(stderr, "Error: invalid answer\n"); // loop again to ask
    }
    return 0;

}





