#ifndef ALIGN_SEQ_H
#define ALIGN_SEQ_H

void combine_seq(char wt_seq[], char patient[], char combined[]);
char *get_outfile_name(char patient_filename[], char outfile[], const char *suffix);
int align(char combined[], char outfile[]);
int report_mutations(char *alignment);

#endif //ALIGN_SEQ_H


