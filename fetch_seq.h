#ifndef FETCH_SEQ_H
#define FETCH_SEQ_H

int get_gene_id(const char *gene_name, char *gene_id);
int get_wt_ntseq(char *gene_id, char cdna_seq[]);
int get_patient_seq(FILE *fasta, char patient_seq[]);

#endif
