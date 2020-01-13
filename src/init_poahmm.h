#ifndef INIT_POAHMM_H
#define INIT_POAHMM_H

extern int poahmm_from_read_structure(struct poahmm** poahmm,struct global_poahmm_param* p, struct read_structure* rs,struct alphabet* a);


//extern int init_nodes_from_read_structure(struct poahmm* poahmm, struct read_structure* rs, struct alphabet* a, int random, int plus_len);


extern int init_nodes_from_read_structure(struct poahmm* poahmm, struct read_structure* rs, struct alphabet* a, int random, int plus_min_len, int plus_max_len);


extern int set_terminal_gap_prob(struct poahmm* poahmm, int seq_len);

extern int set_len_of_unknown_poa(const struct read_structure*rs, int* min_plus_len,int* max_plus_len, int min_seq_len, int max_seq_len);


#endif
