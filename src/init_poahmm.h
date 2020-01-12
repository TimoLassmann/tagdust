#ifndef INIT_POAHMM_H
#define INIT_POAHMM_H

extern int poahmm_from_read_structure(struct poahmm** poahmm,struct global_poahmm_param* p, struct read_structure* rs,struct alphabet* a);


//extern int init_nodes_from_read_structure(struct poahmm* poahmm, struct read_structure* rs, struct alphabet* a, int random, int plus_len);


int init_nodes_from_read_structure(struct poahmm* poahmm, struct read_structure* rs, struct alphabet* a, int random, int plus_min_len, int plus_max_len);

#endif
