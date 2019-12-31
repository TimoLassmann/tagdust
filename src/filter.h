#ifndef FILTER_H
#define FILTER_H


#include "pst.h"
#include "assign_data.h"

#include "tlrng.h"
#include "tlseqio.h"

struct ref{
        struct alphabet* a;
        uint8_t** seq;
        int* hits;
        int* len;
        int num_seq;
};


//extern int run_filteras(struct assign_struct* as, struct ref* ref, int index,int to);
extern int run_filter_exact(struct assign_struct* as, struct ref* ref, int index, int thres);
extern int run_filter_pst(struct assign_struct* as, struct pst* pst, int index, float thres);

//extern int read_reference_sequences(struct ref** r, char* filename,int seed);
extern int init_ref(struct ref** r, struct tl_seq_buffer* sb,struct rng_state* rng);
extern int free_ref(struct ref** r);

#endif
