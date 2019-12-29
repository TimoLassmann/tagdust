#ifndef FILTER_H
#define FILTER_H

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

extern int run_filter(struct assign_struct* as, struct ref* ref, int index, int thres);
extern int read_reference_sequences(struct ref** r, char* filename,int seed);
extern int free_ref(struct ref** r);

#endif
