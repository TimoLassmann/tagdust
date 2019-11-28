#ifndef THREAD_DATA_H
#define THREAD_DATA_H

#include "tldevel.h"

struct thread_data{
        struct arch_bag* ab;
        struct model_bag* mb;
        struct parameters* param;
        struct read_info** ri;
        struct fasta* fasta;
        int numseq; /** Number of sequences.*/
        int start; /** @brief  Starting index of sequences for a particular thread.*/
        int end; /** @brief Endoing index of sequences for a particular thread.*/
};

extern int alloc_thread_data(struct thread_data** td,int num_threads, int num_arch);
extern void free_thread_data(struct thread_data* thread_data);

#endif
