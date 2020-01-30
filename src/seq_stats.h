
#ifndef SEQ_STATS_H
#define SEQ_STATS_H

#include "tlalphabet.h"

#include "read_groups.h"
//#include "arch_lib.h"
//#include "io.h"


/*struct seq_stats{
        struct sequence_stats_info** ssi;
        struct alphabet* a;
        int num;
        };*/

extern int get_sequence_stats(struct read_ensembl* e, struct rng_state* main_rng);
//extern int get_sequence_stats(struct read_groups* rg ,struct rng_state* main_rng);

//extern int get_sequence_stats(struct seq_stats** sequence_stats, struct arch_library* al,char** infiles,int numfiles);

//void free_sequence_stats(struct seq_stats* si);

//extern struct sequence_stats_info* get_sequence_stats(struct parameters* param, struct read_info** ri, int file_num );

#endif
