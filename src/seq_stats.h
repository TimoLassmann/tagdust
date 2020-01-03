
#ifndef SEQ_STATS_H
#define SEQ_STATS_H

#include "tlalphabet.h"

#include "arch_lib.h"
//#include "io.h"
struct sequence_stats_info{
        double background[5];
        double* expected_5_len;
        double* expected_3_len;
        double* mean_5_len;
        double* stdev_5_len;
        double* mean_3_len;
        double* stdev_3_len;
        double average_length;
        double mean_seq_len;
        double stdev_seq_len;
        int max_seq_len;
        int total_num_seq;
};


struct seq_stats{
        struct sequence_stats_info** ssi;
        struct alphabet* a;
        int num;
};

extern int get_sequence_stats(struct seq_stats** sequence_stats, struct arch_library* al,char** infiles,int numfiles,struct rng_state* main_rng);

//extern int get_sequence_stats(struct seq_stats** sequence_stats, struct arch_library* al,char** infiles,int numfiles);

void free_sequence_stats(struct seq_stats* si);

//extern struct sequence_stats_info* get_sequence_stats(struct parameters* param, struct read_info** ri, int file_num );

#endif
