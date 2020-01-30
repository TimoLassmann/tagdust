#ifndef READ_GROUPS_H
#define READ_GROUPS_H

//#include "seq_stats.h"

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
        int min_seq_len;
        int total_num_seq;
};


struct read_ensembl{
        struct sequence_stats_info** ssi;
        struct seq_stats* si;
        char** filenames;
        int* arch_to_read_assignment;
        int num_files;
        struct alphabet* a;

};

struct read_groups{
        struct read_ensembl** e;
        int num_groups;
        int alloc_num_groups;
};


extern int generate_read_groups(struct read_groups** rg, char** infile, int num_infiles);
extern int sort_read_groups_based_on_arch_assign(struct read_groups* rg);
extern void free_read_groups(struct read_groups* rg);
#endif
