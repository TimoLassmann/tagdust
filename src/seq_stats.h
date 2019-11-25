#ifndef SEQ_STATS_H
#define SEQ_STATS_H


#include "io.h"
struct sequence_stats_info{
        double background[5];
        double expected_5_len;
        double expected_3_len;
        double mean_5_len;
        double stdev_5_len;
        double mean_3_len;
        double stdev_3_len;
        double average_length;
        int max_seq_len;
};


extern struct sequence_stats_info* get_sequence_stats(struct parameters* param, struct read_info** ri, int file_num );

#endif
