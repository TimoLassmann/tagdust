
#ifndef INTERFACE_H
#define INTERFACE_H


#if HAVE_CONFIG_H
#include "config.h"
#endif



#include "tldevel.h"
#include <stdio.h>
#include <getopt.h>
#include <stdlib.h>
#include <string.h>


#define OPT_SEG1 1
#define OPT_SEG2 2
#define OPT_SEG3 3
#define OPT_SEG4 4
#define OPT_SEG5 5
#define OPT_SEG6 6
#define OPT_SEG7 7
#define OPT_SEG8 8
#define OPT_SEG9 9
#define OPT_SEG10 10
#define OPT_TRAIN 11
#define OPT_FORMAT 12
#define OPT_START 13
#define OPT_END 14
#define OPT_MINLEN 15
#define OPT_THRESHOLD 16
#define OPT_EXACT5 17
#define OPT_SIM 18
#define OPT_NUMBARCODE 19
#define OPT_FILTER_ERROR 20
#define OPT_FILTER_REFERENCE 21
#define OPT_DUST 22

#define OPT_sim_barlen 23
#define OPT_sim_barnum 24
#define OPT_sim_5seq 25
#define OPT_sim_3seq 26
#define OPT_sim_readlen 27
#define OPT_sim_readlen_mod 28
#define OPT_sim_error_rate 29
#define OPT_sim_InDel_frac 30
#define OPT_sim_numseq 31
#define OPT_sim_random_frac 32
#define OPT_sim_endloss 33

#define OPT_join_paired 34
#define OPT_split 35
#define OPT_archfile 36
#define OPT_seed 37

#define OPT_show_finger_seq 38




struct parameters {
        char** infile;
        //struct arch_library* arch_lib;
        //struct read_structure* read_structure;
        //struct read_structure** architecture_library;
        //struct read_structure** read_structures;
        //struct read_structure* read_structure_R1;
        //struct read_structure* read_structure_R2;
        char** segments;
        char* outfile;
        char* reference_fasta;
        int infiles;
        int quiet_flag;
        int num_query;
        char* format;
        char* filter;
        char* train;
        char* exact5;
        char* messages;
        char* buffer;
        //char errmsg[kslibERRBUFSIZE];
        int gzipped;
        int bzipped;
        int dust;
        int sam;
        int fasta;
        float sequencer_error_rate;
        float indel_frequency;
        int average_read_length;
        int num_threads;
        float confidence_threshold;
        float confidence_threshold_R1;
        float confidence_threshold_R2;
        float* confidence_thresholds;
        float random_prior;
        int matchstart;
        int matchend;
        int minlen;
        int sim;
        int numbarcode;
        int filter_error;
        int print_seq_finger;
        unsigned int seed;
        char* print_artifact;
        char* arch_file;

        char* log;
        int sim_barlen;
        int sim_barnum;
        char* sim_5seq;
        char* sim_3seq;
        int sim_readlen;
        int sim_readlen_mod;
        float sim_error_rate;
        float sim_InDel_frac;
        int sim_numseq;
        float sim_random_frac;
        int sim_end_loss;

        int multiread;

        int join;
        int split;
};

extern int interface(struct parameters** param,int argc, char *argv[]);

int free_param(struct parameters* param);
void usage(void);


#endif
