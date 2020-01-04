
#ifndef INTERFACE_H
#define INTERFACE_H


#include "tldevel.h"
//#include <stdio.h>
#include <getopt.h>
//#include <stdlib.h>
//#include <string.h>






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
        int num_segments;
        int num_infiles;
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
        int bam;
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

extern int free_param(struct parameters* param);



#endif
