
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
        char* filter_fasta;
        char* outfile;
        char* recipe;
        char* book_file;
        int filter_error;
        int num_threads;
        int num_infiles;
        int quiet_flag;
        int seed;
        int bam;

};

extern int interface(struct parameters** param,int argc, char *argv[]);

extern int free_param(struct parameters* param);



#endif
