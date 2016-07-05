/*
 
 Copyright (C) 2013 Timo Lassmann <timolassmann@gmail.com>
 
 This file is part of TagDust.
 
 TagDust is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 TagDust is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with Tagdust.  If not, see <http://www.gnu.org/licenses/>.
 
 */

/*! \file interface.h 
 \brief 
 
 Initializes nucleotide alphabet needed to parse input. Calls parameter parser. Calls functions to process the data. \author Timo Lassmann \bug No known bugs.
 */


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


/**
 * @brief Contains user specified read structure.
 *
 * This structure is used to hold the read architecture as specified by the user. 
 */
struct read_structure{
	char*** sequence_matrix; /**<  @brief Three dimensional character array containing all user specified segments. */
	int* numseq_in_segment;/**< @brief Number of sequences in each segment. */
	char* type; /**<  @brief Type of each segment. */
	int num_segments; /**< @brief Number of segments. */

};

/**
 * @brief All user specified parameters. 
 *
 * This structure is used to hold all command line parameters. 
 */
struct parameters {
	char** infile; /**< @brief Names of input files. */
	struct read_structure* read_structure;
	struct read_structure** read_structures;
	struct read_structure* read_structure_R1;
	struct read_structure* read_structure_R2;
	char* outfile;/**< @brief Output file name. */
	char* reference_fasta;/**< @brief Name of fasta file containing known artifacts to be matched against. . */
	int infiles;/**<  @brief Number of input files. */
	int quiet_flag;
	int num_query;/**< @brief Number of sequences to read at one time. */
	char* format;
	char* filter;
	char* train;
	char* exact5;
	char* messages;
	char* buffer;
	char errmsg[kslibERRBUFSIZE];
	int gzipped;
	int bzipped;
	int dust;
	int sam;
	int fasta;
	float sequencer_error_rate;/**<  @brief Expected error rate of sequencer.  */
	float indel_frequency;/**< @brief Fraction of insertions and deletions among sequencer_error_rate. */
	int average_read_length;/**< @brief Average read length. */
	int num_threads;/**< @brief Number of threads. */
	float confidence_threshold;/**< @brief This threshold is used to determine whether the read matched the HMM. */
	float confidence_threshold_R1;
	float confidence_threshold_R2;
	float* confidence_thresholds; 
	float random_prior;
	int matchstart;
	int matchend;
	int minlen;/**< @brief Minium accepted read length.  */
	int sim;
	int numbarcode;
	int filter_error;
	int print_seq_finger;
	unsigned int seed;
	char* print_artifact;
	char* arch_file;
	
	char* log;/**< @brief Directory where log files are written.  */
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

struct parameters* interface(struct parameters* param,int argc, char *argv[]);
//struct read_structure* assign_segment_sequences(struct read_structure* read_structure, char* tmp, int segment);
int assign_segment_sequences(struct parameters* param, char* tmp, int segment);

int free_param(struct parameters* param);
void usage(void);


struct read_structure* malloc_read_structure(void);
void free_read_structure(struct read_structure* read_structure);
int QC_read_structure(struct parameters* param);


