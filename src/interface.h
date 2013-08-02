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
#include "config.h"

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

/**
 * @brief Contains user specified read structure.
 *
 * This structure is used to hold the read architecture as specified by the user. 
 */
struct read_structure{
	char*** sequence_matrix; /**<  Three dimensional character array containing all user specified segments. */
	int* numseq_in_segment;/**<  Number of sequences in each segment. */
	char* type; /**<  Type of each segment. */
	int num_segments; /**<  Number of segments. */

};

/**
 * @brief All user specified parameters. 
 *
 * This structure is used to hold all command line parameters. 
 */
struct parameters {
	char** infile; /**<  Names of input files. */
	struct read_structure* read_structure; 
	char* outfile;/**<  Output file name. */
	char* reference_fasta;/**<  Name of fasta file containing known artifacts ti be matched against. . */
	int infiles;/**<  Number of input files. */
	int quiet_flag;
	int num_query;/**<  Number of sequences to read at one time. */
	char* format;
	char* filter;
	char* train;
	char* exact5; 
	int gzipped;
	int bzipped;
	int dust;
	int sam;
	int fasta;
	float sequencer_error_rate;/**<  Expected error rate of sequencer.  */
	float indel_frequency;/**<  Fraction of insertions and deletions among sequencer_error_rate. */
	int average_read_length;/**<  Average read length. */
	int num_threads;/**<  Number of threads. */
	float confidence_threshold;
	int matchstart;
	int matchend;
	int minlen;/**<  Minium accepted read length.  */
	int sim;
	int numbarcode;
	int filter_error; 
	
};

struct parameters* interface(struct parameters* param,int argc, char *argv[]);
struct parameters* assign_segment_sequences(struct parameters* param, char* tmp, int segment);
void free_param(struct parameters* param);
void usage(void);



