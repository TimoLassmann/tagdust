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


struct read_structure{
	char*** sequence_matrix;
	int* numseq_in_segment;
	char* type;
	int num_segments;
};

struct parameters {
	char** infile;
	struct read_structure* read_structure;
	char* outfile;
	char* reference_fasta;
	int infiles;
	int quiet_flag;
	int num_query;
	char* format;
	char* filter;
	char* train;
	char* exact5; 
	int gzipped;
	int bzipped;
	int sam;
	int fasta;
	float sequencer_error_rate;
	float indel_frequency;
	int average_read_length;
	int num_threads;
	float confidence_threshold;
	int matchstart;
	int matchend;
	int minlen;
	int sim;
	int numbarcode;
	int filter_error; 
	
};

struct parameters* interface(struct parameters* param,int argc, char *argv[]);
struct parameters* assign_segment_sequences(struct parameters* param, char* tmp, int segment);
void free_param(struct parameters* param);
void usage(void);



