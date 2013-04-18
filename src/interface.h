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
 along with Delve.  If not, see <http://www.gnu.org/licenses/>.
 
 */


#include <stdio.h>
#include <getopt.h>
#include <stdlib.h>
#include <string.h>

#define OPT_UNMAPPED 1
#define OPT_SEG1 2
#define OPT_SEG2 3
#define OPT_SEG3 4
#define OPT_SEG4 5
#define OPT_SEG5 6
#define OPT_SEG6 7
#define OPT_SEG7 8
#define OPT_SEG8 9
#define OPT_SEG9 10
#define OPT_SEG10 11

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
	int kmer_size;
	int infiles;
	int quiet_flag;
	int num_query;
	char* print_unmapped;
	int solexa;
	int print_qual;
	int print_posteriors;
	int k_errors_allowed;
	char*  summary;
	char* format;
	char* filter;
	char* alt_lib_name;
	int gzipped;
	int sam;
	float sequencer_error_rate;
	float indel_frequency;
	int average_read_length;
};

struct parameters* interface(struct parameters* param,int argc, char *argv[]);
struct parameters* assign_segment_sequences(struct parameters* param, char* tmp, int segment);
void free_param(struct parameters* param);
void usage(void);



