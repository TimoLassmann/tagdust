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

#define LIST_STORE_SIZE 1

#include "barcode_hmm.h"

struct read_info{
	char* name;
	char* qual;
	char* seq;
	char* labels;
	unsigned int* strand;
	unsigned int* hits;
	char* cigar;
	char* md;
	//char* xp;
	float mapq;
	float prob;
	float bar_prob; 
	//float* priors;
	//float* identity;
	//int read_start;
	//int read_end;
	int len;
	int errors;
};


FILE* io_handler(FILE* file, int file_num,struct parameters* param);
void print_seq(struct read_info* ri,FILE* out);
int read_sam_chunk(struct read_info** ri,struct parameters* param,FILE* file);
int read_fasta_fastq(struct read_info** ri,struct parameters* param,FILE *file);

int print_trimmed_sequence(struct model_bag* mb, struct parameters* param,  struct read_info* ri,FILE* out);
