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

/*! \file io.h
 \brief functions for reading sequences.
 
 */

/** \def LIST_STORE_SIZE
 \brief Sets maximum number of read mappings. 
 */
#define LIST_STORE_SIZE 1

#define SEEK_START 0
#define SEEK_END 2
/**
 * @brief Deals with fasta files. 
 *
 *  Stores sequences in one big string. 
 */
struct fasta{
	unsigned char** sn;/**<  @brief Sequence names.*/
	unsigned char* string;/**< @brief Holds sequence information.*/
	int* mer_hash;
	int* boost;
	int* s_index;
	int* suffix;
	int numseq;/**< @brief  Number of sequences.*/
	int max_len;
	int string_len;
};

/**
 * @brief Stores reads from SAM/ fastq formatted files.  
 *
 *  
 */
struct read_info{
	char* name;/**<  @brief  Name of read.*/
	char* qual;/**<  @brief Base qualities. */
	char* seq;/**<  @brief Sequence.*/
	char* labels;/**<  @brief Labeling according to HMM.*/
	unsigned int* strand;
	unsigned int* hits;
	char* cigar;
	char* md;
	float mapq;/**<  @brief Mapping Quality.*/
	double prob;/**<  @brief Quality of read.*/
	double bar_prob;/**< @brief Ambiguity */
	int len;/**<  @brief Sequence length.*/
	int errors;
};



/**
 * @brief Used to store info needed to initialize HMM. 
 *
 *
 */
struct sequence_stats_info{
	double background[5];
	int expected_5_len;
	int expected_3_len;
	double mean_5_len;
	double stdev_5_len;
	double mean_3_len;
	double stdev_3_len;
	double average_length;
	int max_seq_len;
	
};

#include "barcode_hmm.h"


FILE* io_handler(FILE* file, int file_num,struct parameters* param);
void print_seq(struct read_info* ri,FILE* out);
int read_sam_chunk(struct read_info** ri,struct parameters* param,FILE* file);
int read_fasta_fastq(struct read_info** ri,struct parameters* param,FILE *file);
void print_sequence(struct read_info* ri,FILE* out);
int print_trimmed_sequence(struct model_bag* mb, struct parameters* param,  struct read_info* ri,FILE* out);
int qsort_ri_prob_compare(const void *a, const void *b);
int qsort_ri_mapq_compare(const void *a, const void *b);

struct sequence_stats_info* get_sequence_stats(struct parameters* param, struct read_info** ri, int file_num );

void concatenate_reads(struct parameters* param,int (*fp)(struct read_info** ,struct parameters*,FILE* ));
void split(struct parameters* param,int (*fp)(struct read_info** ,struct parameters*,FILE* ));
int file_exists (char * name);


struct fasta* read_fasta(struct fasta* f);
struct fasta* get_fasta(struct fasta* p,char *infile);
unsigned char* get_input_into_string(unsigned char* string,char* infile);
void free_fasta(struct fasta*f);




