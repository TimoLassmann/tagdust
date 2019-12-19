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

#ifndef IO_H
#define IO_H


#include "tldevel.h"
#include "interface.h"
#include "assign_data.h"

#define LIST_STORE_SIZE 1

#define SEEK_START 0
#define SEEK_END 2


#ifndef EXTRACTION_OUTCOMES

#define EXTRACTION_OUTCOMES

#define EXTRACT_SUCCESS 0
#define EXTRACT_FAIL_ARCHITECTURE_MISMATCH 1
#define EXTRACT_FAIL_READ_TOO_SHORT 2
#define EXTRACT_FAIL_BAR_FINGER_NOT_FOUND 3

#define EXTRACT_FAIL_MATCHES_ARTIFACTS 5
#define EXTRACT_FAIL_LOW_COMPLEXITY 6





#endif

/**
 * @brief Stores reads from SAM/ fastq formatted files.
 *
 *
 */
struct read_info{
        char* name;/**<  @brief  Name of read.*/
        char* qual;/**<  @brief Base qualities. */
        uint8_t* seq;/**<  @brief Sequence.*/
        char* labels;/**<  @brief Labeling according to HMM.*/
        //unsigned int* strand;
        //unsigned int* hits;
        //char* barcode_string;
        float mapq;/**<  @brief Mapping Quality.*/
        //double prob;/**<  @brief Quality of read.*/
        double bar_prob;/**< @brief Ambiguity */
        int malloc_len;
        int q_len;
        int len;/**<  @brief Sequence length.*/
        int read_type;
        //int barcode;
        //int fingerprint;
};

struct read_info_buffer{
        struct read_info** ri;
        char* left_over_name;
        int num_alloc;
        int num_seq;
        int offset;
};


struct file_handler{
        FILE* f_ptr;
        int gzipped;
        int bzipped;
        int sam;
        int fasta;
};

/**
 * @brief Used to store info needed to initialize HMM.
 *
 */

//#include "barcode_hmm.h"


int io_handler(struct file_handler** fh, char*filename);
//FILE* io_handler(FILE* file, int file_num,struct parameters* param);
void print_seq(struct read_info* ri,FILE* out);
int read_sam_chunk(struct read_info_buffer* rb, struct file_handler* f_handle);
//int read_sam_chunk(struct read_info** ri,struct parameters* param,FILE* file,int* buffer_count);
//int read_sam_chunk(struct read_info** ri,struct parameters* param,FILE* file);

int read_fasta_fastq(struct read_info_buffer* rb, struct file_handler* f_handle);
//int read_fasta_fastq(struct read_info** ri,struct parameters* param,FILE *file,int* buffer_count);
void print_sequence(struct read_info* ri,FILE* out);
//int print_trimmed_sequence(struct model_bag* mb, struct parameters* param,  struct read_info* ri,FILE* out);
//int qsort_ri_mapq_compare(const void *a, const void *b);

int qsort_ri_barcode_compare(const void *a, const void *b);

//struct sequence_stats_info* get_sequence_stats(struct parameters* param, struct read_info** ri, int file_num );

void concatenate_reads(struct parameters* param,int (*fp)(struct read_info** ,struct parameters*,FILE* ));
void split(struct parameters* param,int (*fp)(struct read_info** ,struct parameters*,FILE* ));
int file_exists (char * name);



//void print_split_files(struct parameters* param, struct read_info** ri, int numseq);

int alloc_read_info_buffer(struct read_info_buffer** rb, int size);
void free_read_info_buffer(struct read_info_buffer* rb);

struct read_info** malloc_read_info(struct read_info** ri, int numseq);
struct read_info** clear_read_info(struct read_info** ri, int numseq);
void free_read_info(struct read_info** ri, int numseq);

extern int compare_read_names(char* name1, char* name2);
int check_for_existing_demultiplexed_files(struct parameters* param);

int check_for_existing_demultiplexed_files_multiple(struct parameters* param, int num_reads);

FILE* open_file(struct parameters* param, char* buffer, char* mode);

int print_all(struct read_info*** read_info_container,struct parameters* param, int numseq, char*  read_present);

extern int write_all(const struct assign_struct* as,char* prefix);

#endif
