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


#ifndef _MM_ALIGN16
#ifdef __GNUC__
#define _MM_ALIGN16 __attribute__((aligned (16)))
#endif
#ifdef __MSVC__
#define _MM_ALIGN16 __declspec(align(16))
#endif
#endif




#ifndef tagdust2_barcode_hmm_h
#define tagdust2_barcode_hmm_h

#define MM 0
#define MI 1
#define MD 2
#define II 3
#define IM  4
#define DD 5
#define DM 6

#define MSKIP 7
#define ISKIP 8


#define SELF 0
#define NEXT 1


#define MAX_HMM_SEQ_LEN 150

#define MAX_NUM_SUB_MODELS 64

#define MODE_GET_LABEL 1
#define MODE_TRAIN 2
#define MODE_RUN_RANDOM 3

struct hmm_column{
	float M_foward[MAX_HMM_SEQ_LEN];
	float M_backward[MAX_HMM_SEQ_LEN];
	
	float I_foward[MAX_HMM_SEQ_LEN];
	float I_backward[MAX_HMM_SEQ_LEN];
	
	float D_foward[MAX_HMM_SEQ_LEN];
	float D_backward[MAX_HMM_SEQ_LEN];
	
	float transition[9];
	float transition_e[9];
	
	float m_emit[5];
	float i_emit[5];
	float m_emit_e[5];
	float i_emit_e[5];
	
	int identifier;
}_MM_ALIGN16;

struct hmm{
	struct hmm_column** hmm_column;
	int num_columns;
	
}_MM_ALIGN16;

struct model{
	struct hmm** hmms;
	float background_nuc_frequency[5];
	float* silent_to_M;
	float* silent_to_I;
	float* silent_to_M_e;
	float* silent_to_I_e;
	float silent_forward[MAX_HMM_SEQ_LEN];
	float silent_backward[MAX_HMM_SEQ_LEN];
	float skip;
	float skip_e;
	int average_length;
	int num_hmms;
}_MM_ALIGN16;

struct model_bag{
	struct model** model;
	int num_models;
	float f_score;
	float b_score;
	float bar_score;
	int** path;
	float** dyn_prog_matrix;
	float** transition_matrix;
	int* label;
	int total_hmm_num;
	float model_multiplier;
}_MM_ALIGN16;

struct thread_data{
	struct model_bag* mb;
	struct parameters* param;
	struct read_info** ri;
	int numseq;
	int start;
	int end;
};


struct log_information{
	int probability_distribution[1001];
	
	int total_read;
	int success;
	int prob_failure;
	int len_failure;
	int arch_failure;

};

void hmm_controller(struct parameters* param,int (*fp)(struct read_info** ,struct parameters*,FILE* ),int file_num);

struct model* malloc_model(int main_length, int sub_length, int number_sub_models);

//struct model* init_model(struct model* model);
//struct model* copy_and_malloc_model(struct model* org);
//struct model* add_estimates_to_model(struct model* target, struct model* source);
void free_model(struct model* model);

//struct model* malloc_model_according_to_read_structure(struct read_structure* rs, int key);
struct model* malloc_model_according_to_read_structure(int num_hmm, int length);
struct model* init_model_according_to_read_structure(struct model* model,struct parameters* param , int key, float* background,int assumed_length);
void print_model(struct model* model);


struct model_bag* forward(struct model_bag* mb, char* a, int len);
struct model_bag* backward (struct model_bag* mb, char* a, int len);
struct model_bag* forward_extract_posteriors(struct model_bag* mb, char* a, int len);
//struct model_bag* forward_max_posterior_decoding(struct model_bag* mb, char* a, int len);
struct model_bag* forward_max_posterior_decoding(struct model_bag* mb, struct read_info* ri, char* a, int len);

struct model_bag* copy_model_bag(struct model_bag* org);
struct model_bag* init_model_bag(struct parameters* param,float* back);
struct model* copy_model_parameters(struct model* org, struct model* copy );
struct model* copy_estimated_parameter(struct model* target, struct model* source );
struct model* reestimate(struct model* m, int mode);
void free_model_bag(struct model_bag* mb);

struct model_bag* run_pHMM(struct model_bag* mb,struct read_info** ri,struct parameters* param,int numseq, int mode);
void* do_baum_welch_thread(void *threadarg);
void* do_label_thread(void *threadarg);
void* do_run_random_sequences(void *threadarg);
#endif




