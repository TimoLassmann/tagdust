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

/*! \file barcode_hmm.h
 \brief Functions build and search with user specified HMMs.
 
 Contains all functions for HMM construction, initialization, training and searching. 
 
 \author Timo Lassmann
 \bug No known bugs.
 */




#ifndef tagdust2_barcode_hmm_h
#define tagdust2_barcode_hmm_h

#ifndef MMALLOC
#include "malloc_macro.h"
#endif


#ifndef _MM_ALIGN16
#ifdef __GNUC__
#define _MM_ALIGN16 __attribute__((aligned (16)))
#endif
#ifdef __MSVC__
#define _MM_ALIGN16 __declspec(align(16))
#endif
#endif



//#define COMPARE(a, b) (((a) > (b)) - ((a) < (b)))
/** \def MM
 \brief Index of Match to Match transition.
 
 */

/** \def MI
 \brief Index of Match to Insert transition.
 
 */
/** \def MD
 \brief Index of Match to Delete transition.
 
 */
/** \def II
 \brief Index of Insert to Insert transition.
 
 */
/** \def IM
 \brief Index of Insert to Match transition.
 
 */

/** \def DD
 \brief Index of delete to Delete transition.
 
 */

/** \def DM
 \brief Index of Delete to Match  transition.
 
 */

/** \def MSKIP
 \brief Index of Match to silent transition.
 
 */
/** \def ISKIP
 \brief Index of Insert to silent transition.
 
 */

#define MM 0 
#define MI 1
#define MD 2
#define II 3
#define IM  4
#define DD 5
#define DM 6

#define MSKIP 7
#define ISKIP 8



/** \def MAX_NUM_SUB_MODELS
 \brief Maximum Number of HMM segments. 
 
 */
#define MAX_NUM_SUB_MODELS 64


/** \def MODE_GET_LABEL
 \brief Run sequence labeling in threads. 
 
 */

/** \def MODE_TRAIN
 \brief Run HMM training in threads.
 
 */

/** \def MODE_RUN_RANDOM
 \brief Shuffle Sequences and get probabilities .
 
 */

/** \def MODE_GET_PROB
 \brief Get sequence probabilities but no labels in threads. 
 
 */


#define MODE_GET_LABEL 1
#define MODE_TRAIN 2
#define MODE_RUN_RANDOM 3
#define MODE_GET_PROB 4
#define MODE_ARCH_COMP  5
#define NUM_RANDOM_SCORES 500000



/** \def EXTRACT_SUCCESS
 \brief Index for counting successful extractions. 
 
 */

/** \def EXTRACT_FAIL_BAR_FINGER_NOT_FOUND
 \brief Index for counting reads where barcode / fingerprint was not found.  
 */


/** \def EXTRACT_FAIL_READ_TOO_SHORT
 \brief Index for counting reads which are too short after extraction. 
 */

/** \def EXTRACT_FAIL_AMBIGIOUS_BARCODE
 \brief Index for counting reads can match multiple barcodes. 
 */

/** \def EXTRACT_FAIL_ARCHITECTURE_MISMATCH
 \brief Index for counting reads not matching the given read architecture. 
 */

/** \def EXTRACT_FAIL_MATCHES_ARTIFACTS
 \brief Index for counting reads matching given artifactual reads. 
 */

/** \def EXTRACT_FAIL_LOW_COMPLEXITY
 \brief Index for counting low-complexity reads.
 */

#ifndef EXTRACTION_OUTCOMES

#define EXTRACTION_OUTCOMES


#define EXTRACT_SUCCESS 0
#define EXTRACT_FAIL_BAR_FINGER_NOT_FOUND 1 
#define EXTRACT_FAIL_READ_TOO_SHORT 2
#define EXTRACT_FAIL_AMBIGIOUS_BARCODE 3
#define EXTRACT_FAIL_ARCHITECTURE_MISMATCH 4 
#define EXTRACT_FAIL_MATCHES_ARTIFACTS 5 
#define EXTRACT_FAIL_LOW_COMPLEXITY 6

#endif

/**
 @brief Stores HMM parameters and slices of the dyn. programming matrix.
 Instead of having two-dimentional dynamic programming matrices as usual, I decides to give each HMM column (Match, Delete and Insert state) it's own slice. This makes writing the actual code a litte easier because we can have pointers for current and previous column. 
 
 */
struct hmm_column{
	float* M_foward;//[MAX_HMM_SEQ_LEN]; /**< @brief  Holds forward probabilities for Match states.*/
	float* M_backward;//[MAX_HMM_SEQ_LEN]; /**<@brief  Holds backward probabilities for Match states.*/
	
	float* I_foward;//[MAX_HMM_SEQ_LEN]; /**<@brief  Holds backward probabilities for Insert states.*/
	float* I_backward;//[MAX_HMM_SEQ_LEN]; /**<@brief  Holds backward probabilities for Insert states.*/
	
	float* D_foward;//[MAX_HMM_SEQ_LEN]; /**<@brief  Holds backward probabilities for Delete states.*/
	float* D_backward;//[MAX_HMM_SEQ_LEN]; /**<@brief  Holds backward probabilities for Delete states.*/
	
	float transition[9]; /**<@brief Transition probabilities. */
	float transition_e[9];/**<@brief Estimated transition probabilities. */

	
	float m_emit[5]; /**<@brief Match emision probabilities. */
	float i_emit[5];/**<@brief Insert emision probabilities. */
	float m_emit_e[5];/**<@brief Estimated Match emision probabilities. */
	float i_emit_e[5];/**<@brief Estimated Insert  emision probabilities. */
	
	int identifier;  /**< @brief currently unused. */
}_MM_ALIGN16;



/**
 @brief Collects @ref hmm_column (s) to make a HMM. 
 
 */
struct hmm{
	struct hmm_column** hmm_column;/**<@brief Pointers to @ref hmm_column. */
	int num_columns;/**<@brief Number of columns - HMM length. */
	
}_MM_ALIGN16;



/**
 @brief Collects multiple @ref hmm(s) into a HMM segment. 
 Includes the for silent states connecting multiple segments. 
 
 */
struct model{
	struct hmm** hmms;/**< @brief Pointers to HMMs. */
	float background_nuc_frequency[5];/**<@brief Background nulceotide frequency. */
	float** silent_to_M;/**<@brief  Silent to Match transition probability.  */
	float** silent_to_I;/**<@brief  Silent to Insert transition probability.  */
	float** silent_to_M_e;/**<@brief Estimated Silent to Match transition probability.  */
	float** silent_to_I_e;/**<@brief Estimated Silent to Insert transition probability.  */
	float* silent_forward;//[MAX_HMM_SEQ_LEN]; /**<@brief Dyn. Prog. Matrix for forward silent state. */
	float* silent_backward;//[MAX_HMM_SEQ_LEN]; /**<@brief Dyn. Prog. Matrix for backward silent state. */
	float skip; /**<@brief Probability to skip segment*/
	float skip_e;/**<@brief Estimated probability to skip segment*/
	int average_length; /**<@brief Not used.... */
	int num_hmms;/**<@brief Number of HMMs in segment.*/
}_MM_ALIGN16;

/**
 @brief Collects multiple @ref model (s) into the complete HMM.
 Includes variables to hold the results.  
 */
struct model_bag{
	struct model** model;
	int num_models; /**<@brief Number of segments.*/
	float f_score;/**<@brief Probability of sequence by forward algorithm. */
	float b_score;/**< @brief Probability of sequence by backward algorithm. */
	float r_score;/**<@brief Probability of random model. */
	float bar_score;/**< @brief Max probability of paths going through one barcode HMM. */
	int** path;/**< @brief Matrix to hold viterbi path. */
	float** dyn_prog_matrix;/**<@brief Dyn. Prog. Matrix - used to find consistent max posterior path. */
	float** transition_matrix;/**<@brief Transition scores - used to find consistent max posterior path. */
	int* label; /**<@brief Hold information about HMMs.*/
	//double* random_scores;/**<@brief Holds probabilities of random / shuffled sequences. */
	//int num_random_scores;/**<@brief Number of random probabilities.*/
	
	float* previous_silent;
	float* next_silent;
	
	int average_raw_length;
	int current_dyn_length;
	//float lambda;
	//float mu;
	
	
	int total_hmm_num; /**<@brief Total number of profile HMMs in complete HMM.*/
	float model_multiplier; /**< @brief Number of different profile HMM combinations. */
}_MM_ALIGN16;


struct arch_bag{
	struct model_bag** archs;
	
	float* arch_posterior;
	char** command_line; 
	int num_arch;
};


/**
 @brief Passes data to threads. 
 */
struct thread_data{
	struct arch_bag* ab;
	struct model_bag* mb;
	struct parameters* param;
	struct read_info** ri;
	struct fasta* fasta;
	int numseq; /** Number of sequences.*/
	int start; /** @brief  Starting index of sequences for a particular thread.*/
	int end; /** @brief Endoing index of sequences for a particular thread.*/
};

/**
 @brief Collects Summary Statistics. 
 */
struct log_information{
	int total_read;
	int num_EXTRACT_SUCCESS;
	int num_EXTRACT_FAIL_BAR_FINGER_NOT_FOUND;
	int num_EXTRACT_FAIL_READ_TOO_SHORT;
	int num_EXTRACT_FAIL_AMBIGIOUS_BARCODE;
	int num_EXTRACT_FAIL_ARCHITECTURE_MISMATCH;
	int num_EXTRACT_FAIL_MATCHES_ARTIFACTS;
	int num_EXTRACT_FAIL_LOW_COMPLEXITY;
};

int hmm_controller_multiple(struct parameters* param);

//void hmm_controller_pe(struct parameters* param);
//void hmm_controller(struct parameters* param, int file_num);
//void filter_controller(struct parameters* param, int file_num);

struct model* malloc_model(int main_length, int sub_length, int number_sub_models);

void free_model(struct model* model);

//struct model* malloc_model_according_to_read_structure(struct read_structure* rs, int key);
struct model* malloc_model_according_to_read_structure(int num_hmm, int length,int dyn_length);
struct model* init_model_according_to_read_structure(struct model* model,struct parameters* param , int key, double* background,int assumed_length);
void print_model(struct model* model);


struct model_bag* forward(struct model_bag* mb, char* a, int len);
struct model_bag* backward (struct model_bag* mb, char* a, int len);
struct model_bag* forward_extract_posteriors(struct model_bag* mb, char* a, char* label, int len);
//struct model_bag* forward_extract_posteriors(struct model_bag* mb, char* a, int len);
//struct model_bag* forward_max_posterior_decoding(struct model_bag* mb, char* a, int len);
struct model_bag* forward_max_posterior_decoding(struct model_bag* mb, struct read_info* ri, char* a, int len);

struct model_bag* copy_model_bag(struct model_bag* org);
struct model_bag* init_model_bag(struct parameters* param,struct sequence_stats_info* ssi);
struct model* copy_model_parameters(struct model* org, struct model* copy );
struct model* copy_estimated_parameter(struct model* target, struct model* source );
struct model* reestimate(struct model* m, int mode);
void free_model_bag(struct model_bag* mb);

int run_pHMM(struct arch_bag* ab,struct model_bag* mb,struct read_info** ri,struct parameters* param,struct fasta* reference_fasta ,int numseq, int mode);
int run_rna_dust(struct read_info** ri,struct parameters* param,struct fasta* reference_fasta ,int numseq);

void* do_baum_welch_thread(void *threadarg);
void* do_label_thread(void *threadarg);
void* do_rna_dust(void *threadarg);
void* do_probability_estimation(void *threadarg);
void* do_arch_comparison(void *threadarg);

double pi0_bootstrap(struct read_info** ri, int numseq);
double get_min_pi0(double* x, double* y, int n_points);

int  extract_reads(struct model_bag* mb, struct parameters* param,  struct read_info* ri);
struct read_info*  make_extracted_read(struct model_bag* mb, struct parameters* param,  struct read_info* ri);
struct read_info** match_to_reference(struct thread_data *data);

struct read_info** dust_sequences(struct thread_data *data);

struct model_bag* estimate_length_distribution_of_partial_segments(struct model_bag*mb,struct read_info** ri,struct parameters* param, int numseq);


struct model_bag* estimate_model_from_labels(struct model_bag* mb, struct parameters* param,  struct read_info** ri,int numseq);
struct model_bag* set_model_e_to_laplace(struct model_bag* mb);


int emit_read_sequence(struct model_bag* mb, struct read_info* ri,int average_length,unsigned int* seed );
int emit_random_sequence(struct model_bag* mb, struct read_info* ri,int average_length,unsigned int* seed );

float model_information_content(struct model_bag*mb);

struct hmm* set_hmm_transition_parameters(struct hmm* hmm, int len,double base_error, double indel_freq,  double mean, double stdev);


int estimateQthreshold(struct parameters* param, struct sequence_stats_info* ssi);

int test_architectures(struct parameters* param, int file_num);
char* pretty_print_selected_architecture(char* command_line, char* buffer);

void print_labelled_reads(struct model_bag* mb, struct parameters* param,  struct read_info* ri);

#endif




