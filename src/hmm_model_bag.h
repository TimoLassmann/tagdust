#ifndef HMM_MODEL_BAG_H
#define HMM_MODEL_BAG_H

#include "hmm.h"
#include "seq_stats.h"

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


extern struct model_bag* init_model_bag(struct parameters* param,struct sequence_stats_info* ssi);
extern void free_model_bag(struct model_bag* mb);

#endif
