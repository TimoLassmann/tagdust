#ifndef HMM_H
#define HMM_H

#define MM 0
#define MI 1
#define MD 2
#define II 3
#define IM 4
#define DD 5
#define DM 6

#define MSKIP 7
#define ISKIP 8


#ifndef _MM_ALIGN16
#ifdef __GNUC__
#define _MM_ALIGN16 __attribute__((aligned (16)))
#endif
#ifdef __MSVC__
#define _MM_ALIGN16 __declspec(align(16))
#endif
#endif


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


struct hmm{
        struct hmm_column** hmm_column;/**<@brief Pointers to @ref hmm_column. */
        int num_columns;/**<@brief Number of columns - HMM length. */

}_MM_ALIGN16;



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

#endif
