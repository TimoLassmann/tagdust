#ifndef SIM_SEQ_LIB_H
#define SIM_SEQ_LIB_H


#include "stdint.h"

#include "tlrng.h"

#ifdef TAGDUST_SIM_SEQ_LIB_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN int mutate_seq(char* ref,char* target,int len, float error_rate, struct rng_state* rng);
EXTERN int generate_random_seq(char** seq, int* l, struct rng_state* rng);
EXTERN int insert_seq(char* r, int r_len,char* insert, int i_len, struct rng_state* rng);
EXTERN int seq_to_internal(char* seq, int len, uint8_t** internal, int* i_len);


#undef TAGDUST_SIM_SEQ_LIB_IMPORT
#undef EXTERN

#endif
