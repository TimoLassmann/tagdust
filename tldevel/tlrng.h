#ifndef TLRNG_H
#define TLRNG_H

#include "stdint.h"

#ifdef TLRNG_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif


typedef struct rng_state rngstate;

EXTERN struct rng_state* init_rng(uint64_t seed);
EXTERN struct rng_state* init_rng_from_rng(struct rng_state* rng);

EXTERN void free_rng(struct rng_state* rng);

EXTERN double tl_random_double(struct rng_state* rng);
EXTERN int tl_random_int(struct rng_state* rng,int a);

#undef TLRNG_IMPORT
#undef EXTERN
#endif
