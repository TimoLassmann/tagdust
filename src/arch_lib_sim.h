#ifndef ARCH_LIB_SIM_H
#define ARCH_LIB_SIM_H



#include "tlrng.h"
#include "tlalphabet.h"

#include "arch_lib.h"



int emit_from_rs(const struct read_structure* rs, struct rng_state* rng, struct alphabet* a,  uint8_t** seq, uint8_t** qual, int* len, int sim_len);


#endif
