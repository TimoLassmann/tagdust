#ifndef CALIBRATE_HMM_H
#define CALIBRATE_HMM_H

#include "tldevel.h"
#include "arch_lib.h"

#include "read_groups.h"
#include "tlrng.h"

extern int calibrate_architectures(struct arch_library* al, struct read_ensembl* e,struct rng_state* main_rng);
#endif
