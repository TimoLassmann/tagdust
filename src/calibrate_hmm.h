#ifndef CALIBRATE_HMM_H
#define CALIBRATE_HMM_H

#include "tldevel.h"
#include "arch_lib.h"

#include "seq_stats.h"
#include "tlrng.h"

extern int calibrate_architectures(struct arch_library* al, struct seq_stats* si,struct rng_state* main_rng);
#endif
