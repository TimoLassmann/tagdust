#ifndef CALIBRATE_HMM_H
#define CALIBRATE_HMM_H

#include "tldevel.h"
#include "arch_lib.h"

#include "seq_stats.h"

extern int calibrate_architectures(struct arch_library* al, struct seq_stats* si);

#endif
