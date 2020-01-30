#ifndef EXTRACT_READS_H
#define EXTRACT_READS_H

#include "tldevel.h"

#include "arch_lib.h"
#include "seq_stats.h"
#include "read_groups.h"
#include "interface.h"

extern int extract_reads(struct arch_library* al, struct read_groups* rg,struct parameters* param,struct rng_state* rng);
#endif
