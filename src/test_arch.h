#ifndef TEST_ARCH_H
#define TEST_ARCH_H



#include "tldevel.h"

#include "interface.h"
#include "arch_lib.h"
#include "seq_stats.h"

#include "read_groups.h"
//extern int test_architectures(struct cookbook* cb, struct seq_stats* si, struct parameters* param);
extern int test_architectures(struct cookbook* cb, struct read_groups* rg);

#endif
