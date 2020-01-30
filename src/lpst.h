#ifndef LPST_H
#define LPST_H

#include "pst.h"

#include "seq_stats.h"
#include "arch_lib.h"

/* linear arrangement of pst models  */

extern int lpst_score_read(struct tl_seq_buffer* sb, struct read_structure* rs, struct sequence_stats_info* si, float* score);



#endif
