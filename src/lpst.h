#ifndef LPST_H
#define LPST_H

#include "pst.h"

#include "seq_stats.h"
#include "arch_lib.h"

/* linear arrangement of pst models  */

extern int lpst_score_read( struct read_structure* rs,struct tl_seq_buffer* sb, struct sequence_stats_info* si, float* score);


#endif
