#ifndef CORRECT_H
#define CORRECT_H

/* routines to correct UMI's and error in barcodes  */

#include "khash.h"

#include "tlseqio.h"

#include "base_quality.h"

//#include "assign_data.h"

#ifdef CORRECT_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif


KHASH_MAP_INIT_INT(exact, int)



//EXTERN int ref_correct(khash_t(exact) * h ,struct qsubscore* subm, struct seq_bit* sb, int q_offset);

EXTERN int fill_exact_hash(khash_t(exact) ** hash, struct tl_seq_buffer* sb);


EXTERN uint32_t seq_to_code(char* seq, int len);

/* 4 debugging, */
EXTERN void code_to_seq(uint32_t k,int len);

#undef CORRECT_IMPORT

#undef EXTERN



#endif
