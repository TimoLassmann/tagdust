#ifndef PST_H
#define PST_H

#include "tlseqio.h"

#ifdef PST_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif



struct pst_node;
struct pst;


EXTERN int run_build_pst(struct pst** pst, struct tl_seq_buffer* sb);
EXTERN int scan_read_with_pst(struct pst* pst, char* seq, int len, float* r);
EXTERN void free_pst(struct pst* p);

#undef PST_IMPORT
#undef EXTERN

#endif
