#ifndef PST_H
#define PST_H

#include "tlseqio.h"

#ifdef PST_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif


struct kmer_counts;
struct pst_node;
struct pst;


EXTERN int run_build_pst(struct pst** pst, struct kmer_counts* k);
//EXTERN int run_build_pst(struct pst** pst, struct tl_seq_buffer* sb);
EXTERN int score_pst(struct pst* pst, char* seq, int len, float*r);
EXTERN int score_pst_random(char* seq, int len,float* base_p, float*r);
EXTERN int scan_read_with_pst(struct pst* pst, char* seq, int len, float* r);
EXTERN void free_pst(struct pst* p);

EXTERN int alloc_kmer_counts(struct kmer_counts** k,  int len);
EXTERN int rm_counts(struct kmer_counts* k , struct tl_seq_buffer* sb);
EXTERN int add_counts(struct kmer_counts* k , struct tl_seq_buffer* sb);
EXTERN int test_kmer_counts(struct kmer_counts* k);
EXTERN void free_kmer_counts(struct kmer_counts* k);


#undef PST_IMPORT
#undef EXTERN

#endif
