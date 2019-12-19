#ifndef TLBITVEC_H
#define TLBITVEC_H


#ifdef TLBITVEC_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif


typedef struct bitvec bitvec;


EXTERN int make_bitvector( bitvec** bv, int num_elem);
EXTERN int clear_bitvector(bitvec* bv);
EXTERN int bit_set(bitvec* bv, int i);
EXTERN int bit_clr(bitvec* bv, int i);
EXTERN int bit_test(bitvec* bv, int i, int* ret);
EXTERN int free_bitvector(bitvec** bv);

#undef TLBITVEC_IMPORT
#undef EXTERN


#endif
