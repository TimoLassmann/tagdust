#ifndef MINHASH_H
#define MINHASH_H

#ifdef TLMINHASH_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif

typedef struct Boolean_matrix Boolean_matrix;
typedef struct minhash minhash;

/* allocate and free Boolean matrix */
EXTERN struct Boolean_matrix* init_Bmatrix( int columns,int rows);
EXTERN void free_Boolean_matrix(struct Boolean_matrix* bm);

/* Create minhash from Boolean matrix and free */
EXTERN struct minhash* create_min_hash(struct Boolean_matrix* bm, int num_sig, long int seed);
EXTERN void free_minhash(struct minhash* min_h);

/* calculate jaccard index for set S in based on the minhash and the probability of finding set S in n samples  */

EXTERN int jaccard_sim_min_multihash(struct minhash* min_h , int* S, int n,double* jac_sim, double *p_S_in_X);

#undef TLMINHASH_IMPORT
#undef EXTERN

#endif
