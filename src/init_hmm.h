#ifndef INIT_HMM_H
#define INIT_HMM_H


#include "hmm.h"

#include "tldevel.h"


extern struct model* init_model_according_to_read_structure(struct model* model,struct parameters* param , int key, double* background,int assumed_length);
extern struct model* malloc_model_according_to_read_structure(int num_hmm, int length,int dyn_length);


extern void free_model(struct model* model);
#endif
