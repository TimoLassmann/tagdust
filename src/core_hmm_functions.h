#ifndef CORE_HMM_FUNCTIONS_H
#define CORE_HMM_FUNCTIONS_H


#include "io.h"
#include "tlrng.h"

extern struct model_bag* forward(struct model_bag* mb, const uint8_t* a, int len);
extern int backward(struct model_bag* mb,const uint8_t* a,const int len);
extern int forward_max_posterior_decoding(struct model_bag* mb, struct read_info* ri, const uint8_t* a, int len);

extern struct model_bag* forward_extract_posteriors(struct model_bag* mb,const uint8_t* a, char* label, int len);


extern int  emit_random_sequence(struct model_bag* mb, struct read_info* ri,int average_length, struct rng_state* rng);
extern int emit_read_sequence(struct model_bag* mb, struct read_info* ri,int average_length, struct rng_state* rng);

extern struct model* reestimate(struct model* m, int mode);

extern struct hmm* set_hmm_transition_parameters(struct hmm* hmm, int len,double base_error, double indel_freq,  double mean, double stdev);

#endif
