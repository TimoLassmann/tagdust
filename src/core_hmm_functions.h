#ifndef CORE_HMM_FUNCTIONS_H
#define CORE_HMM_FUNCTIONS_H

#include "barcode_hmm.h"

extern struct model_bag* forward(struct model_bag* mb, char* a, int len);
extern struct model_bag* backward(struct model_bag* mb, char* a, int len);

extern struct model_bag* forward_extract_posteriors(struct model_bag* mb, char* a, char* label, int len);
extern struct model_bag* forward_max_posterior_decoding(struct model_bag* mb, struct read_info* ri, char* a, int len);


extern struct model* reestimate(struct model* m, int mode);

extern struct hmm* set_hmm_transition_parameters(struct hmm* hmm, int len,double base_error, double indel_freq,  double mean, double stdev);

#endif
