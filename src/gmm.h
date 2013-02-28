//
//  gmm.h
//  tagdust2
//
//  Created by lassmann on 2/25/13.
//  Copyright (c) 2013 lassmann. All rights reserved.
//

#ifndef tagdust2_gmm_h
#define tagdust2_gmm_h

struct gmm_model{
	double* p;
	double* mean;
	double* sigma;
	double likelihood;
};

struct gmm_model* run_miniEM( float* x, struct gmm_model*  model, int k, int n);
void run_gmm_on_sequences(struct read_info** ri, int numseq);


#define MAX_NUM_MIXTURES 50

#endif
