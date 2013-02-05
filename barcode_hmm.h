//
//  barcode_hmm.h
//  tagdust2
//
//  Created by lassmann on 2/5/13.
//  Copyright (c) 2013 lassmann. All rights reserved.
//

#ifndef tagdust2_barcode_hmm_h
#define tagdust2_barcode_hmm_h

#define SELF 0
#define NEXT 1
#define EMIT_A 2
#define EMIT_C 3
#define EMIT_G 4
#define EMIT_T 5
#define EMIT_N 6
#define TRANS_START 7

#define MAX_HMM_SEQ_LEN 250


struct hmm_state{
	float prob[64+2+5];
	float e_prob[64+2+5];
	float foward[MAX_HMM_SEQ_LEN];
	float backward[MAX_HMM_SEQ_LEN];
	
	
	int identifier;
};

struct hmm{
	struct hmm_state* states;
	int len;
	
};


#endif
