//
//  pst.h
//  tagdust2
//
//  Created by lassmann on 2/7/13.
//  Copyright (c) 2013 lassmann. All rights reserved.
//

#ifndef tagdust2_pst_h
#define tagdust2_pst_h

#define BITSPERWORD 32
#define SHIFT 5
#define MASK 0x1F



struct pst_node{
	float nuc_probability[5];
	char letter;
};


struct pst {
	struct pst_node* root;
	float* frequencies;
};

#endif
