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

#define MAX_PST_LEN 64


struct pst_node{
	struct pst_node* next[5];
	float nuc_probability[5];
	char* label;
};


struct pst {
	struct pst_node* root;
	char** suffix_array;
	
	float p_min;
	float alpha;
	float lamba;
	float r;
	int L;
	
	float numseq;
	
	int suffix_len;
	int current_suffix_size;
};

void pst_controller(struct parameters* param,int (*fp)(struct read_info** ,struct parameters*,FILE* ),int file_num);

struct pst_node* alloc_node(struct pst_node* n,char* string,int len);
struct pst_node* build_pst(struct pst* pst,struct pst_node* n );
void print_pst(struct pst* pst,struct pst_node* n);
#endif
