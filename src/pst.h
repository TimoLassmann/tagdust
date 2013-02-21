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
	int* bit_occ;
	int occ;
	int last_seen;
	int in_T;
};

struct ranks{
	double value;
	int sample;
};

struct pst {
	struct pst_node* pst_root;
	struct pst_node* ppt_root;
	struct ranks** rank_array;
	char** suffix_array;
	int total_len;
	
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
struct pst_node* build_ppt(struct pst* pst,struct pst_node* n );


struct pst_node* count_patterns(struct read_info** ri,struct pst* pst,struct pst_node* n);
struct read_info**  scan_read_with_pst(struct read_info** ri,struct pst* pst);

struct pst_node*  count_pst_lables(struct pst_node* n, char* string,int target, int pos,int seq_id);
struct pst_node*  count_ppt_lables(struct pst_node* n, char* string,int target, int pos,int seq_id);


float get_pst_prob(struct pst_node* n, char* string,int target, int pos,int seq_id);
float get_ppt_prob(struct pst_node* n, char* string,int target, int pos,int seq_id);
int get_occ(struct pst_node* n, char* string,int target, int pos,int seq_id);

void print_pst(struct pst* pst,struct pst_node* n, struct read_info** ri );

struct pst_node* alloc_bit_occ_pst(struct pst_node* n, int num);

int* bit_set(int*a, int i);
int* bit_clr(int*a, int i);
int bit_test(int*a, int i);


typedef int (*compfn)(const void*, const void*);

int establish_rank(const void *a, const void *b);

#endif
