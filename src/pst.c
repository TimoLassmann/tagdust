//
//  pst.c
//  tagdust2
//
//  Created by lassmann on 2/7/13.
//  Copyright (c) 2013 lassmann. All rights reserved.
//

#include <stdio.h>
#include "tagdust2.h"
#include "interface.h"
#include "io.h"
#include "misc.h"
#include "pst.h"



void pst_controller(struct parameters* param,int (*fp)(struct read_info** ,struct parameters*,FILE* ),int file_num)
{
	char alphabet[] = "ACGTN";
	struct read_info** ri = 0;
	int i,j,c;
	int numseq;
	int total_len;
	FILE* file = 0;
	
	struct pst* pst = 0;
	
	pst = malloc(sizeof(struct pst));
	pst->current_suffix_size = param->num_query* 64;
	pst->suffix_array = malloc(sizeof(char*)* pst->current_suffix_size);
	
	pst->L = MAX_PST_LEN;
	pst->alpha = 0.0f;
	pst->p_min = 0.0001f;
	pst->lamba = 0.001f;
	pst->r = 1.05f;
	pst->root = alloc_node(pst->root,"",0);
		
	ri = malloc(sizeof(struct read_info*) * param->num_query);
	
	for(i = 0; i < param->num_query;i++){
		ri[i] = malloc(sizeof(struct read_info));
		ri[i]->seq = 0;
		ri[i]->name = 0;
		ri[i]->qual = 0;
		ri[i]->len = 0;
		ri[i]->cigar = 0;
		ri[i]->md = 0;
		ri[i]->xp = 0;
		ri[i]->priors = 0;// malloc(sizeof(unsigned int)* (LIST_STORE_SIZE+1));
		ri[i]->strand = malloc(sizeof(unsigned int)* (LIST_STORE_SIZE+1));
		ri[i]->hits = malloc(sizeof(unsigned int)* (LIST_STORE_SIZE+1));
		ri[i]->identity = malloc(sizeof(float)* (LIST_STORE_SIZE+1));
		ri[i]->read_start = -1;
		ri[i]->read_end = -1;
	}
	file =  io_handler(file, file_num,param);
	
	
	while ((numseq = fp(ri, param,file)) != 0){
		total_len = 0;
		/* turn to normal letterss...
		 
		 */
		for(i = 0; i < numseq;i++){
			for(j = 0; j < ri[i]->len;j++){
				ri[i]->seq[j] = alphabet[(int)ri[i]->seq[j]];
			}
		}
	
		
				
		
		if(pst->current_suffix_size < total_len){
			pst->suffix_array = realloc(pst->suffix_array , sizeof(char*)* (total_len+64));
			pst->current_suffix_size =  (total_len+64);
		}
		c = 0;
		for(i = 0; i < numseq;i++){
			for(j = 0; j < ri[i]->len;j++){
				pst->suffix_array[c] = ri[i]->seq +j;
				c++;
			}
			total_len += ri[i]->len;
		}
		pst->suffix_len = c;
		qsort(pst->suffix_array, pst->suffix_len, sizeof(char *), qsort_string_cmp);
		
		pst->numseq = numseq;
		pst->root = build_pst(pst,pst->root );
		print_pst(pst,pst->root);
		
		exit(0);
		//AAAAAAAAAAAAAAA (5).
		/*int a,b;
		for(i = 0; i < 1000000;i++){
		
			a =  count_string("AAAAAAAA",(const char**)pst->suffix_array,pst->suffix_len-1);
			a = binsearch_down("AAAAAAAA",(const char**)pst->suffix_array,pst->suffix_len-1);
			b = binsearch_up("AAAAAAAA",(const char**)pst->suffix_array,pst->suffix_len-1);
		}
		fprintf(stderr,"down:%d\tup:%d	count:%d\n", a,b, b-a);
		for(i = a-1; i < b+1;i++){
			fprintf(stdout,"%d	%s\n",i, pst->suffix_array[i] );
		}*/
		
		
	}
	
	
	for(i = 0; i < param->num_query;i++){
		free(ri[i]->strand);
		free(ri[i]->hits);
		free(ri[i]->identity);
		
		if(ri[i]->cigar){
			free(ri[i]->cigar);
		}
		if(ri[i]->md){
			free(ri[i]->md);
		}
		
		free(ri[i]);
	}
	free(ri);
	//fprintf(stderr,"%p\n",file);
	if(param->sam == 2 || param->sam == 1 || param->gzipped ){
		pclose(file);
	}else{
		//if(file_num != -1){
		fclose(file);
		//}
	}
	
}


struct pst_node* build_pst(struct pst* pst,struct pst_node* n )
{
	char alphabet[] = "ACGTN";
	
	char tmp[MAX_PST_LEN];
	int i;
	int j;
	int c;
	int len = (int) strlen(n->label);
	float sum = 0.0f;
	float sum_s = 0.0f;
	float sum_suf_s = 0.0f;
	
	
	float tmp_counts_s[5];
	
	float tmp_counts_suf_s[5];
	
	
	
	
	if(!len){
	//step 1 count.....
		for(i = 0; i < len;i++){
			tmp[i] = n->label[i];
		}
		for(i = 0;i < 5;i++){
			tmp[len] = alphabet[i];
			tmp[len+1] = 0;//alphabet[i];
			c = count_string(tmp,(const char**)pst->suffix_array,pst->suffix_len-1,len+1);
			n->nuc_probability[i] = c;
			sum+= c;
		}
		for(i = 0;i < 5;i++){
		
			n->nuc_probability[i] /= sum;
			//fprintf(stderr,"%c\t%f\n",alphabet[i], n->nuc_probability[i]);
		}
	}
	
	
	fprintf(stderr,"NODE: %s\n", n->label);
	for(i = 0;i < 5;i++){
		fprintf(stderr,"%s+%c\t%f\n", n->label, alphabet[i],n->nuc_probability[i]);
	}
	
	//step 2 test expansion
	
	//loop though letters at present node
	int go = 1;
	if(len + 1 < 32 ){
		if(len){
			if(count_string(n->label,(const char**)pst->suffix_array,pst->suffix_len-1,len) < pst->numseq * 0.1){
				go = 0;
			}
		}
		if(go){
			for(i = 0; i < len;i++){
				tmp[i] = n->label[i];
			}
			for(i = 0;i < 5;i++){
				if(n->nuc_probability[i] >= pst->p_min){
					
					// do stuff...
					tmp[len] = alphabet[i];
					
					//loop thorugh additions to present node.
					
					////XXXX+ACGT
					
					sum_s = 0.0f;
					sum_suf_s = 0.0f;
					for(j = 0; j < 5;j++){
						tmp[len+1]  = alphabet[j];
						tmp[len+2] = 0;
						c = count_string(tmp,(const char**)pst->suffix_array,pst->suffix_len-1,len+2);
						tmp_counts_s[j] = c;
						sum_s+= c;
						c = count_string(tmp+1,(const char**)pst->suffix_array,pst->suffix_len-1,len+1);
						tmp_counts_suf_s[j] = c;
						sum_suf_s+= c;
						
					}
					c = 0;
					
					if(sum_s){
						for(j = 0; j < 5;j++){
							//if(tmp_counts_s[j] > 1000){ // only continue with strings occuring > 1000 times
							tmp_counts_s[j] /= sum_s;
							tmp_counts_suf_s[j] /= sum_suf_s;
							
							if(tmp_counts_s[j] > 0.95f){ // this looks part of a longer string - extend always
								c++;
							}else if(tmp_counts_s[j] / tmp_counts_suf_s[j] >= pst->r){
								//	fprintf(stderr,"%s	%s	%f	%f\n",tmp, tmp+1, tmp_counts_s[j], tmp_counts_suf_s[j] );
								c++;
							}else if(tmp_counts_s[j] / tmp_counts_suf_s[j] <= (1.0/pst->r)){
								//	fprintf(stderr,"%s	%s	%f	%f\n",tmp, tmp+1, tmp_counts_s[j], tmp_counts_suf_s[j] );
								c++;
							}
							
							//}
							
						}
					}
					
					if(c){
						fprintf(stderr,"ADD: %c\n",alphabet[i] );
						n->nuc_probability[i] = -1;
						tmp[len] = alphabet[i];
						tmp[len+1] = 0;
						
						n->next[i] = alloc_node(pst->root,tmp,len+1);
						for(j = 0; j < 5;j++){
							n->next[i]->nuc_probability[j] = tmp_counts_s[j];
						}
						n->next[i] = build_pst(pst,n->next[i]  );
					}
					//decide to add;
					
					
					
					//alphabet[i];
					
					
					
				}
			}
		}
	}
	
	
	
	return n;
}


struct pst_node* alloc_node(struct pst_node* n,char* string,int len)
{
	int i;
	n = malloc(sizeof(struct pst_node));
	n->label = malloc(sizeof(char) *(len+1));
	for(i = 0; i < len;i++){
		n->label[i] =string[i];
	}
	n->label[len] = 0;
	
	
	for(i =0; i < 5;i++){
		n->next[i] = 0;
		n->nuc_probability[i] = 0.2f;
	}

	return n;
}

void print_pst(struct pst* pst,struct pst_node* n)
{
	int i;
	int internal;
	int len = (int)strlen(n->label);
	char alphabet[] = "ACGTN";
	if(strlen(n->label) > 2){
		internal = 0;
		for(i = 0;i < 5;i++){
			if(n->next[i]){
				internal++;
			}
		}
		if(!internal){
			fprintf(stderr,"%s	%d	%f	%d\n", n->label, count_string(n->label,(const char**)pst->suffix_array,pst->suffix_len-1,len) ,count_string(n->label,(const char**)pst->suffix_array,pst->suffix_len-1,len) / pst->numseq * 100,(int) pst->numseq);
			for(i = 0;i < 5;i++){
				//	fprintf(stderr,"%s+%c\t%f\n", n->label, alphabet[i],n->nuc_probability[i]);
			}
		}
	}
	
	for(i = 0;i < 5;i++){
		if(n->next[i]){
			print_pst(pst,n->next[i]);
		}
	}
}

















