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
#include <time.h>
#include "pst.h"



void pst_controller(struct parameters* param,int (*fp)(struct read_info** ,struct parameters*,FILE* ),int file_num)
{
	char alphabet[] = "ACGTN";
	char tmp[MAX_PST_LEN+5];
	struct read_info** ri = 0;
	int i,j,c;
	int numseq;
	clock_t cStartClock;
	FILE* file = 0;
	double sum;
	
	struct pst* pst = 0;
	
	pst = malloc(sizeof(struct pst));
	pst->current_suffix_size = param->num_query* 64;
	pst->suffix_array = malloc(sizeof(char*)* pst->current_suffix_size);
	
	pst->L = MAX_PST_LEN;
	pst->alpha = 0.0f;
	pst->p_min = 0.01f;
	pst->lamba = 0.001f;
	pst->r = 1.05f;
	pst->total_len = 0;
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
		pst->total_len = 0;
		
		/* turn to normal letterss...
		 
		 */
		for(i = 0; i < numseq;i++){
			for(j = 0; j < ri[i]->len;j++){
				ri[i]->seq[j] = alphabet[(int)ri[i]->seq[j]];
				
			}
			//fprintf(stderr,"%d ",ri[i]->len);
			pst->total_len += ri[i]->len;
		}
	
		
		cStartClock = clock();
		
				
		if(pst->current_suffix_size < pst->total_len){
			pst->suffix_array = realloc(pst->suffix_array , sizeof(char*)* (pst->total_len+64));
			pst->current_suffix_size =  (pst->total_len+64);
		}
		c = 0;
		for(i = 0; i < numseq;i++){
			for(j = 0; j < ri[i]->len;j++){
				pst->suffix_array[c] = ri[i]->seq +j;
				c++;
			}
			
		}
		
		
		
		pst->suffix_len = c;
		
		//fprintf(stderr,"%d\t%d\n",c,pst->total_len);
		//exit(0);
		qsort(pst->suffix_array, pst->suffix_len, sizeof(char *), qsort_string_cmp);
		fprintf(stderr,"built SA in %4.2f seconds\n",(double)( clock() - cStartClock) / (double)CLOCKS_PER_SEC);
		cStartClock = clock();

		pst->numseq = numseq;
		//init root - removes if statement in recursion...
		sum = 0.0;
		for(i = 0;i < 5;i++){
			tmp[0] = alphabet[i];
			tmp[1] = 0;//alphabet[i];
			c = count_string(tmp,(const char**)pst->suffix_array,pst->suffix_len-1,1);
			pst->root->nuc_probability[i] = c;
			sum+= c;
		}
		for(i = 0;i < 5;i++){
			
			pst->root->nuc_probability[i] /= sum;
			//fprintf(stderr,"%c\t%f\n",alphabet[i], n->nuc_probability[i]);
		}
		
		
		pst->root = build_pst(pst,pst->root );
		fprintf(stderr,"built PST in %4.2f seconds\n",(double)( clock() - cStartClock) / (double)CLOCKS_PER_SEC);
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
	
	char tmp[MAX_PST_LEN+5];
	int i;
	int j;
	int c;
	int add;
	int len = (int) strlen(n->label);
	double sum = 0.0f;
	
	
	
	double tmp_counts_s[5];
	
	
	
	fprintf(stderr,"NODE: %s\n", n->label);
	for(i = 0;i < 5;i++){
		fprintf(stderr,"%c+%s\t%f\n",alphabet[i], n->label, n->nuc_probability[i]);
	}
	
	//step 2 test expansion
	
	//loop though letters at present node
	if(len + 1 < 32 ){
		/// search for all strings and record probabilities S+ACGT...
		/// don't search rare strings...
		/// - super mega simple ...
		
		
		for(i = 0; i < 5;i++){
			if(n->nuc_probability[i] >= pst->p_min){ /// string (suffix  + letter is not rare compared to all other letters  ) // e.g. acgC = 0.2 and not acgC = 0.0000001. hence it makes sense to "extend the suffix - test [A,C,G,T] - acgC
								
				
				
				
				
				//init longer suffix
				tmp[0] = alphabet[i];
				for(j = 1; j < len+1;j++){
					tmp[j] = n->label[j-1];
				}
				
				sum = 0.0;
				for(j = 0; j < 5;j++){
					tmp[len+1]  = alphabet[j];
					tmp[len+2] = 0;
					c = count_string(tmp,(const char**)pst->suffix_array,pst->suffix_len-1,len+2);
					tmp_counts_s[j] = c;
					sum+= c;
				}
				
				add = 0;
				for(j = 0; j < 5;j++){
					if(tmp_counts_s[j] > pst->numseq * 0.01){
						add = 1;
						break;
					}
				}
				if(add){
					// here I know that probablility of 'X' is non-neglible AND that there existsa string 'X' - [A,C,G,T,N] which is frequent in the data - hence I add...
					n->next[i] = alloc_node(pst->root,tmp,len+1);
					add = 0;
					for(j = 0; j < 5;j++){
						if((tmp_counts_s[j]/sum) / n->nuc_probability[j] >= pst->r){
							add++;
						}
						
						if((tmp_counts_s[j]/sum) / n->nuc_probability[j] <= 1.0/ pst->r){
							add++;
						}
						
						n->next[i]->nuc_probability[j] = tmp_counts_s[j]/sum;
							
					}
					if(add){
						n->next[i]->in_T = 1;
					}
					n->next[i] = build_pst(pst,n->next[i]  );

				}
			}
		}
		
		
		
		
		/*
		
		if(len){
			if(count_string(n->label,(const char**)pst->suffix_array,pst->suffix_len-1,len) < pst->numseq * 0.01){
				go = 0;
			}
		}
		if(go){
			for(i = 0; i < len;i++){
				tmp[i+1] = n->label[i];
			}
			for(i = 0;i < 5;i++){
				if(n->nuc_probability[i] >= pst->p_min){
					
					// do stuff...
					tmp[0] = alphabet[i];
					
					//loop thorugh additions to present node.
					
					////XXXX+ACGT
					
					sum_s = 0.0f;
					for(j = 0; j < 5;j++){
						tmp[len+1]  = alphabet[j];
						tmp[len+2] = 0;
						c = count_string(tmp,(const char**)pst->suffix_array,pst->suffix_len-1,len+2);
						tmp_counts_s[j] = c;
						sum_s+= c;
												
					}
					c = 0;
					
					if(sum_s){
						for(j = 0; j < 5;j++){
							//if(tmp_counts_s[j] > 1000){ // only continue with strings occuring > 1000 times
							tmp_counts_s[j] /= sum_s;
							
							
						if(tmp_counts_s[j] / n->nuc_probability[j] >= pst->r){
								//	fprintf(stderr,"%s	%s	%f	%f\n",tmp, tmp+1, tmp_counts_s[j], tmp_counts_suf_s[j] );
								c++;
							}else if(tmp_counts_s[j] / n->nuc_probability[j] <= (1.0/pst->r)){
								//	fprintf(stderr,"%s	%s	%f	%f\n",tmp, tmp+1, tmp_counts_s[j], tmp_counts_suf_s[j] );
								c++;
							}
							
							//}
							
						}
					}
					
					if(c){
						fprintf(stderr,"ADD: %c\n",alphabet[i] );
						//n->nuc_probability[i] = -1;
						tmp[0] = alphabet[i];
						for(j = 0; j < len;j++){
							tmp[j+1] = n->label[j];
						}
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
		}*/
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
	n->in_T = 0;
	
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
	//char alphabet[] = "ACGTN";
	if(strlen(n->label) > 2){
		internal = 0;
		for(i = 0;i < 5;i++){
			if(n->next[i]){
				internal++;
			}
		}
		//if(!internal){
			fprintf(stderr,"%s	%d	%d	", n->label,n->in_T, count_string(n->label,(const char**)pst->suffix_array,pst->suffix_len-1,len));
			for(i = 0;i < 5;i++){
				if(n->next[i]){
					fprintf(stderr,"%f G\t",n->nuc_probability[i]);
				}else{
					fprintf(stderr,"%f S\t",n->nuc_probability[i]);
				}
			}
			fprintf(stderr,"\n");
		//}
	}
	
	for(i = 0;i < 5;i++){
		if(n->next[i]){
			print_pst(pst,n->next[i]);
		}
	}
}

















