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
#include "nuc_code.h"
#include <time.h>
#include "pst.h"
#include <xmmintrin.h>
#include "gmm.h"
#include "dbscan.h"
#include <float.h>


void pst_tree(struct parameters* param,int (*fp)(struct read_info** ,struct parameters*,FILE* ),int file_num)
{
	char alphabet[] = "ACGTN";
	struct read_info** ri = 0;
	int i,j,c;
	int numseq;
	int* sample_list = 0;
	
	
	struct pst* pst = 0;
	
	
	FILE* file = 0;
	
	//read in sequences
	ri = malloc(sizeof(struct read_info*) * param->num_query);
	assert(ri != 0);
	
	for(i = 0; i < param->num_query;i++){
		ri[i] = malloc(sizeof(struct read_info));
		ri[i]->seq = 0;
		ri[i]->name = 0;
		ri[i]->qual = 0;
		ri[i]->len = 0;
		ri[i]->cigar = 0;
		ri[i]->md = 0;
		//ri[i]->xp = 0;
		ri[i]->strand = malloc(sizeof(unsigned int)* (LIST_STORE_SIZE+1));
		ri[i]->hits = malloc(sizeof(unsigned int)* (LIST_STORE_SIZE+1));
		//ri[i]->read_start = -1;
		//ri[i]->read_end = -1;
	}
	
	
	file =  io_handler(file, file_num,param);
	
	//struct pst_node** all_patterns = 0;
	
	//init sample_lists
	numseq = fp(ri, param,file);
	
	if(!numseq){
		fprintf(stderr,"ERROR - no sequences could be read....\n");
		exit(-1);
	}
	
	fprintf(stderr,"READ %d sequences\n", numseq);
	
	
	
	
	sample_list = malloc(sizeof(int) * numseq);
	for(i = 0; i < numseq;i++){
		sample_list[i] = 1;
		for(j = 0; j < ri[i]->len;j++){
			ri[i]->seq[j] = alphabet[(int)ri[i]->seq[j]];
		}
	}
	pst = malloc(sizeof(struct pst));
	pst->current_suffix_size = numseq* 64;
	pst->suffix_array = malloc(sizeof(char*)* pst->current_suffix_size);
	
	pst->L = MAX_PST_LEN;
	pst->alpha = 0.0f;
	pst->p_min = 0.0001f;
	pst->lamba = 0.001f;
	pst->r = 1.05f;
	pst->total_len = 0;
	pst->pst_root = alloc_node(pst->pst_root,"",0);
	//pst->ppt_root = alloc_node(pst->ppt_root,"",0);
	pst->rank_array = 0;
	
	for(i = 0; i < numseq;i++){
		pst->total_len += ri[i]->len;
	}
	
	
	//cStartClock = clock();
	if(pst->current_suffix_size < pst->total_len){
		pst->suffix_array = realloc(pst->suffix_array , sizeof(char*)* (pst->total_len+64));
		pst->current_suffix_size =  (pst->total_len+64);
	}
	
	pst->sn = malloc(sizeof(struct suffix_node* ) * pst->total_len);
	for(i = 0; i < pst->total_len;i++){
		pst->sn[i] = malloc(sizeof(struct suffix_node));
		pst->sn[i]->seq_id = -1;
		pst->sn[i]->string = 0;
	}
	
	c = 0;
	pst->mean_length = 0.0;
	for(i = 0; i < numseq;i++){
	
			for(j = 0; j < ri[i]->len;j++){//ri[i]->len;j++){
				pst->suffix_array[c] = ri[i]->seq +j;
				pst->sn[c]->seq_id = i;
				pst->sn[c]->string =  ri[i]->seq +j;
				c++;
			}
			pst->mean_length +=  ri[i]->len;
		
	}
	
	pst->mean_length /= (float)numseq;
	
	pst->suffix_len = c;
	pst->numseq = numseq;
	pst->seq_id_in_suffix = malloc(sizeof(int) * pst->suffix_len );
	
	qsort(pst->sn, pst->suffix_len, sizeof(struct suffix_node *), qsort_suffix_node_string_cmp);
	for(i = 0; i <  pst->suffix_len;i++){
		pst->suffix_array[i] = pst->sn[i]->string;
		pst->seq_id_in_suffix[i] = pst->sn[i]->seq_id;
	}
	
	for(i = 0; i < pst->total_len;i++){
		free(pst->sn[i]);// = malloc(sizeof(struct suffix_node));
	}
	free(pst->sn);

	
	
	// call recursive splitting function...
	
	pst_based_partition(pst,ri,sample_list,numseq,numseq);
	
	
	//end - results are not kept in memory
	
	
	for(i = 0; i < param->num_query;i++){
		free(ri[i]->strand);
		free(ri[i]->hits);
		if(ri[i]->name){
			free(ri[i]->name);
		}
		if(ri[i]->seq){
			free(ri[i]->seq);
		}
		if(ri[i]->qual){
			free(ri[i]->qual);
		}
		
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
	
	//free(sample_list);
	
}

void pst_based_partition(struct pst* pst ,struct read_info** ri, int* samples, int numseq,int active)
{
	//struct pst* pst = 0;
	
	char* seq;
	char alphabet[] = "ACGTN";
	struct pst_node** all_patterns = 0;
	char tmp[MAX_PST_LEN+5];
	double sum;

	int i,j,c;
	int a,b;
	
	int num_patterns;
	float max_reduction = 0.0f;
	
	//pst = malloc(sizeof(struct pst));
	//pst->current_suffix_size = numseq* 64;
	//pst->suffix_array = malloc(sizeof(char*)* pst->current_suffix_size);
	
	//pst->L = MAX_PST_LEN;
	//pst->alpha = 0.0f;
	//pst->p_min = 0.0001f;
	//pst->lamba = 0.001f;
	//pst->r = 1.05f;
	//pst->total_len = 0;
	pst->pst_root = alloc_node(pst->pst_root,"",0);
	//pst->ppt_root = alloc_node(pst->ppt_root,"",0);
	//pst->rank_array = 0;
	
	//for(i = 0; i < numseq;i++){
	//	if(samples[i]){
		//ri[i]->seq[10] = 0;
		//fprintf(stderr,"%d ",ri[i]->len);
	//		pst->total_len += ri[i]->len;
	//	}
	//}
	
	
	//cStartClock = clock();
	//if(pst->current_suffix_size < pst->total_len){
	//	pst->suffix_array = realloc(pst->suffix_array , sizeof(char*)* (pst->total_len+64));
	//	pst->current_suffix_size =  (pst->total_len+64);
	//}
	
	//pst->sn = malloc(sizeof(struct suffix_node* ) * pst->total_len);
	//for(i = 0; i < pst->total_len;i++){
	//	pst->sn[i] = malloc(sizeof(struct suffix_node));
	//	pst->sn[i]->seq_id = -1;
	//	pst->sn[i]->string = 0;
	//}
	
	//c = 0;
	//pst->mean_length = 0.0;
	//for(i = 0; i < numseq;i++){
	//	if(samples[i]){
	//		for(j = 0; j < ri[i]->len;j++){//ri[i]->len;j++){
	//			pst->suffix_array[c] = ri[i]->seq +j;
	//			pst->sn[c]->seq_id = i;
	//			pst->sn[c]->string =  ri[i]->seq +j;
	//			c++;
	//		}
	//		pst->mean_length +=  ri[i]->len;
	//	}
	//}
	
	//pst->mean_length /= (float)numseq;
	
	//pst->suffix_len = c;
	//pst->numseq = numseq;
	//pst->seq_id_in_suffix = malloc(sizeof(int) * pst->suffix_len );
	
	//qsort(pst->sn, pst->suffix_len, sizeof(struct suffix_node *), qsort_suffix_node_string_cmp);
	//for(i = 0; i <  pst->suffix_len;i++){
	//	pst->suffix_array[i] = pst->sn[i]->string;
	//	pst->seq_id_in_suffix[i] = pst->sn[i]->seq_id;
	//}
	
	//for(i = 0; i < pst->total_len;i++){
	//	free(pst->sn[i]);// = malloc(sizeof(struct suffix_node));
	//}
	//free(pst->sn);

	
	//for(i = 0; i < 10000;i++){
	//	fprintf(stderr,"%d	%s\n",pst->sn[i]->seq_id,pst->sn[i]->string  );
	//}
	//exit(0);
	
	//qsort(pst->suffix_array, pst->suffix_len, sizeof(char *), qsort_string_cmp);
	
	//init root - removes if statement in recursion...
	sum = 0.0;
	for(i = 0;i < 5;i++){
		tmp[0] = alphabet[i];
		tmp[1] = 0;//alphabet[i];
		c = count_string(tmp,(const char**)pst->suffix_array,pst->suffix_len-1,1);
		pst->pst_root->nuc_probability[i] = c;
		//pst->ppt_root->nuc_probability[i] = c;
		sum+= c;
	}
	for(i = 0;i < 5;i++){
		
		pst->pst_root->nuc_probability[i] =  pst->pst_root->nuc_probability[i]/ sum;
		//pst->ppt_root->nuc_probability[i] =  pst->ppt_root->nuc_probability[i]/ sum;
		//fprintf(stderr,"%c\t%f\n",alphabet[i], n->nuc_probability[i]);
	}
	
	
	pst->pst_root = build_pst(pst,pst->pst_root );
	
	pst->pst_root = alloc_bit_occ_pst(pst->pst_root , numseq);

	//ri =  scan_read_with_pst( ri, pst);

	for(i = 0; i < pst->numseq;i++){
		if(samples[i]){
			seq = ri[i]->seq;
			for(j = 0; j < ri[i]->len; j++ ){
				pst->pst_root = count_pst_lables(pst->pst_root, seq,  j, i);
			}
		}
	}
	
	num_patterns = 0;
	num_patterns = count_patterns(pst->pst_root, num_patterns);
	//num_patterns = count_patterns(pst->ppt_root,num_patterns);
	//
	//exit(0);
	all_patterns = malloc(sizeof(struct pst_node*)  *num_patterns );
	
	num_patterns = 0;
	num_patterns = add_patterns(all_patterns,pst->pst_root, num_patterns);
	
	
	qsort((void *)  all_patterns, num_patterns, sizeof(struct pst_node* ),(compfn) sort_pst_nodel_according_to_label);
	
	fprintf(stderr,"%d NUMPATTERNS\n", num_patterns );
	
	float left,right,pl1,pl2,pr1,pr2,reduction;
	for(i =0 ; i < num_patterns-1;i++){
		for(j = i+1 ; j < num_patterns;j++){
			//if(i != j){
				//left = (all_patterns[i]->occ) / (float)numseq;
				//right = ((float)numseq - all_patterns[i]->occ)  /(float)numseq ;
				
				
				//1n's left - in region covered by i
				pl1 = intersection(all_patterns[i]->bit_occ, all_patterns[j]->bit_occ ,1+ numseq  / BITSPERWORD);
				
				//if(i ==j){
				//	fprintf(stderr,"%d	%d	%f	%d\n",i,j,pl1, all_patterns[i]->occ );
				//}
				
				//1n's right - remaining ...    - in region covered by i
				pr1 =  all_patterns[j]->occ - pl1;
				
				pl1 /= (float) all_patterns[i]->occ;
				pl2 = 1.0 - pl1;
				
				pr1 /=(float) (active -  all_patterns[i]->occ);
				
				pr2 = 1.0 -pr1;
				
				left = (all_patterns[j]->occ) / (float)active;
				right =1.0 - left;// ((float)numseq - all_patterns[j]->occ)  /(float)numseq ;
				
				reduction = -1.0 * (left * log2(left) + (right) * log2(right));
				
				left = (all_patterns[i]->occ) / (float)active;
				right = ((float)active - all_patterns[i]->occ)  /(float)active ;
				
				
				//fprintf(stderr,"%d		%f	%f	%f	%f	%f	%f	%f	",i,left,right,pl1,pl2,pr1,pr2,reduction);
				
				if(pl1 == 0 || pl2 == 0){
					if(pr1 == 0 || pr2 == 0){
						reduction = reduction;
						//	reduction = -FLT_MAX;
					}else{
						//fprintf`
						reduction = reduction - right*(-1.0*(pr1 * log2(pr1) + (pr2) * log2(pr2)));
						//	reduction = -FLT_MAX;
					}
					
				}else{
					if(pr1 == 0 || pr2 == 0){
						reduction = reduction - left* ( -1.0 * (pl1 * log2(pl1) + (pl2) * log2(pl2)  ));
						
						
						//	reduction = -FLT_MAX;
					}else{
						reduction = reduction - left* ( -1.0 * (pl1 * log2(pl1) + (pl2) * log2(pl2)  )) -   right*(-1.0*(pr1 * log2(pr1) + (pr2) * log2(pr2)));
					}
				}
				
				
				//fprintf(stderr,"%d		%f\n",j, reduction);
				
				
				if(reduction > max_reduction){
					max_reduction = reduction;
					a = i;
					b = j;
				//	fprintf(stderr,"BEST:	%s	%s	%f\n",all_patterns[a]->label,all_patterns[b]->label,max_reduction);
				}
				//if(strlen(patterns[i]->label) > 5 && strlen(patterns[j]->label)  ){
				//	fprintf(stderr,"%d	%d	%s	%s	%f\n",i,j,patterns[i]->label,  patterns[j]->label,  dm[i][j] );
				//}
			}
		//}
	}
	if(max_reduction){
		fprintf(stderr,"BEST:	%s	%f\n",all_patterns[a]->label,max_reduction);
	}
	//i = a;
	//j = b;
	
	if(max_reduction > 0.1){
		
		int* left_samples = 0;
		int* right_samples = 0;
		
		int suffix_len_left = 0;
		
		int suffix_len_right = 0;
		
		int* seq_id_in_suffix_right =  malloc(sizeof(int) * pst->suffix_len );
		char** suffix_right = malloc(sizeof(char*) * pst->suffix_len );
		
		
		left_samples = malloc(sizeof(int)* numseq);
		right_samples = malloc(sizeof(int)* numseq);
		b = 0;
		c = 0;
		for(i = 0; i < numseq;i++){
			if(samples[i]){
				if(bit_test(all_patterns[a]->bit_occ, i)){
					left_samples[i] = 1;
					right_samples[i] = 0;
					b++;
					
				}else{
					left_samples[i] = 0;
					right_samples[i] =1;
					c++;
				}
			}else{
				left_samples[i] = 0;
				right_samples[i] =0;
			}
		}
		
		
		//a = 0;
		//b = 0;
		suffix_len_left = 0;
		suffix_len_right = 0;
		for(i = 0;i < pst->suffix_len;i++){
			if(samples[ pst->seq_id_in_suffix[i]]){
			if(bit_test(all_patterns[a]->bit_occ,  pst->seq_id_in_suffix[i])){
				pst->seq_id_in_suffix[suffix_len_left] = pst->seq_id_in_suffix[i];
				pst->suffix_array[suffix_len_left] = pst->suffix_array[i];
				suffix_len_left++;
			}else{
				seq_id_in_suffix_right[suffix_len_right] =  pst->seq_id_in_suffix[i];
				suffix_right[suffix_len_right] = pst->suffix_array[i];
				suffix_len_right++;
			}
			}
		}
		free(samples);
		free_pst(pst->pst_root);
		//= malloc(sizeof(struct suffix_node* ) * pst->total_len);
		//free(pst->seq_id_in_suffix);
		//free(pst->suffix_array);
		//free(pst);
		free(all_patterns);
		
		fprintf(stderr,"L:%d	%d	%d	%d	%d\n",b,c,suffix_len_left,suffix_len_right , suffix_len_left+suffix_len_right );
		fprintf(stderr,"Going left\n" );
		pst->suffix_len = suffix_len_left;
		
		pst_based_partition(pst,ri,left_samples,  numseq,b);
		fprintf(stderr,"Going right\n" );
		pst->suffix_len = suffix_len_right;
		free(pst->suffix_array);
		free(pst->seq_id_in_suffix);
		pst->suffix_array = suffix_right;
		pst->seq_id_in_suffix = seq_id_in_suffix_right;
		pst_based_partition(pst,ri, right_samples,  numseq,c);
	}else{
		///cluster_reads_based_on_pst_patterns(all_patterns,num_patterns,numseq,ri);
		
		fprintf(stderr,"Reached an END NODE...\n");
		j = 0;
		for(i = 0; i < numseq;i++){
			if(samples[i]){
				fprintf(stderr,"%s\n",ri[i]->seq);
				j++;
			}
			if(j > 50){
				break;
			}
		}

		
		
		//print_pst(pst, pst->pst_root, ri);
		free(samples);
		free_pst(pst->pst_root);
		//= malloc(sizeof(struct suffix_node* ) * pst->total_len);
		//free(pst->seq_id_in_suffix);
		//free(pst->suffix_array);
		//free(pst);
		free(all_patterns);
	}
	
	
	
	//exit(0);

}

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
	int num_patterns = 0;
	struct pst* pst = 0;
	
	pst = malloc(sizeof(struct pst));
	pst->current_suffix_size = param->num_query* 64;
	pst->suffix_array = malloc(sizeof(char*)* pst->current_suffix_size);
	
	pst->L = MAX_PST_LEN;
	pst->alpha = 0.0f;
	pst->p_min = 0.0001f;
	pst->lamba = 0.001f;
	pst->r = 1.05f;
	pst->total_len = 0;
	pst->pst_root = alloc_node(pst->pst_root,"",0);
	pst->ppt_root = alloc_node(pst->ppt_root,"",0);
	pst->rank_array = 0;
		
	ri = malloc(sizeof(struct read_info*) * param->num_query);
	
	for(i = 0; i < param->num_query;i++){
		ri[i] = malloc(sizeof(struct read_info));
		ri[i]->seq = 0;
		ri[i]->name = 0;
		ri[i]->qual = 0;
		ri[i]->len = 0;
		ri[i]->cigar = 0;
		ri[i]->md = 0;
		//ri[i]->xp = 0;
		ri[i]->strand = malloc(sizeof(unsigned int)* (LIST_STORE_SIZE+1));
		ri[i]->hits = malloc(sizeof(unsigned int)* (LIST_STORE_SIZE+1));
		//ri[i]->read_start = -1;
		//ri[i]->read_end = -1;
	}
	file =  io_handler(file, file_num,param);
	
	struct pst_node** all_patterns = 0;
	

	while ((numseq = fp(ri, param,file)) != 0){
		pst->total_len = 0;
		
		/* turn to normal letterss...
		 
		 */
		for(i = 0; i < numseq;i++){
			for(j = 0; j < ri[i]->len;j++){
				ri[i]->seq[j] = alphabet[(int)ri[i]->seq[j]];
			}
			//ri[i]->seq[10] = 0;
			//fprintf(stderr,"%d ",ri[i]->len);
			pst->total_len += ri[i]->len;
		}
	
		
		cStartClock = clock();
		
				
		if(pst->current_suffix_size < pst->total_len){
			pst->suffix_array = realloc(pst->suffix_array , sizeof(char*)* (pst->total_len+64));
			pst->current_suffix_size =  (pst->total_len+64);
		}
		c = 0;
		pst->mean_length = 0.0;
		for(i = 0; i < numseq;i++){
			for(j = 0; j < ri[i]->len;j++){//ri[i]->len;j++){
				pst->suffix_array[c] = ri[i]->seq +j;
				c++;
			}
			pst->mean_length +=  ri[i]->len;
			
		}
		
		pst->mean_length /= (float)numseq;
		
		
		pst->suffix_len = c;
		pst->numseq = numseq;
		pst->rank_array = malloc(sizeof(struct ranks*) * (int)(numseq));
		for(i = 0; i  < numseq;i++){
			pst->rank_array[i] = malloc(sizeof(struct ranks));
		}
		
		
		fprintf(stderr,"%d\t%d\t%d	%f\n",c,pst->total_len, pst->suffix_len*(4 + 4*5),pst->mean_length);
		///exit(0);
		qsort(pst->suffix_array, pst->suffix_len, sizeof(char *), qsort_string_cmp);
		fprintf(stderr,"built SA in %4.2f seconds\n",(double)( clock() - cStartClock) / (double)CLOCKS_PER_SEC);
		
		cStartClock = clock();

		
		//init root - removes if statement in recursion...
		sum = 0.0;
		for(i = 0;i < 5;i++){
			tmp[0] = alphabet[i];
			tmp[1] = 0;//alphabet[i];
			c = count_string(tmp,(const char**)pst->suffix_array,pst->suffix_len-1,1);
			pst->pst_root->nuc_probability[i] = c;
			pst->ppt_root->nuc_probability[i] = c;
			sum+= c;
		}
		for(i = 0;i < 5;i++){
			
			pst->pst_root->nuc_probability[i] =  pst->pst_root->nuc_probability[i]/ sum;
			pst->ppt_root->nuc_probability[i] =  pst->ppt_root->nuc_probability[i]/ sum;
			//fprintf(stderr,"%c\t%f\n",alphabet[i], n->nuc_probability[i]);
		}
		
		
		pst->pst_root = build_pst(pst,pst->pst_root );
		pst->ppt_root = build_ppt(pst,pst->ppt_root );
		
		fprintf(stderr,"built PST in %4.2f seconds\n",(double)( clock() - cStartClock) / (double)CLOCKS_PER_SEC);
		
		
		cStartClock = clock();
		pst->pst_root  = alloc_bit_occ_pst(pst->pst_root , numseq);
		pst->ppt_root = alloc_bit_occ_pst(pst->ppt_root, numseq);
		ri =  scan_read_with_pst( ri, pst);
		fprintf(stderr,"scanned  in %4.2f seconds\n",(double)( clock() - cStartClock) / (double)CLOCKS_PER_SEC);
		
		print_pst(pst,pst->pst_root,ri);
		
		
		num_patterns = 0;
		num_patterns = count_patterns(pst->pst_root, num_patterns);
		num_patterns = count_patterns(pst->ppt_root,num_patterns);
		//
		//exit(0);
		all_patterns = malloc(sizeof(struct pst_node*)  *num_patterns );
		
		num_patterns = 0;
		num_patterns = add_patterns(all_patterns,pst->pst_root, num_patterns);
		num_patterns  = add_patterns(all_patterns,pst->ppt_root,num_patterns);
		
		
		qsort((void *)  all_patterns, num_patterns, sizeof(struct pst_node* ),(compfn) sort_pst_nodel_according_to_label);
		
		fprintf(stderr,"%d numpatterns\n",num_patterns );
		c = 0;
		for(i = 0; i < num_patterns-1;i++){
			if(strcmp(all_patterns[i]->label, all_patterns[i+1]->label)){
				all_patterns[c] = all_patterns[i];
				c++;
			}
			//fprintf(stderr,"%s	%d\n",all_patterns[i]->label,i );
		}
		
		//for(i = 0; i < c;i++){
		//	fprintf(stderr,"%s	%d\n",all_patterns[i]->label,i );
		//}
		// got all strings combined and ready for dbscan.....
		
		cStartClock = clock();
		cluster_reads_based_on_pst_patterns(all_patterns,c,numseq,ri);
		fprintf(stderr,"clustered   in %4.2f seconds\n",(double)( clock() - cStartClock) / (double)CLOCKS_PER_SEC);
		
		//print_pst(pst,pst->pst_root,ri);
		
		exit(0);
		run_gmm_on_sequences(ri,numseq);
		
		
		
		print_pst(pst,pst->pst_root,ri);
		fprintf(stderr,"happily got here...\n");
		print_pst(pst,pst->ppt_root,ri);
		exit(0);
		//cStartClock = clock();
		//ri =  scan_read_with_pst( ri, pst);
		//fprintf(stderr,"scanned  in %4.2f seconds\n",(double)( clock() - cStartClock) / (double)CLOCKS_PER_SEC);
	
		for(i = 0; i  <pst->numseq;i++){
			free(pst->rank_array[i]);// = malloc(sizeof(struct ranks));
		}
		free(pst->rank_array);// = malloc(sizeof(struct ranks*) * (int)(numseq));
		exit(0);
		
		
		cStartClock = clock();
		//pst->pst_root = count_patterns(ri,pst,pst->pst_root);
		fprintf(stderr,"counted  in %4.2f seconds\n",(double)( clock() - cStartClock) / (double)CLOCKS_PER_SEC);
		print_pst(pst,pst->pst_root,ri);
		
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
	
	
	
	//fprintf(stderr,"NODE: %s\n", n->label);
	//for(i = 0;i < 5;i++){
	//	fprintf(stderr,"%c+%s\t%f\n",alphabet[i], n->label, n->nuc_probability[i]);
	//}
	
	//step 2 test expansion
	
	//loop though letters at present node
	if(len + 1 < MAX_PST_LEN ){
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
					if(tmp_counts_s[j] /  (pst->numseq * (pst->mean_length - (len+1))) > 0.001){
						add = 1;
						break;
					}
				}
				if(add){
					// here I know that probablility of 'X' is non-neglible AND that there existsa string 'X' - [A,C,G,T,N] which is frequent in the data - hence I add...
					n->next[i] = alloc_node(n->next[i] ,tmp,len+1);
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
	}
	c= 0;
	for(i = 0; i < 5;i++){
		if(n->next[i]){
			c+= n->next[i]->in_T;
		}
	}
	if(c){
		n->in_T = 1;
	}
	return n;
}



struct pst_node* build_ppt(struct pst* pst,struct pst_node* n )
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
	
	
	
	//fprintf(stderr,"NODE: %s\n", n->label);
	//for(i = 0;i < 5;i++){
	//	fprintf(stderr,"%c+%s\t%f\n",alphabet[i], n->label, n->nuc_probability[i]);
	//}
	
	//step 2 test expansion
	
	//loop though letters at present node
	if(len + 1 < MAX_PST_LEN ){
		/// search for all strings and record probabilities S+ACGT...
		/// don't search rare strings...
		/// - super mega simple ...
		
		
		for(i = 0; i < 5;i++){
			if(n->nuc_probability[i] >= pst->p_min){ /// string (suffix  + letter is not rare compared to all other letters  ) // e.g. acgC = 0.2 and not acgC = 0.0000001. hence it makes sense to "extend the suffix - test [A,C,G,T] - acgC
				
				
				
				
				
				//init longer prefix!!!! 
				
				for(j = 0; j < len;j++){
					tmp[j+1] = n->label[j];
				}
				tmp[len+1] = alphabet[i];
				
				sum = 0.0;
				for(j = 0; j < 5;j++){
					tmp[0]  = alphabet[j];
					tmp[len+2] = 0;
					c = count_string(tmp,(const char**)pst->suffix_array,pst->suffix_len-1,len+2);
					tmp_counts_s[j] = c;
					sum+= c;
				}
				
				add = 0;
				for(j = 0; j < 5;j++){
					if(tmp_counts_s[j] /  (pst->numseq * (pst->mean_length - (len+1))) > 0.001){
						add = 1;
						break;
					}
				}
				if(add){
					// here I know that probablility of 'X' is non-neglible AND that there existsa string 'X' - [A,C,G,T,N] which is frequent in the data - hence I add...
					n->next[i] = alloc_node(n->next[i],tmp+1,len+1);
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
					n->next[i] = build_ppt(pst,n->next[i]);
					
				}
			}
		}
	}
	c= 0;
	for(i = 0; i < 5;i++){
		if(n->next[i]){
			c+= n->next[i]->in_T;
		}
	}
	if(c){
		n->in_T = 1;
	}
	
	return n;
}


struct pst_node* alloc_node(struct pst_node* n,char* string,int len)
{
	int i;
	n = malloc(sizeof(struct pst_node));
	
	assert(n!=0);
	
	n->label = malloc(sizeof(char) *(len+1));
	assert(n->label != 0);
	
	for(i = 0; i < len;i++){
		n->label[i] =string[i];
	}
	n->label[len] = 0;
	n->in_T = 0;
	n->occ = 0;
	n->last_seen = -1;
	//n->bit_occ = 0;
	for(i =0; i < 5;i++){
		n->next[i] = 0;
		n->nuc_probability[i] = 0.2f;
	}

	return n;
}


struct read_info**  scan_read_with_pst(struct read_info** ri,struct pst* pst)
{
	int i,j;
	float P_T;
	//float P_PT;
	float P_R = 0;
	
	float* base_p = pst->pst_root->nuc_probability;
	char* seq;
	char* qual;
	
	float total_T = prob2scaledprob(1.0);
	float total_R = prob2scaledprob(1.0);
	
	float A,B;
	
	//scan to cpount occurances... 
	for(i = 0; i < pst->numseq;i++){
		
		seq = ri[i]->seq;
		for(j = 0; j < ri[i]->len; j++ ){
			pst->pst_root = count_pst_lables(pst->pst_root, seq,  j, i);
			pst->ppt_root = count_ppt_lables(pst->ppt_root, seq, j , i);
			
		}
	}
	
	//fprintf(stderr,"%f\n",pst->numseq);
	
	for(i = 0; i < pst->numseq;i++){
		if(!ri[i]->qual){
			ri[i]->qual = malloc(sizeof(char)* (ri[i]->len+1));
		}
		for(j = 0; j < ri[i]->len; j++ ){
			ri[i]->qual[j] = 48;
		}
		qual = ri[i]->qual;
		seq = ri[i]->seq;
		
		//if(!i){
		
		//}
		
		P_T = prob2scaledprob(1.0);
		P_R = prob2scaledprob(1.0);
		//P_PT =  prob2scaledprob(1.0);
		//fprintf(stdout,"%s\n",seq );
				
		for(j = 0; j < ri[i]->len; j++ ){
			P_R = P_R + prob2scaledprob(base_p[nuc_code5[(int)seq[j]]]);
			
			
			//get_occ
			A = get_pst_prob(pst->pst_root, seq,  nuc_code5[(int)seq[j]], j, i);
			B = get_ppt_prob(pst->ppt_root, seq,  nuc_code5[(int)seq[j]], j, i);
			
			
			P_T = P_T + prob2scaledprob(max(A,B));
			
			//fprintf(stdout,"%d	%f	%f\n", j , A,B);
			ri[i]->qual[j] = 48 + (int) (10.0* (max(A,B)));
		}
		
		total_T = total_T + P_T;
		total_R = total_R + P_R;
		
		ri[i]->mapq = P_T-P_R;
	}
	return ri;
}



int add_patterns(struct pst_node** all_patterns, struct pst_node* n,int num)
{
	int i;
	
	int internal = 0;
	for(i = 0;i < 5;i++){
		if(n->next[i]){
			internal++;
		}
	}
	if (!internal){
	if(strlen(n->label) > 2){
		if(n->in_T){
			all_patterns[num] = n;
			num += 1;
			//fprintf(stderr,"ADDED: %s	%d\n",n->label,num);
		}
	}
		
	}
	
	for(i = 0;i < 5;i++){
		if(n->next[i]){
			//fprintf(stderr,"Going:%d\n",i);
			num = add_patterns(all_patterns,n->next[i],num);
		}
	}
	
	return num;
}

int count_patterns(struct pst_node* n,int num)
{
	int i;
	if(strlen(n->label) > 2){
		if(n->in_T){
			num += 1;
			//fprintf(stderr,"ADDED: %s	%d\n",n->label,num);
		}
	}
	
	for(i = 0;i < 5;i++){
		if(n->next[i]){
			//fprintf(stderr,"Going:%d\n",i);
			num = count_patterns(n->next[i],num);
		}
	}

	return num;
}

float get_pst_prob(struct pst_node* n, char* string,int target, int pos,int seq_id)
{
	int c;
	//if(!seq_id){
	//fprintf(stderr,"%d	%s	%f	(%d)	- OCC:%d\n", target,n->label,n->nuc_probability[target],n->in_T, n->occ);
	//}
	/*if(n->last_seen != seq_id){
		n->occ++;
		n->last_seen = seq_id;
	}*/
	if(pos == 0){
		
		return n->nuc_probability[target];
	}
	pos = pos -1;
	c = nuc_code5[(int)string[pos]];
	if(n->next[c]){
		return get_pst_prob(n->next[c], string, target,pos,seq_id);
	}else{
		
		return n->nuc_probability[target];
	}
}




float get_ppt_prob(struct pst_node* n, char* string,int target, int pos,int seq_id)
{
	int c;
	//if(!seq_id){
	//fprintf(stderr,"%d	%s	%f	(%d)	- OCC:%d\n", target,n->label,n->nuc_probability[target],n->in_T, n->occ);
	//}
	/*if(n->last_seen != seq_id){
	 n->occ++;
	 n->last_seen = seq_id;
	 }*/
	if(string[pos+1] == 0){
		
		return n->nuc_probability[target];
	}
	pos = pos +1;
	c = nuc_code5[(int)string[pos]];
	if(n->next[c]){
		return get_ppt_prob(n->next[c], string, target,pos,seq_id);
	}else{
		
		return n->nuc_probability[target];
	}
}

int get_occ(struct pst_node* n, char* string,int target, int pos,int seq_id)
{
	//if(!seq_id){
	//fprintf(stderr,"%d	%s	%f	(%d)	- OCC:%d\n", target,n->label,n->nuc_probability[target],n->in_T, n->occ);
	//}
	/*if(n->last_seen != seq_id){
	 n->occ++;
	 n->last_seen = seq_id;
	 }*/
	if(pos == 0){
		return n->occ;
	}
	pos = pos -1;
	int c = nuc_code5[(int)string[pos]];
	if(n->next[c]){
		return get_occ(n->next[c], string, target,pos,seq_id);
	}else{
		return n->occ;
	}
}


struct pst_node*  count_pst_lables(struct pst_node* n, char* string, int pos,int seq_id)
{
	//if(!seq_id){
	//fprintf(stderr,"%d	%s	%f	(%d)	- OCC:%d\n", target,n->label,n->nuc_probability[target],n->in_T, n->occ);
	//}
	
	int c = nuc_code5[(int)string[pos]];
	if(n->next[c]){
		if(n->next[c]->last_seen != seq_id){
			if(!bit_test(n->next[c]->bit_occ , seq_id)){
				n->next[c]->occ++;
				n->next[c]->bit_occ = bit_set(n->next[c]->bit_occ , seq_id);
			}
			n->next[c]->last_seen = seq_id;
		}
		pos = pos -1;
		if(pos != -1){
			n->next[c] =  count_pst_lables(n->next[c], string,pos,seq_id);
		}
	}else{
		return n;
	}
	
	/*if(n->last_seen != seq_id){
		n->occ++;
		//n->bit_occ = bit_set(n->bit_occ, seq_id);
		n->last_seen = seq_id;
	}
	//fprintf(stderr,"%s	(L-PST)\n", n->label );
	//fprintf(stderr,"%s\n", string+pos);
	//fprintf(stderr,"%s\n", string);
	
	if(pos == 0){
		return n;
	}
	pos = pos -1;
	
	if(n->next[c]){
		n->next[c] =  count_pst_lables(n->next[c], string,pos,seq_id);
	}else{
		return n;
	}*/
	return n;
}


struct pst_node*  count_ppt_lables(struct pst_node* n, char* string, int pos,int seq_id)
{
	int c = nuc_code5[(int)string[pos]];
	if(!string[pos]){
		return n;
	}
	if(n->next[c]){
		if(n->next[c]->last_seen != seq_id){
			if(!bit_test(n->next[c]->bit_occ , seq_id)){
				n->next[c]->occ++;
				n->next[c]->bit_occ = bit_set(n->next[c]->bit_occ , seq_id);
			}
			n->next[c]->last_seen = seq_id;
		}
		pos = pos +1;
		//n = n->next;
		//if(pos){
			n->next[c] =  count_ppt_lables(n->next[c], string,pos,seq_id);
		//}
	}else{
		return n;
	}

	
	return n;
}




void free_pst(struct pst_node* n)
{
	int i;
	for(i = 0;i < 5;i++){
		if(n->next[i]){
			free_pst(n->next[i]);
		}
	}
	if(n->bit_occ){
		free(n->bit_occ);
	}
	free(n->label);
	free(n);

}

void print_pst(struct pst* pst,struct pst_node* n, struct read_info** ri )
{
	int i;
	int internal;
	
	double p ,e;
	
//	int len = (int)strlen(n->label);
	
//	double N1,N2,U1,U2,R1,R2,Z;
//	int c;
//	struct ranks** rank_array  = pst->rank_array;
	//struct ranks** rank_array = malloc(sizeof(struct ranks*) * (int)(pst->numseq));
	//for(i = 0; i  <pst->numseq;i++){
	//	rank_array[i] = malloc(sizeof(struct ranks));
	//}
	//char alphabet[] = "ACGTN";
	//if(strlen(n->label) > 2){
		internal = 0;
		for(i = 0;i < 5;i++){
			if(n->next[i]){
				internal++;
			}
		}
		if(!internal){
			/*
			//fprintf(stderr,"%p\n",n);
			//fprintf(stderr,"%s	%d	%d\n", n->label,n->in_T, count_string(n->label,(const char**)pst->suffix_array,pst->suffix_len-1,len));
			c = 0;
			N1 = 0;
			N2 = 0;
			//fprintf(stderr,"NUMSEQ:::::%f\n",pst->numseq );
			for(i = 0 ;i  < pst->numseq;i++){
				//fprintf(stderr,"%d ",i);
				if(bit_test(n->bit_occ, i)){
					rank_array[c]->sample = 0;
					rank_array[c]->value = ri[i]->mapq;
					N1++;
				}else{
					rank_array[c]->sample = 1;
					rank_array[c]->value = ri[i]->mapq;
					N2++;
				}
				c++;
				
			}
			
			
			qsort((void *)  rank_array, N1+N2, sizeof(struct ranks* ),(compfn) establish_rank);
			
			U1 = 0.0;
			R1 = 0.0;
			
			U2 = 0.0;
			R2 = 0.0;
			for(i = 0;i < N1+N2;i++){
				
				if(!rank_array[i]->sample){
					R1 += (i+1);
				}else{
					R2 += (i+1);
				}
				//if(i < 10){
				//	fprintf(stderr,"%d\t%d\t%f\n",i, rank_array[i]->sample,rank_array[i]->value);
				//}
			}
			
			U1 = R1 - (N1*(N1 + 1.0))/2.0;
			U2 = R2 - (N2*(N2 + 1.0))/2.0;
			
			if(U2 < U1){
				U1 = U2;
			}
			
			
			Z = (U1- (N1*N2 /2.0) )/  sqrt((N1 * N2 *(N1 + N2 +1) )/  12.0 );
			//	fprintf(stderr,"%f	%f	%f	%f	\n",N1,N2,U1,U2);
			//Z = 1;
			N1 = prob2scaledprob(1.0);
			for(i = 0 ;i < strlen(n->label);i++){
				N1 += prob2scaledprob(pst->pst_root->nuc_probability[nuc_code5[ (int)n->label[i]]] );
			}
			N1 = N1 + prob2scaledprob( pst->mean_length - strlen(n->label)) + prob2scaledprob(pst->numseq) ;
			*/
			
		p =  (double) n->occ / (double)  pst->numseq;
		
		e = p * log2(p);
		
		p = 1.0 -  (double) n->occ / (double)  pst->numseq;
		e+=  p * log2(p);
		e *= -1.0;
			fprintf(stderr,"%s	%d	%d	Entropy:%f	%f\n", n->label,n->in_T,n->occ,  e,pst->numseq);
			//for(i = 0;i < 5;i++){
				//if(n->next[i]){
			//	fprintf(stderr,"%f ",n->nuc_probability[i]);
				//}else{
				//	fprintf(stderr,"%f S\t",n->nuc_probability[i]);
				//}
			//}
			
			
		//	fprintf(stderr,"\n");
		}
	//}
	
	
	for(i = 0;i < 5;i++){
		if(n->next[i]){
			if(n->next[i]->in_T){
			//fprintf(stderr,"Going:%d\n",i);
			print_pst(pst,n->next[i],ri);
			}
		}
	}
}


struct pst_node* alloc_bit_occ_pst(struct pst_node* n, int num)
{
	int i;
	n->bit_occ = _mm_malloc(sizeof(int)* (1+ num / BITSPERWORD) , 16);
	
	assert(n->bit_occ != 0);
	
	for(i = 0; i < (1+ num / BITSPERWORD);i++){
		n->bit_occ[i] = 0;
	}
	
	for(i = 0;i < 5;i++){
		if(n->next[i]){
			n->next[i] = alloc_bit_occ_pst(n->next[i],num);
		}
	}
	return n;
}


int* bit_set(int*a, int i)
{
	a[i >> SHIFT] |= (1 <<(i & MASK));
	return a;
}


int* bit_clr(int*a, int i)
{
	a[i >> SHIFT] &= ~(1 <<(i & MASK));
	return a;
}

int bit_test(int*a, int i)
{
	return a[i >> SHIFT] & (1 <<(i & MASK));
}


int establish_rank(const void *a, const void *b)
{
	struct ranks* const *one = a;
	struct ranks* const *two = b;
	
	if((*one)->value >  (*two)->value){
		return -1;
	}else{
		return 1;
	}
}


int sort_pst_nodel_according_to_label(const void *a, const void *b)
{
	struct pst_node * const *one = a;
	struct pst_node* const *two = b;
	
	return strcmp((*one)->label,(*two)->label) ;
}






int qsort_suffix_node_string_cmp(const void *a, const void *b)
{
	struct suffix_node* const *one  = a;
	struct suffix_node* const *two  = b;
	
	return strcmp((*one)->string,(*two)->string) ;
}






