//
//  barcode_hmm.c
//  tagdust2
//  
//  Created by lassmann on 2/5/13.
//  Copyright (c) 2013 lassmann. All rights reserved.
//

#include <stdio.h>
#include "tagdust2.h"
#include "interface.h"
#include "io.h"
#include "misc.h"
#include "nuc_code.h"

#include "barcode_hmm.h"


void hmm_controller(struct parameters* param,int (*fp)(struct read_info** ,struct parameters*,FILE* ),int file_num)
{
	struct read_info** ri = 0;
	int i,j,g,c;
	int numseq;
	int total_read = 0;
	
	float sum = 0.0;
	
	int segment_length;
	int read_length;
	
	init_logsum();
	
	/*int len,c;
	float max;
	
	for(j = 0; j < 10;j++){
	len = j+1;
	//len = model->hmms[0]->num_columns;
	c = 0;
	max = -1;
	for(i = 0; i< 1000;i++){
		if(binomial_distribution((double)i / 1000.0 , len  ,1)  > max){
			max = binomial_distribution((double)i / 1000.0 ,len ,1 );
			c = i;
		}
	}
	max = (double)c / 1000.0;
		fprintf(stderr,"%d	%f\n",j+1,max	);
	}
	exit(0);*/
	
	
	float* back = 0;
	int average_length = 0;
	back = malloc(sizeof(float)*5);
	for(i = 0; i < 5;i++){
		back[i]= 0.0f;//prob2scaledprob( 0.2);
	}
	
	param->num_query = 100;
	
	/*model = init_model(model , back,24);
	
	struct model* model2 = copy_and_malloc_model(model);
	
	free_model(model);
	free_model(model2);
	free(back);
	
	free_param(param);
	exit(0);
	*/
	 //char command[1000];
	//char  tmp[1000];
	FILE* file = 0;
	FILE* unmapped = 0;
	
	ri = malloc(sizeof(struct read_info*) * param->num_query);
	
	assert(ri !=0);
	
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
	
	if(param->print_unmapped){
		if (!(unmapped = fopen(param->print_unmapped, "w" ))){
			fprintf(stderr,"Cannot open file '%s'\n",param->print_unmapped);
			exit(-1);
		}
	}
	
	/*
	 
	 get backgorund nucleotide distribution - from all reads?
	 
	 */
	
	
	average_length = 0;
	for(i = 0; i < 5;i++){
		back[i] = 0.0;
	}
	total_read = 0;
	
	while ((numseq = fp(ri, param,file)) != 0){
		//fprintf(stderr,"rread: %d\n",numseq);
		for(i = 0; i < numseq;i++){
			average_length += ri[i]->len;
			for(j = 0;j < ri[i]->len;j++){
				back[(int)ri[i]->seq[j]] += 1.0f;
			}
		}
		total_read += numseq;
		if(total_read > 10){
			break;
		}
	}
	
	average_length = average_length / total_read;
	
	param->average_read_length = average_length;
	
	for(i = 0; i < 5;i++){
		sum += back[i];
	}
	
	for(i = 0; i < 5;i++){
		back[i] = prob2scaledprob(back[i]  / sum);
	}
	
	rewind(file);
	
	
	//init models
	
	struct model_bag* mb = malloc(sizeof(struct model_bag));
	
	assert(mb!=0);
	
	mb->model = malloc(sizeof(struct model* ) * param->read_structure->num_segments);
	
	assert(mb->model);
	
	
	mb->f_score = prob2scaledprob(0.0f);
	mb->b_score = prob2scaledprob(0.0f);
	mb->num_models = param->read_structure->num_segments;
	// get read length estimate...
	read_length = average_length;
	for(i = 0; i < mb->num_models;i++){
		//mb->model[i] = malloc_model_according_to_read_structure(param->read_structure,i);
		
		if(param->read_structure->type[i] == 'G'){
			read_length = read_length -2;
			
		}else{
			read_length = read_length - (int)strlen(param->read_structure->sequence_matrix[i][0]);
		}
		
		
		
	}
	fprintf(stderr,"READlength: %d",read_length);
	
	mb->total_hmm_num = 0;
	
	
	for(i = 0; i < mb->num_models;i++){
				mb->model[i] = malloc_model_according_to_read_structure(param->read_structure,i);
		segment_length = 0;
		if(param->read_structure->type[i] == 'G'){
			segment_length = 2;
			
		}
		if(param->read_structure->type[i]  == 'R'){
			segment_length = read_length;
		}
		
		
		
		mb->model[i] = init_model_according_to_read_structure(mb->model[i], param, i,back,segment_length);
		print_model(mb->model[i]);
		mb->total_hmm_num += mb->model[i]->num_hmms;

	}
	
	mb->path = malloc(sizeof(int*) * MAX_SEQ_LEN);
	mb->dyn_prog_matrix = malloc(sizeof(float*) * MAX_SEQ_LEN );

	for (i = 0; i < MAX_SEQ_LEN;i++){
		mb->path[i] = malloc(sizeof(int)* (mb->total_hmm_num +1) );
		mb->dyn_prog_matrix[i] = malloc(sizeof(float) * (mb->total_hmm_num +1) );
	}
	
	mb->transition_matrix = malloc(sizeof(float*) * (mb->total_hmm_num +1));
	mb->label = malloc(sizeof(int) *  (mb->total_hmm_num +1));
	c = 0;
	for(i = 0; i < mb->num_models ;i++){
		for(j = 0; j < mb->model[i]->num_hmms;j++){
			mb->label[c] = (j << 16) | i ;
			if(mb->model[i]->skip != prob2scaledprob(0.0)){
				mb->label[c]  |= 0x80000000;
			}
			fprintf(stderr,"%d %d	%d %d\n",c,mb->label[c],mb->label[c] & 0xFFFF, (mb->label[c] >> 16) & 0x7FFF);
			c++;
			
		}
	}
	
	for(i = 0; i < mb->total_hmm_num+1 ;i++){
		mb->transition_matrix[i] = malloc(sizeof(float) * (mb->total_hmm_num +1));
		for(j = 0; j <  mb->total_hmm_num+1 ;j++){
			mb->transition_matrix[i][j] = 0;
		}
	}
	

	for(i = 0; i < mb->total_hmm_num ;i++){
		//mb->substitution_matrix[i] = malloc(sizeof(float) * (mb->total_hmm_num +1));
		//c = 0; // assume no skipping....
		//if(mb->label[i] & 0x80000000){
		//	c =1;
		//}
		//if(mb->model[i]->skip){
		//	c =1;
		//}
		
		c = 1;
		for(j = i+1; j <  mb->total_hmm_num ;j++){
			mb->transition_matrix[i][j] = 0;
			
			
			
			if(i == j){
				mb->transition_matrix[i][j] = 1;
			}
			
			if((mb->label[i] & 0xFFFF)+1 == ((mb->label[j] & 0xFFFF) ) ){
				mb->transition_matrix[i][j] = 1;
			}
			
			
			if(((mb->label[i] & 0xFFFF) < ((mb->label[j] & 0xFFFF) ) )&& c ){
				mb->transition_matrix[i][j] = 1;
			}
			
			if(!(mb->label[j] & 0x80000000)){
				c =0;
			}
			fprintf(stderr,"%d, %d, %d %d\n ", j,   mb->label[j],mb->label[j] & 0xFFFF, (mb->label[j] >> 16) & 0x7FFF);
			//if(!(mb->label[j] & 0x80000000)){
			//	c =0;
			//}
			//if(!mb->model[mb->label[j] & 0xFFFF]->skip){
			//	c = 0;
			//}
			
			//mb->substitution_matrix[j][i] =  mb->substitution_matrix[i][j] ;
			
		}
		
		// remain in the same state.... 
		mb->transition_matrix[i][i] = 1;
	}
	
	for(i = 0; i < mb->total_hmm_num ;i++){
		for(j = 0; j <  mb->total_hmm_num ;j++){
			fprintf(stderr,"%f ",mb->transition_matrix[i][j] );
		}
		fprintf(stderr,"\n");
	}
	fprintf(stderr,"\n");
	
	
	
	char test[] = "AAAAA";
	
	for(i = 0; i < 5;i++){
		test[i] = nuc_code5[(int)test[i]];
	}
 	
	//mb = forward(mb, test, 5);

	mb = backward(mb, test ,5);
	
	//mb = forward_extract_posteriors(mb, test ,5);
	mb =  forward_max_posterior_decoding(mb, test ,5);
	
	for(i = 0; i < mb->num_models;i++){
		print_model(mb->model[i]);
	}
	
	
	/*char test2[] = "CCCAAAA";
	
	for(i = 0; i < 7;i++){
		test2[i] = nuc_code5[test2[i]];
	}
 	
	//mb = forward(mb, test, 5);
	
	mb = backward(mb, test2 ,7);
	
	mb = forward_extract_posteriors(mb, test2,7);
	
		
	for(i = 0; i < mb->num_models;i++){
		print_model(mb->model[i]);
	}
	*/
	/*
	test[1] = 1;
	test[3] = 1;
	
	mb = backward(mb, test ,5);
	
	mb = forward_extract_posteriors(mb, test ,5);
		
	for(i = 0; i < mb->num_models;i++){
		print_model(mb->model[i]);
	}
	*/
	
	for(i = 0; i < mb->num_models;i++){
		free_model(mb->model[i]);
	}
	free(mb->model);
	free(mb);
	exit(0);
	
	/*
	 
	 got backgorund nucleotide distribution - from all reads?
	 
	 */
	
	/*
	 
	 init model; 
	 
	 */
	
	
	//model = init_model(model);
	/*
	 
	 done init model;
	 
	 */
	
	
	/*
	 
	run Baum Welch  X times using threads
	 
	 */
	
	/*for(i = 0; i < 10;i++){
		
		while ((numseq = fp(ri, param,file)) != 0){
			//fprintf(stderr,"rread: %d\n",numseq);
			for(i = 0; i < numseq;i++){
				model->average_length += ri[i]->len;
				for(j = 0;j < ri[i]->len;j++){
					model->background_nuc_frequency[ri[i]->seq[j]] += 1;
				}
			}
			total_read+= numseq;
			
		}
		
		
		//re-restimate...
		
		
	}*/
	
	/*
	 
	done...  run Baum Welch  X times using threads
	 
	 */
	
	
	if(param->print_unmapped){
		fclose(unmapped);
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


struct model_bag* backward(struct model_bag* mb, char* a, int len)
{
	int i,j;
	int f,g;
	
	int model_len = 0;
	
	struct hmm* hmm = 0;
	struct hmm_column* c_hmm_column = 0;
	struct hmm_column* p_hmm_column = 0;
	
	float previous_silent[MAX_HMM_SEQ_LEN];
	
	float* psilent;
	float* csilent;
	
	
	char* seqa = a -1;
	
	int c;
	
	//init - len+1 set to zero.... 
	
	for(j = 0; j < mb->num_models;j++){
		for(f = 0;f < mb->model[j]->num_hmms;f++){
			model_len = mb->model[j]->hmms[f]->num_columns-1;
			for(g = 0;g < mb->model[j]->hmms[f]->num_columns;g++){
				c_hmm_column = mb->model[j]->hmms[f]->hmm_column[g];
				for(i = 0; i <= len+1;i++){
					c_hmm_column->M_backward[i] = prob2scaledprob(0.0);
					c_hmm_column->I_backward[i] = prob2scaledprob(0.0);
					c_hmm_column->D_backward[i] = prob2scaledprob(0.0);
				}
			}
		}
		for(i = 0; i <= len+1;i++){
			mb->model[j]->silent_backward[i] = prob2scaledprob(0.0f);
		}
	}
	
	for(i = 0; i <= len+1;i++){
		previous_silent[i] = prob2scaledprob(0.0f);
	}
	previous_silent[len+1] = prob2scaledprob(1.0f);
	
	mb->model[mb->num_models-1]->silent_backward[len+1] = prob2scaledprob(1.0) + mb->model[mb->num_models-1]->skip;
	
	for(j = mb->num_models-2 ; j >= 0;j--){
		mb->model[j]->silent_backward[len+1] = mb->model[j+1]->silent_backward[len+1] + mb->model[j]->skip;
	}
	
	
	//start with last segment... 
	for(j = mb->num_models-1 ; j >= 0;j--){
		if(j == mb->num_models-1 ){
			psilent = previous_silent;
		}else{
			psilent = mb->model[j+1]->silent_backward;
		}
		
		
		csilent= mb->model[j]->silent_backward;
		for(f = 0;f < mb->model[j]->num_hmms;f++){
			hmm = mb->model[j]->hmms[f];
			model_len = mb->model[j]->hmms[f]->num_columns-1;
			//previous_silent[len+1] = logsum(previous_silent[len+1],  current_silent[len+1]+ mb->model[j] ->skip );
			//csilent[len+1] =psilent[len+1] + mb->model[j]->skip;
			for(i = len ; i > 0;i-- ){
				
				c = (int)seqa[i+1];
				c_hmm_column = hmm->hmm_column[model_len];
				
				c_hmm_column->M_backward[i] = psilent[i+1] + mb->model[j]->M_to_silent[f] ;
				//fprintf(stderr," Mback at modellen:%f\n", c_hmm_column->M_backward[i] );
				
				c_hmm_column->I_backward[i] =  psilent[i+1]+ mb->model[j]->I_to_silent[f] ;
				
				c_hmm_column->I_backward[i] = logsum(c_hmm_column->I_backward[i] , c_hmm_column->M_backward[i+1] + c_hmm_column->short_transition[IM] + c_hmm_column->m_emit[c]);
				
				c_hmm_column->I_backward[i] = logsum(c_hmm_column->I_backward[i] , c_hmm_column->I_backward[i+1] + c_hmm_column->short_transition[II] + c_hmm_column->i_emit[c]);
				
				
				
				c_hmm_column->D_backward[i] = prob2scaledprob(0.0f);
				for(g = model_len-1;g >= 0;g--){
					c_hmm_column = hmm->hmm_column[g];
					p_hmm_column = hmm->hmm_column[g+1];
					
					c_hmm_column->M_backward[i]  = p_hmm_column->M_backward[i+1] + p_hmm_column->m_emit[c] + c_hmm_column->short_transition[MM];
					
					//insert - emit previous symbol etc. etc.
					c_hmm_column->M_backward[i] = logsum(c_hmm_column->M_backward[i] , c_hmm_column->I_backward[i+1] +c_hmm_column->i_emit[c ]  + c_hmm_column->short_transition[MI]);
					
					//delete - neex to go to previous columns
					
					c_hmm_column->M_backward[i] = logsum(c_hmm_column->M_backward[i],p_hmm_column->D_backward[i] + c_hmm_column->short_transition[MD]);
					
					// insert state..
					// from previous insertion....
					c_hmm_column->I_backward[i] = c_hmm_column->I_backward[i+1] + c_hmm_column->short_transition[II] + c_hmm_column->i_emit[c];
					//from previous match state....
					
					c_hmm_column->I_backward[i] = logsum( c_hmm_column->I_backward[i],p_hmm_column->M_backward[i+1] + c_hmm_column->short_transition[IM] + p_hmm_column->m_emit[c]);
					///GRERRRRRRR 
					
					//delete state
					
					//from previous delection
					c_hmm_column->D_backward[i] = p_hmm_column->D_backward[i] + c_hmm_column->short_transition[DD];
					
					//from previous match (i.e. gap close
					
					c_hmm_column->D_backward[i] = logsum(c_hmm_column->D_backward[i], p_hmm_column->M_backward[i] + p_hmm_column->m_emit[(int) seqa[i]] + c_hmm_column->short_transition[DM]);
					
				}
				c_hmm_column = hmm->hmm_column[0];
				// link j+1 to j... dfor silent;
				csilent[i] = logsum(csilent[i], c_hmm_column->M_backward[i] + mb->model[j]->silent_to_M[f] + c_hmm_column->m_emit[(int)seqa[i]]);
				csilent[i] = logsum(csilent[i], c_hmm_column->I_backward[i] + mb->model[j]->silent_to_I[f] + c_hmm_column->i_emit[(int)seqa[i]]);
				
				//fprintf(stderr,"Looking for Insertyion to silent in segment1: %d	%f\n",f, mb->model[j]->silent_to_I[f]);
				
				//this should come from previous state .....
				csilent[i] = logsum(csilent[i], psilent[i] + mb->model[j]->skip);
			}
		}
	}
	
	mb->b_score = mb->model[0]->silent_backward[1];
	fprintf(stderr,"SCore:%f	%f\n", mb->b_score , scaledprob2prob(mb->b_score) );
	
	fprintf(stderr," BACKWARD:::::::::::\n");
	
	/*for(j = 0; j < mb->num_models;j++){
		for(i = 0; i <= len;i++){
			fprintf(stderr,"%d	%d	%f\n",j,i,mb->model[j]->silent[i]  );
		}
	}
	exit(0);*/
	
	/*for(j = 0; j < mb->num_models;j++){
		
		for(f = 0;f < mb->model[j]->num_hmms;f++){
			for(g = 0;g < mb->model[j]->hmms[f]->num_columns;g++){
				
				c_hmm_column = mb->model[j]->hmms[f]->hmm_column[g];
				for(i = 0; i <= len;i++){
					//c_hmm_column->M_foward[i] = prob2scaledprob(0.0);
					//c_hmm_column->I_foward[i] = prob2scaledprob(0.0);
					//c_hmm_column->D_foward[i] = prob2scaledprob(0.0);
					fprintf(stderr,"segment:%d	HMM:%d	column:%d	i:%d	%f	%f	%f\n",j,f,g,i, scaledprob2prob( c_hmm_column->M_backward[i]) , scaledprob2prob( c_hmm_column->I_backward[i]),scaledprob2prob (c_hmm_column->D_backward[i])   );
				}
			}
		}
	}
	*/
	return mb;
}

struct model_bag* forward(struct model_bag* mb, char* a, int len)
{
	
	int i,j,c;
	int f,g;
	
	struct hmm* hmm = 0;
	struct hmm_column* c_hmm_column = 0;
	struct hmm_column* p_hmm_column = 0;
	
	char* seqa = a -1;
	
	float* psilent;
	float* csilent;
	
	float previous_silent[MAX_HMM_SEQ_LEN];
	//float current_silent[MAX_HMM_SEQ_LEN];
	
	//init
	
	//float silent_start = prob2scaledprob(1.0);
	
	// M state of first set of HMMS.....
	for(j = 0; j < mb->num_models;j++){
		for(f = 0;f < mb->model[j]->num_hmms;f++){
			for(g = 0;g < mb->model[j]->hmms[f]->num_columns;g++){
				c_hmm_column = mb->model[j]->hmms[f]->hmm_column[g];
				//i = 0;
				for(i = 0; i <= len;i++){
					c_hmm_column->M_foward[i] = prob2scaledprob(0.0);
					c_hmm_column->I_foward[i] = prob2scaledprob(0.0);
					c_hmm_column->D_foward[i] = prob2scaledprob(0.0);
					//fprintf(stderr,"segment:%d	HMM:%d	column:%d	i:%d	%f	%f	%f\n",j,f,g,i, c_hmm_column->M_foward[i] ,c_hmm_column->I_foward[i],c_hmm_column->D_foward[i]    );
				}
			}
		}
		for(i = 0; i <= len+1;i++){
			mb->model[j]->silent_forward[i] = prob2scaledprob(0.0f);
		}
	}
	
	//fprintf(stderr,"\n\n\n");
	mb->model[0]->silent_forward[0] = prob2scaledprob(1.0) + mb->model[0]->skip;
	//fprintf(stderr,"Init silent states... \n");
	//fprintf(stderr,"%d	%f\n",0,mb->model[0]->silent_forward[0]   );
	for(j = 1; j < mb->num_models;j++){
		mb->model[j]->silent_forward[0] = mb->model[j-1]->silent_forward[0]  + mb->model[j]->skip ;
		//fprintf(stderr,"%d	%f	%f	%f\n",j,mb->model[j]->silent_forward[0]  ,mb->model[j-1]->silent_forward[0]  , mb->model[j]->skip  );
	}
	
	
 	for(i = 0; i <= len;i++){
		
		
		
		previous_silent[i] = prob2scaledprob(0.0f);
	}
	previous_silent[0] = prob2scaledprob(1.0);

	//loop thorugh the segments
	// in each run the contained HMMS and update silent states;
	for(j = 0; j < mb->num_models;j++){
		if(j == 0){
			psilent = previous_silent;
		}else{
			psilent =  mb->model[j-1]->silent_forward;
		}
		csilent = mb->model[j]->silent_forward;
		for(f = 0;f < mb->model[j]->num_hmms;f++){
			hmm = mb->model[j]->hmms[f];
			for(i = 1; i <= len;i++){
				c = seqa[i];
				
				c_hmm_column = hmm->hmm_column[0];
				// first column  comes from previous state cheekily transferring its pd to M[0[
				c_hmm_column->M_foward[i] = psilent[i-1] + mb->model[j]->silent_to_M[f] + c_hmm_column->m_emit[c];
				
 				c_hmm_column->I_foward[i]    =psilent[i-1] + mb->model[j]->silent_to_I[f] ;
				
				//add transitions to first columns////
				c_hmm_column->I_foward[i] = logsum(c_hmm_column->I_foward[i], c_hmm_column->I_foward[i-1] + c_hmm_column->short_transition[II]);
				
				c_hmm_column->I_foward[i]  = logsum ( c_hmm_column->I_foward[i] ,c_hmm_column->M_foward[i-1] + c_hmm_column->short_transition[MI]);
				
				c_hmm_column->I_foward[i] = c_hmm_column->I_foward[i] + c_hmm_column->i_emit[c];
				
				c_hmm_column->D_foward[i] = prob2scaledprob(0.0f);
				
				for(g = 1;g < hmm->num_columns;g++){
					c_hmm_column = hmm->hmm_column[g];
					p_hmm_column = hmm->hmm_column[g-1];
					
					//Match state
					//transition from previous match state
					c_hmm_column->M_foward[i] = p_hmm_column->M_foward[i-1] + p_hmm_column->short_transition[MM];
					//transition from previous insert state
					c_hmm_column->M_foward[i] = logsum(c_hmm_column->M_foward[i] , p_hmm_column->I_foward[i-1] + p_hmm_column->short_transition[IM] );
					//transition from previous delete state
					c_hmm_column->M_foward[i] = logsum(c_hmm_column->M_foward[i], p_hmm_column->D_foward[i] + p_hmm_column->short_transition[DM]);
					
					// emission promability in curent M state ;
					c_hmm_column->M_foward[i]  = c_hmm_column->M_foward[i]  + c_hmm_column->m_emit[c];
					
					
					// Instertion State ..
					//self loop insertion to insertion
					c_hmm_column->I_foward[i] = c_hmm_column->I_foward[i-1] + c_hmm_column->short_transition[II];
					// start new insertion
					c_hmm_column->I_foward[i]  = logsum ( c_hmm_column->I_foward[i] ,c_hmm_column->M_foward[i-1] + c_hmm_column->short_transition[MI]);
					
					//instertion emission...
					c_hmm_column->I_foward[i]  = c_hmm_column->I_foward[i]  + c_hmm_column->i_emit[c];
					
					
					// deletion state
					//from previous match state.
					c_hmm_column->D_foward[i] = p_hmm_column->M_foward[i] + p_hmm_column->short_transition[MD];
					//from previous delete state
					
					c_hmm_column->D_foward[i] = logsum(c_hmm_column->D_foward[i], p_hmm_column->D_foward[i] + p_hmm_column->short_transition[DD] );
					
				}
				//fprintf(stderr,"%d	%f	%f	%f	%f	%f\n",i,  c_hmm_column->M_foward[i] , mb->model[j]->M_to_silent[f],  c_hmm_column->I_foward[i] , mb->model[j]->I_to_silent[f], psilent[i]);
				csilent[i] =  logsum(csilent[i],  c_hmm_column->M_foward[i] + mb->model[j]->M_to_silent[f]);
				csilent[i] =  logsum(csilent[i],  c_hmm_column->I_foward[i] + mb->model[j]->I_to_silent[f]);
				csilent[i] = logsum(csilent[i], psilent[i] + mb->model[j]->skip);
				
			}
			
		}
	}
	
		
	mb->f_score = mb->model[mb->num_models-1]->silent_forward[len];
	
	fprintf(stderr,"SCORE:%f	%f\n", mb->f_score, scaledprob2prob(mb->f_score));
	
	
	for(j = 0; j < mb->num_models;j++){
		
		for(f = 0;f < mb->model[j]->num_hmms;f++){
			for(g = 0;g < mb->model[j]->hmms[f]->num_columns;g++){
				
				c_hmm_column = mb->model[j]->hmms[f]->hmm_column[g];
				for(i = 0; i <= len;i++){
					//c_hmm_column->M_foward[i] = prob2scaledprob(0.0);
					//c_hmm_column->I_foward[i] = prob2scaledprob(0.0);
					//c_hmm_column->D_foward[i] = prob2scaledprob(0.0);
					fprintf(stderr,"segment:%d	HMM:%d	column:%d	i:%d	%f	%f	%f\n",j,f,g,i,  scaledprob2prob ( c_hmm_column->M_foward[i]) ,scaledprob2prob ( c_hmm_column->I_foward[i]),scaledprob2prob (c_hmm_column->D_foward[i] )   );
				}
			}
		}
	}
	return mb;
}

struct model_bag* forward_extract_posteriors(struct model_bag* mb, char* a, int len)
{
	
	int i,j,c;
	int f,g;
	
	struct hmm* hmm = 0;
	struct hmm_column* c_hmm_column = 0;
	struct hmm_column* p_hmm_column = 0;
	
	char* seqa = a -1;
	
	float* psilent;
	float* csilent;
	float* bsilent;
	
	float previous_silent[MAX_HMM_SEQ_LEN];
	float next_silent[MAX_HMM_SEQ_LEN];
		
	// M state of first set of HMMS.....
	for(j = 0; j < mb->num_models;j++){
		for(f = 0;f < mb->model[j]->num_hmms;f++){
			for(g = 0;g < mb->model[j]->hmms[f]->num_columns;g++){
				c_hmm_column = mb->model[j]->hmms[f]->hmm_column[g];
				//i = 0;
				for(i = 0; i <= len;i++){
					c_hmm_column->M_foward[i] = prob2scaledprob(0.0);
					c_hmm_column->I_foward[i] = prob2scaledprob(0.0);
					c_hmm_column->D_foward[i] = prob2scaledprob(0.0);
					//fprintf(stderr,"segment:%d	HMM:%d	column:%d	i:%d	%f	%f	%f\n",j,f,g,i, c_hmm_column->M_foward[i] ,c_hmm_column->I_foward[i],c_hmm_column->D_foward[i]    );
				}
			}
		}
		for(i = 0; i <= len+1;i++){
			mb->model[j]->silent_forward[i] = prob2scaledprob(0.0f);
		}
	}
	
	//fprintf(stderr,"\n\n\n");
	mb->model[0]->silent_forward[0] = prob2scaledprob(1.0) + mb->model[0]->skip;
	//fprintf(stderr,"Init silent states... \n");
	//fprintf(stderr,"%d	%f\n",0,mb->model[0]->silent_forward[0]   );
	for(j = 1; j < mb->num_models;j++){
		mb->model[j]->silent_forward[0] = mb->model[j-1]->silent_forward[0]  + mb->model[j]->skip ;
		//mb->model[j]->silent_forward[0] = mb->model[j-1]->silent_forward[0]  + mb->model[j-1]->skip ;
	}
	
	
 	for(i = 0; i <= len;i++){
		
		
		next_silent[i] = prob2scaledprob(0.0f);
		previous_silent[i] = prob2scaledprob(0.0f);
	}
	previous_silent[0] = prob2scaledprob(1.0);
	next_silent[len+1] = prob2scaledprob(1.0f);
	
	//loop thorugh the segments
	// in each run the contained HMMS and update silent states;
	for(j = 0; j < mb->num_models;j++){
		if(j == 0){
			psilent = previous_silent;
		}else{
			psilent =  mb->model[j-1]->silent_forward;
		}
		csilent = mb->model[j]->silent_forward;
		if(j +1 != mb->num_models){
			bsilent = mb->model[j+1]->silent_backward;
		}else{
			bsilent = next_silent;
		}
		for(f = 0;f < mb->model[j]->num_hmms;f++){
			hmm = mb->model[j]->hmms[f];
			for(i = 1; i <= len;i++){
				c = seqa[i];
				
				c_hmm_column = hmm->hmm_column[0];
				// first column  comes from previous state cheekily transferring its pd to M[0[
				c_hmm_column->M_foward[i] = psilent[i-1] + mb->model[j]->silent_to_M[f] + c_hmm_column->m_emit[c];
				
				//***************post
				mb->model[j]->silent_to_M_e[f] = logsum(mb->model[j]->silent_to_M_e[f] ,psilent[i-1] + mb->model[j]->silent_to_M[f] + c_hmm_column->m_emit[c] + c_hmm_column->M_backward[i]  -mb->b_score);
				
				
				
				
				c_hmm_column->m_emit_e[c] = logsum(c_hmm_column->m_emit_e[c] , c_hmm_column->M_foward[i]  + c_hmm_column->M_backward[i] -mb->b_score);
				//***************post
				
				
 				c_hmm_column->I_foward[i]    =psilent[i-1] + mb->model[j]->silent_to_I[f] ;
				
				
				
				
				//add transitions to first columns////
				c_hmm_column->I_foward[i] = logsum(c_hmm_column->I_foward[i], c_hmm_column->I_foward[i-1] + c_hmm_column->short_transition[II]);
				
				c_hmm_column->I_foward[i]  = logsum ( c_hmm_column->I_foward[i] ,c_hmm_column->M_foward[i-1] + c_hmm_column->short_transition[MI]);
				
				c_hmm_column->I_foward[i] = c_hmm_column->I_foward[i] + c_hmm_column->i_emit[c];
				
				//***************post
				mb->model[j]->silent_to_I_e[f]  = logsum(mb->model[j]->silent_to_I_e[f] , psilent[i-1] + mb->model[j]->silent_to_I[f]  + c_hmm_column->i_emit[c] +  c_hmm_column->I_backward[i] -mb->b_score);
				
				c_hmm_column->short_transition_e[II] = logsum(c_hmm_column->short_transition_e[II],c_hmm_column->I_foward[i-1] + c_hmm_column->short_transition[II] + c_hmm_column->i_emit[c] +  c_hmm_column->I_backward[i]-mb->b_score);
				
				c_hmm_column->short_transition_e[MI] = logsum(c_hmm_column->short_transition_e[MI] , c_hmm_column->M_foward[i-1] + c_hmm_column->short_transition[MI] +  c_hmm_column->i_emit[c] +  c_hmm_column->I_backward[i]-mb->b_score);
				
				
				c_hmm_column->i_emit_e[c] = logsum( c_hmm_column->i_emit_e[c] , c_hmm_column->I_foward[i] + c_hmm_column->I_backward[i] -mb->b_score);
				
				//***************post
				
				
				c_hmm_column->D_foward[i] = prob2scaledprob(0.0f);
				
				// no post???
				
				//
				
				for(g = 1;g < hmm->num_columns;g++){
					c_hmm_column = hmm->hmm_column[g];
					p_hmm_column = hmm->hmm_column[g-1];
					
					//Match state
					//transition from previous match state
					c_hmm_column->M_foward[i] = p_hmm_column->M_foward[i-1] + p_hmm_column->short_transition[MM];
					//transition from previous insert state
					c_hmm_column->M_foward[i] = logsum(c_hmm_column->M_foward[i] , p_hmm_column->I_foward[i-1] + p_hmm_column->short_transition[IM] );
					//transition from previous delete state
					c_hmm_column->M_foward[i] = logsum(c_hmm_column->M_foward[i], p_hmm_column->D_foward[i] + p_hmm_column->short_transition[DM]);
					
					// emission promability in curent M state ;
					c_hmm_column->M_foward[i]  = c_hmm_column->M_foward[i]  + c_hmm_column->m_emit[c];
					
					
					//***************post
					p_hmm_column->short_transition_e[MM] = logsum(p_hmm_column->short_transition_e[MM] , p_hmm_column->M_foward[i-1] + p_hmm_column->short_transition[MM] +  c_hmm_column->m_emit[c] +  c_hmm_column->M_backward[i] -mb->b_score );
					
					p_hmm_column->short_transition_e[IM] = logsum(p_hmm_column->short_transition_e[IM],p_hmm_column->I_foward[i-1] + p_hmm_column->short_transition[IM] +  c_hmm_column->m_emit[c] +  c_hmm_column->M_backward[i] -mb->b_score);
					
					p_hmm_column->short_transition_e[DM] = logsum(p_hmm_column->short_transition_e[DM],p_hmm_column->D_foward[i] + p_hmm_column->short_transition[DM] + c_hmm_column->m_emit[c] +  c_hmm_column->M_backward[i] -mb->b_score);
					
					c_hmm_column->m_emit_e[c] = logsum(c_hmm_column->m_emit_e[c],  c_hmm_column->M_foward[i] + c_hmm_column->M_backward[i] -mb->b_score );
					//***************post
					
					
					
					
					
					// Instertion State ..
					//self loop insertion to insertion
					c_hmm_column->I_foward[i] = c_hmm_column->I_foward[i-1] + c_hmm_column->short_transition[II];
					// start new insertion
					c_hmm_column->I_foward[i]  = logsum ( c_hmm_column->I_foward[i] ,c_hmm_column->M_foward[i-1] + c_hmm_column->short_transition[MI]);
					
					//instertion emission...
					c_hmm_column->I_foward[i]  = c_hmm_column->I_foward[i]  + c_hmm_column->i_emit[c];
					
					//***************post
					c_hmm_column->short_transition_e[II] = logsum(c_hmm_column->short_transition_e[II] ,  c_hmm_column->I_foward[i-1] + c_hmm_column->short_transition[II] + c_hmm_column->i_emit[c] + c_hmm_column->I_backward[i] - mb->b_score);
					
					c_hmm_column->short_transition_e[MI] = logsum(c_hmm_column->short_transition_e[MI] , c_hmm_column->M_foward[i-1] + c_hmm_column->short_transition[MI] + c_hmm_column->i_emit[c] + c_hmm_column->I_backward[i] - mb->b_score);
					
					c_hmm_column->i_emit_e[c] = logsum(c_hmm_column->i_emit_e[c] , c_hmm_column->I_foward[i]   + c_hmm_column->I_backward[i] - mb->b_score);
					//***************post

					
					
					// deletion state
					//from previous match state.
					c_hmm_column->D_foward[i] = p_hmm_column->M_foward[i] + p_hmm_column->short_transition[MD];
					//from previous delete state
					
					c_hmm_column->D_foward[i] = logsum(c_hmm_column->D_foward[i], p_hmm_column->D_foward[i] + p_hmm_column->short_transition[DD] );
					
					//***************post
					p_hmm_column->short_transition_e[MD] = logsum(p_hmm_column->short_transition_e[MD],  p_hmm_column->M_foward[i] + p_hmm_column->short_transition[MD] + c_hmm_column->D_backward[i] - mb->b_score);
					
					p_hmm_column->short_transition_e[DD] = logsum(p_hmm_column->short_transition_e[DD] ,p_hmm_column->D_foward[i] + p_hmm_column->short_transition[DD] + c_hmm_column->D_backward[i] - mb->b_score);
					//***************post
					
					
				}
				//fprintf(stderr,"%d	%f	%f	%f	%f	%f\n",i,  c_hmm_column->M_foward[i] , mb->model[j]->M_to_silent[f],  c_hmm_column->I_foward[i] , mb->model[j]->I_to_silent[f], psilent[i]);
				csilent[i] =  logsum(csilent[i],  c_hmm_column->M_foward[i] + mb->model[j]->M_to_silent[f]);
				
				
				
				csilent[i] =  logsum(csilent[i],  c_hmm_column->I_foward[i] + mb->model[j]->I_to_silent[f]);
				csilent[i] = logsum(csilent[i], psilent[i] + mb->model[j]->skip);
				
				//***************post
				mb->model[j]->M_to_silent_e[f] = logsum(mb->model[j]->M_to_silent_e[f],  c_hmm_column->M_foward[i] + mb->model[j]->M_to_silent[f] + bsilent[i+1] -mb->b_score);
				
				//fprintf(stderr,"ADDED TO M->S model: %d		%f %f %f %f %f\n", j , scaledprob2prob( mb->model[j]->M_to_silent_e[f]) ,scaledprob2prob(c_hmm_column->M_foward[i]) , scaledprob2prob(mb->model[j]->M_to_silent[f] ),scaledprob2prob( bsilent[i+1]) ,scaledprob2prob( mb->b_score));
				
				mb->model[j]->I_to_silent_e[f] = logsum(mb->model[j]->I_to_silent_e[f] , c_hmm_column->I_foward[i] + mb->model[j]->I_to_silent[f] + bsilent[i+1] -mb->b_score);
				
				mb->model[j]->skip_e =logsum(mb->model[j]->skip_e , psilent[i] + mb->model[j]->skip + bsilent[i+1] -mb->b_score);
 				
				
				
				//***************post
				
			}
			
		}
	}
	
	
	mb->f_score = mb->model[mb->num_models-1]->silent_forward[len];
	
	fprintf(stderr,"SCORE:%f	%f\n", mb->f_score, scaledprob2prob(mb->f_score));
	return mb;
}



struct model_bag* forward_max_posterior_decoding(struct model_bag* mb, char* a, int len)
{
	
	int i,j,c;
	int f,g;
	
	int hmm_counter = 0;
	
	struct hmm* hmm = 0;
	struct hmm_column* c_hmm_column = 0;
	struct hmm_column* p_hmm_column = 0;
	
	char* seqa = a -1;
	
	float* psilent;
	float* csilent;
	float* bsilent;
	
	float previous_silent[MAX_HMM_SEQ_LEN];
	float next_silent[MAX_HMM_SEQ_LEN];
	
	// M state of first set of HMMS.....
	for(j = 0; j < mb->num_models;j++){
		for(f = 0;f < mb->model[j]->num_hmms;f++){
			for(g = 0;g < mb->model[j]->hmms[f]->num_columns;g++){
				c_hmm_column = mb->model[j]->hmms[f]->hmm_column[g];
				//i = 0;
				for(i = 0; i <= len;i++){
					c_hmm_column->M_foward[i] = prob2scaledprob(0.0);
					c_hmm_column->I_foward[i] = prob2scaledprob(0.0);
					c_hmm_column->D_foward[i] = prob2scaledprob(0.0);
					//fprintf(stderr,"segment:%d	HMM:%d	column:%d	i:%d	%f	%f	%f\n",j,f,g,i, c_hmm_column->M_foward[i] ,c_hmm_column->I_foward[i],c_hmm_column->D_foward[i]    );
				}
			}
		}
		for(i = 0; i <= len+1;i++){
			mb->model[j]->silent_forward[i] = prob2scaledprob(0.0f);
		}
	}
	
	//fprintf(stderr,"\n\n\n");
	mb->model[0]->silent_forward[0] = prob2scaledprob(1.0) + mb->model[0]->skip;
	//fprintf(stderr,"Init silent states... \n");
	//fprintf(stderr,"%d	%f\n",0,mb->model[0]->silent_forward[0]   );
	for(j = 1; j < mb->num_models;j++){
		mb->model[j]->silent_forward[0] = mb->model[j-1]->silent_forward[0]  + mb->model[j]->skip ;
		//mb->model[j]->silent_forward[0] = mb->model[j-1]->silent_forward[0]  + mb->model[j-1]->skip ;
	}
	
	for(i = 0; i <= len;i++){
		for(j = 0; j < mb->total_hmm_num;j++){
			mb->dyn_prog_matrix[i][j] = prob2scaledprob(0.0f);
			mb->path[i][j] = -1;
		}
	}

	
	
	
 	for(i = 0; i <= len;i++){
		
		
		next_silent[i] = prob2scaledprob(0.0f);
		previous_silent[i] = prob2scaledprob(0.0f);
	}
	previous_silent[0] = prob2scaledprob(1.0);
	next_silent[len+1] = prob2scaledprob(1.0f);
	
	//loop thorugh the segments
	// in each run the contained HMMS and update silent states;
	for(j = 0; j < mb->num_models;j++){
		if(j == 0){
			psilent = previous_silent;
		}else{
			psilent =  mb->model[j-1]->silent_forward;
		}
		csilent = mb->model[j]->silent_forward;
		if(j +1 != mb->num_models){
			bsilent = mb->model[j+1]->silent_backward;
		}else{
			bsilent = next_silent;
		}
		for(f = 0;f < mb->model[j]->num_hmms;f++){
			
			fprintf(stderr," %d %d %d\n", j , f, hmm_counter);
			hmm = mb->model[j]->hmms[f];
			for(i = 1; i <= len;i++){
				c = seqa[i];
				
				c_hmm_column = hmm->hmm_column[0];
				// first column  comes from previous state cheekily transferring its pd to M[0[
				c_hmm_column->M_foward[i] = psilent[i-1] + mb->model[j]->silent_to_M[f] + c_hmm_column->m_emit[c];
				
				//***************post
				//mb->model[j]->silent_to_M_e[f] = logsum(mb->model[j]->silent_to_M_e[f] ,psilent[i-1] + mb->model[j]->silent_to_M[f] + c_hmm_column->m_emit[c] + c_hmm_column->M_backward[i]  -mb->b_score);
				
				
				mb->dyn_prog_matrix[i][hmm_counter] = logsum(mb->dyn_prog_matrix[i][hmm_counter],  c_hmm_column->M_foward[i]  + c_hmm_column->M_backward[i] -mb->b_score );
				
				//c_hmm_column->m_emit_e[c] = logsum(c_hmm_column->m_emit_e[c] , c_hmm_column->M_foward[i]  + c_hmm_column->M_backward[i] -mb->b_score);
				//***************post
				
				
 				c_hmm_column->I_foward[i]    =psilent[i-1] + mb->model[j]->silent_to_I[f] ;
				
				
				
				
				//add transitions to first columns////
				c_hmm_column->I_foward[i] = logsum(c_hmm_column->I_foward[i], c_hmm_column->I_foward[i-1] + c_hmm_column->short_transition[II]);
				
				c_hmm_column->I_foward[i]  = logsum ( c_hmm_column->I_foward[i] ,c_hmm_column->M_foward[i-1] + c_hmm_column->short_transition[MI]);
				
				c_hmm_column->I_foward[i] = c_hmm_column->I_foward[i] + c_hmm_column->i_emit[c];
				
				
				
				
				//***************post
				
				mb->dyn_prog_matrix[i][hmm_counter] = logsum(mb->dyn_prog_matrix[i][hmm_counter],  c_hmm_column->I_foward[i] + c_hmm_column->I_backward[i] -mb->b_score );
				
				
				//***************post
				
				
				c_hmm_column->D_foward[i] = prob2scaledprob(0.0f);
				
				// no post???
				
				//
				
				for(g = 1;g < hmm->num_columns;g++){
					c_hmm_column = hmm->hmm_column[g];
					p_hmm_column = hmm->hmm_column[g-1];
					
					//Match state
					//transition from previous match state
					c_hmm_column->M_foward[i] = p_hmm_column->M_foward[i-1] + p_hmm_column->short_transition[MM];
					//transition from previous insert state
					c_hmm_column->M_foward[i] = logsum(c_hmm_column->M_foward[i] , p_hmm_column->I_foward[i-1] + p_hmm_column->short_transition[IM] );
					//transition from previous delete state
					c_hmm_column->M_foward[i] = logsum(c_hmm_column->M_foward[i], p_hmm_column->D_foward[i] + p_hmm_column->short_transition[DM]);
					
					// emission promability in curent M state ;
					c_hmm_column->M_foward[i]  = c_hmm_column->M_foward[i]  + c_hmm_column->m_emit[c];
					
					
					//***************post
					mb->dyn_prog_matrix[i][hmm_counter] = logsum(mb->dyn_prog_matrix[i][hmm_counter],  c_hmm_column->M_foward[i]  + c_hmm_column->M_backward[i] -mb->b_score );
					//***************post
					
					
					
					
					
					// Instertion State ..
					//self loop insertion to insertion
					c_hmm_column->I_foward[i] = c_hmm_column->I_foward[i-1] + c_hmm_column->short_transition[II];
					// start new insertion
					c_hmm_column->I_foward[i]  = logsum ( c_hmm_column->I_foward[i] ,c_hmm_column->M_foward[i-1] + c_hmm_column->short_transition[MI]);
					
					//instertion emission...
					c_hmm_column->I_foward[i]  = c_hmm_column->I_foward[i]  + c_hmm_column->i_emit[c];
					
					//***************post
					mb->dyn_prog_matrix[i][hmm_counter] = logsum(mb->dyn_prog_matrix[i][hmm_counter],  c_hmm_column->I_foward[i] + c_hmm_column->I_backward[i] -mb->b_score );
					//***************post
					
					
					
					// deletion state
					//from previous match state.
					c_hmm_column->D_foward[i] = p_hmm_column->M_foward[i] + p_hmm_column->short_transition[MD];
					//from previous delete state
					
					c_hmm_column->D_foward[i] = logsum(c_hmm_column->D_foward[i], p_hmm_column->D_foward[i] + p_hmm_column->short_transition[DD] );
					
					//***************post
					p_hmm_column->short_transition_e[MD] = logsum(p_hmm_column->short_transition_e[MD],  p_hmm_column->M_foward[i] + p_hmm_column->short_transition[MD] + c_hmm_column->D_backward[i] - mb->b_score);
					
					p_hmm_column->short_transition_e[DD] = logsum(p_hmm_column->short_transition_e[DD] ,p_hmm_column->D_foward[i] + p_hmm_column->short_transition[DD] + c_hmm_column->D_backward[i] - mb->b_score);
					//***************post
					
					
				}
				//fprintf(stderr,"%d	%f	%f	%f	%f	%f\n",i,  c_hmm_column->M_foward[i] , mb->model[j]->M_to_silent[f],  c_hmm_column->I_foward[i] , mb->model[j]->I_to_silent[f], psilent[i]);
				csilent[i] =  logsum(csilent[i],  c_hmm_column->M_foward[i] + mb->model[j]->M_to_silent[f]);
				
				
				
				csilent[i] =  logsum(csilent[i],  c_hmm_column->I_foward[i] + mb->model[j]->I_to_silent[f]);
				csilent[i] = logsum(csilent[i], psilent[i] + mb->model[j]->skip);
			}
			hmm_counter++;
			
		}
		//hmm_counter++;
	}
	
	
	mb->f_score = mb->model[mb->num_models-1]->silent_forward[len];
	
	for(i = 0; i <= len;i++){
		fprintf(stderr,"%d ",i);
		for(j = 0; j < mb->total_hmm_num;j++){
			fprintf(stderr,"%0.3f ",scaledprob2prob( mb->dyn_prog_matrix[i][j]));
			mb->dyn_prog_matrix[i][j] = scaledprob2prob( mb->dyn_prog_matrix[i][j]);
		}
		fprintf(stderr,"\n");
	}
	
	float max = 0;
	float tmp;
	int move = -1;
	
	for(i = 1;i <= len;i++){
		for(j = 0; j < mb->total_hmm_num;j++){
			max = -1;
			for(c = 0 ;c <= j ;c++){
				tmp =  mb->dyn_prog_matrix[i-1][c] * mb->transition_matrix[c][j];
				
				if(tmp > max){
					move = c;
					max = tmp;
				}
			
			//	fprintf(stderr,"%0.3f ",scaledprob2prob( mb->dyn_prog_matrix[i][j]));
			}
			mb->dyn_prog_matrix[i][j]+= max;
			mb->path[i][j] = move;
		}
	}
	fprintf(stderr,"MATRIX:\n");
	for(i = 0; i <= len;i++){
		fprintf(stderr,"%d ",i);
		for(j = 0; j < mb->total_hmm_num;j++){
			fprintf(stderr,"%0.3f ", mb->dyn_prog_matrix[i][j]);
			//mb->dyn_prog_matrix[i][j] = scaledprob2prob( mb->dyn_prog_matrix[i][j]);
		}
		fprintf(stderr,"\n");
	}
	fprintf(stderr,"PATH:\n");
	for(i = 0; i <= len;i++){
		fprintf(stderr,"%d ",i);
		for(j = 0; j < mb->total_hmm_num;j++){
			fprintf(stderr,"%d ", mb->path[i][j]);
			//mb->dyn_prog_matrix[i][j] = scaledprob2prob( mb->dyn_prog_matrix[i][j]);
		}
		fprintf(stderr,"\n");
	}
	
	char path[100];
	
	i = len;
	max = -1;
	for(j = 0; j < mb->total_hmm_num;j++){
		if(mb->dyn_prog_matrix[i][j] > max){
			max = mb->dyn_prog_matrix[i][j];
			move = j;
		}
	}
	
	for(i = 0; i <= len;i++){
		path[i] = 0;
		//fprintf(stderr,"%d %d\n",i, path[i]);
	}
	
	path[len] = move;
	
	
	for(i = len ;i > 0;i--){
		move = mb->path[i][move];
		path[i-1] = move;
	}
	
	for(i = 0; i <= len;i++){
		fprintf(stderr,"%d %d\n",i, (int)path[i]);
	}
	
	
	
	
	
	
	
	fprintf(stderr,"SCORE:%f	%f\n", mb->f_score, scaledprob2prob(mb->f_score));
	return mb;
}



struct model* malloc_model(int main_length, int sub_length, int number_sub_models)
{
	struct model* model = NULL;
	int i = 0;
	int j = 0;
	
	
	assert(number_sub_models  <=MAX_NUM_SUB_MODELS );
	
	model = malloc(sizeof(struct model));
	assert(model != 0);
	
	model->num_hmms =  (1+ number_sub_models);
	model->hmms = malloc(sizeof(struct hmm*) * (1+ number_sub_models));
	assert(model->hmms !=0);
	for(i = 0; i < model->num_hmms;i++){
		model->hmms[i] = malloc(sizeof(struct hmm) );
		assert(model->hmms[i]  != 0);
	}
	
	
		
	
	
	model->hmms[0]->num_columns = main_length;
	model->hmms[0]->hmm_column = malloc(sizeof(struct hmm_column*) * main_length);
	assert(model->hmms[0]->hmm_column !=0);
	
	for(j = 0; j < main_length;j++){
		model->hmms[0]->hmm_column[j] = malloc(sizeof(struct hmm_column));
		assert(model->hmms[0]->hmm_column[j] !=0);
	}

	for(i = 1; i < model->num_hmms;i++){
		model->hmms[i]->num_columns = sub_length;
		model->hmms[i]->hmm_column = malloc(sizeof(struct hmm_column*) * sub_length);
		assert(model->hmms[i]->hmm_column !=0);
		for(j = 0; j < sub_length;j++){
			model->hmms[i]->hmm_column[j] = malloc(sizeof(struct hmm_column));
			assert(model->hmms[i]->hmm_column[j] !=0);
			//model->hmms[i]->hmm_column[j]->identifier = -1;
		}
	}
	
	
	
	return model;
}


struct model* malloc_model_according_to_read_structure(struct read_structure* rs, int key)
{
	struct model* model = NULL;
	int i = 0;
	int j = 0;
	int len = 0;

	model = malloc(sizeof(struct model));
	assert(model != 0);
	
	model->num_hmms =  (rs->numseq_in_segment[key]);
	model->hmms = malloc(sizeof(struct hmm*) * (rs->numseq_in_segment[key]));
	assert(model->hmms !=0);
	
	for(i = 0; i < model->num_hmms;i++){
		model->hmms[i] = malloc(sizeof(struct hmm) );
		assert(model->hmms[i]  != 0);
	}
	model->M_to_silent = malloc(sizeof(float) * model->num_hmms);
	model->silent_to_M = malloc(sizeof(float) * model->num_hmms);
	
	model->I_to_silent = malloc(sizeof(float) * model->num_hmms);
	model->silent_to_I = malloc(sizeof(float) * model->num_hmms);
	
	
	model->M_to_silent_e = malloc(sizeof(float) * model->num_hmms);
	model->silent_to_M_e = malloc(sizeof(float) * model->num_hmms);
	
	model->I_to_silent_e = malloc(sizeof(float) * model->num_hmms);
	model->silent_to_I_e = malloc(sizeof(float) * model->num_hmms);
	
	for(i = 0 ;i  < model->num_hmms;i++){
		model->M_to_silent[i] = 0.0f;
		model->silent_to_M[i] = 0.0f;
		model->I_to_silent[i] = 0.0f;
		model->silent_to_I[i] = 0.0f;
	}
	
	len = (int)strlen(rs->sequence_matrix[key][0]);
	for(i = 0; i < model->num_hmms;i++){
		model->hmms[i]->num_columns = len;
		model->hmms[i]->hmm_column = malloc(sizeof(struct hmm_column*) * len);
		assert(model->hmms[i]->hmm_column != 0);
		for(j = 0; j < len;j++){
			model->hmms[i]->hmm_column[j] = malloc(sizeof(struct hmm_column));
			assert(model->hmms[i]->hmm_column[j] != 0);
		}
	}
	return model;
}

struct model* init_model_according_to_read_structure(struct model* model,struct parameters* param , int key, float* background,int assumed_length)
{

	struct read_structure* rs = param->read_structure;
	float base_error = param->sequencer_error_rate;
	float indel_freq = param->indel_frequency;
	struct hmm_column* col =0;
	int i,j,c,len;
	int current_nuc;
	char* tmp = 0;
	
	for(i= 0;i < 5;i++){
		model->background_nuc_frequency[i]= background[i];
		//fprintf(stderr,"%f\n",background[i]);
	}
	
	for(i = 0; i < model->num_hmms;i++){
		len = model->hmms[i]->num_columns;
		tmp = rs->sequence_matrix[key][i];
		for(j = 0; j < len;j++){
			col = model->hmms[i]->hmm_column[j];
			
			current_nuc = nuc_code5[(int) tmp[j]];
			col->identifier = -1;
			if(current_nuc != 4){
				
				for(c = 0; c < 5;c++){
					if(c == current_nuc){
						col->m_emit[c] = prob2scaledprob(1.0 - base_error* (1.0- indel_freq));
					}else{
						col->m_emit[c] =  prob2scaledprob( base_error* (1.0- indel_freq)/ 4.0);
					}
					col->i_emit[c] = background[c];
					col->i_emit_e[c] =  prob2scaledprob(0.0f);
					col->m_emit_e[c] =  prob2scaledprob(0.0f);
				}
			}else{
				for(c = 0; c < 5;c++){
	
					col->m_emit[c] =  background[c];
					col->i_emit[c] =  background[c];
					col->i_emit_e[c] =  prob2scaledprob(0.0f);
					col->m_emit_e[c] =  prob2scaledprob(0.0f);
				}
			}
			
			if(j == len-1){
				col->short_transition[MM] = prob2scaledprob(0.0f );// 1.0 - base_error * indel_freq);
				col->short_transition[MI] = prob2scaledprob(0.0f );//(base_error * indel_freq) +  prob2scaledprob(0.5);
				col->short_transition[MD] = prob2scaledprob(0.0f );//(base_error * indel_freq) +  prob2scaledprob(0.5);
				
				col->short_transition[II] = prob2scaledprob(1.0 - 0.999);
				col->short_transition[IM] = prob2scaledprob(0.999);
				
				col->short_transition[DD] = prob2scaledprob(0.0f );//(1.0 - 0.999);
				col->short_transition[DM] = prob2scaledprob(0.0f );//0.999);
			}else{
			
			col->short_transition[MM] = prob2scaledprob( 1.0 - base_error * indel_freq);
			col->short_transition[MI] = prob2scaledprob(base_error * indel_freq) +  prob2scaledprob(0.5);
			col->short_transition[MD] = prob2scaledprob(base_error * indel_freq) +  prob2scaledprob(0.5);
			
			col->short_transition[II] = prob2scaledprob(1.0 - 0.999);
			col->short_transition[IM] = prob2scaledprob(0.999);
			
			col->short_transition[DD] = prob2scaledprob(1.0 - 0.999);
			col->short_transition[DM] = prob2scaledprob(0.999);
			}
			
			
			col->short_transition_e[MM] =  prob2scaledprob(0.0);
			col->short_transition_e[MI] =  prob2scaledprob(0.0);
			col->short_transition_e[MD] =  prob2scaledprob(0.0);
			
			col->short_transition_e[II] =  prob2scaledprob(0.0);
			col->short_transition_e[IM] =  prob2scaledprob(0.0);
			
			col->short_transition_e[DD] =  prob2scaledprob(0.0);
			col->short_transition_e[DM] =  prob2scaledprob(0.0);
		}
		
	}
	
	// init all probs to 0
	
	for(i = 0 ; i < model->num_hmms;i++){
		model->silent_to_M[i] = prob2scaledprob(0.0f);
		model->M_to_silent[i] = prob2scaledprob(0.0f);
		
		model->silent_to_I[i] = prob2scaledprob(0.0f);
		model->I_to_silent[i] = prob2scaledprob(0.0f);
		
		model->silent_to_M_e[i] = prob2scaledprob(0.0f);
		model->M_to_silent_e[i] = prob2scaledprob(0.0f);
		
		model->silent_to_I_e[i] = prob2scaledprob(0.0f);
		model->I_to_silent_e[i] = prob2scaledprob(0.0f);

		
	}
	model->skip = prob2scaledprob(0.0f);
	model->skip_e = prob2scaledprob(0.0f);
	
	
	model->random_next = prob2scaledprob(0.0);
	model->random_self = prob2scaledprob(0.0);
	
	if(rs->type[key] == 'B'){// barcodes all have same length & equal prior probability... 
		for(i = 0 ; i < model->num_hmms;i++){
			model->silent_to_M[i] = prob2scaledprob(1.0 / (float) model->num_hmms);
			model->M_to_silent[i] = prob2scaledprob(1.0);
			
			model->silent_to_I[i] = prob2scaledprob(0.0f);
			model->I_to_silent[i] = prob2scaledprob(0.0f);
			
			
		}
		model->skip = prob2scaledprob(0.0);
	}
	
	if(rs->type[key] == 'F'){// fingerprint all have same length & equal prior probability... (of course we specify 1 with NNNNNNNN
		for(i = 0 ; i < model->num_hmms;i++){
			model->silent_to_M[i] = prob2scaledprob(1.0 / (float) model->num_hmms);
			model->M_to_silent[i] = prob2scaledprob(1.0);
		}
		model->skip = prob2scaledprob(0.0);
	}
		
	if(rs->type[key] == 'O'){ // optional - like a G, GG or GGG priot probability set to 0.5  - assume length 2 for now,
		len = model->hmms[0]->num_columns;
		for(i = 0 ; i < model->num_hmms;i++){
			//model->silent_to_M[i] = prob2scaledprob(1.0 / (float) model->num_hmms) + prob2scaledprob(0.5);
			//model->M_to_silent[i] = prob2scaledprob(1.0);
			
			model->silent_to_I[i] = prob2scaledprob(1.0 / (float) model->num_hmms) + prob2scaledprob(0.5);
			model->I_to_silent[i] = prob2scaledprob(1.0 / (float) (len+1));
			
			//len = model->hmms[i]->num_columns;
			//tmp = rs->sequence_matrix[key][i];
			for(j = 0; j < len;j++){
				col = model->hmms[i]->hmm_column[j];
				for(c = 0; c < 5;c++){
					col->i_emit[c] = col->m_emit[c];
					col->m_emit[c] = prob2scaledprob(0.0);
				}
			}
		}
		model->skip = prob2scaledprob(0.5);
		col = model->hmms[0]->hmm_column[0];
		col->short_transition[MM] = prob2scaledprob( 0.0 );
		col->short_transition[MI] = prob2scaledprob(0.0);
		col->short_transition[MD] = prob2scaledprob(0.0);
		
		//col->short_transition[MQUIT] = prob2scaledprob(1.0 / (float) 2);
		
		col->short_transition[II] = prob2scaledprob(1.0 - 1.0 / (float)(len+1) );
		col->short_transition[IM] = prob2scaledprob(0.0);
		
		col->short_transition[DD] = prob2scaledprob(0.0);
		col->short_transition[DM] = prob2scaledprob(0.0);
		
		
		col->short_transition_e[MM] =  prob2scaledprob(0.0);
		col->short_transition_e[MI] =  prob2scaledprob(0.0);
		col->short_transition_e[MD] =  prob2scaledprob(0.0);
		
		col->short_transition_e[II] =  prob2scaledprob(0.0);
		col->short_transition_e[IM] =  prob2scaledprob(0.0);
		
		col->short_transition_e[DD] =  prob2scaledprob(0.0);
		col->short_transition_e[DM] =  prob2scaledprob(0.0);
	}
	
	if(rs->type[key] == 'R'){// read - skip impossible; 
		for(i = 0 ; i < model->num_hmms;i++){
			model->silent_to_I[i] = prob2scaledprob(1.0 / (float) model->num_hmms);
			model->I_to_silent[i] = prob2scaledprob(1.0 / (float) assumed_length);
		}
		col = model->hmms[0]->hmm_column[0];
		for(c = 0; c < 5;c++){
			
			col->m_emit[c] =background[c];
			
			col->i_emit[c] = background[c];
			col->i_emit_e[c] =  prob2scaledprob(0.0f);
			col->m_emit_e[c] =  prob2scaledprob(0.0f);
		}
		col->short_transition[MM] = prob2scaledprob( 0.0);
		col->short_transition[MI] = prob2scaledprob(0.0);
		col->short_transition[MD] = prob2scaledprob(0.0);
		
		//col->short_transition[MQUIT] = prob2scaledprob(1.0 / (float) assumed_length);
		
		col->short_transition[II] = prob2scaledprob(1.0 - 1.0 / (float) assumed_length );
		col->short_transition[IM] = prob2scaledprob(0.0);
		
		col->short_transition[DD] = prob2scaledprob(0.0);
		col->short_transition[DM] = prob2scaledprob(0.0);
		
		
		col->short_transition_e[MM] =  prob2scaledprob(0.0);
		col->short_transition_e[MI] =  prob2scaledprob(0.0);
		col->short_transition_e[MD] =  prob2scaledprob(0.0);
		
		col->short_transition_e[II] =  prob2scaledprob(0.0);
		col->short_transition_e[IM] =  prob2scaledprob(0.0);
		
		col->short_transition_e[DD] =  prob2scaledprob(0.0);
		col->short_transition_e[DM] =  prob2scaledprob(0.0);
		
		model->skip = prob2scaledprob(0.0);
		
	}
	
	
	return model;
}

void print_model(struct model* model)
{
	int i,j,c;
	int len;
	float sum = 0;
	struct hmm_column* col =0;
	fprintf(stderr,"Skip:%f Self:%f Next:%f\n",scaledprob2prob(model->skip), scaledprob2prob(model->random_self) , scaledprob2prob(model->random_next));
	for(i = 0; i < model->num_hmms;i++){
		fprintf(stderr,"HMM%d:		(get there:%f)\n",i, scaledprob2prob(model->silent_to_M[i]));
		
		
		
		len = model->hmms[i]->num_columns;
		//tmp = rs->sequence_matrix[key][i];
		for(j = 0; j < len;j++){
			col = model->hmms[i]->hmm_column[j];
			fprintf(stderr,"\t%d\t",j);
			sum = 0.0;
			for(c = 0; c < 5;c++){
				sum += scaledprob2prob(col->m_emit[c]);
				fprintf(stderr,"%0.4f ", scaledprob2prob(col->m_emit[c]));
				
			}
			fprintf(stderr,"%0.1fc ",sum);
			
			sum = 0.0;
			for(c = 0; c < 5;c++){
				sum += scaledprob2prob(col->i_emit[c]);
				fprintf(stderr,"%0.4f ", scaledprob2prob(col->i_emit[c]));
				
			}
			fprintf(stderr,"%0.1fc ",sum);
			
			for(c = 0; c < 7;c++){
				fprintf(stderr,"%0.4f ", scaledprob2prob(col->short_transition[c]));

			}
			fprintf(stderr,"%0.3f %0.3f %0.3f\n",scaledprob2prob(col->short_transition[MM])+scaledprob2prob(col->short_transition[MI])+scaledprob2prob(col->short_transition[MD])  ,scaledprob2prob(col->short_transition[II])+scaledprob2prob(col->short_transition[IM]),  scaledprob2prob(col->short_transition[DD]) + scaledprob2prob(col->short_transition[DM]) );
		}
	}
	fprintf(stderr," ESTIMATED::::: \n");
	for(i = 0; i < model->num_hmms;i++){
		fprintf(stderr,"HMM%d:		(get there:%f)\n",i, scaledprob2prob(model->silent_to_M[i]));
		
		
		
		len = model->hmms[i]->num_columns;
		//tmp = rs->sequence_matrix[key][i];
		for(j = 0; j < len;j++){
			col = model->hmms[i]->hmm_column[j];
			fprintf(stderr,"\t%d\t",j);
			sum = 0.0;
			for(c = 0; c < 5;c++){
				sum += scaledprob2prob(col->m_emit_e[c]);
				fprintf(stderr,"%0.4f ", scaledprob2prob(col->m_emit_e[c]));
				
			}
			fprintf(stderr,"%0.1fc ",sum);
			
			sum = 0.0;
			for(c = 0; c < 5;c++){
				sum += scaledprob2prob(col->i_emit_e[c]);
				fprintf(stderr,"%0.4f ", scaledprob2prob(col->i_emit_e[c]));
				
			}
			fprintf(stderr,"%0.1fc ",sum);
			
			for(c = 0; c < 7;c++){
				fprintf(stderr,"%0.4f ", scaledprob2prob(col->short_transition_e[c]));
				
			}
			fprintf(stderr,"%0.3f %0.3f %0.3f\n",scaledprob2prob(col->short_transition_e[MM])+scaledprob2prob(col->short_transition_e[MI])+scaledprob2prob(col->short_transition_e[MD])  ,scaledprob2prob(col->short_transition_e[II])+scaledprob2prob(col->short_transition_e[IM]),  scaledprob2prob(col->short_transition_e[DD]) + scaledprob2prob(col->short_transition_e[DM]) );
		}
		
		
		
		
	}
	fprintf(stderr,"Links:silent to\n");
	
	for(i = 0; i < model->num_hmms;i++){
		fprintf(stderr,"%d	%f	%f	%f	%f\n",i, scaledprob2prob(  model->silent_to_M[i]), scaledprob2prob(  model->silent_to_I[i]),scaledprob2prob(   model->silent_to_M_e[i]),scaledprob2prob(  model->silent_to_I_e[i]));
	}
	fprintf(stderr,"Links:to silent \n");
	
	for(i = 0; i < model->num_hmms;i++){
		fprintf(stderr,"%d	%f	%f	%f	%f\n",i, scaledprob2prob(  model->M_to_silent[i]), scaledprob2prob(  model->I_to_silent[i]),scaledprob2prob(   model->M_to_silent_e[i]),scaledprob2prob(  model->I_to_silent_e[i]));
	}
	
	fprintf(stderr,"SKIP:\n");
	fprintf(stderr,"%f	%f\n", scaledprob2prob(model->skip) , scaledprob2prob(model->skip_e));
	

	
}




void free_model(struct model* model)
{
	int i = 0;
	int j = 0;

		
	for(i = 0; i < model->num_hmms;i++){
		
		for(j = 0; j < model->hmms[i]->num_columns;j++){
			free(model->hmms[i]->hmm_column[j]);// = malloc(sizeof(struct hmm_column));
			//assert(model->hmms[i]->hmm_column[j] !=0);
			//model->hmms[i]->hmm_column[j]->identifier = -1;
			
		}
		//model->hmms[i]->num_columns = sub_length;
		free(model->hmms[i]->hmm_column);// = malloc(sizeof(struct hmm_column*) * sub_length);
		//assert(model->hmms[i]->hmm_column !=0);
	}
	//for(j = 0; j < 	model->hmms[0]->num_columns;j++){
	//	free(model->hmms[0]->hmm_column[j]);// = malloc(sizeof(struct hmm_column));
		//assert(model->hmms[0]->hmm_column[j] !=0);
	//}

	//model->hmms[0]->num_columns = main_length;
	//free(model->hmms[0]->hmm_column);// = malloc(sizeof(struct hmm_column*) * main_length);
	//assert(model->hmms[0]->hmm_column !=0);
	
	///assert(model->hmms !=0);
	

	for(i = 0; i < model->num_hmms;i++){
		free(model->hmms[i]);// = malloc(sizeof(struct hmm) );
		//assert(model->hmms[i]  != 0);
	}
	free(model->hmms);// = malloc(sizeof(struct hmm*) * (1+ number_sub_models));
	if(model->silent_to_M){
		free(model->silent_to_M);
		free(model->silent_to_M_e);
	}
	if(model->M_to_silent){
		free(model->M_to_silent );
		free(model->M_to_silent_e );
	}
	if(model->silent_to_I){
		free(model->silent_to_I);
		free(model->silent_to_I_e);
	}
	if(model->I_to_silent){
		free(model->I_to_silent );
		free(model->I_to_silent_e );
	}
	
	free(model);// = malloc(sizeof(struct model));
	//assert(model != 0);

	
	//return model;
}





struct model* copy_and_malloc_model(struct model* org)
{
	struct model* model = NULL;
	int i = 0;
	int j = 0;
	int len = 0;
	int c = 0;

	struct hmm_column* col =0;

	
	
	model = malloc(sizeof(struct model));
	assert(model != 0);
	
	model->num_hmms = org->num_hmms;
	model->hmms = malloc(sizeof(struct hmm*) * org->num_hmms);
	assert(model->hmms !=0);
	for(i = 0; i < model->num_hmms;i++){
		model->hmms[i] = malloc(sizeof(struct hmm) );
		assert(model->hmms[i]  != 0);
	}
	
	
	
	
	
	model->hmms[0]->num_columns = org->hmms[0]->num_columns;
	model->hmms[0]->hmm_column = malloc(sizeof(struct hmm_column*) * model->hmms[0]->num_columns );
	assert(model->hmms[0]->hmm_column !=0);
	
	for(j = 0; j < model->hmms[0]->num_columns ;j++){
		model->hmms[0]->hmm_column[j] = malloc(sizeof(struct hmm_column));
		assert(model->hmms[0]->hmm_column[j] !=0);
	}
	
	for(i = 1; i < model->num_hmms;i++){
		model->hmms[i]->num_columns = org->hmms[i]->num_columns;
		model->hmms[i]->hmm_column = malloc(sizeof(struct hmm_column*) * model->hmms[i]->num_columns);
		assert(model->hmms[i]->hmm_column !=0);
		for(j = 0; j < model->hmms[i]->num_columns;j++){
			model->hmms[i]->hmm_column[j] = malloc(sizeof(struct hmm_column));
			assert(model->hmms[i]->hmm_column[j] !=0);
			//model->hmms[i]->hmm_column[j]->identifier = -1;
		}
	}
	
	len = model->hmms[0]->num_columns;
	
	for(i = 0; i < len;i++){
		col = model->hmms[0]->hmm_column[i];
		
		//col->identifier = -1;
		for(j = 0; j < 5;j++){
			col->i_emit[j] =  org->hmms[0]->hmm_column[i]->i_emit[j];// background[j];
			col->m_emit[j] =  org->hmms[0]->hmm_column[i]->m_emit[j];
			col->i_emit_e[j] = org->hmms[0]->hmm_column[i]->i_emit_e[j];
			col->m_emit_e[j] = org->hmms[0]->hmm_column[i]->m_emit_e[j];
		}
		col->short_transition[NEXT] = org->hmms[0]->hmm_column[i]->short_transition[NEXT]; //  prob2scaledprob( max);
		col->short_transition[SELF] =  org->hmms[0]->hmm_column[i]->short_transition[SELF];// prob2scaledprob(1.0 - max) +  prob2scaledprob(0.5);
		
		col->short_transition_e[NEXT] = org->hmms[0]->hmm_column[i]->short_transition_e[NEXT]; //  prob2scaledprob( max);
		col->short_transition_e[SELF] =  org->hmms[0]->hmm_column[i]->short_transition_e[SELF];// prob2scaledprob(1.0 - max) +  prob2scaledprob(0.5);
		
		
		for(j = 0; j< model->num_hmms-1;j++){
			col->long_transition[j] = org->hmms[0]->hmm_column[i]->long_transition[j];// prob2scaledprob(1.0 - max) +  prob2scaledprob(0.5  /  (float) (model->num_hmms-1));
			col->long_transition_e[j] = org->hmms[0]->hmm_column[i]->long_transition_e[j];
		}
		
		
	}
	
	/*
	 
	 sub models model....
	 
	 */
	len = model->hmms[1]->num_columns;
	for(i = 1; i< model->num_hmms;i++){
		for(j = 0; j < len;j++){
			col = model->hmms[i]->hmm_column[j];
			///col->identifier = -1;
			for(c = 0; c < 5;c++){
				col->i_emit[c] = org->hmms[i]->hmm_column[j]->i_emit[c];
				col->m_emit[c] = org->hmms[i]->hmm_column[j]->m_emit[c];
				
				col->i_emit_e[c] = org->hmms[i]->hmm_column[j]->i_emit_e[c];
				col->m_emit_e[c] =  org->hmms[i]->hmm_column[j]->m_emit_e[c];
			}
			
			
			
			
			col->short_transition[MM] =   org->hmms[i]->hmm_column[j]->short_transition[MM];//  prob2scaledprob( 0.999);
			col->short_transition[MI] = org->hmms[i]->hmm_column[j]->short_transition[MI];// prob2scaledprob(1.0 - 0.999) +  prob2scaledprob(0.5);
			col->short_transition[MD] = org->hmms[i]->hmm_column[j]->short_transition[MD];// prob2scaledprob(1.0 - 0.999) +  prob2scaledprob(0.5);
			
			col->short_transition[II] = org->hmms[i]->hmm_column[j]->short_transition[II];// prob2scaledprob(1.0 - 0.999);
			col->short_transition[IM] = org->hmms[i]->hmm_column[j]->short_transition[IM];// prob2scaledprob(0.999);
			
			col->short_transition[DD] = org->hmms[i]->hmm_column[j]->short_transition[DD];// prob2scaledprob(1.0 - 0.999);
			col->short_transition[DM] = org->hmms[i]->hmm_column[j]->short_transition[DM];// prob2scaledprob(0.999);
			
			col->short_transition_e[MM] =   org->hmms[i]->hmm_column[j]->short_transition_e[MM];//  prob2scaledprob( 0.999);
			col->short_transition_e[MI] = org->hmms[i]->hmm_column[j]->short_transition_e[MI];// prob2scaledprob(1.0 - 0.999) +  prob2scaledprob(0.5);
			col->short_transition_e[MD] = org->hmms[i]->hmm_column[j]->short_transition_e[MD];// prob2scaledprob(1.0 - 0.999) +  prob2scaledprob(0.5);
			
			col->short_transition_e[II] = org->hmms[i]->hmm_column[j]->short_transition_e[II];// prob2scaledprob(1.0 - 0.999);
			col->short_transition_e[IM] = org->hmms[i]->hmm_column[j]->short_transition_e[IM];// prob2scaledprob(0.999);
			
			col->short_transition_e[DD] = org->hmms[i]->hmm_column[j]->short_transition_e[DD];// prob2scaledprob(1.0 - 0.999);
			col->short_transition_e[DM] = org->hmms[i]->hmm_column[j]->short_transition_e[DM];// prob2scaledprob(0.999);
			
			
			
			for(c = 0;c <  (model->num_hmms-1);c++){
				col->long_transition[c] = org->hmms[i]->hmm_column[j]->long_transition[c];
				col->long_transition_e[c] = org->hmms[i]->hmm_column[j]->long_transition_e[c];
			}
			
			
			
		}
		
	}
	
	return model;
}


struct model* add_estimates_to_model(struct model* target, struct model* source)
{
	
	int i = 0;
	int j = 0;
	int len = 0;
	int c = 0;
	
	struct hmm_column* col_target =0;
	struct hmm_column* col_source =0;
	len = target->hmms[0]->num_columns;
	
	for(i = 0; i < len;i++){
		col_target = target->hmms[0]->hmm_column[i];
		col_source = source->hmms[0]->hmm_column[i];
		
		//col->identifier = -1;
		for(j = 0; j < 5;j++){
			col_target->i_emit_e[j] = logsum(col_target->m_emit_e[j] ,col_source->m_emit[j] ); /// org->hmms[0]->hmm_column[i]->i_emit_e[j];
			col_target->i_emit_e[j] = logsum(col_target->i_emit_e[j] ,col_source->i_emit[j] ); ///  org->hmms[0]->hmm_column[i]->m_emit_e[j];
		}
		col_target->short_transition_e[NEXT] = logsum(col_target->short_transition_e[NEXT] ,col_source->short_transition_e[NEXT] );
		col_target->short_transition_e[SELF] = logsum(col_target->short_transition_e[SELF] ,col_source->short_transition_e[SELF] );

		for(j = 0; j < target->num_hmms-1;j++){
			//col_target->long_transition[j] = org->hmms[0]->hmm_column[i]->long_transition[j];// prob2scaledprob(1.0 - max) +  prob2scaledprob(0.5  /  (float) (model->num_hmms-1));
			col_target->long_transition_e[j] = logsum(col_target->long_transition_e[j], col_source->long_transition_e[j]);// org->hmms[0]->hmm_column[i]->long_transition_e[j];
		}
		
		
	}
	
	/*
	 
	 sub models model....
	 
	 */
	len = target->hmms[1]->num_columns;
	for(i = 1; i< target->num_hmms;i++){
		for(j = 0; j < len;j++){
			col_target = target->hmms[i]->hmm_column[i];
			col_source = source->hmms[i]->hmm_column[i];
			//col = model->hmms[i]->hmm_column[j];
			//col->identifier = -1;
			for(c = 0; c < 5;c++){
				//col->i_emit[c] = org->hmms[i]->hmm_column[j]->i_emit[c];
				//col->m_emit[c] = org->hmms[i]->hmm_column[j]->m_emit[c];
				
				col_target->i_emit_e[c] = logsum(col_target->i_emit_e[c] , col_source->i_emit_e[c] ); //  org->hmms[i]->hmm_column[j]->i_emit_e[c];
				col_target->m_emit_e[c] =  logsum(col_target->m_emit_e[c] , col_source->m_emit_e[c] );/// org->hmms[i]->hmm_column[j]->m_emit_e[c];
			}
			
			
			
			
						
			col_target->short_transition_e[MM] = logsum(col_target->short_transition_e[MM] ,col_source->short_transition_e[MM] );//
			col_target->short_transition_e[MI] = logsum(col_target->short_transition_e[MI] ,col_source->short_transition_e[MI] );//
			col_target->short_transition_e[MD] = logsum(col_target->short_transition_e[MD] ,col_source->short_transition_e[MD] );//
			
			col_target->short_transition_e[II] = logsum(col_target->short_transition_e[II] ,col_source->short_transition_e[II] );//
			col_target->short_transition_e[IM] = logsum(col_target->short_transition_e[IM] ,col_source->short_transition_e[IM] );//
			
			col_target->short_transition_e[DD] = logsum(col_target->short_transition_e[DD] ,col_source->short_transition_e[DD] );//
			col_target->short_transition_e[DM] = logsum(col_target->short_transition_e[DM] ,col_source->short_transition_e[DM] );//
		
			
			
			
			for(c = 0;c <  (target->num_hmms-1);c++){

				col_target->long_transition_e[c] = logsum(col_target->long_transition_e[c],col_source->long_transition_e[c] ); //  org->hmms[i]->hmm_column[j]->long_transition_e[c];
			}
			
			
			
		}
		
	}
	return target;
}



struct model* init_model(struct model* model)
{
	
	int i,j,c;
	double max;
	int len = 0;
	
	int kmer_len = log(model->num_hmms-1) / log(4);
	int x;
	
	float start_prob = 1.0;
	
	struct hmm_column* col =0;
	assert(model != 0);
	len = model->hmms[0]->num_columns;
	c = 0;
	max = -1;
	for(i = 0; i< 1000;i++){
		if(binomial_distribution((double)i / 1000.0 , model->average_length  ,len)  > max){
			max = binomial_distribution((double)i / 1000.0 ,model->average_length ,len );
			c = i;
		}
	}
	max = (double)c / 1000.0;
	
	
	model->random_next =  prob2scaledprob( max);
	model->random_self = prob2scaledprob(1.0-max);
	
/*
 
 Init main model....
 
 */
	
	
	for(i = 0; i < len;i++){
		col = model->hmms[0]->hmm_column[i];
		
		col->identifier = -1;
		for(j = 0; j < 5;j++){
			col->i_emit[j] = model->background_nuc_frequency[j];
			col->m_emit[j] = model->background_nuc_frequency[j];
			col->i_emit_e[j] =  prob2scaledprob(0.0f);
			col->m_emit_e[j] =  prob2scaledprob(0.0f);
		}
		col->short_transition[NEXT] = prob2scaledprob( max);
		col->short_transition[SELF] = prob2scaledprob((1.0 - max) * 0.5);
		
		col->short_transition_e[NEXT] = prob2scaledprob( 0.0);
		col->short_transition_e[SELF] = prob2scaledprob(0.0);
		
		
		for(j = 0; j< model->num_hmms-1;j++){
			col->long_transition[j] = prob2scaledprob(((1.0 - max) * 0.5)  /  (float) (model->num_hmms-1));
			col->long_transition_e[j] = prob2scaledprob(0.0);
		}
	}
	
	/*
	 
	 sub models model....
	 
	 */
	len = model->hmms[1]->num_columns;
	for(i = 1; i< model->num_hmms;i++){
		for(j = 0; j < len;j++){
			col = model->hmms[i]->hmm_column[j];
			col->identifier = -1;
			for(c = 0; c < 5;c++){
				col->i_emit[c] = model->background_nuc_frequency[c];
				
				col->i_emit_e[c] =  prob2scaledprob(0.0f);
				col->m_emit_e[c] =  prob2scaledprob(0.0f);
			}
			
			/*
			
			 This horrible piece of code checks if we are at the start of the sub - model - if so it assigns a high probability to nucleotides of the #kmer given by i-1 -> the current sub-hmm (i.e. first hmm is AAAXXXXX, second AACXXXXX and so forth....
			 
			*/
			if(j < kmer_len){
				for(c = 0; c < 5;c++){
					col->m_emit[c] =prob2scaledprob(0.0001);
				}
				x = (kmer_len - 1 - j)  <<1  ;
				
				col->m_emit[  ((i-1) >> x) & 0x3] = prob2scaledprob(1.0 - (0.0001*4.0) );
			}else{
				for(c = 0; c < 5;c++){
					col->m_emit[c] = model->background_nuc_frequency[c];
				}
			}
			
			//for(c = 0; c < 5;c++){
			//	fprintf(stderr,"%d	%d	%f	%f	%f	%f\n",i,j, scaledprob2prob(col->m_emit[0]),  scaledprob2prob(col->m_emit[1]),  scaledprob2prob(col->m_emit[2]),  scaledprob2prob(col->m_emit[3]));
			//}
			
			
			if(j < kmer_len){
				start_prob = 1.0;
				for(c = 0;c <  (model->num_hmms-1);c++){
					col->long_transition[c] = prob2scaledprob(0.0);
					col->long_transition_e[c] = prob2scaledprob(0.0);
				}
				
				
				
			}else{
				start_prob = 0.9;
				for(c = 0;c <  (model->num_hmms-1);c++){
					col->long_transition[c] = prob2scaledprob(0.1);
					col->long_transition_e[c] = prob2scaledprob(0.0);
				}

			}
			
			col->short_transition[MM] = prob2scaledprob( start_prob * 0.99);
			col->short_transition[MI] = prob2scaledprob(1.0 - start_prob * 0.99) +  prob2scaledprob(0.5);
			col->short_transition[MD] = prob2scaledprob(1.0 - start_prob * 0.99) +  prob2scaledprob(0.5);
			
			col->short_transition[II] = prob2scaledprob(1.0 - 0.99);
			col->short_transition[IM] = prob2scaledprob(0.99);
			
			col->short_transition[DD] = prob2scaledprob(1.0 - 0.99);
			col->short_transition[DM] = prob2scaledprob(0.99);
			
			
			col->short_transition_e[MM] =  prob2scaledprob(0.0);
			col->short_transition_e[MI] =  prob2scaledprob(0.0);
			col->short_transition_e[MD] =  prob2scaledprob(0.0);
			
			col->short_transition_e[II] =  prob2scaledprob(0.0);
			col->short_transition_e[IM] =  prob2scaledprob(0.0);
			
			col->short_transition_e[DD] =  prob2scaledprob(0.0);
			col->short_transition_e[DM] =  prob2scaledprob(0.0);
		}
	}
	return model;
}



















