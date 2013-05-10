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
#include <pthread.h>

#include "barcode_hmm.h"


void hmm_controller(struct parameters* param,int (*fp)(struct read_info** ,struct parameters*,FILE* ),int file_num)
{
	struct read_info** ri = 0;
	
	FILE* outfile;
	int i,j,c;
	int numseq;
	int total_read = 0;
	float sum = 0.0;
	
	init_logsum();
	
	float* back = 0;
	int average_length = 0;
	back = malloc(sizeof(float)*5);
	for(i = 0; i < 5;i++){
		back[i]= 0.0f;//prob2scaledprob( 0.2);
	}
	
	param->num_query = 5000000;
	
	FILE* file = 0;
	
	ri = malloc(sizeof(struct read_info*) * param->num_query);
	
	assert(ri !=0);
	
	for(i = 0; i < param->num_query;i++){
		ri[i] = malloc(sizeof(struct read_info));
		ri[i]->seq = 0;
		ri[i]->name = 0;
		ri[i]->qual = 0;
		ri[i]->labels = 0;
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
		
		if(param->matchstart!= -1 || param->matchend !=-1){
			for(i = 0; i < numseq;i++){
				average_length += param->matchend - param->matchstart;
				for(j = param->matchstart;j < param->matchend;j++){
					back[(int)ri[i]->seq[j]] += 1.0f;
				}
			}
		}else{
			for(i = 0; i < numseq;i++){
				average_length += ri[i]->len;
				for(j = 0;j < ri[i]->len;j++){
					back[(int)ri[i]->seq[j]] += 1.0f;
				}
			}
		}
		total_read += numseq;
		if(total_read > 5000000){
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
	
	pclose(file);
	file =  io_handler(file, file_num,param);
	
	
	//rewind(file);
	
	struct model_bag* mb = init_model_bag(param, back);
	
	if(!param->train ){
	
	}else if( !strcmp( param->train , "full")){
		for(i = 0; i < 10;i++){
			fprintf(stderr,"Iteration %d\n",i);
			while ((numseq = fp(ri, param,file)) != 0){
				//	numseq = fp(ri, param,file);
				mb =  run_pHMM(mb,ri,param,numseq,MODE_TRAIN);
				
				
				//fprintf(stderr,"\n\n\n\n\n\n");
				
				
				
			}
			pclose(file);
			file =  io_handler(file, file_num,param);
			//rewind(file);
			for(j = 0; j < mb->num_models;j++){
				//print_model(mb->model[i]);
				mb->model[j] = reestimate(mb->model[j], 0);
				//	print_model(mb->model[i]);
			}
		}
		
	}else if (!strcmp(param->train, "half" )){
		for(i = 0; i < 10;i++){
			fprintf(stderr,"Iteration %d\n",i);
			while ((numseq = fp(ri, param,file)) != 0){
				//	numseq = fp(ri, param,file);
				mb =  run_pHMM(mb,ri,param,numseq,MODE_TRAIN);
				
				
				//fprintf(stderr,"\n\n\n\n\n\n");
				
				
				
			}
			pclose(file);
			file =  io_handler(file, file_num,param);
			//rewind(file);
			for(j = 0; j < mb->num_models;j++){
				//print_model(mb->model[i]);
				mb->model[j] = reestimate(mb->model[j], 2);
				//	print_model(mb->model[i]);
			}
		}
	}
	
	pclose(file);
	file =  io_handler(file, file_num,param);
	
	if(param->outfile){
		if ((outfile = fopen( param->outfile, "w")) == NULL){
			fprintf(stderr,"can't open output\n");
			exit(-1);
		}
	}else{
		outfile= stdout;
	}
	
	struct log_information* li = 0;
	
	li = malloc(sizeof(struct log_information));
	
	li->total_read = 0;
	li->success = 0;
	li->len_failure = 0;
	li->prob_failure = 0;
	li->arch_failure = 0;
	for(i = 0;i < 1001;i++){
		li->probability_distribution[i] = 0;
	}
	//int c1,c2,c3,key,bar,mem,fingerlen,required_finger_len,ret;
	//char alpha[5] = "ACGTN";
	
	total_read = 0;
	while ((numseq = fp(ri, param,file)) != 0){
		//	numseq = fp(ri, param,file);
		mb =  run_pHMM(mb,ri,param,numseq,MODE_GET_LABEL);
		li->total_read += numseq;
		
		for(i = 0; i < numseq;i++){
			ri[i]->prob = expf( ri[i]->prob) / (1.0f + expf(ri[i]->prob ));
			li->probability_distribution[  (int) floor(ri[i]->prob* 1000.0     + 0.5)]+= 1;
			c = print_trimmed_sequence(mb, param,  ri[i],outfile);
			
			switch (c) {
				case 1:
					li->success++;
					break;
				case -1:
					li->prob_failure++;
					break;
				case -2:
					li->len_failure++;
					break;
				case -3:
					li->arch_failure++;
					break;
				default:
					break;
			}
			
			/*key = 0;
			bar = -1;
			mem = -1;
			//if( expf( ri[i]->prob) / (1.0 + expf(ri[i]->prob )) < 0.5){
			fprintf(stderr,"%f	%f\n",  expf( ri[i]->prob) / (1.0 + expf(ri[i]->prob )) ,ri[i]->prob);
			
			
			fprintf(stderr,"%s\n", ri[i]->name);
			for(j = 0; j < ri[i]->len;j++){
				c1 = mb->label[(int)ri[i]->labels[j+1]];
				c2 = c1 & 0xFFFF;
				c3 = (c1 >> 16) & 0x7FFF;
				fprintf(stderr,"%c",  alpha[(int)ri[i]->seq[j]] );
			}
			fprintf(stderr,"\n");
			for(j = 0; j < ri[i]->len;j++){
				c1 = mb->label[(int)ri[i]->labels[j+1]];
				c2 = c1 & 0xFFFF;
				c3 = (c1 >> 16) & 0x7FFF;
				fprintf(stderr,"%c",   param->read_structure->type[c2] );
				if(param->read_structure->type[c2] == 'F'){
					key = (key << 2 )|  (ri[i]->seq[j] & 0x3);
				}
				
			}
			fprintf(stderr,"	key = %d\n",key);
			
			for(j = 0; j < ri[i]->len;j++){
				c1 = mb->label[(int)ri[i]->labels[j+1]];
				c2 = c1 & 0xFFFF;
				c3 = (c1 >> 16) & 0x7FFF;
				fprintf(stderr,"%d", c3   );
				if(param->read_structure->type[c2] == 'B'){
					bar = c3;
					mem = c2;
				}
			}
				fprintf(stderr,"\n");
			//	bar = 1;
			//fprintf(stderr,"	bar= %d (%s)\n",bar,   param->read_structure->sequence_matrix[mem][bar] );
			//}
			if(i == 1000){
				exit(0);
			}*/
		}
		//fprintf(stderr,"%0.1f%% extracted (%d read so far...)\n",  (float) success / (float) total_read  *100.0f,total_read);
		//fprintf(stderr,"%d	successfully extracted\n" , success);
		//fprintf(stderr,"%d	low probability\n" , prob_failure);
		//fprintf(stderr,"%d	too short\n" , len_failure);
		//fprintf(stderr,"%d	problems with architecture\n" , arch_failure);
	}
	
	//for(i = 0; i < mb->num_models;i++){
	//	print_model(mb->model[i]);
	
	//}
	
	fprintf(stderr,"%d\n", li->total_read);
	fprintf(stderr,"%d	successfully extracted\n" ,li->success);
	fprintf(stderr,"%d	low probability\n" , li->prob_failure);
	fprintf(stderr,"%d	too short\n" , li->len_failure);
	fprintf(stderr,"%d	problems with architecture\n" , li->arch_failure);
	
	fprintf(stderr,"%0.1f%% extracted\n",  (float) li->success / (float) li->total_read  *100.0f);
	
	for(i = 0;i < 1001;i++){
		fprintf(stderr,"%f	%d\n",(float) i / 1000.0f,li->probability_distribution[i]);
	}
	
	free_model_bag(mb);
	free(li);
	for(i = 0; i < param->num_query;i++){
		free(ri[i]->strand);
		free(ri[i]->hits);
		
		if(ri[i]->cigar){
			free(ri[i]->cigar);
		}
		if(ri[i]->md){
			free(ri[i]->md);
		}
		if(ri[i]->name){
			free(ri[i]->name);
		}
		if(ri[i]->seq){
			free(ri[i]->seq);
		}
		if(ri[i]->qual){
			free(ri[i]->qual );
		}
		
		free(ri[i]);
	}
	free(back);
	free(ri);
	if(param->sam == 2 || param->sam == 1 || param->gzipped ){
		pclose(file);
	}else{
		fclose(file);
	}
	if(param->outfile){
		fclose(outfile);
	}
	
}


struct model_bag* run_pHMM(struct model_bag* mb,struct read_info** ri,struct parameters* param,int numseq, int mode)
{
	struct thread_data* thread_data = 0;
	thread_data = malloc(sizeof(struct thread_data)* param->num_threads);
	pthread_t threads[param->num_threads];
	pthread_attr_t attr;
	int i,t;
	int interval = 0;
	int rc;
	
	interval =  (int)((double)numseq /(double)param->num_threads);
	
	for(t = 0;t < param->num_threads ;t++) {
		thread_data[t].ri = ri;
		thread_data[t].mb = copy_model_bag(mb);
		thread_data[t].start = t*interval;
		thread_data[t].end = t*interval + interval;
		thread_data[t].param = param;
	//m	fprintf(stderr,"%d %d %d %d\n",t,thread_data[t].start,thread_data[t].end,numseq );
	}
	thread_data[param->num_threads-1].end = numseq;
	
	
	
	//for(i = 0; i < mb->num_models;i++){
	//	print_model(thread_data[0].mb->model[i]);
	//}
	
	
	//fprintf(stderr,"got here...\n");
	//exit(0);
	rc = pthread_attr_init(&attr);
	if(rc){
		fprintf(stderr,"ERROR; return code from pthread_attr_init() is %d\n", rc);
		exit(-1);
	}
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	
	for(t = 0;t < param->num_threads;t++) {
		switch (mode) {
			case MODE_GET_LABEL:
				rc = pthread_create(&threads[t], &attr, do_label_thread, (void *) &thread_data[t]);
				break;
			case MODE_TRAIN:
				rc = pthread_create(&threads[t], &attr, do_baum_welch_thread, (void *) &thread_data[t]);
				break;
		}
		
		
		
		if (rc) {
			fprintf(stderr,"ERROR; return code from pthread_create() is %d\n", rc);
			exit(-1);
		}
	}
	
	pthread_attr_destroy(&attr);
	
	for (t = 0;t < param->num_threads;t++){
		rc = pthread_join(threads[t], NULL);
		if (rc){
			fprintf(stderr,"ERROR; return code from pthread_join()is %d\n", rc);
			exit(-1);
		}
	}
	
	for (t = 0;t < param->num_threads;t++){
		for(i = 0; i < mb->num_models;i++){

			mb->model[i] = copy_estimated_parameter(mb->model[i], thread_data[t].mb->model[i]);
		}
	}
	//fprintf(stderr,"got here...\n");
	
	for(t = 0;t < param->num_threads;t++) {
		free_model_bag(thread_data[t].mb);
	}
	
	free(thread_data);
	return mb;
}


void* do_label_thread(void *threadarg)
{
	struct thread_data *data;
	data = (struct thread_data *) threadarg;
	
	struct read_info** ri  = data->ri;
	struct model_bag* mb = data->mb;
	
	int matchstart = data->param->matchstart;
	int matchend = data->param->matchend;
	
	int start = data->start;
	int end = data->end;
	int i;
	int tmp = 0;
	
	if(matchstart != -1 || matchend != -1){
		for(i = start; i < end;i++){
			tmp = matchend - matchstart ;
			
			
			
			mb = backward(mb, ri[i]->seq + matchstart , tmp);
			//mb = forward(mb,  ri[i]->seq ,ri[i]->len);
			mb = forward_max_posterior_decoding(mb, ri[i] , ri[i]->seq+matchstart ,tmp );
			
		}
	}else{
		for(i = start; i < end;i++){
			mb = backward(mb, ri[i]->seq ,ri[i]->len);
			//mb = forward(mb,  ri[i]->seq ,ri[i]->len);
			mb = forward_max_posterior_decoding(mb, ri[i], ri[i]->seq ,ri[i]->len);
		
		}
	}
	pthread_exit((void *) 0);
}


void* do_baum_welch_thread(void *threadarg)
{
	struct thread_data *data;
	data = (struct thread_data *) threadarg;
	
	struct read_info** ri  = data->ri;
	struct model_bag* mb = data->mb;
	
	int start = data->start;
	int end = data->end;
	int i;
	
	for(i = start; i < end;i++){
		mb = backward(mb, ri[i]->seq ,ri[i]->len);
		mb = forward_extract_posteriors(mb, ri[i]->seq ,ri[i]->len);
	}
	pthread_exit((void *) 0);
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
				
				//c_hmm_column->M_backward[i] = psilent[i+1] + mb->model[j]->M_to_silent[f] ;
				
				c_hmm_column->M_backward[i] = psilent[i+1] + c_hmm_column->transition[MSKIP];// mb->model[j]->M_to_silent[f] ;
				
				
				//fprintf(stderr," Mback at modellen:%d %f %f\n",i, c_hmm_column->M_backward[i] ,c_hmm_column->transition[MSKIP]);
				
				//c_hmm_column->I_backward[i] =  psilent[i+1]+ mb->model[j]->I_to_silent[f] ;
				
				c_hmm_column->I_backward[i] =  psilent[i+1] + c_hmm_column->transition[ISKIP];//  mb->model[j]->I_to_silent[f] ;
				
				c_hmm_column->I_backward[i] = logsum(c_hmm_column->I_backward[i] , c_hmm_column->M_backward[i+1] + c_hmm_column->transition[IM] + c_hmm_column->m_emit[c]);
				
				c_hmm_column->I_backward[i] = logsum(c_hmm_column->I_backward[i] , c_hmm_column->I_backward[i+1] + c_hmm_column->transition[II] + c_hmm_column->i_emit[c]);
				
				
				
				c_hmm_column->D_backward[i] = prob2scaledprob(0.0f);
				for(g = model_len-1;g >= 0;g--){
					c_hmm_column = hmm->hmm_column[g];
					p_hmm_column = hmm->hmm_column[g+1];
					
					c_hmm_column->M_backward[i]  = p_hmm_column->M_backward[i+1] + p_hmm_column->m_emit[c] + c_hmm_column->transition[MM];
					
					
					c_hmm_column->M_backward[i] = logsum (c_hmm_column->M_backward[i],psilent[i+1] + c_hmm_column->transition[MSKIP]);
					
					//insert - emit previous symbol etc. etc.
					c_hmm_column->M_backward[i] = logsum(c_hmm_column->M_backward[i] , c_hmm_column->I_backward[i+1] +c_hmm_column->i_emit[c ]  + c_hmm_column->transition[MI]);
					
					//delete - neex to go to previous columns
					
					c_hmm_column->M_backward[i] = logsum(c_hmm_column->M_backward[i],p_hmm_column->D_backward[i] + c_hmm_column->transition[MD]);
					
					// insert state..
					// from previous insertion....
					c_hmm_column->I_backward[i] = c_hmm_column->I_backward[i+1] + c_hmm_column->transition[II] + c_hmm_column->i_emit[c];
					
					c_hmm_column->I_backward[i] = logsum(c_hmm_column->I_backward[i], psilent[i+1] + c_hmm_column->transition[ISKIP]);
					//from previous match state....
					
					c_hmm_column->I_backward[i] = logsum( c_hmm_column->I_backward[i],p_hmm_column->M_backward[i+1] + c_hmm_column->transition[IM] + p_hmm_column->m_emit[c]);
					///GRERRRRRRR 
					
					//delete state
					
					//from previous delection
					c_hmm_column->D_backward[i] = p_hmm_column->D_backward[i] + c_hmm_column->transition[DD];
					
					//from previous match (i.e. gap close
					
					c_hmm_column->D_backward[i] = logsum(c_hmm_column->D_backward[i], p_hmm_column->M_backward[i] + p_hmm_column->m_emit[(int) seqa[i]] + c_hmm_column->transition[DM]);
					
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
	//fprintf(stderr,"SCore:%f	%f\n", mb->b_score , scaledprob2prob(mb->b_score) );
	
	//fprintf(stderr," BACKWARD:::::::::::\n");
	
	/*for(j = 0; j < mb->num_models;j++){
		for(i = 0; i <= len;i++){
			fprintf(stderr,"%d	%d	%f\n",j,i,mb->model[j]->silent[i]  );
		}
	}
	exit(0);*/
	/*
	for(j = 0; j < mb->num_models;j++){
		
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
	//exit(0);
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
				c_hmm_column->I_foward[i] = logsum(c_hmm_column->I_foward[i], c_hmm_column->I_foward[i-1] + c_hmm_column->transition[II]);
				
				c_hmm_column->I_foward[i]  = logsum ( c_hmm_column->I_foward[i] ,c_hmm_column->M_foward[i-1] + c_hmm_column->transition[MI]);
				
				c_hmm_column->I_foward[i] = c_hmm_column->I_foward[i] + c_hmm_column->i_emit[c];
				
				c_hmm_column->D_foward[i] = prob2scaledprob(0.0f);
				
				csilent[i] =  logsum(csilent[i],  c_hmm_column->M_foward[i] + c_hmm_column->transition[MSKIP]);// mb->model[j]->M_to_silent[f]);
				csilent[i] =  logsum(csilent[i],  c_hmm_column->I_foward[i] + c_hmm_column->transition[ISKIP]);
				
				for(g = 1;g < hmm->num_columns;g++){
					c_hmm_column = hmm->hmm_column[g];
					p_hmm_column = hmm->hmm_column[g-1];
					
					//Match state
					//transition from previous match state
					c_hmm_column->M_foward[i] = p_hmm_column->M_foward[i-1] + p_hmm_column->transition[MM];
					//transition from previous insert state
					c_hmm_column->M_foward[i] = logsum(c_hmm_column->M_foward[i] , p_hmm_column->I_foward[i-1] + p_hmm_column->transition[IM] );
					//transition from previous delete state
					c_hmm_column->M_foward[i] = logsum(c_hmm_column->M_foward[i], p_hmm_column->D_foward[i] + p_hmm_column->transition[DM]);
					
					// emission promability in curent M state ;
					c_hmm_column->M_foward[i]  = c_hmm_column->M_foward[i]  + c_hmm_column->m_emit[c];
					
					
					// Instertion State ..
					//self loop insertion to insertion
					c_hmm_column->I_foward[i] = c_hmm_column->I_foward[i-1] + c_hmm_column->transition[II];
					// start new insertion
					c_hmm_column->I_foward[i]  = logsum ( c_hmm_column->I_foward[i] ,c_hmm_column->M_foward[i-1] + c_hmm_column->transition[MI]);
					
					//instertion emission...
					c_hmm_column->I_foward[i]  = c_hmm_column->I_foward[i]  + c_hmm_column->i_emit[c];
					
					
					// deletion state
					//from previous match state.
					c_hmm_column->D_foward[i] = p_hmm_column->M_foward[i] + p_hmm_column->transition[MD];
					//from previous delete state
					
					c_hmm_column->D_foward[i] = logsum(c_hmm_column->D_foward[i], p_hmm_column->D_foward[i] + p_hmm_column->transition[DD] );
					
					csilent[i] =  logsum(csilent[i],  c_hmm_column->M_foward[i] + c_hmm_column->transition[MSKIP]);// mb->model[j]->M_to_silent[f]);
					csilent[i] =  logsum(csilent[i],  c_hmm_column->I_foward[i] + c_hmm_column->transition[ISKIP]);
					
				}
				//fprintf(stderr,"%d	%f	%f	%f	%f	%f\n",i,  c_hmm_column->M_foward[i] , mb->model[j]->M_to_silent[f],  c_hmm_column->I_foward[i] , mb->model[j]->I_to_silent[f], psilent[i]);
				//csilent[i] =  logsum(csilent[i],  c_hmm_column->M_foward[i] + mb->model[j]->M_to_silent[f]);
				//csilent[i] =  logsum(csilent[i],  c_hmm_column->I_foward[i] + mb->model[j]->I_to_silent[f]);
				csilent[i] = logsum(csilent[i], psilent[i] + mb->model[j]->skip);
				
			}
			
		}
	}
	
		
	mb->f_score = mb->model[mb->num_models-1]->silent_forward[len];
	
	//fprintf(stderr,"SCORE:%f	%f\n", mb->f_score, scaledprob2prob(mb->f_score));
	
	/*
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
	}*/
	//exit(0);
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
				c_hmm_column->I_foward[i] = logsum(c_hmm_column->I_foward[i], c_hmm_column->I_foward[i-1] + c_hmm_column->transition[II]);
				
				c_hmm_column->I_foward[i]  = logsum ( c_hmm_column->I_foward[i] ,c_hmm_column->M_foward[i-1] + c_hmm_column->transition[MI]);
				
				c_hmm_column->I_foward[i] = c_hmm_column->I_foward[i] + c_hmm_column->i_emit[c];
				
				//***************post
				mb->model[j]->silent_to_I_e[f]  = logsum(mb->model[j]->silent_to_I_e[f] , psilent[i-1] + mb->model[j]->silent_to_I[f]  + c_hmm_column->i_emit[c] +  c_hmm_column->I_backward[i] -mb->b_score);
				
				c_hmm_column->transition_e[II] = logsum(c_hmm_column->transition_e[II],c_hmm_column->I_foward[i-1] + c_hmm_column->transition[II] + c_hmm_column->i_emit[c] +  c_hmm_column->I_backward[i]-mb->b_score);
				
				c_hmm_column->transition_e[MI] = logsum(c_hmm_column->transition_e[MI] , c_hmm_column->M_foward[i-1] + c_hmm_column->transition[MI] +  c_hmm_column->i_emit[c] +  c_hmm_column->I_backward[i]-mb->b_score);
				
				
				c_hmm_column->i_emit_e[c] = logsum( c_hmm_column->i_emit_e[c] , c_hmm_column->I_foward[i] + c_hmm_column->I_backward[i] -mb->b_score);
				
				//***************post
				
				
				c_hmm_column->D_foward[i] = prob2scaledprob(0.0f);
				
				// no post???
				
				//
				
				csilent[i] =  logsum(csilent[i],  c_hmm_column->M_foward[i] +c_hmm_column->transition[MSKIP]);
				csilent[i] =  logsum(csilent[i],  c_hmm_column->I_foward[i] +c_hmm_column->transition[ISKIP]);
				
				//***************post
				c_hmm_column->transition_e[MSKIP] = logsum(c_hmm_column->transition_e[MSKIP], c_hmm_column->M_foward[i] +c_hmm_column->transition[MSKIP] +  bsilent[i+1] -mb->b_score);
				
				c_hmm_column->transition_e[ISKIP] = logsum(c_hmm_column->transition_e[ISKIP], c_hmm_column->I_foward[i] +c_hmm_column->transition[ISKIP] +  bsilent[i+1] -mb->b_score);
				
				
				
				//***************post
				
				
				for(g = 1;g < hmm->num_columns;g++){
					c_hmm_column = hmm->hmm_column[g];
					p_hmm_column = hmm->hmm_column[g-1];
					
					//Match state
					//transition from previous match state
					c_hmm_column->M_foward[i] = p_hmm_column->M_foward[i-1] + p_hmm_column->transition[MM];
					//transition from previous insert state
					c_hmm_column->M_foward[i] = logsum(c_hmm_column->M_foward[i] , p_hmm_column->I_foward[i-1] + p_hmm_column->transition[IM] );
					//transition from previous delete state
					c_hmm_column->M_foward[i] = logsum(c_hmm_column->M_foward[i], p_hmm_column->D_foward[i] + p_hmm_column->transition[DM]);
					
					// emission promability in curent M state ;
					c_hmm_column->M_foward[i]  = c_hmm_column->M_foward[i]  + c_hmm_column->m_emit[c];
					
					
					//***************post
					p_hmm_column->transition_e[MM] = logsum(p_hmm_column->transition_e[MM] , p_hmm_column->M_foward[i-1] + p_hmm_column->transition[MM] +  c_hmm_column->m_emit[c] +  c_hmm_column->M_backward[i] -mb->b_score );
					
					p_hmm_column->transition_e[IM] = logsum(p_hmm_column->transition_e[IM],p_hmm_column->I_foward[i-1] + p_hmm_column->transition[IM] +  c_hmm_column->m_emit[c] +  c_hmm_column->M_backward[i] -mb->b_score);
					
					p_hmm_column->transition_e[DM] = logsum(p_hmm_column->transition_e[DM],p_hmm_column->D_foward[i] + p_hmm_column->transition[DM] + c_hmm_column->m_emit[c] +  c_hmm_column->M_backward[i] -mb->b_score);
					
					c_hmm_column->m_emit_e[c] = logsum(c_hmm_column->m_emit_e[c],  c_hmm_column->M_foward[i] + c_hmm_column->M_backward[i] -mb->b_score );
					//***************post
					
					
					
					
					
					// Instertion State ..
					//self loop insertion to insertion
					c_hmm_column->I_foward[i] = c_hmm_column->I_foward[i-1] + c_hmm_column->transition[II];
					// start new insertion
					c_hmm_column->I_foward[i]  = logsum ( c_hmm_column->I_foward[i] ,c_hmm_column->M_foward[i-1] + c_hmm_column->transition[MI]);
					
					//instertion emission...
					c_hmm_column->I_foward[i]  = c_hmm_column->I_foward[i]  + c_hmm_column->i_emit[c];
					
					//***************post
					c_hmm_column->transition_e[II] = logsum(c_hmm_column->transition_e[II] ,  c_hmm_column->I_foward[i-1] + c_hmm_column->transition[II] + c_hmm_column->i_emit[c] + c_hmm_column->I_backward[i] - mb->b_score);
					
					c_hmm_column->transition_e[MI] = logsum(c_hmm_column->transition_e[MI] , c_hmm_column->M_foward[i-1] + c_hmm_column->transition[MI] + c_hmm_column->i_emit[c] + c_hmm_column->I_backward[i] - mb->b_score);
					
					c_hmm_column->i_emit_e[c] = logsum(c_hmm_column->i_emit_e[c] , c_hmm_column->I_foward[i]   + c_hmm_column->I_backward[i] - mb->b_score);
					//***************post

					
					
					// deletion state
					//from previous match state.
					c_hmm_column->D_foward[i] = p_hmm_column->M_foward[i] + p_hmm_column->transition[MD];
					//from previous delete state
					
					c_hmm_column->D_foward[i] = logsum(c_hmm_column->D_foward[i], p_hmm_column->D_foward[i] + p_hmm_column->transition[DD] );
					
					//***************post
					p_hmm_column->transition_e[MD] = logsum(p_hmm_column->transition_e[MD],  p_hmm_column->M_foward[i] + p_hmm_column->transition[MD] + c_hmm_column->D_backward[i] - mb->b_score);
					
					p_hmm_column->transition_e[DD] = logsum(p_hmm_column->transition_e[DD] ,p_hmm_column->D_foward[i] + p_hmm_column->transition[DD] + c_hmm_column->D_backward[i] - mb->b_score);
					//***************post
					csilent[i] =  logsum(csilent[i],  c_hmm_column->M_foward[i] +c_hmm_column->transition[MSKIP]);
					csilent[i] =  logsum(csilent[i],  c_hmm_column->I_foward[i] +c_hmm_column->transition[ISKIP]);
					
					//***************post
					c_hmm_column->transition_e[MSKIP] = logsum(c_hmm_column->transition_e[MSKIP], c_hmm_column->M_foward[i] +c_hmm_column->transition[MSKIP] +  bsilent[i+1] -mb->b_score);
					
					c_hmm_column->transition_e[ISKIP] = logsum(c_hmm_column->transition_e[ISKIP], c_hmm_column->I_foward[i] +c_hmm_column->transition[ISKIP] +  bsilent[i+1] -mb->b_score);
					
					
					
					//***************post
					
					
					
				}
				//fprintf(stderr,"%d	%f	%f	%f	%f	%f\n",i,  c_hmm_column->M_foward[i] , mb->model[j]->M_to_silent[f],  c_hmm_column->I_foward[i] , mb->model[j]->I_to_silent[f], psilent[i]);
				//csilent[i] =  logsum(csilent[i],  c_hmm_column->M_foward[i] + mb->model[j]->M_to_silent[f]);
				
				
				
				//csilent[i] =  logsum(csilent[i],  c_hmm_column->I_foward[i] + mb->model[j]->I_to_silent[f]);
				csilent[i] = logsum(csilent[i], psilent[i] + mb->model[j]->skip);
				
				//***************post
				//mb->model[j]->M_to_silent_e[f] = logsum(mb->model[j]->M_to_silent_e[f],  c_hmm_column->M_foward[i] + mb->model[j]->M_to_silent[f] + bsilent[i+1] -mb->b_score);
				
				//fprintf(stderr,"ADDED TO M->S model: %d		%f %f %f %f %f\n", j , scaledprob2prob( mb->model[j]->M_to_silent_e[f]) ,scaledprob2prob(c_hmm_column->M_foward[i]) , scaledprob2prob(mb->model[j]->M_to_silent[f] ),scaledprob2prob( bsilent[i+1]) ,scaledprob2prob( mb->b_score));
				
				//mb->model[j]->I_to_silent_e[f] = logsum(mb->model[j]->I_to_silent_e[f] , c_hmm_column->I_foward[i] + mb->model[j]->I_to_silent[f] + bsilent[i+1] -mb->b_score);
				
				mb->model[j]->skip_e =logsum(mb->model[j]->skip_e , psilent[i-1] + mb->model[j]->skip + bsilent[i] -mb->b_score);
 				
				
				
				//***************post
				
			}
			
		}
	}
	
	
	mb->f_score = mb->model[mb->num_models-1]->silent_forward[len];
	
	//fprintf(stderr,"SCORE:%f	%f\n", mb->f_score, scaledprob2prob(mb->f_score));
	return mb;
}



struct model_bag* forward_max_posterior_decoding(struct model_bag* mb, struct read_info* ri, char* a, int len)
{
	
	//char* a = ri->seq;
	//int len = ri->len;
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
				}
			}
		}
		for(i = 0; i <= len+1;i++){
			mb->model[j]->silent_forward[i] = prob2scaledprob(0.0f);
		}
	}
	
	mb->model[0]->silent_forward[0] = prob2scaledprob(1.0) + mb->model[0]->skip;
	
	for(j = 1; j < mb->num_models;j++){
		mb->model[j]->silent_forward[0] = mb->model[j-1]->silent_forward[0]  + mb->model[j]->skip ;
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
			
			//fprintf(stderr," %d %d %d\n", j , f, hmm_counter);
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
				
				
 				c_hmm_column->I_foward[i] = psilent[i-1] + mb->model[j]->silent_to_I[f] ;
				
				
				
				
				//add transitions to first columns////
				c_hmm_column->I_foward[i] = logsum(c_hmm_column->I_foward[i], c_hmm_column->I_foward[i-1] + c_hmm_column->transition[II]);
				
				c_hmm_column->I_foward[i] = logsum ( c_hmm_column->I_foward[i] ,c_hmm_column->M_foward[i-1] + c_hmm_column->transition[MI]);
				
				c_hmm_column->I_foward[i] = c_hmm_column->I_foward[i] + c_hmm_column->i_emit[c];
				
				
				
				
				//***************post
				
				mb->dyn_prog_matrix[i][hmm_counter] = logsum(mb->dyn_prog_matrix[i][hmm_counter],  c_hmm_column->I_foward[i] + c_hmm_column->I_backward[i] -mb->b_score );
				
				
				//***************post
				
				
				c_hmm_column->D_foward[i] = prob2scaledprob(0.0f);
				
				// no post???
				
				//
				
				csilent[i] =  logsum(csilent[i],  c_hmm_column->M_foward[i] + c_hmm_column->transition[MSKIP]);
				csilent[i] =  logsum(csilent[i],  c_hmm_column->I_foward[i] + c_hmm_column->transition[ISKIP]);

				
				for(g = 1;g < hmm->num_columns;g++){
					c_hmm_column = hmm->hmm_column[g];
					p_hmm_column = hmm->hmm_column[g-1];
					
					//Match state
					//transition from previous match state
					c_hmm_column->M_foward[i] = p_hmm_column->M_foward[i-1] + p_hmm_column->transition[MM];
					//transition from previous insert state
					c_hmm_column->M_foward[i] = logsum(c_hmm_column->M_foward[i] , p_hmm_column->I_foward[i-1] + p_hmm_column->transition[IM] );
					//transition from previous delete state
					c_hmm_column->M_foward[i] = logsum(c_hmm_column->M_foward[i], p_hmm_column->D_foward[i] + p_hmm_column->transition[DM]);
					
					// emission promability in curent M state ;
					c_hmm_column->M_foward[i]  = c_hmm_column->M_foward[i]  + c_hmm_column->m_emit[c];
					
					
					//***************post
					mb->dyn_prog_matrix[i][hmm_counter] = logsum(mb->dyn_prog_matrix[i][hmm_counter],  c_hmm_column->M_foward[i]  + c_hmm_column->M_backward[i] -mb->b_score );
					//***************post
					
					
					
					
					
					// Instertion State ..
					//self loop insertion to insertion
					c_hmm_column->I_foward[i] = c_hmm_column->I_foward[i-1] + c_hmm_column->transition[II];
					// start new insertion
					c_hmm_column->I_foward[i]  = logsum ( c_hmm_column->I_foward[i] ,c_hmm_column->M_foward[i-1] + c_hmm_column->transition[MI]);
					
					//instertion emission...
					c_hmm_column->I_foward[i]  = c_hmm_column->I_foward[i]  + c_hmm_column->i_emit[c];
					
					//***************post
					mb->dyn_prog_matrix[i][hmm_counter] = logsum(mb->dyn_prog_matrix[i][hmm_counter],  c_hmm_column->I_foward[i] + c_hmm_column->I_backward[i] -mb->b_score );
					//***************post
					
					
					
					// deletion state
					//from previous match state.
					c_hmm_column->D_foward[i] = p_hmm_column->M_foward[i] + p_hmm_column->transition[MD];
					//from previous delete state
					
					c_hmm_column->D_foward[i] = logsum(c_hmm_column->D_foward[i], p_hmm_column->D_foward[i] + p_hmm_column->transition[DD] );
					
					//***************post
					p_hmm_column->transition_e[MD] = logsum(p_hmm_column->transition_e[MD],  p_hmm_column->M_foward[i] + p_hmm_column->transition[MD] + c_hmm_column->D_backward[i] - mb->b_score);
					
					p_hmm_column->transition_e[DD] = logsum(p_hmm_column->transition_e[DD] ,p_hmm_column->D_foward[i] + p_hmm_column->transition[DD] + c_hmm_column->D_backward[i] - mb->b_score);
					//***************post
					
					csilent[i] =  logsum(csilent[i],  c_hmm_column->M_foward[i] +c_hmm_column->transition[MSKIP]);
					csilent[i] =  logsum(csilent[i],  c_hmm_column->I_foward[i] +c_hmm_column->transition[ISKIP]);

					
				}
				//fprintf(stderr,"%d	%f	%f	%f	%f	%f\n",i,  c_hmm_column->M_foward[i] , mb->model[j]->M_to_silent[f],  c_hmm_column->I_foward[i] , mb->model[j]->I_to_silent[f], psilent[i]);
				//csilent[i] =  logsum(csilent[i],  c_hmm_column->M_foward[i] + mb->model[j]->M_to_silent[f]);
				
				
				
				//csilent[i] =  logsum(csilent[i],  c_hmm_column->I_foward[i] + mb->model[j]->I_to_silent[f]);
				csilent[i] = logsum(csilent[i], psilent[i] + mb->model[j]->skip);
			}
			hmm_counter++;
			
		}
		//hmm_counter++;
	}
	
	
	mb->f_score = mb->model[mb->num_models-1]->silent_forward[len];
	
	for(i = 0; i <= len;i++){
	//	fprintf(stderr,"%d ",i);
		for(j = 0; j < mb->total_hmm_num;j++){
	//		fprintf(stderr,"%0.3f ",scaledprob2prob( mb->dyn_prog_matrix[i][j]));
			mb->dyn_prog_matrix[i][j] = scaledprob2prob( mb->dyn_prog_matrix[i][j]);
		}
	//	fprintf(stderr,"\n");
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
				if(tmp == max && c == j){
					move = c;
					max = tmp;
				}
				
			//	fprintf(stderr,"%0.3f ",scaledprob2prob( mb->dyn_prog_matrix[i][j]));
			}
			
			mb->dyn_prog_matrix[i][j]+= max;
			mb->path[i][j] = move;
		}
	}
	/*fprintf(stderr,"MATRIX:\n");
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
	*/
	//char path[100];
	
	i = len;
	max = -1;
	for(j = 0; j < mb->total_hmm_num;j++){
		if(mb->dyn_prog_matrix[i][j] > max){
			max = mb->dyn_prog_matrix[i][j];
			move = j;
		}
	}
	
	for(i = 0; i <= len;i++){
		ri->labels[i] = 0;
	}
	
	//path[len] = move;
	ri->labels[len] = move;
	
	for(i = len ;i > 0;i--){
		move = mb->path[i][move];
	//	path[i-1] = move;
		ri->labels[i-1] = move;
	}

	next_silent[0] = prob2scaledprob(1.0);
	
	for(i = 1; i <= len;i++){
		c = seqa[i];
		next_silent[0] = next_silent[0] + mb->model[0]->background_nuc_frequency[c] + prob2scaledprob(1.0 - (1.0 / (float)len));
	}
	next_silent[0] += prob2scaledprob(1.0 / (float)len);
	
	
	ri->prob  = mb->f_score - next_silent[0];
  	
	//fprintf(stderr,"F:%f B:%f\n", mb->f_score,mb->b_score );
	//fprintf(stderr,"SCORE:%f	%f	%f\n", mb->f_score,next_silent[0], scaledprob2prob(next_silent[0]) );
	///exit(0);
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


struct model* malloc_model_according_to_read_structure(int num_hmm, int length)
{
	struct model* model = NULL;
	int i = 0;
	int j = 0;
	int len = 0;

	model = malloc(sizeof(struct model));
	assert(model != 0);
	
	model->num_hmms = num_hmm;// (rs->numseq_in_segment[key]);
	model->hmms = malloc(sizeof(struct hmm*) * model->num_hmms  );//(rs->numseq_in_segment[key]));
	assert(model->hmms !=0);
	
	for(i = 0; i < model->num_hmms;i++){
		model->hmms[i] = malloc(sizeof(struct hmm) );
		assert(model->hmms[i]  != 0);
	}
	//model->M_to_silent = malloc(sizeof(float) * model->num_hmms);
	model->silent_to_M = malloc(sizeof(float) * model->num_hmms);
	
	//model->I_to_silent = malloc(sizeof(float) * model->num_hmms);
	model->silent_to_I = malloc(sizeof(float) * model->num_hmms);
	
	
	//model->M_to_silent_e = malloc(sizeof(float) * model->num_hmms);
	model->silent_to_M_e = malloc(sizeof(float) * model->num_hmms);
	
	//model->I_to_silent_e = malloc(sizeof(float) * model->num_hmms);
	model->silent_to_I_e = malloc(sizeof(float) * model->num_hmms);
	
	for(i = 0 ;i  < model->num_hmms;i++){
	//	model->M_to_silent[i] = 0.0f;
		model->silent_to_M[i] = 0.0f;
	//	model->I_to_silent[i] = 0.0f;
		model->silent_to_I[i] = 0.0f;
	}
	
	len = length;// (int)strlen(rs->sequence_matrix[key][0]);
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
				col->transition[MM] = prob2scaledprob(0.0f );// 1.0 - base_error * indel_freq);
				col->transition[MI] = prob2scaledprob(0.0f );//(base_error * indel_freq) +  prob2scaledprob(0.5);
				col->transition[MD] = prob2scaledprob(0.0f );//(base_error * indel_freq) +  prob2scaledprob(0.5);
				col->transition[MSKIP] = prob2scaledprob(1.0);
				
				col->transition[II] = prob2scaledprob(0.00);
				col->transition[IM] = prob2scaledprob(0.0);
				col->transition[ISKIP] = prob2scaledprob(1.0f);
				
				col->transition[DD] = prob2scaledprob(0.0f );//(1.0 - 0.999);
				col->transition[DM] = prob2scaledprob(0.0f );//0.999);
			}else{
			
				col->transition[MM] = prob2scaledprob( 1.0 - base_error * indel_freq);
				col->transition[MI] = prob2scaledprob(base_error * indel_freq) +  prob2scaledprob(0.5);
				col->transition[MD] = prob2scaledprob(base_error * indel_freq) +  prob2scaledprob(0.5);
				col->transition[MSKIP] = prob2scaledprob(0.0);
				
				col->transition[II] = prob2scaledprob(1.0 - 0.999);
				col->transition[IM] = prob2scaledprob(0.999);
				col->transition[ISKIP] = prob2scaledprob(0.0f);
				
				col->transition[DD] = prob2scaledprob(1.0 - 0.999);
				col->transition[DM] = prob2scaledprob(0.999);
			}
			
			
			col->transition_e[MM] =  prob2scaledprob(0.0);
			col->transition_e[MI] =  prob2scaledprob(0.0);
			col->transition_e[MD] =  prob2scaledprob(0.0);
			col->transition_e[MSKIP] =  prob2scaledprob(0.0);
			
			col->transition_e[II] =  prob2scaledprob(0.0);
			col->transition_e[IM] =  prob2scaledprob(0.0);
			col->transition_e[ISKIP] =  prob2scaledprob(0.0);
			
			col->transition_e[DD] =  prob2scaledprob(0.0);
			col->transition_e[DM] =  prob2scaledprob(0.0);
		}
		
	}
	
	// init all probs to 0
	
	for(i = 0 ; i < model->num_hmms;i++){
		model->silent_to_M[i] = prob2scaledprob(0.0f);
		//model->M_to_silent[i] = prob2scaledprob(0.0f);
		
		model->silent_to_I[i] = prob2scaledprob(0.0f);
		//model->I_to_silent[i] = prob2scaledprob(0.0f);
		
		model->silent_to_M_e[i] = prob2scaledprob(0.0f);
		//model->M_to_silent_e[i] = prob2scaledprob(0.0f);
		
		model->silent_to_I_e[i] = prob2scaledprob(0.0f);
		//model->I_to_silent_e[i] = prob2scaledprob(0.0f);

		
	}
	model->skip = prob2scaledprob(0.0f);
	model->skip_e = prob2scaledprob(0.0f);
	
	if(rs->type[key] == 'B'){// barcodes all have same length & equal prior probability... 
		for(i = 0 ; i < model->num_hmms;i++){
			model->silent_to_M[i] = prob2scaledprob(1.0 / (float) model->num_hmms);
			//model->M_to_silent[i] = prob2scaledprob(1.0);
			
			model->silent_to_I[i] = prob2scaledprob(0.0f);
			//model->I_to_silent[i] = prob2scaledprob(0.0f);
			
			
		}
		model->skip = prob2scaledprob(0.0);
	}
	
	if(rs->type[key] == 'F'){// fingerprint all have same length & equal prior probability... (of course we specify 1 with NNNNNNNN
		for(i = 0 ; i < model->num_hmms;i++){
			model->silent_to_M[i] = prob2scaledprob(1.0 / (float) model->num_hmms);
			//model->M_to_silent[i] = prob2scaledprob(1.0);
		}
		model->skip = prob2scaledprob(0.0);
	}
	
	
	if(rs->type[key] == 'S'){// fingerprint all have same length & equal prior probability... (of course we specify 1 with NNNNNNNN
		for(i = 0 ; i < model->num_hmms;i++){
			model->silent_to_M[i] = prob2scaledprob(1.0 / (float) model->num_hmms);
			//model->M_to_silent[i] = prob2scaledprob(1.0);
		}
		model->skip = prob2scaledprob(0.0);
	}
	
	if(rs->type[key] == 'P'){// Partial - can skip and exit at every M / I state.... 
		len = model->hmms[0]->num_columns;
		for(i = 0 ; i < model->num_hmms;i++){
			model->silent_to_M[i] = prob2scaledprob(1.0 / (float) model->num_hmms) + prob2scaledprob(1.0 - 0.01);
			//model->M_to_silent[i] = prob2scaledprob(1.0);
			
			for(j = 0; j < len;j++){
				col = model->hmms[i]->hmm_column[j];
				col->transition[MM] = prob2scaledprob( 1.0 - base_error * indel_freq ) + prob2scaledprob(0.99f);
				col->transition[MI] = prob2scaledprob(base_error * indel_freq) +  prob2scaledprob(0.5)+ prob2scaledprob(0.99f);
				col->transition[MD] = prob2scaledprob(base_error * indel_freq) +  prob2scaledprob(0.5)+ prob2scaledprob(0.99f);
				col->transition[MSKIP] = prob2scaledprob(0.01f);
				
				col->transition[II] = prob2scaledprob(1.0 - 0.999)+ prob2scaledprob(0.99f);
				col->transition[IM] = prob2scaledprob(0.999)+ prob2scaledprob(0.99f);
				col->transition[ISKIP] = prob2scaledprob(0.01f);
			
			}
			
		}
		
		model->skip = prob2scaledprob(0.01);
	}


		
	if(rs->type[key] == 'O'){ // optional - like a G, GG or GGG priot probability set to 0.5  - assume length 2 for now,
		len = model->hmms[0]->num_columns;
		for(i = 0 ; i < model->num_hmms;i++){
			//model->silent_to_M[i] = prob2scaledprob(1.0 / (float) model->num_hmms) + prob2scaledprob(0.5);
			//model->M_to_silent[i] = prob2scaledprob(1.0);
			
			model->silent_to_I[i] = prob2scaledprob(1.0 / (float) model->num_hmms) + prob2scaledprob(0.5);
			//model->I_to_silent[i] = prob2scaledprob(1.0 / (float) (len+1));
			
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
		col->transition[MM] = prob2scaledprob( 0.0 );
		col->transition[MI] = prob2scaledprob(0.0);
		col->transition[MD] = prob2scaledprob(0.0);
		col->transition[MSKIP] = prob2scaledprob(0.0);
		
		//col->transition[MQUIT] = prob2scaledprob(1.0 / (float) 2);
		
		col->transition[II] = prob2scaledprob(1.0 - 1.0 / (float)(len+1) );
		col->transition[IM] = prob2scaledprob(0.0);
		col->transition[ISKIP] =  prob2scaledprob(1.0 / (float) (len+1));
		
		
		col->transition[DD] = prob2scaledprob(0.0);
		col->transition[DM] = prob2scaledprob(0.0);
		
		
		col->transition_e[MM] =  prob2scaledprob(0.0);
		col->transition_e[MI] =  prob2scaledprob(0.0);
		col->transition_e[MD] =  prob2scaledprob(0.0);
		
		col->transition_e[II] =  prob2scaledprob(0.0);
		col->transition_e[IM] =  prob2scaledprob(0.0);
		
		col->transition_e[DD] =  prob2scaledprob(0.0);
		col->transition_e[DM] =  prob2scaledprob(0.0);
		
		
		
	}
	
	if(rs->type[key] == 'G'){ // optional - like a G, GG or GGG priot probability set to 0.5  - assume length 2 for now,
		len = model->hmms[0]->num_columns;
		for(i = 0 ; i < model->num_hmms;i++){
			//model->silent_to_M[i] = prob2scaledprob(1.0 / (float) model->num_hmms) + prob2scaledprob(0.5);
			//model->M_to_silent[i] = prob2scaledprob(1.0);
			
			model->silent_to_I[i] = prob2scaledprob(0.8935878);
			//model->I_to_silent[i] = prob2scaledprob(1.0 - 0.195);
			
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
		model->skip = prob2scaledprob(1.0 - 0.8935878);
		col = model->hmms[0]->hmm_column[0];
		col->transition[MM] = prob2scaledprob( 0.0 );
		col->transition[MI] = prob2scaledprob(0.0);
		col->transition[MD] = prob2scaledprob(0.0);
		
		//col->transition[MQUIT] = prob2scaledprob(1.0 / (float) 2);
		
		col->transition[II] = prob2scaledprob(0.195);
		col->transition[IM] = prob2scaledprob(0.0);
		
		col->transition[DD] = prob2scaledprob(0.0);
		col->transition[DM] = prob2scaledprob(0.0);
		
		
		col->transition_e[MM] =  prob2scaledprob(0.0);
		col->transition_e[MI] =  prob2scaledprob(0.0);
		col->transition_e[MD] =  prob2scaledprob(0.0);
		
		col->transition_e[II] =  prob2scaledprob(0.0);
		col->transition_e[IM] =  prob2scaledprob(0.0);
		
		col->transition_e[DD] =  prob2scaledprob(0.0);
		col->transition_e[DM] =  prob2scaledprob(0.0);
	}
	
	if(rs->type[key] == 'R'){// read - skip impossible; 
		for(i = 0 ; i < model->num_hmms;i++){
			model->silent_to_I[i] = prob2scaledprob(1.0 / (float) model->num_hmms);
	//		model->I_to_silent[i] = prob2scaledprob(1.0 / (float) assumed_length);
		}
		col = model->hmms[0]->hmm_column[0];
		for(c = 0; c < 5;c++){
			
			col->m_emit[c] =background[c];
			
			col->i_emit[c] = background[c];
			col->i_emit_e[c] =  prob2scaledprob(0.0f);
			col->m_emit_e[c] =  prob2scaledprob(0.0f);
		}
		col->transition[MM] = prob2scaledprob( 0.0);
		col->transition[MI] = prob2scaledprob(0.0);
		col->transition[MD] = prob2scaledprob(0.0);
		col->transition[MSKIP] = prob2scaledprob(0.0);
		
		//col->transition[MQUIT] = prob2scaledprob(1.0 / (float) assumed_length);
		
		col->transition[II] = prob2scaledprob(1.0 - 1.0 / (float) assumed_length );
		col->transition[IM] = prob2scaledprob(0.0);
		col->transition[ISKIP] = prob2scaledprob(1.0 / (float) assumed_length);
		
		col->transition[DD] = prob2scaledprob(0.0);
		col->transition[DM] = prob2scaledprob(0.0);
		
		
		col->transition_e[MM] =  prob2scaledprob(0.0);
		col->transition_e[MI] =  prob2scaledprob(0.0);
		col->transition_e[MD] =  prob2scaledprob(0.0);
		
		col->transition_e[II] =  prob2scaledprob(0.0);
		col->transition_e[IM] =  prob2scaledprob(0.0);
		
		col->transition_e[DD] =  prob2scaledprob(0.0);
		col->transition_e[DM] =  prob2scaledprob(0.0);
		
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
	//fprintf(stderr,"Skip:%f Self:%f Next:%f\n",scaledprob2prob(model->skip), scaledprob2prob(model->random_self) , scaledprob2prob(model->random_next));
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
			
			for(c = 0; c < 9;c++){
				fprintf(stderr,"%0.4f ", scaledprob2prob(col->transition[c]));

			}
			fprintf(stderr,"%0.3f %0.3f %0.3f\n",scaledprob2prob(col->transition[MM])+scaledprob2prob(col->transition[MI])+scaledprob2prob(col->transition[MD])  ,scaledprob2prob(col->transition[II])+scaledprob2prob(col->transition[IM]),  scaledprob2prob(col->transition[DD]) + scaledprob2prob(col->transition[DM]) );
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
			
			for(c = 0; c < 9;c++){
				fprintf(stderr,"%0.4f ", scaledprob2prob(col->transition_e[c]));
				
			}
			fprintf(stderr,"%0.3f %0.3f %0.3f\n",scaledprob2prob(col->transition_e[MM])+scaledprob2prob(col->transition_e[MI])+scaledprob2prob(col->transition_e[MD])  ,scaledprob2prob(col->transition_e[II])+scaledprob2prob(col->transition_e[IM]),  scaledprob2prob(col->transition_e[DD]) + scaledprob2prob(col->transition_e[DM]) );
		}
		
		
		
		
	}
	fprintf(stderr,"Links:silent to\n");
	
	for(i = 0; i < model->num_hmms;i++){
		fprintf(stderr,"%d	%f	%f	%f	%f\n",i, scaledprob2prob(  model->silent_to_M[i]), scaledprob2prob(  model->silent_to_I[i]),scaledprob2prob(   model->silent_to_M_e[i]),scaledprob2prob(  model->silent_to_I_e[i]));
	}
	fprintf(stderr,"Links:to silent \n");
	
	
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
	
	if(model->silent_to_I){
		free(model->silent_to_I);
		free(model->silent_to_I_e);
	}
	
	
	free(model);// = malloc(sizeof(struct model));
	//assert(model != 0);

	
	//return model;
}


/*


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
		col->transition[NEXT] = org->hmms[0]->hmm_column[i]->transition[NEXT]; //  prob2scaledprob( max);
		col->transition[SELF] =  org->hmms[0]->hmm_column[i]->transition[SELF];// prob2scaledprob(1.0 - max) +  prob2scaledprob(0.5);
		
		col->transition_e[NEXT] = org->hmms[0]->hmm_column[i]->transition_e[NEXT]; //  prob2scaledprob( max);
		col->transition_e[SELF] =  org->hmms[0]->hmm_column[i]->transition_e[SELF];// prob2scaledprob(1.0 - max) +  prob2scaledprob(0.5);
		
		
		for(j = 0; j< model->num_hmms-1;j++){
			col->long_transition[j] = org->hmms[0]->hmm_column[i]->long_transition[j];// prob2scaledprob(1.0 - max) +  prob2scaledprob(0.5  /  (float) (model->num_hmms-1));
			col->long_transition_e[j] = org->hmms[0]->hmm_column[i]->long_transition_e[j];
		}
		
		
	}
	

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
			
			
			
			
			col->transition[MM] =   org->hmms[i]->hmm_column[j]->transition[MM];//  prob2scaledprob( 0.999);
			col->transition[MI] = org->hmms[i]->hmm_column[j]->transition[MI];// prob2scaledprob(1.0 - 0.999) +  prob2scaledprob(0.5);
			col->transition[MD] = org->hmms[i]->hmm_column[j]->transition[MD];// prob2scaledprob(1.0 - 0.999) +  prob2scaledprob(0.5);
			
			col->transition[II] = org->hmms[i]->hmm_column[j]->transition[II];// prob2scaledprob(1.0 - 0.999);
			col->transition[IM] = org->hmms[i]->hmm_column[j]->transition[IM];// prob2scaledprob(0.999);
			
			col->transition[DD] = org->hmms[i]->hmm_column[j]->transition[DD];// prob2scaledprob(1.0 - 0.999);
			col->transition[DM] = org->hmms[i]->hmm_column[j]->transition[DM];// prob2scaledprob(0.999);
			
			col->transition_e[MM] =   org->hmms[i]->hmm_column[j]->transition_e[MM];//  prob2scaledprob( 0.999);
			col->transition_e[MI] = org->hmms[i]->hmm_column[j]->transition_e[MI];// prob2scaledprob(1.0 - 0.999) +  prob2scaledprob(0.5);
			col->transition_e[MD] = org->hmms[i]->hmm_column[j]->transition_e[MD];// prob2scaledprob(1.0 - 0.999) +  prob2scaledprob(0.5);
			
			col->transition_e[II] = org->hmms[i]->hmm_column[j]->transition_e[II];// prob2scaledprob(1.0 - 0.999);
			col->transition_e[IM] = org->hmms[i]->hmm_column[j]->transition_e[IM];// prob2scaledprob(0.999);
			
			col->transition_e[DD] = org->hmms[i]->hmm_column[j]->transition_e[DD];// prob2scaledprob(1.0 - 0.999);
			col->transition_e[DM] = org->hmms[i]->hmm_column[j]->transition_e[DM];// prob2scaledprob(0.999);
			
			
			
			for(c = 0;c <  (model->num_hmms-1);c++){
				col->long_transition[c] = org->hmms[i]->hmm_column[j]->long_transition[c];
				col->long_transition_e[c] = org->hmms[i]->hmm_column[j]->long_transition_e[c];
			}
			
			
			
		}
		
	}
	
	return model;
}
*/
/*
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
		col_target->transition_e[NEXT] = logsum(col_target->transition_e[NEXT] ,col_source->transition_e[NEXT] );
		col_target->transition_e[SELF] = logsum(col_target->transition_e[SELF] ,col_source->transition_e[SELF] );

		for(j = 0; j < target->num_hmms-1;j++){
			//col_target->long_transition[j] = org->hmms[0]->hmm_column[i]->long_transition[j];// prob2scaledprob(1.0 - max) +  prob2scaledprob(0.5  /  (float) (model->num_hmms-1));
			col_target->long_transition_e[j] = logsum(col_target->long_transition_e[j], col_source->long_transition_e[j]);// org->hmms[0]->hmm_column[i]->long_transition_e[j];
		}
		
		
	}
	
	
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
			
			
			
			
						
			col_target->transition_e[MM] = logsum(col_target->transition_e[MM] ,col_source->transition_e[MM] );//
			col_target->transition_e[MI] = logsum(col_target->transition_e[MI] ,col_source->transition_e[MI] );//
			col_target->transition_e[MD] = logsum(col_target->transition_e[MD] ,col_source->transition_e[MD] );//
			
			col_target->transition_e[II] = logsum(col_target->transition_e[II] ,col_source->transition_e[II] );//
			col_target->transition_e[IM] = logsum(col_target->transition_e[IM] ,col_source->transition_e[IM] );//
			
			col_target->transition_e[DD] = logsum(col_target->transition_e[DD] ,col_source->transition_e[DD] );//
			col_target->transition_e[DM] = logsum(col_target->transition_e[DM] ,col_source->transition_e[DM] );//
		
			
			
			
			for(c = 0;c <  (target->num_hmms-1);c++){

				col_target->long_transition_e[c] = logsum(col_target->long_transition_e[c],col_source->long_transition_e[c] ); //  org->hmms[i]->hmm_column[j]->long_transition_e[c];
			}
			
			
			
		}
		
	}
	return target;
}

*/

struct model_bag* copy_model_bag(struct model_bag* org)
{
	struct model_bag* copy = 0;
	int i,j;
	copy =  malloc(sizeof(struct model_bag));
	
	assert(copy!=0);
	
	copy->model = malloc(sizeof(struct model* ) * org->num_models);//   param->read_structure->num_segments);
	
	assert(copy->model);
	copy->num_models  = org->num_models;
	copy->total_hmm_num = org->total_hmm_num;
	for(i = 0; i < org->num_models;i++){
		copy->model[i] = malloc_model_according_to_read_structure(org->model[i]->num_hmms,  org->model[i]->hmms[0]->num_columns);
		
		copy->model[i]  = copy_model_parameters(org->model[i],copy->model[i]) ;
	}
	
	copy->path = malloc(sizeof(int*) * MAX_SEQ_LEN);
	copy->dyn_prog_matrix = malloc(sizeof(float*) * MAX_SEQ_LEN );
	
	for (i = 0; i < MAX_SEQ_LEN;i++){
		copy->path[i] = malloc(sizeof(int)* (copy->total_hmm_num +1) );
		copy->dyn_prog_matrix[i] = malloc(sizeof(float) * (copy->total_hmm_num +1) );
	}
	
	copy->transition_matrix = malloc(sizeof(float*) * (copy->total_hmm_num +1));
	copy->label = malloc(sizeof(int) *  (copy->total_hmm_num +1));
		
	for(i = 0; i < copy->total_hmm_num +1;i++){
		copy->label[i] = org->label[i];
	}
	
	for(i = 0; i < copy->total_hmm_num+1 ;i++){
		copy->transition_matrix[i] = malloc(sizeof(float) * (copy->total_hmm_num +1));
		for(j = 0; j <  copy->total_hmm_num+1 ;j++){
			copy->transition_matrix[i][j] = org->transition_matrix[i][j];
		}
	}
	// hmm parameters....
	
	
	
	
	
	
	
	
	return copy;
}

struct model* copy_model_parameters(struct model* org, struct model* copy )
{
	int i,j,c;
	
	struct hmm_column* org_col = 0;
	struct hmm_column* copy_col = 0;
	for(i = 0; i < 5;i++){
		copy->background_nuc_frequency[i] = org->background_nuc_frequency[i];
	}
	
	for(i = 0; i < org->num_hmms;i++){
		copy->silent_to_I[i] = org->silent_to_I[i];
		copy->silent_to_I_e[i] = org->silent_to_I_e[i];
		copy->silent_to_M[i] = org->silent_to_M[i];
		copy->silent_to_M_e[i] = org->silent_to_M_e[i];
		
		for(j = 0; j < org->hmms[0]->num_columns;j++){
			org_col = org->hmms[i]->hmm_column[j];
			copy_col = copy->hmms[i]->hmm_column[j];
			for(c = 0; c< 5;c++){
				copy_col->i_emit[c] = org_col->i_emit[c];
				copy_col->i_emit_e[c] = org_col->i_emit_e[c];
				
				copy_col->m_emit[c] = org_col->m_emit[c];
				copy_col->m_emit_e[c] = org_col->m_emit_e[c];
			}
			
			for(c = 0; c < 9;c++){
				copy_col->transition[c] = org_col->transition[c];
				copy_col->transition_e[c] = org_col->transition_e[c];
			}
		}
	}
	copy->skip = org->skip;
	copy->skip_e = org->skip_e;
	copy->num_hmms = org->num_hmms;
	return copy;
}



struct model* reestimate(struct model* m, int mode)
{
	int i,j,c;
	
	struct hmm_column* m_col = 0;
	//struct hmm_column* copy_col = 0;
	
	float sum = 0.0f;
	
	// silent to M /I ....
	// add pseudocount of 1;
	
	// mode 0
	// train everything...
	
	//mode 1
	// train everything apart from ssilent to & skip....
	
	//mode2
	// only train emission probabilities. ....
	
	
	
	if(mode < 1){
		sum = prob2scaledprob(0.0);
		for(i = 0; i < m->num_hmms;i++){
			sum = logsum(sum, logsum(m->silent_to_I_e[i],prob2scaledprob(1.0)));
			sum = logsum(sum, logsum(m->silent_to_M_e[i] , prob2scaledprob(1.0)) );
			//fprintf(stderr," silent to I: %f",m->silent_to_I_e[i]);
			//fprintf(stderr," silent to M: %f",m->silent_to_M_e[i]);
		}
		sum = logsum(sum,logsum( m->skip_e , prob2scaledprob(1.0)));
		//fprintf(stderr,"estimated skip: %f\n", m->skip_e);
		
		for(i = 0; i < m->num_hmms;i++){
			m->silent_to_I[i]  =  logsum(m->silent_to_I_e[i] ,prob2scaledprob(1.0)) - sum;
			m->silent_to_M[i]  = logsum(m->silent_to_M_e[i],prob2scaledprob(1.0)) - sum;
		}
		
		m->skip = logsum(m->skip_e ,prob2scaledprob(1.0)) - sum;
		
		//fprintf(stderr,"SKIP: %f\n", m->skip );
		
		//clear counts....
		for(i = 0; i < m->num_hmms;i++){
			m->silent_to_I_e[i] = prob2scaledprob(0.0);
			m->silent_to_M_e[i] = prob2scaledprob(0.0);
		}
		m->skip_e = prob2scaledprob(0.0);
	}
	
	for(i = 0; i < m->num_hmms;i++){
		//copy->silent_to_I[i] = org->silent_to_I[i];
		//copy->silent_to_I_e[i] = org->silent_to_I_e[i];
		//copy->silent_to_M[i] = org->silent_to_M[i];
		//copy->silent_to_M_e[i] = org->silent_to_M_e[i];
		
		//copy->I_to_silent[i] = org->I_to_silent[i];
		//copy->I_to_silent_e[i] =org->I_to_silent_e[i];
		//copy->M_to_silent[i] = org->M_to_silent[i];
		//copy->M_to_silent_e[i] = org->M_to_silent_e[i];
		
		for(j = 0; j < m->hmms[0]->num_columns;j++){
			m_col = m->hmms[i]->hmm_column[j];
			//copy_col = copy->hmms[i]->hmm_column[j];
			sum = prob2scaledprob(0.0f);
			
			for(c = 0; c< 5;c++){
				sum = logsum(sum, logsum(m_col->i_emit_e[c] , prob2scaledprob(1.0f)));
			}
			
			for(c = 0; c< 5;c++){
				m_col->i_emit[c] = logsum(m_col->i_emit_e[c] , prob2scaledprob(1.0f)) - sum;
			}
			
			for(c = 0; c< 5;c++){
				m_col->i_emit_e[c] =  prob2scaledprob(0.0f);
			}
			
			
			sum = prob2scaledprob(0.0f);
			
			for(c = 0; c< 5;c++){
				sum = logsum(sum, logsum(m_col->m_emit_e[c] , prob2scaledprob(1.0f)));
			}
			
			for(c = 0; c< 5;c++){
				m_col->m_emit[c] = logsum(m_col->m_emit_e[c] , prob2scaledprob(1.0f)) - sum;
			}
			
			for(c = 0; c< 5;c++){
				m_col->m_emit_e[c] =  prob2scaledprob(0.0f);
			}

			if(mode < 2){
				// internal hmm states...
				if(j != m->hmms[0]->num_columns-1){
					sum = prob2scaledprob(0.0f);
					
					sum = logsum(sum, logsum(m_col->transition_e[MM] , prob2scaledprob(1.0)));
					sum = logsum(sum, logsum(m_col->transition_e[MI] , prob2scaledprob(1.0)));
					sum = logsum(sum, logsum(m_col->transition_e[MD] , prob2scaledprob(1.0)));
					if(m_col->transition[MSKIP] != prob2scaledprob(0.0)){
						sum = logsum(sum, logsum(m_col->transition_e[MSKIP] , prob2scaledprob(1.0)));
					}
					
					
					m_col->transition[MM] =  logsum(m_col->transition_e[MM] , prob2scaledprob(1.0)) - sum;
					m_col->transition[MI] =  logsum(m_col->transition_e[MI] , prob2scaledprob(1.0)) - sum;
					m_col->transition[MD] =  logsum(m_col->transition_e[MD], prob2scaledprob(1.0)) - sum;
					if(m_col->transition[MSKIP] != prob2scaledprob(0.0)){
						m_col->transition[MSKIP] =  logsum(m_col->transition_e[MSKIP], prob2scaledprob(1.0)) - sum;
					}
					
					
					sum = prob2scaledprob(0.0f);
					
					sum = logsum(sum, logsum(m_col->transition_e[II] ,prob2scaledprob(1.0)));
					sum = logsum(sum, logsum(m_col->transition_e[IM] , prob2scaledprob(1.0)));
					if(m_col->transition[ISKIP] != prob2scaledprob(0.0)){
						sum = logsum(sum, logsum(m_col->transition_e[ISKIP] , prob2scaledprob(1.0)));
					}
					
					m_col->transition[II] =  logsum(m_col->transition_e[II] , prob2scaledprob(1.0)) - sum;
					m_col->transition[IM] =  logsum(m_col->transition_e[IM] , prob2scaledprob(1.0)) - sum;
					if(m_col->transition[ISKIP] != prob2scaledprob(0.0)){
						m_col->transition[ISKIP] =  logsum(m_col->transition_e[ISKIP] , prob2scaledprob(1.0)) - sum;
					}
					sum = prob2scaledprob(0.0f);
					
					sum = logsum(sum, logsum(m_col->transition_e[DD] , prob2scaledprob(1.0)));
					sum = logsum(sum, logsum(m_col->transition_e[DM] , prob2scaledprob(1.0)));
					
					m_col->transition[DD] =  logsum(m_col->transition_e[DD] , prob2scaledprob(1.0)) - sum;
					m_col->transition[DM] =  logsum(m_col->transition_e[DM] , prob2scaledprob(1.0)) - sum;
					
					
					
				}else{ // last hmm column...
					// no transitions from M possible....
					m_col->transition[MM] =  prob2scaledprob(0.0);
					m_col->transition[MI] =  prob2scaledprob(0.0);
					m_col->transition[MD] =  prob2scaledprob(0.0);
					m_col->transition[MSKIP] = prob2scaledprob(1.0);
					
					//either continue i or goto silent state....
					sum = prob2scaledprob(0.0f);
					
					sum = logsum(sum, logsum(m_col->transition_e[II] , prob2scaledprob(1.0)));
					sum = logsum(sum, logsum(m_col->transition_e[ISKIP] , prob2scaledprob(1.0)));
					
					//sum = logsum(sum, logsum(m->I_to_silent_e[i] , prob2scaledprob(1.0)));
					
					m_col->transition[II] =  logsum(m_col->transition_e[II] , prob2scaledprob(1.0)) - sum;
					m_col->transition[ISKIP] =  logsum(m_col->transition_e[ISKIP] , prob2scaledprob(1.0)) - sum;
					//m->I_to_silent[i] =  logsum(m->I_to_silent_e[i] , prob2scaledprob(1.0)) - sum;
					
					
					//no transtition from D possible.
					
					m_col->transition[DD] = prob2scaledprob(0.0);// m_col->transition_e[DD] + prob2scaledprob(1.0) - sum;
					m_col->transition[DM] = prob2scaledprob(0.0);//  m_col->transition_e[DM] + prob2scaledprob(1.0) - sum;
					
					
				}
				
				//m->I_to_silent_e[i] = prob2scaledprob(0.0);
				//m->M_to_silent_e[i] = prob2scaledprob(0.0);
				
			}
			
			for(c = 0; c< 9;c++){
				m_col->transition_e[c] =  prob2scaledprob(0.0f);
			}
		}
	}
	//copy->skip = org->skip;
	//copy->skip_e = copy->skip_e;
	return m;
}


struct model* copy_estimated_parameter(struct model* target, struct model* source )
{
	int i,j,c;
	
	struct hmm_column* target_col = 0;
	struct hmm_column* source_col = 0;
	
	for(i = 0; i < target->num_hmms;i++){
		//copy->silent_to_I[i] = org->silent_to_I[i];
		target->silent_to_I_e[i] = logsum(target->silent_to_I_e[i], source->silent_to_I_e[i]);
		//copy->silent_to_M[i] = org->silent_to_M[i];
		target->silent_to_M_e[i] = logsum(target->silent_to_M_e[i] ,source->silent_to_M_e[i]);//  org->silent_to_M_e[i];
		
		//copy->I_to_silent[i] = org->I_to_silent[i];
		//target->I_to_silent_e[i] = logsum(target->I_to_silent_e[i] , source->I_to_silent_e[i]);//   org->I_to_silent_e[i];
		//copy->M_to_silent[i] = org->M_to_silent[i];
		//target->M_to_silent_e[i] = logsum(target->M_to_silent_e[i],source->M_to_silent_e[i]); // org->M_to_silent_e[i];
		
		for(j = 0; j < target->hmms[0]->num_columns;j++){
			target_col = target->hmms[i]->hmm_column[j];
			source_col = source->hmms[i]->hmm_column[j];
			for(c = 0; c< 5;c++){
				//copy_col->i_emit[c] = org_col->i_emit[c];
				target_col->i_emit_e[c] = logsum(target_col->i_emit_e[c] , source_col->i_emit_e[c] );  //org_col->i_emit_e[c];
				
				//copy_col->m_emit[c] = org_col->m_emit[c];
				target_col->m_emit_e[c] = logsum(target_col->m_emit_e[c], source_col->m_emit_e[c] );// org_col->m_emit_e[c];
				
				
			}
			
			for(c = 0; c < 9;c++){
				//copy_col->transition[c] = org_col->transition[c];
				target_col->transition_e[c] = logsum (target_col->transition_e[c], source_col->transition_e[c] ); // org_col->transition_e[c];
				
			}
		}
	}
	//copy->skip = org->skip;
	target->skip_e = logsum(target->skip_e, source->skip_e );// copy->skip_e;
	
	
	
	return target;
}



struct model_bag* init_model_bag(struct parameters* param,float* back)
{
	int i,j,c;
	//int average_length = 12;
	int read_length = 1;
	int segment_length;
	
	struct model_bag* mb = 0;
	mb = malloc(sizeof(struct model_bag));
	
	assert(mb!=0);
	
	mb->model = malloc(sizeof(struct model* ) * param->read_structure->num_segments);
	
	assert(mb->model);
	
	
	mb->f_score = prob2scaledprob(0.0f);
	mb->b_score = prob2scaledprob(0.0f);
	mb->num_models = param->read_structure->num_segments;
	// get read length estimate...
	read_length = param->average_read_length;
	for(i = 0; i < mb->num_models;i++){
		//mb->model[i] = malloc_model_according_to_read_structure(param->read_structure,i);
		//fprintf(stderr," %d\n",read_length );
		if(param->read_structure->type[i] == 'G'){
			read_length = read_length -2;
			
		}else{
			read_length = read_length - (int)strlen(param->read_structure->sequence_matrix[i][0]);
		}
		
		
		
	}
	fprintf(stderr,"READlength: %d\n",read_length);
	
	mb->total_hmm_num = 0;
	
	
	for(i = 0; i < mb->num_models;i++){
		mb->model[i] = malloc_model_according_to_read_structure(param->read_structure->numseq_in_segment[i],(int)strlen(param->read_structure->sequence_matrix[i][0]));
		segment_length = 0;
		if(param->read_structure->type[i] == 'G'){
			segment_length = 2;
			
		}
		if(param->read_structure->type[i]  == 'R'){
			segment_length = read_length;
		}
		
		
		
		mb->model[i] = init_model_according_to_read_structure(mb->model[i], param, i,back,segment_length);
		//print_model(mb->model[i]);
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
			//fprintf(stderr,"%d %d	%d %d\n",c,mb->label[c],mb->label[c] & 0xFFFF, (mb->label[c] >> 16) & 0x7FFF);
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
			//fprintf(stderr,"%d, %d, %d %d\n ", j,   mb->label[j],mb->label[j] & 0xFFFF, (mb->label[j] >> 16) & 0x7FFF);
		}
		
		// remain in the same state....
		mb->transition_matrix[i][i] = 1;
	}
	/*
	for(i = 0; i < mb->total_hmm_num ;i++){
		for(j = 0; j <  mb->total_hmm_num ;j++){
			fprintf(stderr,"%f ",mb->transition_matrix[i][j] );
		}
		fprintf(stderr,"\n");
	}
	fprintf(stderr,"\n");
	*/
	
	return mb;
}


void free_model_bag(struct model_bag* mb)
{
	int i;
	
	
	//mb->transition_matrix = malloc(sizeof(float*) * (mb->total_hmm_num +1));
	//mb->label = malloc(sizeof(int) *  (mb->total_hmm_num +1));
	
	for (i = 0; i < MAX_SEQ_LEN;i++){
		free(mb->path[i]);// = malloc(sizeof(int)* (mb->total_hmm_num +1) );
		free(mb->dyn_prog_matrix[i]);// = malloc(sizeof(float) * (mb->total_hmm_num +1) );
	}
	
	free(mb->path);// = malloc(sizeof(int*) * MAX_SEQ_LEN);
	free(mb->dyn_prog_matrix);// = malloc(sizeof(float*) * MAX_SEQ_LEN );
	
	
	for(i = 0; i < mb->total_hmm_num+1 ;i++){
		free(mb->transition_matrix[i]);//  = malloc(sizeof(float) * (mb->total_hmm_num +1));
		
	}
	free(mb->transition_matrix);
	free(mb->label);
	
	
	for(i = 0; i < mb->num_models;i++){
		free_model(mb->model[i]);
	}

	
	
	free(mb->model);// = malloc(sizeof(struct model* ) * param->read_structure->num_segments);
	
	
	free(mb);// = malloc(sizeof(struct model_bag));

	
}


















