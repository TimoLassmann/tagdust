#if HAVE_CONFIG_H
#include "config.h"
#endif

#include "kslib.h"
#include "tagdust2.h"
#include "interface.h"
#include "io.h"
#include "barcode_hmm.h"


#include <stdio.h>
#include "misc.h"



int estimateQthreshold(struct parameters* param, struct sequence_stats_info* ssi)
{
	int i,j;
	int status;
	
	//Warning - I set the seed to number 42 because the results are completely robust with 400000 sequences...
	// if we reduce the latter we need to test for consistency between runs...
	//unsigned int seed = (unsigned int) (time(NULL) * ( 42));
	unsigned int seed = 0;
	
	if(param->seed){
		seed = param->seed;
	}else{
		seed = (unsigned int) (time(NULL) * ( 42));
	}
	
	srand(seed);
	
	int binsize = 0;
	
#if DEBUG
#if RTEST
	int num_test_sequences = 4000;
#else
	int num_test_sequences = 400;
#endif
#else
#if RTEST
	int num_test_sequences = 4000;
#else
	int num_test_sequences = 400000;
#endif
#endif
	struct model_bag* mb = 0;
	
	struct read_info** ri = 0;
	
	ri = malloc_read_info(ri,num_test_sequences);
	double TP,FP,TN,FN;
	TP = 0.0;
	FP = 0.0;
	TN = 0.0;
	FN = 0.0;
	
	binsize = (int) (num_test_sequences / 4);
	
	int readnum = 0;
	// simulate reads from model with increased fluffyness...
	param->sequencer_error_rate = 0.05;//sim_errors[c];
	
	mb = init_model_bag(param, ssi);
	
	//fprintf(stderr,"Testing the model %e\n", sim_errors[c]);
	for(i = 0; i < mb->num_models;i++){
		if(param->read_structure->type[i] == 'B'){
			//hasbarcode = i;
			for(j = 0 ; j < mb->model[i]->num_hmms-1;j++){
				mb->model[i]->silent_to_M[j][0] = prob2scaledprob(1.0 / (float)( mb->model[i]->num_hmms-1));
			}
			mb->model[i]->silent_to_M[mb->model[i]->num_hmms-1][0] = prob2scaledprob(0.0);
		}
		if(param->read_structure->type[i] == 'S'){
			//hasbarcode = i;
			for(j = 0 ; j < mb->model[i]->num_hmms-1;j++){
				mb->model[i]->silent_to_M[j][0] = prob2scaledprob(1.0 / (float)( mb->model[i]->num_hmms-1));
			}
			mb->model[i]->silent_to_M[mb->model[i]->num_hmms-1][0] = prob2scaledprob(0.0);
		}
		//print_model(mb->model[i]);
	}
	KSL_DPRINTF2(("Will emit reads\n" ));
	for(i = 0; i < binsize*2;i++){
		
		
		if((status = emit_read_sequence(mb, ri[readnum],ssi->average_length,&seed)) != kslOK) KSLIB_XFAIL(kslFAIL,param->errmsg,"Failed to emit a sequence from read HMM.\n");
		//ri[readnum] = emit_read_sequence(mb, ri[readnum],ssi->average_length,&seed);
		
		
		
		
		ri[readnum]->read_type = 0;
		FN++;
		readnum++;
	}
	KSL_DPRINTF2(("Will emit random\n" ));
	for(i = 0; i <  binsize+binsize ;i++){
		if((status = emit_random_sequence(mb, ri[readnum],ssi->average_length,&seed)) != kslOK) KSLIB_XFAIL(kslFAIL,param->errmsg,"Failed to emit a random sequence.\n");
		
		//ri[readnum] = emit_random_sequence(mb, ri[readnum],ssi->average_length,&seed);
		ri[readnum]->read_type = 1;
		TN++;
		
		readnum++;
		if(readnum == num_test_sequences){
			break;
		}
	}
	
	free_model_bag(mb);
	
	param->sequencer_error_rate = 0.05;
	
	mb = init_model_bag(param, ssi);

	j = 0;
	for(i = 0; i < readnum;i++){
		if(ri[i]->len >=ssi->max_seq_len){
			ssi->max_seq_len = ri[i]->len;
			j = 1;
		}
	}
	
	if(j){
		sprintf(param->buffer,"Long sequence found. Need to realloc model...\n");
		param->messages = append_message(param->messages, param->buffer);
		free_model_bag(mb);
		mb = init_model_bag(param, ssi);
	}
	KSL_DPRINTF2(("Got here\n" ));
	if((status =run_pHMM(0,mb,ri,param,0, readnum,MODE_GET_PROB)) != kslOK) KSLIB_XFAIL(kslFAIL,param->errmsg,"run_pHMM failed\n");
	
	//mb =  run_pHMM(0,mb,ri,param,0, readnum,MODE_GET_PROB);
#if DEBUG
	//print_model(mb->model[0]);
	fprintf(stderr,"LENGTH: %f\n",(float)ssi->average_length);
	fprintf(stderr,"in random:  %f\n",1.0 - (1.0 / (float)ssi->average_length));
#endif
	
	free_model_bag(mb);
	qsort(ri,readnum, sizeof(struct read_info*), qsort_ri_mapq_compare);
	
	
	
	double kappa = 0.0;
	double tmp = 0.0;
	
	double P_e = 0.0;
	
	double P_o = 0.0;

	
	
	float thres[6];
	thres[0] = 1000.0;
	thres[1] = 1000.0;
	thres[2] = 1000.0;
	thres[3] = 0.0;
	thres[4] = 1000.0;
	thres[5] = 1000.0;
	
	float sensitivity, specificity;
	for(i = 0; i < readnum;i++){
		if(ri[i]->read_type){
			FP += 1.0;
			TN -= 1.0;
		}else{
			TP += 1.0;
			FN -= 1.0;
		}
		//fprintf(stderr,"%d	%f	%d\n",i,ri[i]->mapq, ri[i]->read_type);
		
		sensitivity = TP/( TP + FN );
		specificity =  TN / ( TN + FP);
		
		if(FP /(FP+ TP) < 0.01){
			thres[0] =  ri[i]->mapq;
		}else if(FP /(FP+ TP) < 0.05){
			thres[1] = ri[i]->mapq;
		}else if(FP /(FP+ TP) < 0.1){
			thres[2] = ri[i]->mapq;
		}
		
		if(sensitivity + specificity > thres[3]){
			thres[3] = specificity + sensitivity;
			thres[4] = ri[i]->mapq;
		}
		P_e = ((TP+FN) / (double)readnum) * ((TP+FP) / (double)readnum) +  ( ((FP+TN) / (double)readnum  ) * ((FN+TN) / (double)readnum));
		P_o =(TP+TN)/(double)readnum ;
		
		tmp = (P_o - P_e) / (1.0 - P_e);
		
		
		if(tmp > kappa ){
			
			kappa = tmp;
			thres[5] = ri[i]->mapq;
			//fprintf(stderr,"%d	KAPPA:%f	%f\n",i,kappa,ri[i]->mapq);
		}
		
	}
	
	if(thres[4] < 20){
		param->confidence_threshold  =thres[4];
	}else{
		param->confidence_threshold  = 20;
	}


	
#if DEBUG
	fprintf(stderr,"FDR:0.01: %f\n", thres[0]);
	fprintf(stderr,"FDR:0.05: %f\n", thres[1]);
	fprintf(stderr,"FDR:0.1: %f\n", thres[2]);
	fprintf(stderr,"Sen+spe: %f\n", thres[4]);
	
	fprintf(stderr,"Kappa: %f at %f\n",  kappa,thres[5]   );
	
	fprintf(stderr,"Selected Threshold: %f\n", param->confidence_threshold );
	//fprintf(stderr,"Sensitivity: %f\n", sensitivity);
	//fprintf(stderr,"Specificity: %f\n", specificity);
#endif
	sprintf(param->buffer,"Selected Threshold:: %f\n", param->confidence_threshold );

	param->messages = append_message(param->messages, param->buffer);
	free_read_info(ri,  num_test_sequences);
	return kslOK;
ERROR:
	return kslFAIL;
}


