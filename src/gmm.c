//
//  gmm.c
//  tagdust2
//
//  Created by lassmann on 2/25/13.
//  Copyright (c) 2013 lassmann. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <float.h>
#include <limits.h>
#include <math.h>

#include "interface.h"
#include "io.h"

#include "misc.h"
#include "gmm.h"



void run_gmm_on_sequences(struct read_info** ri, int numseq)
{
	float* x = 0;
	int i;
	
	x = malloc(sizeof(float)* numseq );
	
	for(i = 0; i < numseq;i++){
		x[i] = ri[i]->mapq;
	}
	
	struct gmm_model* model = 0;
	model  = malloc(sizeof(struct gmm_model));
	model->mean = malloc(sizeof(double) * MAX_NUM_MIXTURES);
	model->sigma = malloc(sizeof(double) * MAX_NUM_MIXTURES);
	model->p = malloc(sizeof(double) * MAX_NUM_MIXTURES);
	
	
	for(i = 0; i < 4;i++){
		model = run_miniEM( x, model , i+1 , numseq);
		fprintf(stderr,"\n\n");
		//exit(0);

	}	
	
	free(model->mean);// = malloc(sizeof(double) * MAX_NUM_MIXTURES);
	free(model->sigma);// = malloc(sizeof(double) * MAX_NUM_MIXTURES);
	free(model->p);// = malloc(sizeof(double) * MAX_NUM_MIXTURES);
	
	free(model);//  = malloc(sizeof(struct gmm_model));
	
	
	
}


struct gmm_model* run_miniEM(float* x, struct gmm_model*  model, int k, int n)
{
	int i,j,g;
	int changed = INT_MAX;
	int max_try = 10;
	int try = 0;
	double tmp = 0;
	double likelihood = 0;
	double old_likelihood = 1e-10f;
	double* sum = 0;
	//double* raw = 0;
	//double t = 0.0;
	double likelihood_threshold = 1e-6f;
	
	double* mean = 0;
	double* sigma = 0;
	double* e_mean = 0;
	double* e_sigma = 0;
	
	double* p = 0;
	double* e_p = 0;
	
	
	double S0,S1,S2,stdev,org_mean;
	
	S0 = 0.0;
	S1 = 0.0;
	S2 = 0.0;
	for(i = 0; i < n;i++){
		S0 += 1;
		S1 += x[i];
		S2 += x[i] *x[i];
	}
	
	//fprintf(stderr,"MEAN:  %f	STDEV:%f\n",  S1 / S0,sqrt(  (S0 * S2 - pow(S1,2.0))   /  (  S0 *(S0-1.0) )) );
	org_mean = S1 / S0;
	stdev = sqrt(  (S0 * S2 - pow(S1,2.0))   /  (  S0 *(S0-1.0) )) ;
	
	int* assign = malloc(sizeof(int) * n);
	
	mean = malloc(sizeof(double) * k);
	sum = malloc(sizeof(double) * k);
	sigma = malloc(sizeof(double) * k);
	p = malloc(sizeof(double) * k);
	
	e_mean = malloc(sizeof(double) * k);
	e_sigma = malloc(sizeof(double) * k);
	e_p = malloc(sizeof(double) * k);
	
	sum = malloc(sizeof(double) * k);
	
	double** tmp_p_matrix = 0;
	
	tmp_p_matrix = malloc(sizeof(double*) * n);
	for(i = 0; i < n;i++){
		tmp_p_matrix[i] = malloc(sizeof(double)*k);
	}
	old_likelihood = FLT_MAX;
	
	for(try = 0;try < max_try;try++){
		
		//for(i = 0; i < n;i++){
		//	assign[i] = (int)( (double) rand_r(seed) / (double) RAND_MAX * (double) k) ;
		
		//}
			
		for(i = 0 ; i < k;i++){
	//		e_mean[i] =   e_mean[i] / sum[i];
			
			e_mean[i] =  org_mean +  stdev  / (k+1)  *(i+1)  -  (stdev / 2 );
		//	fprintf(stderr,"Starting:%d	%f\n",i, e_mean[i]);
		//	e_mean[i] = ( (double) rand_r(seed) / (double) RAND_MAX * (double) n) ;
		}
		
		for(i = 0; i < n;i++){
			tmp = FLT_MAX;
			g = -1;
			for(j = 0; j < k;j++){
				if(fabs((float)x[i] -e_mean[j]) < tmp){
					tmp = fabs((float)x[i] -e_mean[j]) ;
					g = j;
				}
			}
			assign[i]  = g;
		}
		
		
		changed = INT_MAX;
		
		while (changed){
			changed = 0;
			likelihood = 0.0;
			for(i = 0; i < n;i++){
				tmp = FLT_MAX;
				g = -1;
				for(j = 0; j < k;j++){
					if(fabs((float)x[i] -e_mean[j]) < tmp){
						tmp = fabs((float)x[i] -e_mean[j]) ;
						g = j;
					}
				}
				likelihood += tmp;
				if(g != assign[i]  ){
					assign[i] = g;
					changed++;
				}
				
			}
			for(i = 0 ; i < k;i++){
				e_mean[i] = 0;
				sum[i] = 0.0;
			}
			for(i = 0; i < n;i++){
				e_mean[assign[i]] += x[i];
				sum[assign[i]] += 1.0;
			}
			for(i = 0 ; i < k;i++){
				
				if(!sum[i]){
					e_mean[i]  = -FLT_MAX;
				}else{
					e_mean[i] = e_mean[i] / sum[i];
				}
			//	fprintf(stderr,"Mean:%f\n",e_mean[i]);
			}
			
			//fprintf(stderr,"changed:%d\n",changed);
		}
		//fprintf(stderr,"TRY:%d	%f\n",try,likelihood);
		
		if(likelihood < old_likelihood){
			//fprintf(stderr,"Found better:%d	%f\n",try,likelihood);
			for(i = 0 ; i < k;i++){
				mean[i]  = e_mean[i];
				//	fprintf(stderr,"%f	",mean[i]);
			}
			//fprintf(stderr,"\n");
			old_likelihood = likelihood;
		}
		
	}
	fprintf(stderr,"Model:\n");
	for(i = 0 ; i < k;i++){
		//mean[i]  = e_mean[i];
			fprintf(stderr,"%d	%f\n",i,mean[i]);
	}
	//fprintf(stderr,"\n");

	//return model;
	
	
	for(i = 0 ; i < k;i++){
		tmp = -FLT_MAX;
		for(j = 0; j < n;j++){
			sigma[i] = (x[j] - mean[i])* (x[j] - mean[i]);
			sum[i] += 1.0;
		}
	}
	for(i = 0; i < k;i++){
		if(sum[i]){
			sigma[i] = sqrtf( sigma[i] / sum[i]);
			fprintf(stderr,"SIGMA:%f\n", sigma[i]);
		}else{
			sigma[i] = MIN_STDEV;
		}
		//if(sigma[i] < MIN_STDEV){
		//	sigma[i] = MIN_STDEV;//FLT_MIN;
		//}
	}
	j = 0;
	
	for(i  = 0;i < k;i++){
		if(mean[i] !=  -FLT_MAX){
			j++;
		}
	}
	//for(i = 0; i < n;i++){
	
	//	fprintf(stderr,"%f ",x[i]);
	//}
	//fprintf(stderr,"\n");
	
	for(i  = 0;i < k;i++){
		if(mean[i] ==  -FLT_MAX){
			p[i] = prob2scaledprob(0.0);// fprintf(stderr,"GAGA\n");
			mean[i] = 0;
			sigma[i] = 0;
		}else{
			p[i] = prob2scaledprob( 1.0 / (float)j);
		}
		fprintf(stderr,"INIT:	%d	%f	%f	%f	%d\n",i,p[i],mean[i],sigma[i],j);
	}
	old_likelihood = 1e-10f;
	
	//return model;
	
	for(g = 0; g < 100;g++){
		likelihood = prob2scaledprob(1.0);
		for(j = 0; j < k;j++){
			sum[j] =prob2scaledprob( 0.0);
			e_mean[j] = prob2scaledprob(0.0);
			e_p[j] = prob2scaledprob(0.0);
			e_sigma[j] = 0.0;
		}
		
		//estimate new means;
		
		for(i = 0; i < n;i++){
			//t = (double)i / (double)n - (double) n / 2.0;
 			tmp = prob2scaledprob(0.0);
			for(j = 0; j < k;j++){
				//raw[j]  = p[j] * pdf((float)i, mean[j],sigma[j]);
				if(p[j] != prob2scaledprob(0.0)){
					tmp_p_matrix[i][j]  = p[j] + log_pdf(x[i], mean[j],sigma[j]);
					
					//tmp_p_matrix[i][j]  = p[j] + log_truncated_pdf(i, mean[j],sigma[j],0,n);
					
					//fprintf(stderr,"%d	%d	%f\n",i,j,p[j] + log_pdf(i, mean[j],sigma[j]));
					tmp = logsum(tmp  ,tmp_p_matrix[i][j]) ;// raw[j];
				}
				//fprintf(stderr,"%f	\n",raw[j]);
			}
			if(tmp == prob2scaledprob(0.0)){
				fprintf(stderr,"GAGA:%d\n",i);
				for(j = 0; j < k;j++){
					//raw[j]  = p[j] * pdf((float)i, mean[j],sigma[j]);
					if(p[j] != prob2scaledprob(0.0)){
						//tmp_p_matrix[i][j]  = p[j] + log_pdf(i, mean[j],sigma[j]);
						fprintf(stderr,"%f	%f	%f\n",p[j] ,log_truncated_pdf(x[i], mean[j],sigma[j],0,n), log_pdf(x[i], mean[j],sigma[j]));
						
					}
					//fprintf(stderr,"%f	\n",raw[j]);
				}
				
				
				tmp = 10e-7;
			}
			likelihood = likelihood + tmp;// prob2scaledprob(tmp);
			
			//fprintf(stderr,"len:%d	%f	%f	%f\n",i,x[i],likelihood, tmp );
			
			for(j = 0; j < k;j++){
				if(p[j] != prob2scaledprob(0.0)){
					tmp_p_matrix[i][j]  = tmp_p_matrix[i][j]  - tmp;
					//raw[j]  = raw[j] / tmp;
					e_p[j]  = logsum(e_p[j], tmp_p_matrix[i][j]);// + prob2scaledprob( 1.0x[i]));//raw[j] *x[i];
					
					sum[j] = logsum(sum[j], prob2scaledprob(1.0));// + prob2scaledprob(  x[i]));//raw[j]*x[i];
					e_mean[j] = logsum(e_mean[j],  tmp_p_matrix[i][j]+ prob2scaledprob(x[i]));
				}
			}
			//for(j = 0; j < k;j++){
			
			//}
		}
		tmp = prob2scaledprob(0.0);
		for(j = 0; j < k;j++){
			e_mean[j] = scaledprob2prob( e_mean[j] - sum[j]);
			fprintf(stderr,"NEW mean:%d	%d	%f\n",g,j, e_mean[j]);
			tmp = logsum(tmp, e_p[j]);
			//e_p[j] =  e_p[j]/ sum[j];//1.0/(float)n_observations * e_p[j];
			//fprintf(stderr,"%f	%f\n", 1.0/(float)n_observations * e_p[j], e_p[j]/ sum[j]);
		}
		for(j = 0; j < k;j++){
			e_p[j] =  e_p[j] -  tmp;//1.0/(float)n_observations * e_p[j];
			//fprintf(stderr,"%f	%f\n", 1.0/(float)n_observations * e_p[j], e_p[j]);
		}
		
		
		//estimate new sigma;
		for(i = 0; i < n;i++){
			for(j = 0; j < k;j++){
				if(p[j] != prob2scaledprob(0.0)){
					e_sigma[j] +=  scaledprob2prob( tmp_p_matrix[i][j] )* ( x[i] - e_mean[j]  )  * (x[i] - e_mean[j]  );
				}
			}
		}
		
		for(j = 0; j < k;j++){
			e_sigma[j]  =  sqrtf(    e_sigma[j] / scaledprob2prob(sum[j]));
			//if(e_sigma[j] < MIN_STDEV  ){
			//	e_sigma[j] = MIN_STDEV  ;//FLT_MIN;
			//}
		}
		for(j = 0; j < k;j++){
			sigma[j] = e_sigma[j] ;
			mean[j] = e_mean[j];
			p[j] = e_p[j];
			
		}
		
		if(fabs((likelihood/old_likelihood)-1) < likelihood_threshold ){
			break;
		}
		
		old_likelihood = likelihood;
		
	}
	
	for(j = 0; j < k;j++){
		model->p[j] = p[j];
		model->mean[j] = mean[j];
		model->sigma[j] = sigma[j];
		fprintf(stderr,"%d	%f	%f	%f\n",j,p[j],mean[j],sigma[j]);
	}
	model->likelihood = (double) likelihood ;
	
	
	for(i = 0; i < n;i++){
		free(tmp_p_matrix[i]);// = malloc(sizeof(float)*k);
	}
	free(tmp_p_matrix);// = malloc(sizeof(float*) * n);
	
	free(mean);
	free(sigma);
	free(p);
	free(e_p);
	free(e_mean);
	free(e_sigma);
	
	free(assign);
	
	//free(raw);
	free(sum);
	
	return model;
}
