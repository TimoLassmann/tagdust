/*
 
 Copyright (C) 2013 Timo Lassmann <timolassmann@gmail.com>
 
 This file is part of TagDust.
 
 TagDust is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 TagDust is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with Tagdust.  If not, see <http://www.gnu.org/licenses/>.
 
 */


/*! \file barcode_hmm.c
 \brief Functions build and search with user specified HMMs.
 
 Contains all functions for HMM construction, initialization, training and searching.
 
 \author Timo Lassmann
 \bug No known bugs.
 */

#include <stdio.h>
#include "tagdust2.h"
#include "interface.h"
#include "io.h"
#include "misc.h"
#include "nuc_code.h"
#include <pthread.h>
#include <assert.h>
#include <float.h>
#include <time.h>
#include "barcode_hmm.h"


/** \fn void hmm_controller(struct parameters* param,int (*fp)(struct read_info** ,struct parameters*,FILE* ),int file_num)
 \brief Runs all functions.
 
 Constructs HMM, runs it on sequences.
 
 \param param @ref parameters.
 \param fp Pointer to function used to read sequences (either SAM or fastq).
 \param filenum Number of input files.
 
 */
void hmm_controller(struct parameters* param,int (*fp)(struct read_info** ,struct parameters*,FILE* ),int file_num)
{
	struct read_info** ri = 0;
	
	FILE* outfile;
	int i,j,c,g;
	int numseq;
	int total_read = 0;
	int barcode_length;
	int min_distance; 
	double sum = 0.0;
	
	double max_bar;
	double max_prob;
		
	init_logsum();
	
	double* back = 0;
	int average_length = 0;
	
	back = malloc(sizeof(double)*5);
	
#if DEBUG
	//printf("Debug\n");
	param->num_query = 501;
#else
	//printf("No Debug\n");
	param->num_query = 1000000;
#endif

	//param->num_query = 500000;
	
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
		ri[i]->bar_prob = 0;
		ri[i]->md = 0;
		ri[i]->mapq = -1.0;
		ri[i]->strand = malloc(sizeof(unsigned int)* (LIST_STORE_SIZE+1));
		ri[i]->hits = malloc(sizeof(unsigned int)* (LIST_STORE_SIZE+1));
	}
	file =  io_handler(file, file_num,param);
	
	/*
	 
	 get backgorund nucleotide distribution - from all reads?
	 
	 */
	
	average_length = 0;
	for(i = 0; i < 5;i++){
		back[i] = 1.0; // pseudocount..... -- caused crashes on tiny datasets.... (even N does not occur in input, the M emit states of the HMM would have a N probability becasue of the sequencing error. But the random model was initializes straight from the background...)
	}
	total_read = 0;
	
	while ((numseq = fp(ri, param,file)) != 0){
		for(i = 0; i < numseq;i++){
			average_length += ri[i]->len;
			for(j = 0;j < ri[i]->len;j++){
				//fprintf(stderr,"%d ",(int)ri[i]->seq[j] );
				back[(int)ri[i]->seq[j]] += 1.0f;
			}
		}
		
		total_read += numseq;
#if DEBUG
		if(total_read > 501){
			break;
		}
#else
		if(total_read > 1000000){
			break;
		}
#endif
	}
	if(param->matchstart!= -1 || param->matchend !=-1){
		average_length = (param->matchend - param->matchstart )* total_read;
	}
	
	average_length = average_length / total_read;
	
	param->average_read_length = average_length;
	
	sum = 0.0;
	for(i = 0; i < 5;i++){
	//	fprintf(stderr,"%f\n",(back[i])  );
		sum += back[i];
	}
	
	for(i = 0; i < 5;i++){
		//back[i] = prob2scaledprob(0.25);
		back[i] = prob2scaledprob(back[i]  / sum);
	//	fprintf(stderr,"%f\n",scaledprob2prob(back[i])  );
	}
	//exit(0);
	pclose(file);
	
	
	// Inits model.
	
	struct model_bag* mb = init_model_bag(param, back);
		
	// Let's check the hamming distance between barcodes...
	
	g = 0;
	min_distance = 1000; // hamming distance.
	for(i = 0; i < mb->num_models;i++){
		if(param->read_structure->type[i] == 'B'){
			barcode_length = (int)strlen(param->read_structure->sequence_matrix[i][0]);
			for(j = 0; j < mb->model[i]->num_hmms;j++){
				for(c = j +1;c <  mb->model[i]->num_hmms ;c++){
					numseq = bpm(param->read_structure->sequence_matrix[i][j] ,param->read_structure->sequence_matrix[i][c], barcode_length,barcode_length);
					
					if(numseq < min_distance){
						min_distance = numseq;
						g = 1;
					//	fprintf(stderr,"%s\n%s\t%d\n",param->read_structure->sequence_matrix[i][j],param->read_structure->sequence_matrix[i][c] ,numseq);
					}else if(numseq == min_distance ){
						g++;
					}
				}
			}
		}
	}
	if(min_distance != 1000){
		fprintf(stderr,"Minumum edit distance among barcodes: %d, %d pairs\n", min_distance,g);
	}
	
	
	// Estimates lengths of partial segments by using exact matching + normal distribution. 
	
	file =  io_handler(file, file_num,param);
	
	numseq = fp(ri, param,file);
	mb = estimate_length_distribution_of_partial_segments(mb,ri, param,  numseq);

	pclose(file);
	//fprintf(stderr,"%f	%f\n", model_information_content(mb),   exp(model_information_content(mb) ) / (1 + exp(model_information_content(mb))) * 10 + 5) ;
	/*
	for(i = 1; i < 1000;i++){
		sum = (double)i / 100.0;
		fprintf(stderr,"%f	%f\n",sum,exp(sum ) / (1.0 + exp(sum)) * 10 + 5) ;
	}
	
	exit(0);*/
	//sum = model_information_content(mb) ;
	//param->random_prior = 0.1 * 0.5  + 0.9 *( 1- exp(-1.0 * sum));
	
	//fprintf(stderr,"%f	%f\n",param->random_prior, 1- exp(-1.0 * sum));
	
	//exit(0);
	/*if(!param->confidence_threshold){
		sum = model_information_content(mb) ;
		param->confidence_threshold =  exp(sum ) / (1.0 + exp(sum)) * 20.0;
		fprintf(stderr,"Setting Threshold to %f\n",param->confidence_threshold );
	}*/
	// HMM training - not used in this version..
	
	
	
	/*for(i = 0; i < mb->num_models;i++){
		print_model(mb->model[i]);
	}
	exit(0);
	*/
	/*
	file =  io_handler(file, file_num,param);
	numseq = fp(ri, param,file);
	mb =  run_pHMM(mb,ri,param,numseq,MODE_GET_PROB);
	mb =  run_pHMM(mb,ri,param,numseq,MODE_RUN_RANDOM);
	fprintf(stderr,"read: %d	%d\n",numseq,param->num_query );
	param->confidence_threshold =  set_Q_threshold(mb, ri,  numseq);
	pclose(file);
	*/
	
	file =  io_handler(file, file_num,param);

	if(!param->train ){
	
	}else if( !strcmp( param->train , "full")){
		for(i = 0; i < 10;i++){
			fprintf(stderr,"Iteration %d\n",i);
			mb = set_model_e_to_laplace(mb);
			while ((numseq = fp(ri, param,file)) != 0){
				//	numseq = fp(ri, param,file);
				mb =  run_pHMM(mb,ri,param,numseq,MODE_TRAIN);
			}
			pclose(file);
			file =  io_handler(file, file_num,param);
			//rewind(file);
			for(j = 0; j < mb->num_models;j++){
				print_model(mb->model[j]);
				mb->model[j] = reestimate(mb->model[j], 0);
				print_model(mb->model[j]);
			}
		}
		
	}else if (!strcmp(param->train, "half" )){
		for(i = 0; i < 10;i++){
			fprintf(stderr,"Iteration %d\n",i);
			while ((numseq = fp(ri, param,file)) != 0){
				//	numseq = fp(ri, param,file);
				mb =  run_pHMM(mb,ri,param,numseq,MODE_TRAIN);
			}
			pclose(file);
			file =  io_handler(file, file_num,param);
			//rewind(file);
			for(j = 0; j < mb->num_models;j++){
				print_model(mb->model[i]);
				mb->model[j] = reestimate(mb->model[j], 2);
				print_model(mb->model[i]);
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

	li->num_EXTRACT_SUCCESS = 0;
	li->num_EXTRACT_FAIL_BAR_FINGER_NOT_FOUND = 0;
	li->num_EXTRACT_FAIL_READ_TOO_SHORT = 0;
	li->num_EXTRACT_FAIL_AMBIGIOUS_BARCODE = 0;
	li->num_EXTRACT_FAIL_ARCHITECTURE_MISMATCH = 0;
	li->num_EXTRACT_FAIL_MATCHES_ARTIFACTS = 0;
	li->num_EXTRACT_FAIL_LOW_COMPLEXITY= 0;
		
	total_read = 0;
	while ((numseq = fp(ri, param,file)) != 0){
		
		
		
		//	numseq = fp(ri, param,file);
		mb =  run_pHMM(mb,ri,param,numseq,MODE_GET_LABEL);
		li->total_read += numseq;
		
		max_bar = -FLT_MAX;
		max_prob = -FLT_MAX;
		
		for(i = 0; i < numseq;i++){
			switch ((int) ri[i]->prob) {
					
				case EXTRACT_SUCCESS:
					print_sequence(ri[i],outfile);
					li->num_EXTRACT_SUCCESS++;
					break;
				case EXTRACT_FAIL_BAR_FINGER_NOT_FOUND:
					li->num_EXTRACT_FAIL_BAR_FINGER_NOT_FOUND++;
					break;
				case  EXTRACT_FAIL_READ_TOO_SHORT:
					li->num_EXTRACT_FAIL_READ_TOO_SHORT++;
					break;
				case  EXTRACT_FAIL_AMBIGIOUS_BARCODE:
					li->num_EXTRACT_FAIL_AMBIGIOUS_BARCODE++;
					break;
				case  EXTRACT_FAIL_ARCHITECTURE_MISMATCH:
					li->num_EXTRACT_FAIL_ARCHITECTURE_MISMATCH++;
					break;
				case  EXTRACT_FAIL_MATCHES_ARTIFACTS:
					li->num_EXTRACT_FAIL_MATCHES_ARTIFACTS++;
				case  EXTRACT_FAIL_LOW_COMPLEXITY:
					li->num_EXTRACT_FAIL_LOW_COMPLEXITY++;
					break;
			}
		}
	}
	
	
	if(param->outfile){
		fclose(outfile);
	}
	
	struct tm *ptr;
	int hour;
	char am_or_pm;
	char logfile[100];

	time_t current = time(NULL);
	ptr = localtime(&current);
	hour = ptr->tm_hour;
	if (hour <= 11)
		am_or_pm = 'a';
	else {
		hour -= 12;
		am_or_pm = 'p';
	}
	if (hour == 0){
		hour = 12;
	}
	
	
	sprintf (logfile, "%s_tagdust_log.txt",shorten_pathname(param->infile[file_num]));

	if ((outfile = fopen( logfile, "w")) == NULL){
		fprintf(stderr,"can't open output\n");
		exit(-1);
	}
	
	fprintf(outfile,"%s	%.2d-%.2d-%d	%2d:%.2d%cm\n",param->infile[file_num],ptr->tm_mon + 1,ptr->tm_mday, ptr->tm_year + 1900,hour,ptr->tm_min, am_or_pm );
	fprintf(outfile,"%d	total input reads\n", li->total_read);
	fprintf(outfile,"%f	selected threshold\n", param->confidence_threshold);

	
	
	fprintf(outfile,"%d	successfully extracted\n" ,li->num_EXTRACT_SUCCESS);
	fprintf(outfile,"%d	barcode / UMI not found\n" ,li->num_EXTRACT_FAIL_BAR_FINGER_NOT_FOUND);
	fprintf(outfile,"%d	too short\n" , li->num_EXTRACT_FAIL_READ_TOO_SHORT);
	fprintf(outfile,"%d	ambiguous barcode\n" , li->num_EXTRACT_FAIL_AMBIGIOUS_BARCODE);
	fprintf(outfile,"%d	problems with architecture\n" , li->num_EXTRACT_FAIL_ARCHITECTURE_MISMATCH);
	fprintf(outfile,"%d	matches artifacts\n" , li->num_EXTRACT_FAIL_MATCHES_ARTIFACTS);
	fprintf(outfile,"%d	low complexity\n" , li->num_EXTRACT_FAIL_LOW_COMPLEXITY);
	
	fprintf(outfile,"%0.1f%%	extracted\n",  (float) li->num_EXTRACT_SUCCESS / (float) li->total_read  *100.0f);
	
	fclose(outfile);
	
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
}


		   
		
		   
/** \fn float model_information_content(struct model_bag*mb)
 
 \brief Calculated the information content of a model based on the match states.  
 
 \f[
 
 IC =  - \sum_{i,j} p_{i,j} \log \left( p_{i,j} b_j \right)   
 \f]
 
 
 
 \param mb  @ref model_bag - contains the HMM model.

 */


float model_information_content(struct model_bag*mb)
{
	float IC =0.0;
	int i,j,c,f;
	struct hmm_column* col = 0;
	float test_prob = 0;
	for(i = 0; i < mb->num_models;i++){
		
		for(j = 0;j < mb->model[i]->hmms[0]->num_columns ;j++){
			for(c= 0; c < 5;c++){
				test_prob = 0.0;
				for(f = 0;f < mb->model[i]->num_hmms;f++){
					col = mb->model[i]->hmms[f]->hmm_column[j];
					test_prob += scaledprob2prob( col->m_emit[c])  * 1.0 / (float)mb->model[i]->num_hmms;
					
				}
				IC +=test_prob* log2 (   test_prob  /  scaledprob2prob( mb->model[0]->background_nuc_frequency[c]));
				
			}
		}
	}
	return  IC;
}



/** \fn double set_Q_threshold(struct model_bag* mb, struct read_info** ri, int numseq)
 
 \brief Sets the confidence threshold based on comparison of Q scores obtained freom real and reversed sequences.
 
 \param mb  @ref model_bag - contains the HMM model.
\param ri  @ref read_info - contains the reads.
\param numseq - number of sequences.
 */


double set_Q_threshold(struct model_bag* mb, struct read_info** ri, int numseq)
{
	double estimated_threshold = -1.0;
	double fdr = 0.05;
	int i,j,c,g;
	int extracted;
	int real[5000];
	int fake[5000];
	for(i =0; i < 5000;i++){
		real[i] = 0;
		fake[i] = 0;
	}
	
	for(i =0; i < numseq;i++){
		
		//fprintf(stderr,"%d	\n",i);
		//fprintf(stderr,"%f rand	\n", mb->random_scores[i]);
		//fprintf(stderr,"%f	 mapq\n",  ri[i]->mapq);
		
		if((int) (ri[i]->mapq*100.0) > 4999){
			real[4999]++;
		}else{
			real[(int) (ri[i]->mapq*100.0)]++;
		}
		
		if((int) (mb->random_scores[i]*100) > 4999){
			fake[4999]++;
		}else{
			fake[(int) (mb->random_scores[i]*100.0)]++;
		}
	}
	
	//while(estimated_threshold == -1.0){
		
		c = 0; g = 0;
		
		j = 0;
		extracted = 0;
		
		for(i = 4999; i >= 0;i--){
			c +=real[i];
			g += fake[i];
			j = c+g;
			if((float)g / (float)(c+g)  <= fdr && j >= extracted){
				extracted = j;
				estimated_threshold = (double)i/100.0;
			}
			if(real[i]+fake[i]){
			fprintf(stderr,"%d	%d	%d	%f	%d	%f\n", i,real[i],fake[i],  (float)g / (float)(c+g) , c+g,estimated_threshold);
			}
		}
		if(estimated_threshold == -1.0){
			fprintf(stderr,"warning - cannot find a good separation with FDR: %f\n",fdr);
			fdr += 0.01;
			return 1000;
		}else{
			fprintf(stderr,"Good separation found at FDR: %f\n",fdr);

		}
		/*if(fdr >= 0.05){
			fprintf(stderr,"Setting threshold to 50 - cannot extract reads\n");
			return 50;
		}*/

	//}
	
	//if(estimated_threshold < 10.0){
	//	estimated_threshold = 10.0;
	//}
	

	return estimated_threshold;
}


/** \fn struct model_bag* estimate_length_distribution_of_partial_segments(struct model_bag*mb,struct read_info** ri,struct parameters* param, int numseq)
 
 \brief Initialized length related probabilities for partial segments.
 
 This function matches partial sequences to reads exactly. Based on the exact hits HMM parameters are derived using a normal distribution.
 \param mb  @ref model_bag - contains the HMM model.
 \param ri @ref read_info - contains the sequences.
 \param param @ref parameters - contain the parameters. 
 \param numseq number of sequences.
 
 \bug Seems to set some probabilities to dead states at the end of HMMs... 
 */

struct model_bag* estimate_length_distribution_of_partial_segments(struct model_bag*mb,struct read_info** ri,struct parameters* param, int numseq)
{
	int i,j,c;
	char* test_sequence = 0;
	
	struct model* model = 0;
	struct hmm_column* col = 0;
	
	double mean;
	double stdev;
	
	double sum_prob = 0;
	double s0,s1,s2;
	int len = 0;
	
	float base_error = param->sequencer_error_rate;
	float indel_freq = param->indel_frequency;
	
	//5'
	if(param->read_structure->type[0] == 'P'){
		test_sequence = param->read_structure->sequence_matrix[0][0];
		len = (int) strlen(test_sequence);
		
		for(i = 0; i < len;i++){
			test_sequence[i] = nuc_code[(int) test_sequence[i]];
		}
		//for(c = 0;c < len;c++){
		//	fprintf(stderr,"%c",*(test_sequence +c) + 65);
		//
		//}
		//fprintf(stderr,"\n");
		
		mean = 0;
		s0 = 0;
		s1 = 0;
		s2 = 0;
		
		
		for(i = 0; i < numseq;i++){
			for(j = 0;j <= len ;j++){
				for(c = 0;c < len-j;c++){
					if(ri[i]->seq[c] != test_sequence[j +c]){
						break;
					}
				}
				if(c == len-j && c > 3 ){
					s0++;
					s1 += len -j;
					s2 += (len-j) * (len-j);
					break;
				}
				
			}
		}
		mean = s1 / s0;
		stdev = sqrt(  (s0 * s2 - pow(s1,2.0))   /  (  s0 *(s0-1.0) )) ;
		if(stdev < 1){
			stdev = 1;
		}
		fprintf(stderr,"5: %f %f	%f\n", mean,  stdev,s0);

		if(mean <= 1){
			fprintf(stderr,"WARNING: 5' partial segment seems not to be present in the data (length < 1).\n");
		}

		sum_prob = 0;
		
		for(i = 0; i <  len;i++){
		//	fprintf(stderr,"%f ",sum_prob );
			sum_prob +=gaussian_pdf(i , mean ,stdev);
		}
		
		//fprintf(stderr,"\n");
		
		//Init model ....
		model = mb->model[0];
		model->skip = prob2scaledprob(  gaussian_pdf(0 , mean ,stdev) / sum_prob    );
		
		//fprintf(stderr,"5SKIP:%f\t%f\n", model->skip, scaledprob2prob(model->skip));
		s1 = prob2scaledprob(0.0);
		//fprintf(stderr,"SUMPROB:%f\n", s1);
		s1 = logsum(s1, model->skip);
		//fprintf(stderr,"SUMPROB:%f\n", s1);
		//if(rs->type[key] == 'P'){// Partial - can skip and exit at every M / I state....
		len = model->hmms[0]->num_columns;
		for(i = 0 ; i < model->num_hmms;i++){
			
			//model->M_to_silent[i] = prob2scaledprob(1.0);
			
			for(j = 0; j < len;j++){
				col = model->hmms[i]->hmm_column[j];
				model->silent_to_M[i][j]  = prob2scaledprob(1.0 / (float) model->num_hmms) + prob2scaledprob(   gaussian_pdf(len-j , mean ,stdev) / sum_prob);
				
				s1 = logsum(s1, model->silent_to_M[i][j]);
				
			}
			
			model->hmms[i] = set_hmm_transition_parameters(model->hmms[i],len, base_error, indel_freq, -1.0, -1.0);
		}
		
		model->skip = model->skip - s1;
		s2 = model->skip;
		for(i = 0 ; i < model->num_hmms;i++){
			
			//model->M_to_silent[i] = prob2scaledprob(1.0);
			
			for(j = 0; j < len;j++){
				model->silent_to_M[i][j]  = model->silent_to_M[i][j]  - s1;
			}
		}
		
		
			
		/*fprintf(stderr,"%f	skip\n", scaledprob2prob(model->skip));
		for(j = 0; j < len;j++){
			fprintf(stderr,"%f	%d len \n", scaledprob2prob(model->silent_to_M[0][j] ),j);
		}*/
			
		//}
		

	
		
//#endif
		
	}
	
	for(c = 1; c < mb->num_models-1;c++){
		if(param->read_structure->type[c] == 'P'){
			
			model = mb->model[c];
			len = model->hmms[0]->num_columns;
			for(i = 0 ; i < model->num_hmms;i++){
				
				j = 0;
				model->hmms[i] = set_hmm_transition_parameters(model->hmms[i],len, base_error, indel_freq, 0.1, -1.0); // 0.1 is exit parameter for partial internal sequence (I don't know how to set this.....)
			}
		}
	}
	
	
	//3'
	if(param->read_structure->type[mb->num_models-1] == 'P'){
		test_sequence = param->read_structure->sequence_matrix[ mb->num_models-1][0];
		len = (int) strlen(test_sequence);
		for(i = 0; i < len;i++){
			test_sequence[i] = nuc_code[(int) test_sequence[i]];
		}
		
		mean = 0;
		s0 = 0;
		s1 = 0;
		s2 = 0;
		
		for(i = 0; i < numseq;i++){
			for(j = 0;j <= len ;j++){
				
				for(c = 0;c < len-j;c++){
					if(ri[i]->seq[ri[i]->len - (len-j -c)] != test_sequence[c]){
						break;
					}
				}
				if(c == len-j  && c > 3){
					s0++;
					s1 += len -j;
					s2 += (len-j) * (len-j);
					break;
				}
							}
		}
		mean = s1 / s0;
		stdev = sqrt(  (s0 * s2 - pow(s1,2.0))   /  (  s0 *(s0-1.0) )) ;
		if(stdev < 1){
			stdev = 1;
		}
		
		fprintf(stderr,"3: %f %f\n", mean,  stdev);
		if(mean <= 1){
			fprintf(stderr,"WARNING: 3' partial segment seems not to be present in the data (length < 1).\n");
		}
		
		sum_prob = 0;
		
		for(i = 0; i <  len;i++){
			sum_prob +=gaussian_pdf(i , mean ,stdev);
		}

		//Init model ....
		model = mb->model[mb->num_models-1];
		model->skip = prob2scaledprob(  gaussian_pdf(0 , mean ,stdev) / sum_prob    );
		
		
		s1 = model->skip;
		
		//fprintf(stderr,"3SKIP:%f\t%f\n", model->skip, scaledprob2prob(model->skip));
		//if(rs->type[key] == 'P'){// Partial - can skip and exit at every M / I state....
		//len = model->hmms[mb->num_models-1]->num_columns;
		for(i = 0 ; i < model->num_hmms;i++){
			
			//model->M_to_silent[i] = prob2scaledprob(1.0);
			model->silent_to_M[i][0] = prob2scaledprob(1.0 / (float) model->num_hmms) + prob2scaledprob(1.0 -  gaussian_pdf(0 , mean ,stdev) / sum_prob );
			//fprintf(stderr,"move:%d:%f\t%f\n", i,model->silent_to_M[i][0], scaledprob2prob(model->silent_to_M[i][0]));
			
			model->hmms[i] = set_hmm_transition_parameters(model->hmms[i],len, base_error, indel_freq, mean, stdev);
		}
		
	}
	return mb;
}



/** \fnstruct hmm* set_hmm_transition_parameters(struct hmm* hmm, int len,double base_error, double indel_freq,  double mean, double stdev)
 
 \brief Initialized transition probabilities of HMM segments.
 
 This function sets the transition probabilities within HMMs.
 
 If mean and stdev = -1.0 MSKIP transitions are set to 0 - the sequence has to match the entire HMM.
 If mean is > -1.0 and stdev = -1.0 MSKIP is set to a constant throughout the model. 
 If mean and stdev are not -1.0 MSKIP is set according to a gaussian distribiution - i.e. leaving the model towards the end is more probable than at the beginning...
 
 
 
 
 \param hmm  @ref hmm - contains the HMM model.
 \param len - length of the HMM.
 \param base_error - seuqencer error rate.
 \param indel_freq - InDel frequency.
 
 \param mean - mean length of modelled 3' sequence
 \param stdev - standard deviation of read length.
 
 */




struct hmm* set_hmm_transition_parameters(struct hmm* hmm, int len,double base_error, double indel_freq,  double mean, double stdev)
{
	//cases:
	//1) no MSKIP (apart from last columns
	//2) MSKIP determined by mean & stdev.
	//3 constant MSKIP...
	
	
	int i;
	struct hmm_column* col = 0;
	
	
	double sum_prob = 0.0;
	if(mean > 0.0 && stdev > 0.0){
		for(i = 0; i <=  len;i++){
			sum_prob +=gaussian_pdf(i , mean ,stdev);
		}
	}
	
	
	if(len == 1){
		//single state - only silent to / from M everything else set to zero....
		col = hmm->hmm_column[0];
		col->transition[MM] = prob2scaledprob(0.0f);
		col->transition[MI] = prob2scaledprob(0.0f);
		col->transition[MD] = prob2scaledprob(0.0f);
		col->transition[MSKIP] = prob2scaledprob(1.0f);
		
		col->transition[II] = prob2scaledprob(0.0f);
		col->transition[IM] = prob2scaledprob(0.0f);
		col->transition[ISKIP] = prob2scaledprob(0.0f);
		
		col->transition[DD] = prob2scaledprob(0.0);
		col->transition[DM] = prob2scaledprob(0.0);
	}else if(len == 2){
		
		//first column
		col = hmm->hmm_column[0];
		
		if(mean == -1.0 && stdev == -1.0){
			col->transition[MSKIP] = prob2scaledprob(0.0);
		}else if(mean > -1.0 && stdev == -1.0){
			col->transition[MSKIP] = prob2scaledprob(mean / (float)(len-1));
		}else{
			col->transition[MSKIP] = prob2scaledprob(   gaussian_pdf(1 , mean ,stdev) / sum_prob);
		}
		
		col->transition[MM] = prob2scaledprob( 1.0 - base_error * indel_freq ) + prob2scaledprob (1.0 - scaledprob2prob(col->transition[MSKIP]  ));
		col->transition[MI] = prob2scaledprob(base_error * indel_freq)+ prob2scaledprob (1.0 - scaledprob2prob(col->transition[MSKIP]  ));
		
		col->transition[MD] = prob2scaledprob(base_error * indel_freq*0.0)+ prob2scaledprob (1.0 - scaledprob2prob(col->transition[MSKIP]  ));
		
		
		col->transition[II] = prob2scaledprob(1.0 - 0.999);
		col->transition[IM] = prob2scaledprob(0.999);
		col->transition[ISKIP] = prob2scaledprob(0.0f);
		
		col->transition[DD] = prob2scaledprob(0.0);
		col->transition[DM] = prob2scaledprob(0.0);
		
		//second column
		col = hmm->hmm_column[1];
		col->transition[MM] = prob2scaledprob(0.0f );// 1.0 - base_error * indel_freq);
		col->transition[MI] = prob2scaledprob(0.0f );//(base_error * indel_freq) +  prob2scaledprob(0.5);
		col->transition[MD] = prob2scaledprob(0.0f );//(base_error * indel_freq) +  prob2scaledprob(0.5);
		col->transition[MSKIP] = prob2scaledprob(1.0);
		
		col->transition[II] = prob2scaledprob(0.00);
		col->transition[IM] = prob2scaledprob(0.0);
		col->transition[ISKIP] = prob2scaledprob(0.0f);
		
		col->transition[DD] = prob2scaledprob(0.0f );//(1.0 - 0.999);
		col->transition[DM] = prob2scaledprob(0.0f );//0.999);
		
	}else{
		//first column....
		col = hmm->hmm_column[0];
		
		if(mean == -1.0 && stdev == -1.0){
			col->transition[MSKIP] = prob2scaledprob(0.0);
		}else if(mean > -1.0 && stdev == -1.0){
			col->transition[MSKIP] = prob2scaledprob(mean / (float)(len-1));
		}else{
			col->transition[MSKIP] = prob2scaledprob(   gaussian_pdf(1 , mean ,stdev) / sum_prob);
		}
		
		col->transition[MM] = prob2scaledprob( 1.0 - base_error * indel_freq)+ prob2scaledprob (1.0 - scaledprob2prob(col->transition[MSKIP]  ));
		col->transition[MI] = prob2scaledprob(base_error * indel_freq*0.5)+ prob2scaledprob (1.0 - scaledprob2prob(col->transition[MSKIP]  ));
		
		col->transition[MD] = prob2scaledprob(base_error * indel_freq*0.5)+ prob2scaledprob (1.0 - scaledprob2prob(col->transition[MSKIP]  ));
		
		
		
		col->transition[II] = prob2scaledprob(1.0 - 0.999);
		col->transition[IM] = prob2scaledprob(0.999);
		col->transition[ISKIP] = prob2scaledprob(0.0f);
		
		col->transition[DD] = prob2scaledprob(0.0);
		col->transition[DM] = prob2scaledprob(0.0);
		
		//middle columns...
		for(i = 1; i < len-2;i++){
			col = hmm->hmm_column[i];
			
			if(mean == -1.0 && stdev == -1.0){
				col->transition[MSKIP] = prob2scaledprob(0.0);
			}else if(mean > -1.0 && stdev == -1.0){
				col->transition[MSKIP] = prob2scaledprob(mean / (float)(len-1));
			}else{
				col->transition[MSKIP] = prob2scaledprob(   gaussian_pdf(i +1.0 , mean ,stdev) / sum_prob);
			}
			
			col->transition[MM] = prob2scaledprob( 1.0 - base_error * indel_freq)+ prob2scaledprob (1.0 - scaledprob2prob(col->transition[MSKIP]  ));
			col->transition[MI] = prob2scaledprob(base_error * indel_freq*0.5)+ prob2scaledprob (1.0 - scaledprob2prob(col->transition[MSKIP]  ));
			
			col->transition[MD] = prob2scaledprob(base_error * indel_freq*0.5)+ prob2scaledprob (1.0 - scaledprob2prob(col->transition[MSKIP]  ));
			
			
						
			col->transition[II] = prob2scaledprob(1.0 - 0.999);
			col->transition[IM] = prob2scaledprob(0.999);
			col->transition[ISKIP] = prob2scaledprob(0.0f);
			
			
			col->transition[DD] = prob2scaledprob(1.0 - 0.999);
			col->transition[DM] = prob2scaledprob(0.999);
		}
		
		//second last...
		col = hmm->hmm_column[len -2];
		
		
		
		if(mean == -1.0 && stdev == -1.0){
			col->transition[MSKIP] = prob2scaledprob(0.0);
		}else if(mean > -1.0 && stdev == -1.0){
			col->transition[MSKIP] = prob2scaledprob(mean / (float)(len-1));
		}else{
			col->transition[MSKIP] = prob2scaledprob(   gaussian_pdf( len-1.0 , mean ,stdev) / sum_prob);
		}
		
		col->transition[MM] = prob2scaledprob( 1.0 - base_error * indel_freq)+ prob2scaledprob (1.0 - scaledprob2prob(col->transition[MSKIP]  ));
		col->transition[MI] = prob2scaledprob(base_error * indel_freq*1.0)+ prob2scaledprob (1.0 - scaledprob2prob(col->transition[MSKIP]  ));
		
		col->transition[MD] = prob2scaledprob(base_error * indel_freq*0.0)+ prob2scaledprob (1.0 - scaledprob2prob(col->transition[MSKIP]  ));
		
				
		col->transition[II] = prob2scaledprob(1.0 - 0.999);
		col->transition[IM] = prob2scaledprob(0.999);
		col->transition[ISKIP] = prob2scaledprob(0.0f);
		
		col->transition[DD] = prob2scaledprob(0.0);
		col->transition[DM] = prob2scaledprob(1.0);
		//col->transition[DD] = prob2scaledprob(0.0);
		//col->transition[DM] = prob2scaledprob(0.0);
		
		col = hmm->hmm_column[len -1];
		
		col->transition[MM] = prob2scaledprob(0.0f );// 1.0 - base_error * indel_freq);
		col->transition[MI] = prob2scaledprob(0.0f );//(base_error * indel_freq) +  prob2scaledprob(0.5);
		col->transition[MD] = prob2scaledprob(0.0f );//(base_error * indel_freq) +  prob2scaledprob(0.5);
		col->transition[MSKIP] = prob2scaledprob(1.0);
		
		col->transition[II] = prob2scaledprob(0.00);
		col->transition[IM] = prob2scaledprob(0.0);
		col->transition[ISKIP] = prob2scaledprob(0.0f);
		
		col->transition[DD] = prob2scaledprob(0.0f );//(1.0 - 0.999);
		col->transition[DM] = prob2scaledprob(0.0f );//0.999);
	}
	return hmm;
}


/** \fn struct model_bag* run_pHMM(struct model_bag* mb,struct read_info** ri,struct parameters* param,int numseq, int mode)
 \brief Starts threads to run HMM functions on subsets of sequences.
 Depending on the specifies mode, this function startes threads with different function pointers. 
 
 
 \param ri @ref read_info.
 \param param @ref parameters.
 \param numseq Number of sequences. 
 \param mode Determines which function to run.
 
 */
struct model_bag* run_pHMM(struct model_bag* mb,struct read_info** ri,struct parameters* param,int numseq, int mode)
{
	struct thread_data* thread_data = 0;
	
	
	thread_data = malloc(sizeof(struct thread_data)* param->num_threads);
	pthread_t threads[param->num_threads];
	pthread_attr_t attr;
	int i,t;
	
	int interval = 0;
	int rc;
	
	
	struct fasta* reference_fasta = 0;
	
	if(param->reference_fasta){
		reference_fasta = get_fasta(reference_fasta,param->reference_fasta);
	}

	
	interval =  (int)((double)numseq /(double)param->num_threads);
	
	for(t = 0;t < param->num_threads ;t++) {
		thread_data[t].fasta = reference_fasta;
		thread_data[t].ri = ri;
		thread_data[t].mb = copy_model_bag(mb);
		thread_data[t].start = t*interval;
		thread_data[t].end = t*interval + interval;
		thread_data[t].param = param;
	}
	thread_data[param->num_threads-1].end = numseq;
	//unsigned int seed = (unsigned int) (time(NULL) * ( 42));
	
	rc = pthread_attr_init(&attr);
	if(rc){
		fprintf(stderr,"ERROR; return code from pthread_attr_init() is %d\n", rc);
		exit(-1);
	}
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	
	for(t = 0;t < param->num_threads;t++) {
		switch (mode) {
			case MODE_GET_PROB:
				rc = pthread_create(&threads[t], &attr, do_probability_estimation, (void *) &thread_data[t]);
				break;
			case MODE_GET_LABEL:
				rc = pthread_create(&threads[t], &attr, do_label_thread, (void *) &thread_data[t]);
				break;
			case MODE_TRAIN:
				rc = pthread_create(&threads[t], &attr, do_baum_welch_thread, (void *) &thread_data[t]);
				break;
				
			case MODE_RUN_RANDOM:
				rc = pthread_create(&threads[t], &attr, do_run_random_sequences, (void *) &thread_data[t]);
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
	if(mode == MODE_RUN_RANDOM){
		for (t = 0;t < param->num_threads;t++){
			for(i = thread_data[t].start; i <  thread_data[t].end;i++){
				mb->random_scores[i] = thread_data[t].mb->random_scores[i];
			}
		}
	}
	
	for(t = 0;t < param->num_threads;t++) {
		free_model_bag(thread_data[t].mb);
	}
	
	free(thread_data);
	
	if(reference_fasta){
		free_fasta(reference_fasta);
	}
	
	return mb;
}


/** \fn void* do_probability_estimation(void *threadarg)
 \brief Estimates the probability of a sequence but without labelling it.
 
 This function calculates the equivalent of a mapping Q value to reflext the confidence we have in a match to the user specified read architecture. 
 

 \f[
 
 P_{wrong} = 1 - \frac{P(x|M_i) }{P(x|M) + P(x | R)}
 \f]

 where \f$P(x|M_i)\f$ is the probability of selecing one particular barcode, \f$P(x|M)\f$ is the total probability of the model and \f$P(x|R)\f$ is the probability of the random model. The resulting probability is converted into a phred scaled quality value:
 
 \f[
 Q = -10 * log_{10}( P_{wrong})
 \f]
 
 \param threadarg  A @ref thread_data object used to pass data / parameters to function.
 */
void* do_probability_estimation(void *threadarg)
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
	
	float pbest = 0;
	float Q = 0;
	if(matchstart != -1 || matchend != -1){
		for(i = start; i < end;i++){
			tmp = matchend - matchstart ;
			mb = backward(mb, ri[i]->seq + matchstart , tmp);
			mb = forward_max_posterior_decoding(mb, ri[i] , ri[i]->seq+matchstart ,tmp );
			pbest = 1.0 - scaledprob2prob(  (ri[i]->bar_prob + mb->f_score) - logsum(mb->f_score,  mb->r_score));
			if(!pbest){
				Q = 40.0;
			}else if(pbest == 1.0){
				Q = 0.0;
				
			}else{
				Q = -10.0 * log10(pbest) ;
			}
			ri[i]->mapq = Q;
		}
	}else{
		for(i = start; i < end;i++){
			
			//fprintf(stderr,"%d	%d	%d	%d\n",i,ri[i]->len,start,end);
			mb = backward(mb, ri[i]->seq ,ri[i]->len);			
			mb = forward_max_posterior_decoding(mb, ri[i], ri[i]->seq ,ri[i]->len);
			
			pbest = 1.0 - scaledprob2prob(  (ri[i]->bar_prob + mb->f_score) - logsum(mb->f_score,  mb->r_score));
			
			if(pbest < 0.0){
				exit(-1);
			}
			if(!pbest){
				Q = 40.0;
			}else if(pbest == 1.0){
				Q = 0.0;
	
			}else{
				Q = -10.0 * log10(pbest) ;
			}
			//fprintf(stderr,"%d	score: 1 - %f+%f  / %f = %f    => Q:%f\n", i, ri[i]->bar_prob , mb->f_score , logsum(mb->f_score,  mb->r_score),pbest,Q);
			
			ri[i]->mapq = Q;
		}
	}
		
	pthread_exit((void *) 0);
}



/** \fn void* do_label_thread(void *threadarg)
 \brief Estimates the probability of a sequence and labels it. 

 This function rund the forward and backward algorithm on each sequence, obtails posterior label probabilities for each residue and calculates the optimal sequence labeling based on the posteriors. In addition the quality value to reflext the confidence we have in a matchis calculated:
 
 \f[
 
 P_{wrong} = 1 - \frac{P(x|M_i) }{P(x|M) + P(x | R)}
 \f]
 
 where \f$P(x|M_i)\f$ is the probability of selecing one particular barcode, \f$P(x|M)\f$ is the total probability of the model and \f$P(x|R)\f$ is the probability of the random model. The resulting probability is converted into a phred scaled quality value:
 
 \f[
 Q = -10 * log_{10}( P_{wrong})
 \f]
 
 After labelling each sequence is compared against user defined artifact sequences and low comlexity sequences are filtered out. 
 
 \param threadarg  A @ref thread_data object used to pass data / parameters to function.
 */
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
	float pbest,Q;
	//if(ri[0]->mapq == -1){ // skip sequence scoring if Q values were already calculated....
	//add new random model if barcodes are present..,.
	int bar = -1;
	int r = -1;
	tmp = 100000000;
	//switch read and next shortest segment. ...
	struct model* tmp_model = 0;
	for(i = 0; i < mb->num_models;i++){
		if(data->param->read_structure->type[i] == 'R'){
			r = i;
		}
		if(data->param->read_structure->type[i] == 'B'){
			//if(mb->model[i]->hmms[0]->num_columns < tmp){
				//tmp =mb->model[i]->hmms[0]->num_columns;
				bar= i;
			//}
		}
		//mb->model[i]->hmms[0]->num_columns
	}
	if(bar != -1){
		tmp_model = mb->model[r];
		mb->model[r] = mb->model[bar];
		mb->model[bar] = tmp_model;
		
		if(matchstart != -1 || matchend != -1){
			for(i = start; i < end;i++){
				tmp = matchend - matchstart ;
				mb = backward(mb, ri[i]->seq + matchstart , tmp);
				ri[i]->mapq = mb->b_score;
				/*
				mb = forward_max_posterior_decoding(mb, ri[i] , ri[i]->seq+matchstart ,tmp );
				
				pbest = logsum(mb->f_score, mb->r_score);
				pbest = 1.0 - scaledprob2prob(  (ri[i]->bar_prob + mb->f_score) - pbest);
				if(!pbest){
					Q = 40.0;
				}else if(pbest == 1.0){
					Q = 0.0;
					
				}else{
					Q = -10.0 * log10(pbest) ;
				}
				ri[i]->mapq = Q;*/
				
			}
		}else{
			for(i = start; i < end;i++){
				mb = backward(mb, ri[i]->seq ,ri[i]->len);
				
				ri[i]->mapq = mb->b_score;
				
				/*mb = forward_max_posterior_decoding(mb, ri[i], ri[i]->seq ,ri[i]->len);
				pbest = logsum(mb->f_score, mb->r_score);
				
				pbest = 1.0 - scaledprob2prob(  (ri[i]->bar_prob + mb->f_score ) - pbest);
				if(!pbest){
					Q = 40.0;
				}else if(pbest == 1.0){
					Q = 0.0;
					
				}else{
					Q = -10.0 * log10(pbest) ;
				}
				ri[i]->mapq = Q;*/
				
				
			}
		}
		tmp_model = mb->model[r];
		mb->model[r] = mb->model[bar];
		mb->model[bar] = tmp_model;
	}else{
		for(i = start; i < end;i++){
			ri[i]->mapq = prob2scaledprob(0.0);
		}
	}
	
	
	
	if(matchstart != -1 || matchend != -1){
		for(i = start; i < end;i++){
			tmp = matchend - matchstart ;
			mb = backward(mb, ri[i]->seq + matchstart , tmp);
			mb = forward_max_posterior_decoding(mb, ri[i] , ri[i]->seq+matchstart ,tmp );
		
			
			
			pbest =ri[i]->mapq;
			
			pbest = logsum(pbest, mb->f_score);
			pbest = logsum(pbest, mb->r_score);
			
			pbest = 1.0 - scaledprob2prob(  (ri[i]->bar_prob + mb->f_score) - pbest);
			if(!pbest){
				Q = 40.0;
			}else if(pbest == 1.0){
				Q = 0.0;
				
			}else{
				Q = -10.0 * log10(pbest) ;
			}
			ri[i]->mapq = Q;
		}
	}else{
		for(i = start; i < end;i++){
			mb = backward(mb, ri[i]->seq ,ri[i]->len);
			mb = forward_max_posterior_decoding(mb, ri[i], ri[i]->seq ,ri[i]->len);
			
			
			pbest =ri[i]->mapq;
			
			pbest = logsum(pbest, mb->f_score);
			pbest = logsum(pbest, mb->r_score);
			
			//fprintf(stderr,"%f	%f	%f\n",mb->f_score+ri[i]->bar_prob,ri[i]->mapq,mb->r_score);
			
			pbest = 1.0 - scaledprob2prob(  (ri[i]->bar_prob + mb->f_score ) - pbest);
			if(!pbest){
				Q = 40.0;
			}else if(pbest == 1.0){
				Q = 0.0;
				
			}else{
				Q = -10.0 * log10(pbest) ;
			}
			fprintf(stderr,"%f	%f	%f:Q: %f\n",mb->f_score+ri[i]->bar_prob,ri[i]->mapq,mb->r_score,Q);
			ri[i]->mapq = Q;
		}
	}
	
	for(i = start; i < end;i++){
		ri[i]->bar_prob = 100;
		ri[i] = extract_reads(mb,data->param,ri[i]);
	}
	
	if(data->param->reference_fasta){
		ri = match_to_reference(data);
	}
	
	if(data->param->dust){
		ri = dust_sequences(data);
	}
	/*
#if DEBUG	
	for(i = start; i < end;i++){
		
		char alpha[5] = "ACGTN";
		
		switch ((int) ri[i]->prob) {
				
			case EXTRACT_SUCCESS:
				fprintf(stderr,"Success!!!\n");
				break;
			case EXTRACT_FAIL_BAR_FINGER_NOT_FOUND:
				fprintf(stderr,"FAIL: barcode not found!!!\n");
				break;
			case  EXTRACT_FAIL_READ_TOO_SHORT:
				fprintf(stderr,"FAIL: read too short !!!\n");
				break;
			case  EXTRACT_FAIL_AMBIGIOUS_BARCODE:
				fprintf(stderr,"FAIL: ambigious barcode  !!!\n");
				break;
			case  EXTRACT_FAIL_ARCHITECTURE_MISMATCH:
				fprintf(stderr,"FAIL: architecture does not match  !!!\n");
				break;
			case  EXTRACT_FAIL_MATCHES_ARTIFACTS:
				fprintf(stderr,"FAIL: matches user supplied artifact  !!!\n");
				break;
			case  EXTRACT_FAIL_LOW_COMPLEXITY:
				fprintf(stderr,"FAIL: matches user supplied artifact  !!!\n");
				break;
		}
		
		fprintf(stderr,"%s\n",ri[i]->name);
		int j;
		for(j = 0; j < ri[i]->len;j++){
			fprintf(stderr,"%c", alpha[(int) ri[i]->seq[j]]);
		}
		fprintf(stderr,"\n+\n%s\n" ,ri[i]->qual);
	}
	
#endif
	*/
	pthread_exit((void *) 0);
}



/** \fn struct read_info** dust_sequences(struct thread_data *data)
 \brief Removel low complexity sequences.
 
 Runs the DUST algorithm on the first 64 nucleotides of a read. 
 \param threadarg  A @ref thread_data object used to pass data / parameters to function.
 */

 struct read_info** dust_sequences(struct thread_data *data)
{
	struct read_info** ri = data->ri;
	const int start = data->start;
	const int end = data->end;
	
	int dust_cut = data->param->dust ;
	
	int i,j;
	int key = 0;
	double triplet[64];
	double s = 0.0;
	int len;
	for(j = 0;j < 64;j++){
		triplet[j] = 0.0;
	}
//	

	
	for(i = start; i < end;i++){
		key = ((ri[i]->seq[0] & 0x3 ) << 2 )|  (ri[i]->seq[1] & 0x3 );
		
		len = ri[i]->len;
		if(len > 64){
			len = 64;
		}
		
		for(j = 2;j < len;j++){
			key = key << 2 | (ri[i]->seq[j] & 0x3 );
			triplet[key & 0x3F]++;
		}
		s = 0.0;
		for(j = 0;j < 64;j++){
			
			s+= triplet[j] * (triplet[j] -1.0) / 2.0;
			triplet[j] = 0.0;
		}
		s = s / (double)(len-3) *10.0; //should be number of triplets -2 : i.e. len -2 triples -1 = len -3;
	
		if(s > dust_cut){
#if DEBUG
			char alphabet[] = "ACGTN";
			for(j = 0;j < len;j++){
				fprintf(stderr,"%c",alphabet[(int)ri[i]->seq[j]]);
			}
			fprintf(stderr,"\tLOW:%f\n",s);
#endif
			ri[i]->prob = EXTRACT_FAIL_LOW_COMPLEXITY;
		}
	}
	return ri;
}


/** \fn struct read_info** match_to_reference(struct thread_data *data)
 \brief Matches reads to artifactual sequences.
 
Exhaustively compares reads to a fasta file of known artifact sequences. Uses a SSE version of the Myers bit-parallel dynamic programming algorithm.
 
 \param threadarg  A @ref thread_data object used to pass data / parameters to function.
 */

 struct read_info** match_to_reference(struct thread_data *data)
{
	struct read_info** ri = data->ri;
	struct fasta* reference = data->fasta ;
	const int start = data->start;
	const int end = data->end;
	
	int error_cut = data->param->filter_error ;
	
	int i,j,c;
	int test = 1;
	int reverse = 0;
	unsigned char* seq[4];
	
	int _MM_ALIGN16 lengths[4];
	int _MM_ALIGN16 errors[4];
	
	for(i = start; i <= end-4;i+=4){
		test = 1;
		reverse = 0;
		for(c = 0;c < 4;c++){
			errors[c] = 100000;
			
		}
		for(j =0; j < reference->numseq;j++){
			seq[0] = (unsigned char* ) ri[i]->seq;
			seq[1] = (unsigned char* ) ri[i+1]->seq;
			seq[2] = (unsigned char* ) ri[i+2]->seq;
			seq[3] = (unsigned char* ) ri[i+3]->seq;
			lengths[0] =  ri[i]->len;
			lengths[1] =  ri[i+1]->len;
			lengths[2] =  ri[i+2]->len;
			lengths[3] =  ri[i+3]->len;
			validate_bpm_sse(seq,lengths,reference->string +  reference->s_index[j],reference->s_index[j+1] - reference->s_index[j],4);
			for(c = 0;c < 4;c++){
				if(lengths[c] < errors[c]){
					errors[c] = lengths[c];
				}
			}
			
			
			seq[0] = reverse_complement((unsigned char* ) ri[i]->seq,ri[i]->len);
			seq[1] = reverse_complement((unsigned char* ) ri[i+1]->seq,ri[i+1]->len);
			seq[2] = reverse_complement((unsigned char* ) ri[i+2]->seq,ri[i+2]->len);
			seq[3] = reverse_complement((unsigned char* ) ri[i+3]->seq,ri[i+3]->len);
			lengths[0] =  ri[i]->len;
			lengths[1] =  ri[i+1]->len;
			lengths[2] =  ri[i+2]->len;
			lengths[3] =  ri[i+3]->len;
			validate_bpm_sse(seq,lengths,reference->string +  reference->s_index[j],reference->s_index[j+1] - reference->s_index[j],4);
			for(c = 0;c < 4;c++){
				if(lengths[c] < errors[c]){
					errors[c] = lengths[c];
				}
			}
			
			seq[0] = reverse_complement((unsigned char* ) ri[i]->seq,ri[i]->len);
			seq[1] = reverse_complement((unsigned char* ) ri[i+1]->seq,ri[i+1]->len);
			seq[2] = reverse_complement((unsigned char* ) ri[i+2]->seq,ri[i+2]->len);
			seq[3] = reverse_complement((unsigned char* ) ri[i+3]->seq,ri[i+3]->len);
		}
		for(c = 0;c < 4;c++){
			if(errors[c] <= error_cut){
				if(ri[i+c]->prob == EXTRACT_SUCCESS){
					ri[i+c]->prob  =  EXTRACT_FAIL_MATCHES_ARTIFACTS;
				}
			}
		}
	}
	
	while(i < end){
		//fprintf(stderr,"Looking at %d	%d	%d\n",i,start,end);
		test = 1;
		reverse = 0;
		for(j =0; j < reference->numseq;j++){
			c = bpm_check_error(reference->string +  reference->s_index[j], (unsigned char* )ri[i]->seq,reference->s_index[j+1] - reference->s_index[j] , ri[i]->len);
			if(c <= error_cut){
				test = 0;
				break;
			}
			ri[i]->seq = (char* )reverse_complement((unsigned char* ) ri[i]->seq,   ri[i]->len);
			c = bpm_check_error(reference->string +  reference->s_index[j], (unsigned char* )ri[i]->seq,reference->s_index[j+1] - reference->s_index[j] , ri[i]->len);
			reverse = 1;
			if(c <= error_cut){
				ri[i]->seq =(char* ) reverse_complement( (unsigned char* )ri[i]->seq,   ri[i]->len);
				test = 0;
				break;
			}
			ri[i]->seq = (char* )reverse_complement((unsigned char* ) ri[i]->seq,   ri[i]->len);
		}
		if(!test){
			if(ri[i]->prob == EXTRACT_SUCCESS){
				ri[i]->prob  =  EXTRACT_FAIL_MATCHES_ARTIFACTS;
			}
		}
		i++;
	}

	
	return ri;
}


/** \fn struct read_info* emit_random_sequence(struct model_bag* mb, struct read_info* ri,int average_length,unsigned int* seed )
 \brief Emits sequences from random HMM model. 
 
 \param mb  The model @ref model_bag . 
 \param ri  @ref read_info - emitted sequences are written here. 
 \param average_length Average length of the sequences. 
 \param seed Seed used for randomization.
 
 */


struct read_info* emit_random_sequence(struct model_bag* mb, struct read_info* ri,int average_length,unsigned int* seed )
{
	
	int current_length = 0;
	int allocated_length = 100;
	double r = (float)rand_r(seed)/(float)RAND_MAX;
	double sum = prob2scaledprob(0.0f);
	//char alpha[] = "ACGTN";
	int nuc;
	
	free(ri->seq);
	free(ri->name);
	free(ri->qual);
	free(ri->labels);
	ri->seq = 0;
	ri->name = 0;
	ri->qual = 0;
	ri->labels = 0;
	ri->len = 0;
	ri->md = 0;
	ri->prob = 0;
	//ri[i]->xp = 0;
	ri->cigar = 0;
	ri->errors = 0;
	ri->seq = malloc(sizeof(char) * allocated_length);
	
	while(current_length < average_length){
		
		while(1){
			//emission
			sum = prob2scaledprob(0.0f);
			for(nuc = 0;nuc < 5;nuc++){
				sum = logsum(sum, mb->model[0]->background_nuc_frequency[nuc] );
				if(r <  scaledprob2prob(sum)){
					ri->seq[current_length] = nuc;
					//fprintf(stderr,"%c",alpha[(int)nuc]);
					current_length++;
					break;
				}
			}
			if(current_length == allocated_length){
				allocated_length = allocated_length*2;
				ri->seq = realloc(ri->seq, sizeof(char) * allocated_length );
			}
			//transition
			r = (float)rand_r(seed)/(float)RAND_MAX;
			// prob2scaledprob(1.0 - (1.0 / (float)len));
			if(r > 1.0 - (1.0 / (float)average_length)){
				break;
			}
		}
		//fprintf(stderr,"	%d\n",current_length);
		if(current_length < average_length){
			current_length = 0;
		}
		
		if(current_length >= MAX_HMM_SEQ_LEN){
			current_length = 0;
		}
	}
	
	ri->seq = realloc(ri->seq, sizeof(char) * (current_length+1));
	ri->seq[current_length] = 0;
	ri->len = current_length;
	
	ri->name = malloc(sizeof(char) *2);
	ri->name[0] = 'N';
	ri->name[1] = 0;
	
	ri->labels = malloc(sizeof(char) * (current_length+1));
	//ri->labels[0] = 'N';
	//ri->labels[1] = 0;
	
	ri->qual = malloc(sizeof(char) *2);
	ri->qual[0] = 'N';
	ri->qual[1] = 0;
	
	return ri;
}


/** \fn struct read_info* emit_read_sequence(struct model_bag* mb, struct read_info* ri,int average_length,unsigned int* seed )
 \brief Emits sequences from read HMM model.
 
 \param mb  The model @ref model_bag .
 \param ri  @ref read_info - emitted sequences are written here.
 \param average_length Average length of the sequences.
 \param seed Seed used for randomization.
 
 */


struct read_info* emit_read_sequence(struct model_bag* mb, struct read_info* ri,int average_length,unsigned int* seed )
{
	
	int i,j,nuc;
	int state = 0; //0 silent ; 1 M , 2 I , 3 D
	int column = 0;
	int hmm = 0;
	int segment= 0;
	//nt region = 0;
	//int start = 1;
	int len;//mb->model[segment]->hmms[0]->num_columns;
	//char alpha[] = "ACGTN";
	//int parashute = 0;
	
	double r = (float)rand_r(seed)/(float)RAND_MAX;
	double sum = prob2scaledprob(0.0f);
	
	double prob = prob2scaledprob(1.0f);
	
	int current_length = 0;
	int allocated_length = 100;
	free(ri->seq);
	free(ri->name);
	free(ri->qual);
	free(ri->labels);
	ri->seq = 0;
	ri->name = 0;
	ri->qual = 0;
	ri->labels = 0;
	ri->len = 0;
	ri->md = 0;
	ri->prob = 0;
	//ri[i]->xp = 0;
	ri->cigar = 0;
	ri->errors = 0;
	
	
	ri->seq = malloc(sizeof(char) * allocated_length);
	
	while(current_length < average_length){
		state = 0; //0 silent ; 1 M , 2 I , 3 D
		column = 0;
		hmm = 0;
		segment= 0;
		
		while(1){
			
			//transition
			r = (float)rand_r(seed)/(float)RAND_MAX;
			sum = prob2scaledprob(0.0f);
			switch (state) {
				case 0:
					//fprintf(stderr,"AM in silent... %f\n",r);
					len = mb->model[segment]->hmms[0]->num_columns;
					
					for(i = 0; i < mb->model[segment]->num_hmms;i++){
						for(j = 0;  j < len;j++){
							sum = logsum(sum,mb->model[segment]->silent_to_M[i][j]);
							//fprintf(stderr,"Trying: hmm:%d Mstate:%d	%f\n", i,j, scaledprob2prob(sum) );
							if(r <  scaledprob2prob(sum) ){
								prob += mb->model[segment]->silent_to_M[i][j];
								state = 1;
								column = j;
								hmm = i;
								//if(!start){
								//	segment++;
								//}
								i = 0x7FFFFFF;
								j =  0x7FFFFFF;
								//start = 0;
								break;
							}
							sum = logsum(sum,mb->model[segment]->silent_to_I[i][j]);
							//fprintf(stderr,"Trying: hmm:%d Istate:%d	%f\n", i,j, scaledprob2prob(sum) );
							if(r <  scaledprob2prob(sum) ){
								prob += mb->model[segment]->silent_to_I[i][j];
								state = 2;
								column = j;
								hmm = i;
								//if(!start){
								//	segment++;
								//}
								i = 0x7FFFFFF;
								j =  0x7FFFFFF;
								//start = 0;
								break;
							}
							
							
						}
					}
					
					
					
					
					
					
					break;
				case 1:
					// MM
					sum = logsum(sum, mb->model[segment]->hmms[hmm]->hmm_column[column]->transition[MM] );
					
					if(r <  scaledprob2prob(sum)){
						prob += mb->model[segment]->hmms[hmm]->hmm_column[column]->transition[MM];
						state = 1;
						column++;
						//hmm = hmm;
						break;
					}
					
					// MI
					sum = logsum(sum, mb->model[segment]->hmms[hmm]->hmm_column[column]->transition[MI] );
					
					if(r <  scaledprob2prob(sum)){
						prob += mb->model[segment]->hmms[hmm]->hmm_column[column]->transition[MI];
						state = 2;
						//column;
						//hmm = hmm;
						break;
					}
					
					// MD
					sum = logsum(sum, mb->model[segment]->hmms[hmm]->hmm_column[column]->transition[MD] );
					
					if(r <  scaledprob2prob(sum)){
						prob += mb->model[segment]->hmms[hmm]->hmm_column[column]->transition[MD];
						state = 2;
						column++;
						//hmm = hmm;
						break;
					}
					
					// MSKIP;
					sum = logsum(sum, mb->model[segment]->hmms[hmm]->hmm_column[column]->transition[MSKIP] );
					prob += mb->model[segment]->hmms[hmm]->hmm_column[column]->transition[MSKIP];
					//if(r <  scaledprob2prob(sum)){
					state = 0;
					segment++;
					column = 0;
					hmm = 0;
					//column++;
					//hmm = hmm;
					//	break;
					//}
					
					
					
					break;
					
				case 2:
					// II
					sum = logsum(sum, mb->model[segment]->hmms[hmm]->hmm_column[column]->transition[II] );
					
					if(r <  scaledprob2prob(sum)){
						prob += mb->model[segment]->hmms[hmm]->hmm_column[column]->transition[II];
						state = 2;
						//column++;
						//hmm = hmm;
						break;
					}
					
					// IM
					sum = logsum(sum, mb->model[segment]->hmms[hmm]->hmm_column[column]->transition[IM] );
					
					if(r <  scaledprob2prob(sum)){
						prob += mb->model[segment]->hmms[hmm]->hmm_column[column]->transition[IM];
						state = 1;
						column++;
						//hmm = hmm;
						break;
					}
					
					// ISKIP
					sum = logsum(sum, mb->model[segment]->hmms[hmm]->hmm_column[column]->transition[ISKIP] );
					prob += mb->model[segment]->hmms[hmm]->hmm_column[column]->transition[ISKIP];
					
					//if(r <  scaledprob2prob(sum)){
					state = 0;
					segment++;
					column = 0;
					hmm = 0;
					//column++;
					//hmm = hmm;
					//	break;
					//}
					
					
					break;
				case 3:
					// DD
					sum = logsum(sum, mb->model[segment]->hmms[hmm]->hmm_column[column]->transition[DD] );
					
					if(r <  scaledprob2prob(sum)){
						prob += mb->model[segment]->hmms[hmm]->hmm_column[column]->transition[DD];
						state = 3;
						column++;
						//hmm = hmm;
						break;
					}
					
					// DM
					sum = logsum(sum, mb->model[segment]->hmms[hmm]->hmm_column[column]->transition[DM] );
					
					//if(r <  scaledprob2prob(sum)){
					prob += mb->model[segment]->hmms[hmm]->hmm_column[column]->transition[DM];
					state = 1;
					column++;
					//hmm = hmm;
					//	break;
					//}
					
					
					break;
					
					
				default:
					break;
			}
			
			//fprintf(stderr,"%d segment \n",segment );
			//fprintf(stderr,"%d hmm \n",hmm );
			//fprintf(stderr,"%d column \n",column );
			//fprintf(stderr,"%d state \n",state );
			//emit...
			r = (float)rand_r(seed)/(float)RAND_MAX;
			sum = prob2scaledprob(0.0f);
			
			if(state == 1){
				for(nuc = 0;nuc < 5;nuc++){
					sum = logsum(sum, mb->model[segment]->hmms[hmm]->hmm_column[column]->m_emit[nuc]);
					if(r <  scaledprob2prob(sum)){
						prob += mb->model[segment]->hmms[hmm]->hmm_column[column]->m_emit[nuc];
						ri->seq[current_length] = nuc;
						
						//fprintf(stderr,"%c",alpha[(int)nuc]);
						current_length++;
						//fprintf(stderr,"Letter: %d	Segment:%d	hmm:%d	column:%d	state:%d\n",nuc, segment,hmm,column,state );
						break;
					}
				}
				if(nuc == 4){
					//fprintf(stderr,"R:%f ",r);
					sum = prob2scaledprob(0.0f);
					for(nuc = 0;nuc < 5;nuc++){
						
						sum = logsum(sum, mb->model[segment]->hmms[hmm]->hmm_column[column]->m_emit[nuc]);
					//	fprintf(stderr,"%f ",scaledprob2prob(sum));
					}
					//fprintf(stderr,"\n");
				}
			}
			
			if(state == 2){
				for(nuc = 0;nuc < 5;nuc++){
					sum = logsum(sum, mb->model[segment]->hmms[hmm]->hmm_column[column]->i_emit[nuc]);
					prob += mb->model[segment]->hmms[hmm]->hmm_column[column]->i_emit[nuc];
					if(r <  scaledprob2prob(sum)){
						ri->seq[current_length] = nuc;
						//fprintf(stderr,"%c",alpha[(int)nuc]);
						current_length++;
						//fprintf(stderr,"Letter: %d	Segment:%d	hmm:%d	column:%d	state:%d\n",nuc, segment,hmm,column,state );
						break;
					}
				}
			}
			
			if(current_length == allocated_length){
				allocated_length = allocated_length*2;
				ri->seq = realloc(ri->seq, sizeof(char) * allocated_length );
			}
			
			//fprintf(stderr,"segement: %d %d\n", segment,mb->num_models);
			if(segment == mb->num_models){
				break;
			}
			
			
			
		}
		if(current_length < average_length){
			current_length = 0;
		}
		
		if(current_length >= MAX_HMM_SEQ_LEN){
			current_length = 0;
		}
		
	}
	//fprintf(stderr,"	%f\n", prob);
	
	
	ri->seq = realloc(ri->seq, sizeof(char) * (current_length+1));
	ri->seq[current_length] = 0;
	ri->len = current_length;
	
	ri->name = malloc(sizeof(char) *2);
	ri->name[0] = 'P';
	ri->name[1] = 0;
	
	ri->labels = malloc(sizeof(char) * (current_length+1));
	
	
	ri->qual = malloc(sizeof(char) *2);
	ri->qual[0] = 'P';
	ri->qual[1] = 0;
	
	return ri;
}


/** \fn struct model_bag* estimate_model_from_labels(struct model_bag* mb, struct parameters* param,  struct read_info** ri,int numseq)
 \brief Estimate model based on labelled sequences 

 \bug not complete - is very buggy. 
\deprecated not used.
 */

struct model_bag* estimate_model_from_labels(struct model_bag* mb, struct parameters* param,  struct read_info** ri,int numseq)
{
	int i,j,c1,c2,c3,g;//,bar,mem;
	char alpha[5] = "ACGTN";
	//set all counts to 1;
	
	mb = set_model_e_to_laplace(mb);
	
	
	char seq[100];
	
	int current_position = 0;
	int current_hmm = 0;
	int current_segment = -1;
	for(i =0; i < numseq;i++){
		if(ri[i]->prob == EXTRACT_SUCCESS){
			current_position = 0;
			current_hmm = -1;
			current_segment = -1;
			// Do stuff with labels & sequences...
			for(j = 0; j < ri[i]->len;j++){
				
				c1 = mb->label[(int)ri[i]->labels[j]];
				c2 = c1 & 0xFFFF; //which segment
				c3 = (c1 >> 16) & 0x7FFF; // which HMM in a particular segment....
				
				if(c2 != current_segment){ // new segment!
					if(current_segment != -1){
						fprintf(stderr,"%c segment:\n", param->read_structure->type[current_segment] );
						for(g = 0;g < current_position;g++){
							fprintf(stderr,"%d ", seq[g]);
							
						}
						fprintf(stderr,"\n");
						
						fprintf(stderr,"%s\n",param->read_structure->sequence_matrix[current_segment][current_hmm]);
					}
					
					current_position = 0;
					current_segment = c2;
					current_hmm = c3;
					
					//update M->silent
					//update I->silent...
					//fprintf(stderr,"%d	%d\n", current_segment,current_hmm);
					
					
					
				}else{ // move into next HMM column.....
					seq[current_position] = ri[i]->seq[j];
					
					current_position++;
				}
				
				
				
				
				
			}
			int key = 0;
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
	//				bar = c3;
	//				mem = c2;
				}
			}
	//		fprintf(stderr,"	bar= %d (%s)\n",bar,   param->read_structure->sequence_matrix[mem][bar] );			
		}
		
		
		
	}
	
	
	
	return mb;
}



/** \fn  struct read_info*  extract_reads(struct model_bag* mb, struct parameters* param,  struct read_info* ri)
 \brief Extracts reads from labelled raw sequences.
 
This function extracts the mappable reads from the raw sequences. Barcodes and Fingerprint sequences are decoded and appended to re read names using the BC and FP tag.  
 
\param mb The HMM model @ref model_bag.
 \param param @ref parameters .
 \param ri The reads.
 
 */


 struct read_info*  extract_reads(struct model_bag* mb, struct parameters* param,  struct read_info* ri)
{
	int j,c1,c2,c3,key,bar,mem,fingerlen,required_finger_len,ret;
	char buffer[MAX_HMM_SEQ_LEN];
	int s_pos = 0;
	key = 0;
	bar = -1;
	mem = -1;
	ret = 0;
	int offset = 0;
	int len;
	int hmm_has_barcode = 0;
	int read_start = -1;
	
	len = ri->len;
	if(param->matchstart != -1 || param->matchend != -1){
		offset = param->matchstart;
		len = param->matchend - param->matchstart;
	}
	required_finger_len = 0;
	for(j = 0; j < param->read_structure->num_segments;j++){
		if(param->read_structure->type[j] == 'F'){
			required_finger_len += (int) strlen(param->read_structure->sequence_matrix[j][0]);
		}
	}
	
	//ri[i]->prob = expf( ri[i]->prob) / (1.0f + expf(ri[i]->prob ));
	
	if(param->confidence_threshold <=  ri->mapq){
		
		if(0.5 <=  ri->bar_prob){
			fingerlen = 0;
			//required_finger_len = 0;
			
			for(j = 0; j < len;j++){
				c1 = mb->label[(int)ri->labels[j+1]];
				c2 = c1 & 0xFFFF; //which segment
				c3 = (c1 >> 16) & 0x7FFF; // which HMM in a particular segment.... 
				//fprintf(stderr,"%c",   param->read_structure->type[c2] );
				if(param->read_structure->type[c2] == 'F'){
					//	required_finger_len += (int) strlen(param->read_structure->sequence_matrix[c2][0]);
					fingerlen++;
					key = (key << 2 )|  (ri->seq[j+offset] & 0x3);
				}
				if(param->read_structure->type[c2] == 'B'){
					hmm_has_barcode = 1;
					bar = c3;
					mem = c2;
				}
				if(param->read_structure->type[c2] == 'R'){
					if(read_start == -1){
						read_start = j+offset;
					}
					s_pos++;
				}
			}
			for(j = len; j < ri->len;j++){
				s_pos++;
			}
			
			if(s_pos >= param->minlen){
				
				if(hmm_has_barcode && required_finger_len){
					if(fingerlen == required_finger_len && bar != -1){
						buffer[0] = 0;
						sprintf (buffer, "@%s;BC:%s;FP:%d;",ri->name,param->read_structure->sequence_matrix[mem][bar],key);
						//strcat (buffer, tmp);
						ri->name = realloc(ri->name, sizeof(char) * (strlen(buffer) + 1) );
						
						strcpy(ri->name, buffer);
						for(j = 0; j < s_pos;j++){
							ri->seq[j] = ri->seq[read_start+j];
							ri->qual[j] = ri->qual[read_start+j];
						}
						ri->len = s_pos;
						ri->prob = EXTRACT_SUCCESS;
						//ret = 1;
						//fprintf(out,"@%s;BC:%s;FP:%d\n",ri->name,param->read_structure->sequence_matrix[mem][bar],key);
						//fprintf(out,"%s\n+\n%s\n", out_seq,out_qual);
					}else{
						ri->prob  = EXTRACT_FAIL_BAR_FINGER_NOT_FOUND; // something wrong with the architecture
					}
				}else if(hmm_has_barcode){
					if(bar != -1){
						
						buffer[0] = 0;
						sprintf (buffer, "@%s;BC:%s;",ri->name,param->read_structure->sequence_matrix[mem][bar]);
						//strcat (buffer, tmp);
						ri->name = realloc(ri->name, sizeof(char) * (strlen(buffer) + 1) );
						
						strcpy(ri->name, buffer);
						for(j = 0; j < s_pos;j++){
							ri->seq[j] = ri->seq[read_start+j];
							ri->qual[j] = ri->qual[read_start+j];
						}
						ri->len = s_pos;
						ri->prob = EXTRACT_SUCCESS;
					}else{
						ri->prob  = EXTRACT_FAIL_BAR_FINGER_NOT_FOUND;
					}
					
				}else if(required_finger_len){
					if(fingerlen == required_finger_len){
						buffer[0] = 0;
						sprintf (buffer, "@%s;FP:%d",ri->name,key);
						//strcat (buffer, tmp);
						ri->name = realloc(ri->name, sizeof(char) * (strlen(buffer) + 1) );
						
						strcpy(ri->name, buffer);
						for(j = 0; j < s_pos;j++){
							ri->seq[j] = ri->seq[read_start+j];
							ri->qual[j] = ri->qual[read_start+j];
						}
						ri->len = s_pos;
						ri->prob = EXTRACT_SUCCESS;
					}else{
						ri->prob  = EXTRACT_FAIL_BAR_FINGER_NOT_FOUND;
					}
				}else{
					buffer[0] = 0;
					sprintf (buffer, "@%s",ri->name);
					//strcat (buffer, tmp);
					ri->name = realloc(ri->name, sizeof(char) * (strlen(buffer) + 1) );
					
					strcpy(ri->name, buffer);
					for(j = 0; j < s_pos;j++){
						ri->seq[j] = ri->seq[read_start+j];
						ri->qual[j] = ri->qual[read_start+j];
					}
					ri->len = s_pos;
					ri->prob = EXTRACT_SUCCESS;
				}
			}else{
				ri->prob = EXTRACT_FAIL_READ_TOO_SHORT;
			}
		}else{
			ri->prob = EXTRACT_FAIL_AMBIGIOUS_BARCODE;
		}
	}else{
		ri->prob = EXTRACT_FAIL_ARCHITECTURE_MISMATCH;
	}
	
	ri->qual[ri->len] = 0;
	
	return ri;
}



/** \fn void* do_baum_welch_thread(void *threadarg)
 \brief Runs Baum-Welch procedure. 
 
We implemented the HMM training procedure to verify that the forward and backward recursion work as expected. In other words this is only used for testing. 
 
 \param threadarg  A @ref thread_data object used to pass data / parameters to function.
 */
void* do_baum_welch_thread(void *threadarg)
{
	struct thread_data *data;
	data = (struct thread_data *) threadarg;
	
	struct read_info** ri  = data->ri;
	struct model_bag* mb = data->mb;
	
	int start = data->start;
	int end = data->end;
	int i;
	float pbest;
	float Q;
	for(i = start; i < end;i++){
		mb = backward(mb, ri[i]->seq ,ri[i]->len);
		mb = forward_max_posterior_decoding(mb, ri[i], ri[i]->seq ,ri[i]->len);
		
		pbest = 1.0 - scaledprob2prob(  (ri[i]->bar_prob + mb->f_score) - logsum(mb->f_score,  mb->r_score));
		if(!pbest){
			Q = 40.0;
		}else if(pbest == 1.0){
			Q = 0.0;
			
		}else{
			Q = -10.0 * log10(pbest) ;
		}
		if(Q >= 10){
			mb = forward_extract_posteriors(mb, ri[i]->seq,ri[i]->labels,ri[i]->len);
		}
	}
	pthread_exit((void *) 0);
}


/** \fn void* do_run_random_sequences(void *threadarg)
 \brief Runs Forward and BAckward oin random sequences. 
 The read probabilities are stored in @ref model_bag. 
 
 
 
 \param threadarg  A @ref thread_data object used to pass data / parameters to function.
 \deprecated We do not need this anymore. 
 
 */

void* do_run_random_sequences(void *threadarg)
{
	struct thread_data *data;
	data = (struct thread_data *) threadarg;
	
	struct read_info** ri  = data->ri;
	struct model_bag* mb = data->mb;
	
	int start = data->start;
	int end = data->end;
	int i,tmp;
	int r,bar;
	
	int matchstart = data->param->matchstart;
	int matchend = data->param->matchend;
	//int c;
	//int j;
	//f/loat r;
	//char random[MAX_HMM_SEQ_LEN];
	
	//unsigned int seed = (unsigned int) (time(NULL) * ( 42));
	
	
	//fprintf(stderr," %d - %d\n", start,end);

	/*float a,c,g;
	
	a = scaledprob2prob( mb->model[0]->background_nuc_frequency[0]);
	
	c = a +  scaledprob2prob(mb->model[0]->background_nuc_frequency[1]);
	
	g = c +  scaledprob2prob(mb->model[0]->background_nuc_frequency[2]);
	
	*/
	bar = -1;
	r = -1;
	//switch read and barcode...
	struct model* tmp_model = 0;
	for(i = 0; i < mb->num_models;i++){
		if(data->param->read_structure->type[i] == 'R'){
			r = i;
		}
		if(data->param->read_structure->type[i] == 'B'){
			bar= i;
		}
	}
	if(bar != -1){
		tmp_model = mb->model[r];
		mb->model[r] = mb->model[bar];
		mb->model[bar] = tmp_model;
	}
	
	
	double pbest,Q;
	if(matchstart != -1 || matchend != -1){
		for(i = start; i < end;i++){
			tmp = matchend - matchstart ;
			if(bar == -1){
				reverse_sequence(ri[i]->seq,  ri[i]->len);
			}
			mb = backward(mb, ri[i]->seq + matchstart , tmp);
			mb = forward_max_posterior_decoding(mb, ri[i] , ri[i]->seq+matchstart ,tmp );
			if(bar == -1){
				reverse_sequence(ri[i]->seq,  ri[i]->len);
			}
			
			pbest = logsum(mb->f_score, mb->r_score);
			pbest = 1.0 - scaledprob2prob(  (ri[i]->bar_prob + mb->f_score) - pbest);
			if(!pbest){
				Q = 40.0;
			}else if(pbest == 1.0){
				Q = 0.0;
				
			}else{
				Q = -10.0 * log10(pbest) ;
			}
			mb->random_scores[i] = Q;
		}
	}else{
		for(i = start; i < end;i++){
		
			
			if(bar == -1){
				reverse_sequence(ri[i]->seq,  ri[i]->len);
			}
			
			mb = backward(mb, ri[i]->seq,ri[i]->len);
			mb = forward_max_posterior_decoding(mb, ri[i], ri[i]->seq ,ri[i]->len);
			if(bar == -1){
				reverse_sequence(ri[i]->seq,  ri[i]->len);
			}
			fprintf(stderr,"%f	%f	%f	%f\n",mb->f_score,mb->b_score,mb->r_score,mb->f_score-mb->b_score  );
			pbest = logsum(mb->f_score, mb->r_score);
			
			pbest = 1.0 - scaledprob2prob(  (ri[i]->bar_prob + mb->f_score) - pbest);
			if(!pbest){
				Q = 40.0;
			}else if(pbest == 1.0){
				Q = 0.0;
				
			}else{
				Q = -10.0 * log10(pbest) ;
			}
			//ri[i]->prob = Q;
			
			
			
			//fprintf(stdout,"%f\n", ri[i]->prob);
			
			mb->random_scores[i] = Q;
		}
		//mb = forward_extract_posteriors(mb, ri[i]->seq ,ri[i]->len);
	}
	
	if(bar != -1){
		tmp_model = mb->model[r];
		mb->model[r] = mb->model[bar];
		mb->model[bar] = tmp_model;
	}
	
	
	pthread_exit((void *) 0);
}




/** \fn struct model_bag* backward(struct model_bag* mb, char* a, int len)
 \brief Runs the backard algorithm .
 
 
 
 
 \param mb The HMM model. 
 \param a The seuqence.
 \param len Sequence length. 
 */


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
				
				//for(i = len-1 ; i >= 0;i--){ /// DONT MESS WITH THIS - the code takes the first letter to calculate the silent states - EVERYHTHIN is ok..
				// in the last column an I state will emit a letter not present in the seuqence BUT this will not make it when joined to the forward... messy but works.  
				//	c = (int) seqa[i+1];
				
				c = (int)seqa[i+1];
								
				c_hmm_column = hmm->hmm_column[model_len];
				
				//c_hmm_column->M_backward[i] = psilent[i+1] + mb->model[j]->M_to_silent[f] ;
				
				c_hmm_column->M_backward[i] = psilent[i+1] + c_hmm_column->transition[MSKIP];// mb->model[j]->M_to_silent[f] ;
				
				
				//fprintf(stderr," Mback at modellen:%d %f %f\n",i, c_hmm_column->M_backward[i] ,c_hmm_column->transition[MSKIP]);
				
				//c_hmm_column->I_backward[i] =  psilent[i+1]+ mb->model[j]->I_to_silent[f] ;
				
				c_hmm_column->I_backward[i] =  psilent[i+1] + c_hmm_column->transition[ISKIP];//  mb->model[j]->I_to_silent[f] ;
				
				
				c_hmm_column->I_backward[i] = logsum(c_hmm_column->I_backward[i] , c_hmm_column->M_backward[i+1] + c_hmm_column->transition[IM] + c_hmm_column->m_emit[c]);
				
				c_hmm_column->I_backward[i] = logsum(c_hmm_column->I_backward[i] , c_hmm_column->I_backward[i+1] + c_hmm_column->transition[II] + c_hmm_column->i_emit[c]);
				
				//##################################
				csilent[i] = logsum(csilent[i], c_hmm_column->M_backward[i] + mb->model[j]->silent_to_M[f][model_len] + c_hmm_column->m_emit[(int)seqa[i]]);
				csilent[i] = logsum(csilent[i], c_hmm_column->I_backward[i]+ mb->model[j]->silent_to_I[f][model_len] + c_hmm_column->i_emit[(int)seqa[i]]);
				//##################################
				//fprintf(stderr," M:%f I:%f \n", c_hmm_column->M_backward[i] ,c_hmm_column->I_backward[i] );
				
				
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
					//GRERRRRRRR
					
					//delete state
					
					//from previous delection
					c_hmm_column->D_backward[i] = p_hmm_column->D_backward[i] + c_hmm_column->transition[DD];
					
					//from previous match (i.e. gap close
					
					c_hmm_column->D_backward[i] = logsum(c_hmm_column->D_backward[i], p_hmm_column->M_backward[i] + p_hmm_column->m_emit[(int) seqa[i]] + c_hmm_column->transition[DM]);
					
					
					//##################################
					csilent[i] = logsum(csilent[i], c_hmm_column->M_backward[i] + mb->model[j]->silent_to_M[f][g] + c_hmm_column->m_emit[(int)seqa[i]]);
					csilent[i] = logsum(csilent[i], c_hmm_column->I_backward[i]+ mb->model[j]->silent_to_I[f][g] + c_hmm_column->i_emit[(int)seqa[i]]);
					//##################################
					
					
				}
				
				//##################################
				
				//c_hmm_column = hmm->hmm_column[0];
				// link j+1 to j... dfor silent;
				//csilent[i] = logsum(csilent[i], c_hmm_column->M_backward[i] + mb->model[j]->silent_to_M[f][0] + c_hmm_column->m_emit[(int)seqa[i]]);
				//csilent[i] = logsum(csilent[i], c_hmm_column->I_backward[i]+ mb->model[j]->silent_to_I[f][0] + c_hmm_column->i_emit[(int)seqa[i]]);
				
				//##################################
				
				//fprintf(stderr,"Looking for Insertyion to silent in segment1: %d	%f\n",f, mb->model[j]->silent_to_I[f]);
				
				//this should come from previous state .....
				csilent[i] = logsum(csilent[i], psilent[i] + mb->model[j]->skip);
				
				
				
				//exit(0);
				
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



/** \fn struct model_bag* backward(struct model_bag* mb, char* a, int len)
 \brief Runs the forward algorithm .
 
 
 
 
 \param mb The HMM model.
 \param a The seuqence.
 \param len Sequence length.
 */

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
				c_hmm_column->M_foward[i] = psilent[i-1] + mb->model[j]->silent_to_M[f][0] + c_hmm_column->m_emit[c];
				
 				c_hmm_column->I_foward[i]    =psilent[i-1] + mb->model[j]->silent_to_I[f][0] ;
				
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
					
					c_hmm_column->M_foward[i] = psilent[i-1] + mb->model[j]->silent_to_M[f][g];
					
					//transition from previous match state
					c_hmm_column->M_foward[i] = logsum( c_hmm_column->M_foward[i] , p_hmm_column->M_foward[i-1] + p_hmm_column->transition[MM]);
					//transition from previous insert state
					c_hmm_column->M_foward[i] = logsum(c_hmm_column->M_foward[i] , p_hmm_column->I_foward[i-1] + p_hmm_column->transition[IM] );
					//transition from previous delete state
					c_hmm_column->M_foward[i] = logsum(c_hmm_column->M_foward[i], p_hmm_column->D_foward[i] + p_hmm_column->transition[DM]);
					
					// emission promability in curent M state ;
					c_hmm_column->M_foward[i]  = c_hmm_column->M_foward[i]  + c_hmm_column->m_emit[c];
					
					
					// Instertion State ..
					
					c_hmm_column->I_foward[i]    =psilent[i-1] + mb->model[j]->silent_to_I[f][g] ;
					
					//self loop insertion to insertion
					c_hmm_column->I_foward[i] =  logsum ( c_hmm_column->I_foward[i] ,c_hmm_column->I_foward[i-1] + c_hmm_column->transition[II]);
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



/** \fn struct model_bag* backward(struct model_bag* mb, char* a, int len)
 \brief Runs the forward algorithm and extracts estimated transition and emission probabilities.
 \warning Need to run @ref backward first!
 
 
 
 \param mb The HMM model.
 \param label The label of the sequence.
 \param a The seuqence.
 \param len Sequence length.
 */

struct model_bag* forward_extract_posteriors(struct model_bag* mb, char* a, char* label, int len)
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
				c_hmm_column->M_foward[i] = psilent[i-1] + mb->model[j]->silent_to_M[f][0] + c_hmm_column->m_emit[c];
				
				//***************post
				
				//if(label[i] == hmm_counter){
				
				mb->model[j]->silent_to_M_e[f][0] = logsum(mb->model[j]->silent_to_M_e[f][0] ,psilent[i-1] + mb->model[j]->silent_to_M[f][0] + c_hmm_column->m_emit[c] + c_hmm_column->M_backward[i]  -mb->b_score);
				
				
				c_hmm_column->m_emit_e[c] = logsum(c_hmm_column->m_emit_e[c] , c_hmm_column->M_foward[i]  + c_hmm_column->M_backward[i] -mb->b_score);
				//}
				//***************post
				
				
 				c_hmm_column->I_foward[i]    =psilent[i-1] + mb->model[j]->silent_to_I[f][0] ;
				
				
				
				
				//add transitions to first columns////
				c_hmm_column->I_foward[i] = logsum(c_hmm_column->I_foward[i], c_hmm_column->I_foward[i-1] + c_hmm_column->transition[II]);
				
				c_hmm_column->I_foward[i]  = logsum ( c_hmm_column->I_foward[i] ,c_hmm_column->M_foward[i-1] + c_hmm_column->transition[MI]);
				
				c_hmm_column->I_foward[i] = c_hmm_column->I_foward[i] + c_hmm_column->i_emit[c];
				
				//***************post
				//if(label[i] == hmm_counter){
				mb->model[j]->silent_to_I_e[f][0]  = logsum(mb->model[j]->silent_to_I_e[f][0] , psilent[i-1] + mb->model[j]->silent_to_I[f][0]  + c_hmm_column->i_emit[c] +  c_hmm_column->I_backward[i] -mb->b_score);
				
				c_hmm_column->transition_e[II] = logsum(c_hmm_column->transition_e[II],c_hmm_column->I_foward[i-1] + c_hmm_column->transition[II] + c_hmm_column->i_emit[c] +  c_hmm_column->I_backward[i]-mb->b_score);
				
				c_hmm_column->transition_e[MI] = logsum(c_hmm_column->transition_e[MI] , c_hmm_column->M_foward[i-1] + c_hmm_column->transition[MI] +  c_hmm_column->i_emit[c] +  c_hmm_column->I_backward[i]-mb->b_score);
				
				
				c_hmm_column->i_emit_e[c] = logsum( c_hmm_column->i_emit_e[c] , c_hmm_column->I_foward[i] + c_hmm_column->I_backward[i] -mb->b_score);
				//}
				//***************post
				
				
				c_hmm_column->D_foward[i] = prob2scaledprob(0.0f);
				
				// no post???
				
				//
				
				csilent[i] =  logsum(csilent[i],  c_hmm_column->M_foward[i] +c_hmm_column->transition[MSKIP]);
				csilent[i] =  logsum(csilent[i],  c_hmm_column->I_foward[i] +c_hmm_column->transition[ISKIP]);
				
				//***************post
				//if(label[i] == hmm_counter){
				c_hmm_column->transition_e[MSKIP] = logsum(c_hmm_column->transition_e[MSKIP], c_hmm_column->M_foward[i] +c_hmm_column->transition[MSKIP] +  bsilent[i+1] -mb->b_score);
				
				c_hmm_column->transition_e[ISKIP] = logsum(c_hmm_column->transition_e[ISKIP], c_hmm_column->I_foward[i] +c_hmm_column->transition[ISKIP] +  bsilent[i+1] -mb->b_score);
				//}
				
				
				//***************post
				
				
				for(g = 1;g < hmm->num_columns;g++){
					c_hmm_column = hmm->hmm_column[g];
					p_hmm_column = hmm->hmm_column[g-1];
					
					//Match state
					
					
					c_hmm_column->M_foward[i] = psilent[i-1] + mb->model[j]->silent_to_M[f][g];
					
					//transition from previous match state
					c_hmm_column->M_foward[i] =  logsum(c_hmm_column->M_foward[i] , p_hmm_column->M_foward[i-1] + p_hmm_column->transition[MM]);
					//transition from previous insert state
					c_hmm_column->M_foward[i] = logsum(c_hmm_column->M_foward[i] , p_hmm_column->I_foward[i-1] + p_hmm_column->transition[IM] );
					//transition from previous delete state
					c_hmm_column->M_foward[i] = logsum(c_hmm_column->M_foward[i], p_hmm_column->D_foward[i] + p_hmm_column->transition[DM]);
					
					// emission promability in curent M state ;
					c_hmm_column->M_foward[i]  = c_hmm_column->M_foward[i]  + c_hmm_column->m_emit[c];
					
					
					//***************post
					//if(label[i] == hmm_counter){
					
					//I HOPE THIS IS CORRECT....
					mb->model[j]->silent_to_M_e[f][g] = logsum(mb->model[j]->silent_to_M_e[f][g] ,psilent[i-1] + mb->model[j]->silent_to_M[f][g] + c_hmm_column->m_emit[c] + c_hmm_column->M_backward[i]  -mb->b_score);
					
					//I HOPE THIS IS CORRECT.... 
					
					
					p_hmm_column->transition_e[MM] = logsum(p_hmm_column->transition_e[MM] , p_hmm_column->M_foward[i-1] + p_hmm_column->transition[MM] +  c_hmm_column->m_emit[c] +  c_hmm_column->M_backward[i] -mb->b_score );
					
					p_hmm_column->transition_e[IM] = logsum(p_hmm_column->transition_e[IM],p_hmm_column->I_foward[i-1] + p_hmm_column->transition[IM] +  c_hmm_column->m_emit[c] +  c_hmm_column->M_backward[i] -mb->b_score);
					
					p_hmm_column->transition_e[DM] = logsum(p_hmm_column->transition_e[DM],p_hmm_column->D_foward[i] + p_hmm_column->transition[DM] + c_hmm_column->m_emit[c] +  c_hmm_column->M_backward[i] -mb->b_score);
					
					c_hmm_column->m_emit_e[c] = logsum(c_hmm_column->m_emit_e[c],  c_hmm_column->M_foward[i] + c_hmm_column->M_backward[i] -mb->b_score );
					//}
					//***************post
					
					
					
					
					
					// Instertion State ..
					
					c_hmm_column->I_foward[i]    =psilent[i-1] + mb->model[j]->silent_to_I[f][g] ;
					
					
					//self loop insertion to insertion
					c_hmm_column->I_foward[i] = logsum ( c_hmm_column->I_foward[i] ,c_hmm_column->I_foward[i-1] + c_hmm_column->transition[II]);
					// start new insertion
					c_hmm_column->I_foward[i]  = logsum ( c_hmm_column->I_foward[i] ,c_hmm_column->M_foward[i-1] + c_hmm_column->transition[MI]);
					
					//instertion emission...
					c_hmm_column->I_foward[i]  = c_hmm_column->I_foward[i]  + c_hmm_column->i_emit[c];
					
					//***************post
					//if(label[i] == hmm_counter){
					//I HOPE THIS IS CORRECT....

					mb->model[j]->silent_to_I_e[f][g]  = logsum(mb->model[j]->silent_to_I_e[f][g] , psilent[i-1] + mb->model[j]->silent_to_I[f][g]  + c_hmm_column->i_emit[c] +  c_hmm_column->I_backward[i] -mb->b_score);
					
					//I HOPE THIS IS CORRECT....

					
					c_hmm_column->transition_e[II] = logsum(c_hmm_column->transition_e[II] ,  c_hmm_column->I_foward[i-1] + c_hmm_column->transition[II] + c_hmm_column->i_emit[c] + c_hmm_column->I_backward[i] - mb->b_score);
					
					c_hmm_column->transition_e[MI] = logsum(c_hmm_column->transition_e[MI] , c_hmm_column->M_foward[i-1] + c_hmm_column->transition[MI] + c_hmm_column->i_emit[c] + c_hmm_column->I_backward[i] - mb->b_score);
					
					c_hmm_column->i_emit_e[c] = logsum(c_hmm_column->i_emit_e[c] , c_hmm_column->I_foward[i]   + c_hmm_column->I_backward[i] - mb->b_score);
					//}
					//***************post

					
					
					// deletion state
					//from previous match state.
					c_hmm_column->D_foward[i] = p_hmm_column->M_foward[i] + p_hmm_column->transition[MD];
					//from previous delete state
					
					c_hmm_column->D_foward[i] = logsum(c_hmm_column->D_foward[i], p_hmm_column->D_foward[i] + p_hmm_column->transition[DD] );
					
					//***************post
					//if(label[i] == hmm_counter){
					p_hmm_column->transition_e[MD] = logsum(p_hmm_column->transition_e[MD],  p_hmm_column->M_foward[i] + p_hmm_column->transition[MD] + c_hmm_column->D_backward[i] - mb->b_score);
					
					p_hmm_column->transition_e[DD] = logsum(p_hmm_column->transition_e[DD] ,p_hmm_column->D_foward[i] + p_hmm_column->transition[DD] + c_hmm_column->D_backward[i] - mb->b_score);
					//}
					//***************post
					csilent[i] =  logsum(csilent[i],  c_hmm_column->M_foward[i] +c_hmm_column->transition[MSKIP]);
					csilent[i] =  logsum(csilent[i],  c_hmm_column->I_foward[i] +c_hmm_column->transition[ISKIP]);
					
					//***************post
					//if(label[i] == hmm_counter){
					c_hmm_column->transition_e[MSKIP] = logsum(c_hmm_column->transition_e[MSKIP], c_hmm_column->M_foward[i] +c_hmm_column->transition[MSKIP] +  bsilent[i+1] -mb->b_score);
					
					c_hmm_column->transition_e[ISKIP] = logsum(c_hmm_column->transition_e[ISKIP], c_hmm_column->I_foward[i] +c_hmm_column->transition[ISKIP] +  bsilent[i+1] -mb->b_score);
					//}
					
					
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
				//if(label[i] == hmm_counter){
				mb->model[j]->skip_e =logsum(mb->model[j]->skip_e , psilent[i-1] + mb->model[j]->skip + bsilent[i] -mb->b_score);
 				
				//}
				
				//***************post
				
			}
			hmm_counter++;
		}
	}
	
	
	mb->f_score = mb->model[mb->num_models-1]->silent_forward[len];
	
	//fprintf(stderr,"SCORE:%f	%f\n", mb->f_score, scaledprob2prob(mb->f_score));
	return mb;
}


/** \fn struct model_bag* forward_max_posterior_decoding(struct model_bag* mb, struct read_info* ri, char* a, int len)
 \brief Runs the forward algorithm and labels sequences.
 
 Main function in tagdust2. Runs the forward algorithm and estimates posterior label probabilities of all residures. Then runs normal dynamic programming algorithm to obtain the most likely path given the posteriors. Gap penalties ensure that the resulting labeling matches the user specified architecture.
 
 \warning Need to run @ref backward first!
 
 
 
 \param mb The HMM model.
 \param label The label of the sequence.
 \param a The seuqence.
 \param len Sequence length.
 */


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

	
	float total_prob[100];
	for(j = 0; j < mb->total_hmm_num;j++){
		total_prob[j] = prob2scaledprob(0.0);
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
				c_hmm_column->M_foward[i] = psilent[i-1] + mb->model[j]->silent_to_M[f][0] + c_hmm_column->m_emit[c];
				
				//***************post
				//mb->model[j]->silent_to_M_e[f] = logsum(mb->model[j]->silent_to_M_e[f] ,psilent[i-1] + mb->model[j]->silent_to_M[f] + c_hmm_column->m_emit[c] + c_hmm_column->M_backward[i]  -mb->b_score);
				
				
				
				total_prob[hmm_counter] = logsum(total_prob[hmm_counter], c_hmm_column->M_foward[i]  +  c_hmm_column->M_backward[i] -mb->b_score );
				mb->dyn_prog_matrix[i][hmm_counter] = logsum(mb->dyn_prog_matrix[i][hmm_counter],  c_hmm_column->M_foward[i]  + c_hmm_column->M_backward[i] -mb->b_score );
				
				//c_hmm_column->m_emit_e[c] = logsum(c_hmm_column->m_emit_e[c] , c_hmm_column->M_foward[i]  + c_hmm_column->M_backward[i] -mb->b_score);
				//***************post
				
				
 				c_hmm_column->I_foward[i] = psilent[i-1] + mb->model[j]->silent_to_I[f][0] ;
				
				
				
				
				//add transitions to first columns////
				c_hmm_column->I_foward[i] = logsum(c_hmm_column->I_foward[i], c_hmm_column->I_foward[i-1] + c_hmm_column->transition[II]);
				
				c_hmm_column->I_foward[i] = logsum ( c_hmm_column->I_foward[i] ,c_hmm_column->M_foward[i-1] + c_hmm_column->transition[MI]);
				
				c_hmm_column->I_foward[i] = c_hmm_column->I_foward[i] + c_hmm_column->i_emit[c];
				
				
				
				
				//***************post
				
				
				total_prob[hmm_counter] = logsum(total_prob[hmm_counter], psilent[i-1] + mb->model[j]->silent_to_I[f][0]  + c_hmm_column->i_emit[c] + c_hmm_column->I_backward[i] -mb->b_score );
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
					
					//going from silent to internal HMM columns.
					c_hmm_column->M_foward[i] = psilent[i-1] + mb->model[j]->silent_to_M[f][g];// + c_hmm_column->m_emit[c];
					
					//transition from previous match state
					c_hmm_column->M_foward[i] = logsum(c_hmm_column->M_foward[i] , p_hmm_column->M_foward[i-1] + p_hmm_column->transition[MM]);
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
					
					c_hmm_column->I_foward[i] = psilent[i-1] + mb->model[j]->silent_to_I[f][g] ;
					
					//self loop insertion to insertion
					c_hmm_column->I_foward[i] = logsum (c_hmm_column->I_foward[i] ,  c_hmm_column->I_foward[i-1] + c_hmm_column->transition[II]);
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
	
	// get barcode score....
	hmm_counter = 0;
	g = 1;
	next_silent[0] = prob2scaledprob(0.0);
	next_silent[1] = prob2scaledprob(0.0);
	for(j = 0; j < mb->num_models;j++){
		if(mb->model[j]->num_hmms > 1){
			g = 0;
			next_silent[1] = prob2scaledprob(0.0);
			for(f = 0;f < mb->model[j]->num_hmms;f++){
				if(total_prob[hmm_counter] > next_silent[0]){
					next_silent[0] = total_prob[hmm_counter];
					next_silent[1] = logsum(next_silent[1] , total_prob[hmm_counter]);
				}
				//fprintf(stderr,"%d %f	%f\n",f,total_prob[hmm_counter],scaledprob2prob( total_prob[hmm_counter]));
				hmm_counter++;
			}
			next_silent[0] = next_silent[0] - next_silent[1]; // this ensures that barprob is never > 1 (happens due to numerical inaccuracy... )
			//fprintf(stderr,"SUM:%f	%f\n\n", next_silent[1] , scaledprob2prob(next_silent[1] ));
			
			
			
		}else{
			hmm_counter+= mb->model[j]->num_hmms;
		}
	}
	
	if(g){
		ri->bar_prob = prob2scaledprob(1.0);
	}else{
		ri->bar_prob  = next_silent[0];
	}
	
	for(i = 0; i <= len;i++){
	//	fprintf(stderr,"%d ",i);
		for(j = 0; j < mb->total_hmm_num;j++){
	//		fprintf(stderr,"%0.3f ",scaledprob2prob( mb->dyn_prog_matrix[i][j]));
			mb->dyn_prog_matrix[i][j] = scaledprob2prob( mb->dyn_prog_matrix[i][j]);
		//	total_prob[j] = logsum(total_prob[j] ,  mb->dyn_prog_matrix[i][j]);
		}
		
	//	fprintf(stderr,"\n");
	}
	/*fprintf(stderr,"totalprob: \n");
	for(j = 0; j < mb->total_hmm_num;j++){
		fprintf(stderr,"%d	%d	%f	%f\n", j, mb->label[j], total_prob[j], scaledprob2prob(total_prob[j]));
	}*/
	
	
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

	mb->r_score  = prob2scaledprob(1.0);
	
	for(i = 1; i <= len;i++){
		c = seqa[i];
		mb->r_score  = mb->r_score  + mb->model[0]->background_nuc_frequency[c] + prob2scaledprob(1.0 - (1.0 / (float)len));
		
		//fprintf(stderr,"%d,%f	%e	%f	%f	%f\n",   i,scaledprob2prob(next_silent[0]),   scaledprob2prob(next_silent[0]),next_silent[0] , scaledprob2prob(mb->model[0]->background_nuc_frequency[c] ) , 1.0 - (1.0 / (float)len) );
	}
	mb->r_score  += prob2scaledprob(1.0 / (float)len);
	return mb;
}




/** \fn struct model* malloc_model(int main_length, int sub_length, int number_sub_models)
 \brief Allocates a HMM segment. 
 
 \param main_length Length of the first HMM. 
 \param sub_length Length of second... N HMM.
 \param number_sub_models Number of HMMs. 
 \deprecated Not used anymore. 
 */


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



/** \fn struct model* malloc_model_according_to_read_structure(int num_hmm, int length)
 \brief Allocates a HMM segment.
 
 \param length Length of all HMMs.
 \param num_hmm Number of HMMs.
 */

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

	model->silent_to_M = malloc(sizeof(float*) * model->num_hmms);
	model->silent_to_I = malloc(sizeof(float*) * model->num_hmms);
	model->silent_to_M_e = malloc(sizeof(float*) * model->num_hmms);
	model->silent_to_I_e = malloc(sizeof(float*) * model->num_hmms);
	
	
	len = length;// (int)strlen(rs->sequence_matrix[key][0]);
	
	for(i = 0 ;i  < model->num_hmms;i++){
		model->silent_to_M[i] =  malloc(sizeof(float) * len);
		model->silent_to_I[i] =  malloc(sizeof(float) * len);
		
		model->silent_to_M_e[i] =  malloc(sizeof(float) * len);
		model->silent_to_I_e[i] =  malloc(sizeof(float) * len);

		for(j = 0; j < len;j++){
			model->silent_to_M[i][j] =  0.0f;
			model->silent_to_I[i][j] = 0.0f;
			model->silent_to_M_e[i][j] =  0.0f;
			model->silent_to_I_e[i][j] = 0.0f;
		}
	}
	
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


/** \fn struct model* init_model_according_to_read_structure(struct model* model,struct parameters* param , int key, double* background,int assumed_length)

 \brief Initialized whole HMM.
 
 \param model The model to be initialized.
 \param param @ref parameters.
 \param key The key to the HMM type. 
 \param background Background nucleotide probabilities.
 \param assumed_length Length of the HMM.
 */


struct model* init_model_according_to_read_structure(struct model* model,struct parameters* param , int key, double* background,int assumed_length)
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
		
		//sets emission probabilities...
		for(j = 0; j < len;j++){
			col = model->hmms[i]->hmm_column[j];
			col->transition_e[MM] =  prob2scaledprob(0.0);
			col->transition_e[MI] =  prob2scaledprob(0.0);
			col->transition_e[MD] =  prob2scaledprob(0.0);
			col->transition_e[MSKIP] =  prob2scaledprob(0.0);
			
			col->transition_e[II] =  prob2scaledprob(0.0);
			col->transition_e[IM] =  prob2scaledprob(0.0);
			col->transition_e[ISKIP] =  prob2scaledprob(0.0);
			
			col->transition_e[DD] =  prob2scaledprob(0.0);
			col->transition_e[DM] =  prob2scaledprob(0.0);
			current_nuc = nuc_code[(int) tmp[j]];
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
		}
		//sets transition probabilities....
		model->hmms[i] = set_hmm_transition_parameters(model->hmms[i],len, base_error, indel_freq, -1.0, -1.0);
		/*if(len == 1){
			//single state - only silent to / from M everything else set to zero.... 
			col = model->hmms[i]->hmm_column[0];
			col->transition[MM] = prob2scaledprob(0.0f);
			col->transition[MI] = prob2scaledprob(0.0f);
			col->transition[MD] = prob2scaledprob(0.0f);
			col->transition[MSKIP] = prob2scaledprob(1.0f);
			
			col->transition[II] = prob2scaledprob(0.0f);
			col->transition[IM] = prob2scaledprob(0.0f);
			col->transition[ISKIP] = prob2scaledprob(0.0f);
			
			col->transition[DD] = prob2scaledprob(0.0);
			col->transition[DM] = prob2scaledprob(0.0);
		}else if(len == 2){
			
			//first column
			col = model->hmms[i]->hmm_column[0];
			col->transition[MM] = prob2scaledprob( 1.0 - base_error * indel_freq);
			col->transition[MI] = prob2scaledprob(base_error * indel_freq);
			
			col->transition[MD] =  prob2scaledprob(0.0);
			col->transition[MSKIP] = prob2scaledprob(0.0);
			
			col->transition[II] = prob2scaledprob(1.0 - 0.999);
			col->transition[IM] = prob2scaledprob(0.999);
			col->transition[ISKIP] = prob2scaledprob(0.0f);
			
			col->transition[DD] = prob2scaledprob(0.0);
			col->transition[DM] = prob2scaledprob(0.0);
			
			//second column
			col = model->hmms[i]->hmm_column[1];
			col->transition[MM] = prob2scaledprob(0.0f );// 1.0 - base_error * indel_freq);
			col->transition[MI] = prob2scaledprob(0.0f );//(base_error * indel_freq) +  prob2scaledprob(0.5);
			col->transition[MD] = prob2scaledprob(0.0f );//(base_error * indel_freq) +  prob2scaledprob(0.5);
			col->transition[MSKIP] = prob2scaledprob(1.0);
			
			col->transition[II] = prob2scaledprob(0.00);
			col->transition[IM] = prob2scaledprob(0.0);
			col->transition[ISKIP] = prob2scaledprob(0.0f);
			
			col->transition[DD] = prob2scaledprob(0.0f );//(1.0 - 0.999);
			col->transition[DM] = prob2scaledprob(0.0f );//0.999);
			
		}else{
			//first column....
			col = model->hmms[i]->hmm_column[0];
			col->transition[MM] = prob2scaledprob( 1.0 - base_error * indel_freq);
			col->transition[MI] = prob2scaledprob(base_error * indel_freq) +  prob2scaledprob(0.5);
			col->transition[MD] = prob2scaledprob(base_error * indel_freq) +  prob2scaledprob(0.5);
			col->transition[MSKIP] = prob2scaledprob(0.0);
			
			col->transition[II] = prob2scaledprob(1.0 - 0.999);
			col->transition[IM] = prob2scaledprob(0.999);
			col->transition[ISKIP] = prob2scaledprob(0.0f);
			
			col->transition[DD] = prob2scaledprob(0.0);
			col->transition[DM] = prob2scaledprob(0.0);
		
			//middle columns...
			for(j = 1; j < len-2;j++){
				col = model->hmms[i]->hmm_column[j];

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
			
			//second last...
			col = model->hmms[i]->hmm_column[len -2];
			col->transition[MM] = prob2scaledprob( 1.0 - base_error * indel_freq);
			col->transition[MI] = prob2scaledprob(base_error * indel_freq);
			
			col->transition[MD] =  prob2scaledprob(0.0);

			col->transition[MSKIP] = prob2scaledprob(0.0);
			
			col->transition[II] = prob2scaledprob(1.0 - 0.999);
			col->transition[IM] = prob2scaledprob(0.999);
			col->transition[ISKIP] = prob2scaledprob(0.0f);
			
			col->transition[DD] = prob2scaledprob(0.0);
			col->transition[DM] = prob2scaledprob(1.0);
			//col->transition[DD] = prob2scaledprob(0.0);
			//col->transition[DM] = prob2scaledprob(0.0);

			col = model->hmms[i]->hmm_column[len -1];

			col->transition[MM] = prob2scaledprob(0.0f );// 1.0 - base_error * indel_freq);
			col->transition[MI] = prob2scaledprob(0.0f );//(base_error * indel_freq) +  prob2scaledprob(0.5);
			col->transition[MD] = prob2scaledprob(0.0f );//(base_error * indel_freq) +  prob2scaledprob(0.5);
			col->transition[MSKIP] = prob2scaledprob(1.0);
			
			col->transition[II] = prob2scaledprob(0.00);
			col->transition[IM] = prob2scaledprob(0.0);
			col->transition[ISKIP] = prob2scaledprob(0.0f);
			
			col->transition[DD] = prob2scaledprob(0.0f );//(1.0 - 0.999);
			col->transition[DM] = prob2scaledprob(0.0f );//0.999);
			
						
			
			
			
		}*/
		
	}
	
	// init all probs to 0
	
	for(i = 0 ; i < model->num_hmms;i++){
		len = model->hmms[i]->num_columns;
		for(j = 0; j < len;j++){
			
			model->silent_to_M[i][j] = prob2scaledprob(0.0f);
			model->silent_to_I[i][j] = prob2scaledprob(0.0f);
			model->silent_to_M_e[i][j] = prob2scaledprob(0.0f);
			model->silent_to_I_e[i][j] = prob2scaledprob(0.0f);
		}

		
	}
	model->skip = prob2scaledprob(0.0f);
	model->skip_e = prob2scaledprob(0.0f);
	
	if(rs->type[key] == 'B'){// barcodes all have same length & equal prior probability... 
		for(i = 0 ; i < model->num_hmms;i++){
			model->silent_to_M[i][0] = prob2scaledprob(1.0 / (float) model->num_hmms);// + prob2scaledprob(0.9);
			//model->M_to_silent[i] = prob2scaledprob(1.0);
			
			model->silent_to_I[i][0] = prob2scaledprob(0.0f);
			//model->I_to_silent[i] = prob2scaledprob(0.0f);
			
			
		}
		model->skip = prob2scaledprob(0.0);
	}
	
	if(rs->type[key] == 'F'){// fingerprint all have same length & equal prior probability... (of course we specify 1 with NNNNNNNN
		for(i = 0 ; i < model->num_hmms;i++){
			model->silent_to_M[i][0] = prob2scaledprob(1.0 / (float) model->num_hmms);
			//model->M_to_silent[i] = prob2scaledprob(1.0);
		}
		model->skip = prob2scaledprob(0.0);
	}
	
	
	if(rs->type[key] == 'S'){// fingerprint all have same length & equal prior probability... (of course we specify 1 with NNNNNNNN
		for(i = 0 ; i < model->num_hmms;i++){
			model->silent_to_M[i][0] = prob2scaledprob(1.0 / (float) model->num_hmms);
			//model->M_to_silent[i] = prob2scaledprob(1.0);
		}
		model->skip = prob2scaledprob(0.0);
	}
	
	if(rs->type[key] == 'P'){// Partial - can skip and exit at every M / I state.... 
		len = model->hmms[0]->num_columns;
		for(i = 0 ; i < model->num_hmms;i++){
			model->silent_to_M[i][0] = prob2scaledprob(1.0 / (float) model->num_hmms) + prob2scaledprob(1.0 - 0.01);
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
			
			model->silent_to_I[i][0] = prob2scaledprob(1.0 / (float) model->num_hmms) + prob2scaledprob(0.5);
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
			
			model->silent_to_I[i][0] = prob2scaledprob(0.8935878);
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
			model->silent_to_I[i][0] = prob2scaledprob(1.0 / (float) model->num_hmms);
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


/** \fn void print_model(struct model* model)
 
 \brief Prints out all parameters of a HMM. 
 
 \param model The model to be printed out.

 */


void print_model(struct model* model)
{
	int i,j,c;
	int len;
	float sum = 0;
	struct hmm_column* col =0;
	//fprintf(stderr,"Skip:%f Self:%f Next:%f\n",scaledprob2prob(model->skip), scaledprob2prob(model->random_self) , scaledprob2prob(model->random_next));
	for(i = 0; i < model->num_hmms;i++){
		fprintf(stderr,"HMM%d:		(get there:%f)\n",i, scaledprob2prob(model->silent_to_M[i][0]));
		
		
		
		len = model->hmms[i]->num_columns;
		//tmp = rs->sequence_matrix[key][i];
		for(j = 0; j < len;j++){
			col = model->hmms[i]->hmm_column[j];
			fprintf(stderr,"\t%d\t",j);
			sum = 0.0;
			for(c = 0; c < 5;c++){
				sum += scaledprob2prob(col->m_emit[c]);
			//	fprintf(stderr,"%0.2f ", scaledprob2prob(col->m_emit[c]));
				
			}
			fprintf(stderr,"%0.1fc ",sum);
			
			sum = 0.0;
			for(c = 0; c < 5;c++){
				sum += scaledprob2prob(col->i_emit[c]);
			//	fprintf(stderr,"%0.2f ", scaledprob2prob(col->i_emit[c]));
				
			}
			fprintf(stderr,"%0.1fc ",sum);
			
			for(c = 0; c < 9;c++){
				fprintf(stderr,"%0.4f ", scaledprob2prob(col->transition[c]));

			}
			fprintf(stderr,"%0.4f %0.4f %0.4f\n",scaledprob2prob(col->transition[MSKIP])+scaledprob2prob(col->transition[MM])+scaledprob2prob(col->transition[MI])+scaledprob2prob(col->transition[MD])  ,scaledprob2prob(col->transition[ISKIP])+scaledprob2prob(col->transition[II])+scaledprob2prob(col->transition[IM]),  scaledprob2prob(col->transition[DD]) + scaledprob2prob(col->transition[DM]) );
		}
	}
	/*fprintf(stderr," ESTIMATED::::: \n");
	for(i = 0; i < model->num_hmms;i++){
		fprintf(stderr,"HMM%d:		(get there:%f)\n",i, scaledprob2prob(model->silent_to_M[i][0]));
		
		
		
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
	}*/
	fprintf(stderr,"Links:silent to\n");
	double sumM = prob2scaledprob( 0.0);
	double sumI =  prob2scaledprob( 0.0);

	for(i = 0; i < model->num_hmms;i++){
		len = model->hmms[i]->num_columns;
		//tmp = rs->sequence_matrix[key][i];
		for(j = 0; j < len;j++){
			fprintf(stderr,"%d	%f	%f\n",i, scaledprob2prob(  model->silent_to_M[i][j]), scaledprob2prob(  model->silent_to_I[i][j]));
			sumM = logsum(sumM, model->silent_to_M[i][j]);
			sumI = logsum(sumI, model->silent_to_I[i][j]);
		}
		fprintf(stderr,"SANITY: %f\t%f\n",scaledprob2prob(sumM),scaledprob2prob(sumI));
		//fprintf(stderr,"%d	%f	%f	%f	%f\n",i, scaledprob2prob(  model->silent_to_M[i][0]), scaledprob2prob(  model->silent_to_I[i][0]),scaledprob2prob(   model->silent_to_M_e[i][0]),scaledprob2prob(  model->silent_to_I_e[i][0]));
	}
	fprintf(stderr,"Links:to silent \n");
	
	
	fprintf(stderr,"SKIP:\n");
	fprintf(stderr,"%f	%f\n", scaledprob2prob(model->skip) , scaledprob2prob(model->skip_e));
	

	
}


/** \fn void free_model(struct model* model)
 
 \brief Frees model.
 
 \param model The model to be freed .
 
 */



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
	
	//assert(model->hmms !=0);
	

	for(i = 0; i < model->num_hmms;i++){
		free(model->hmms[i]);// = malloc(sizeof(struct hmm) );
		//assert(model->hmms[i]  != 0);
	}
	free(model->hmms);// = malloc(sizeof(struct hmm*) * (1+ number_sub_models));
	if(model->silent_to_M){
		for(i = 0; i < model->num_hmms;i++){
			free(model->silent_to_M[i]);
			free(model->silent_to_M_e[i]);
			free(model->silent_to_I[i]);
			free(model->silent_to_I_e[i]);
		}
		free(model->silent_to_M);
		free(model->silent_to_M_e);
		free(model->silent_to_I);
		free(model->silent_to_I_e);
	}
	
	
	free(model);// = malloc(sizeof(struct model));
	//assert(model != 0);

	
	//return model;
}




/** \fn struct model_bag* copy_model_bag(struct model_bag* org)
 
 \brief Copies HMM into new HMM. 
 Used to copy HMM for each thread. 
 
 \param org  The @ref model_bag to be copied.
 
 */



struct model_bag* copy_model_bag(struct model_bag* org)
{
	struct model_bag* copy = 0;
	int i,j;
	copy =  malloc(sizeof(struct model_bag));
	
	assert(copy!=0);
	
	copy->model = malloc(sizeof(struct model* ) * org->num_models);//   param->read_structure->num_segments);
	
	
	
	
	assert(copy->model);
	
	copy->random_scores = malloc(sizeof(double) * org->num_random_scores);
	assert(copy-> random_scores);

	for(i = 0; i < org->num_random_scores;i++){
		copy->random_scores[i] = org->random_scores[i];
	}
	
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


/** \fn struct model* copy_model_parameters(struct model* org, struct model* copy )
 
 \brief Copies HMM parametes into new HMM.

 \param org  The original @ref model.
  \param copy  The copy @ref model to be copied.
 
 */

struct model* copy_model_parameters(struct model* org, struct model* copy )
{
	int i,j,c;
	
	struct hmm_column* org_col = 0;
	struct hmm_column* copy_col = 0;
	for(i = 0; i < 5;i++){
		copy->background_nuc_frequency[i] = org->background_nuc_frequency[i];
	}
	
	for(i = 0; i < org->num_hmms;i++){
		
		
		for(j = 0; j < org->hmms[i]->num_columns;j++){
			org_col = org->hmms[i]->hmm_column[j];
			copy_col = copy->hmms[i]->hmm_column[j];
				      
				      copy->silent_to_I[i][j] = org->silent_to_I[i][j];
				      copy->silent_to_I_e[i][j] = org->silent_to_I_e[i][j];
				      copy->silent_to_M[i][j] = org->silent_to_M[i][j];
				      copy->silent_to_M_e[i][j] = org->silent_to_M_e[i][j];
				      
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




/** \fn struct model_bag* set_model_e_to_laplace(struct model_bag* mb)
 
 \brief Sets all estiamted probabilities to 1.
 
 \param mb  The model.
 
 */

struct model_bag* set_model_e_to_laplace(struct model_bag* mb)
{
	int g,i,j,c;
	
	struct model* m = 0;
	struct hmm_column* hmm_column = 0;
	
	
	for(g = 0;g < mb->num_models;g++){
		m = mb->model[g];
	
		
		for(i = 0; i < m->num_hmms;i++){
			
			
			for(j = 0; j < m->hmms[i]->num_columns;j++){
				hmm_column = m->hmms[i]->hmm_column[j];
				//copy->silent_to_I[i][j] = org->silent_to_I[i][j];
				if(m->silent_to_I[i][j] != prob2scaledprob(0.0)){
					m->silent_to_I_e[i][j] = prob2scaledprob( 1.0f); //org->silent_to_I_e[i][j];
				}
				//copy->silent_to_M[i][j] = org->silent_to_M[i][j];
				if(m->silent_to_M[i][j] != prob2scaledprob(0.0)){
					m->silent_to_M_e[i][j] = prob2scaledprob( 1.0f); //org->silent_to_M_e[i][j];
				}
				
				for(c = 0; c< 5;c++){
					//copy_col->i_emit[c] = org_col->i_emit[c];
					hmm_column->i_emit_e[c] = prob2scaledprob(1.0f);
					
					//copy_col->m_emit[c] = org_col->m_emit[c];
					hmm_column->m_emit_e[c] = prob2scaledprob(1.0f);
				}
				
				for(c = 0; c < 9;c++){
					//c/opy_col->transition[c] = org_col->transition[c];
					if(hmm_column->transition[c] != prob2scaledprob(0.0)){
						hmm_column->transition_e[c] = prob2scaledprob(1.0f);
					}
				}
			}
		}
		//copy->skip = org->skip;
		if(m->skip != prob2scaledprob(0.0)){
			m->skip_e = prob2scaledprob( 1.0f);//?//?org->skip_e;
		}
		//copy->num_hmms = org->num_hmms;
	}
	return mb;
}




/** \fn struct model* reestimate(struct model* m, int mode)
 
 \brief Sets new model parameters based on estiamted probabilities from data.
 
 \param mb  The model.
  \param mode Sets which parameters should be optimized. 
 
 */


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
			for(j = 0; j < m->hmms[i]->num_columns;j++){
				      sum = logsum(sum, logsum(m->silent_to_I_e[i][j],prob2scaledprob(1.0)));
				      sum = logsum(sum, logsum(m->silent_to_M_e[i][j] , prob2scaledprob(1.0)) );
			}
			//fprintf(stderr," silent to I: %f",m->silent_to_I_e[i]);
			//fprintf(stderr," silent to M: %f",m->silent_to_M_e[i]);
		}
		sum = logsum(sum,logsum( m->skip_e , prob2scaledprob(1.0)));
		//fprintf(stderr,"estimated skip: %f\n", m->skip_e);
		
		for(i = 0; i < m->num_hmms;i++){
			for(j = 0; j < m->hmms[i]->num_columns;j++){
				m->silent_to_I[i][j]  =  logsum(m->silent_to_I_e[i][j] ,prob2scaledprob(1.0)) - sum;
				m->silent_to_M[i][j]  = logsum(m->silent_to_M_e[i][j],prob2scaledprob(1.0)) - sum;
			}
		}
		
		m->skip = logsum(m->skip_e ,prob2scaledprob(1.0)) - sum;
		
		//fprintf(stderr,"SKIP: %f\n", m->skip );
		
		//clear counts....
		for(i = 0; i < m->num_hmms;i++){
			for(j = 0; j < m->hmms[i]->num_columns;j++){
				      m->silent_to_I_e[i][j] = prob2scaledprob(0.0);
				      m->silent_to_M_e[i][j] = prob2scaledprob(0.0);
			}
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
				sum = logsum(sum, m_col->i_emit_e[c]);
			}
			
			for(c = 0; c< 5;c++){
				m_col->i_emit[c] = m_col->i_emit_e[c] - sum;
			}
			
			for(c = 0; c< 5;c++){
				m_col->i_emit_e[c] =  prob2scaledprob(0.0f);
			}
			
			
			sum = prob2scaledprob(0.0f);
			
			for(c = 0; c< 5;c++){
				sum = logsum(sum, m_col->m_emit_e[c] );
			}
			
			for(c = 0; c< 5;c++){
				m_col->m_emit[c] = m_col->m_emit_e[c] - sum;
			}
			
			for(c = 0; c< 5;c++){
				m_col->m_emit_e[c] =  prob2scaledprob(0.0f);
			}

			if(mode < 2){
				// internal hmm states...
				if(j != m->hmms[0]->num_columns-1){
					sum = prob2scaledprob(0.0f);
					
					//detect &assign 
					if(m_col->transition[MM] != prob2scaledprob(0.0)){
						sum = logsum(sum, m_col->transition_e[MM]);
					}
					if(m_col->transition[MI] != prob2scaledprob(0.0)){
						sum = logsum(sum, m_col->transition_e[MI]);
					}
					if(m_col->transition[MD] != prob2scaledprob(0.0)){
						sum = logsum(sum, m_col->transition_e[MD]);
					}
					if(m_col->transition[MSKIP] != prob2scaledprob(0.0)){
						sum = logsum(sum, m_col->transition_e[MSKIP]);
					}
					
					//set new prob... 
					if(m_col->transition[MM] != prob2scaledprob(0.0)){
						m_col->transition[MM] = m_col->transition_e[MM]  - sum;
					}
					if(m_col->transition[MI] != prob2scaledprob(0.0)){
						m_col->transition[MI] = m_col->transition_e[MI]  - sum;

						//sum = logsum(sum, m_col->transition_e[MI]);
					}
					if(m_col->transition[MD] != prob2scaledprob(0.0)){
						m_col->transition[MD] = m_col->transition_e[MD]  - sum;
						//sum = logsum(sum, m_col->transition_e[MD]);
					}
					if(m_col->transition[MSKIP] != prob2scaledprob(0.0)){
						m_col->transition[MSKIP] = m_col->transition_e[MSKIP]  - sum;
						//sum = logsum(sum, m_col->transition_e[MSKIP]);
					}
					
					
					
					
					/*sum = logsum(sum, logsum(m_col->transition_e[MM] , prob2scaledprob(1.0)));
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
					*/
					
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




/** \fn struct model* copy_estimated_parameter(struct model* target, struct model* source )
 
 \brief Sums eatimated parameters from different models. 
 Used to merge results from multiple threads.
 
 \param target  Model to hold the sums. 
 \param source Model from which to copy parameters.
 
 */

struct model* copy_estimated_parameter(struct model* target, struct model* source )
{
	int i,j,c;
	
	struct hmm_column* target_col = 0;
	struct hmm_column* source_col = 0;
	
	
	
	
	for(i = 0; i < target->num_hmms;i++){
				
		//copy->I_to_silent[i] = org->I_to_silent[i];
		//target->I_to_silent_e[i] = logsum(target->I_to_silent_e[i] , source->I_to_silent_e[i]);//   org->I_to_silent_e[i];
		//copy->M_to_silent[i] = org->M_to_silent[i];
		//target->M_to_silent_e[i] = logsum(target->M_to_silent_e[i],source->M_to_silent_e[i]); // org->M_to_silent_e[i];
		
		for(j = 0; j < target->hmms[0]->num_columns;j++){
				      
			//copy->silent_to_I[i] = org->silent_to_I[i];
			target->silent_to_I_e[i][j] = logsum(target->silent_to_I_e[i][j], source->silent_to_I_e[i][j]);
			//copy->silent_to_M[i] = org->silent_to_M[i];
			target->silent_to_M_e[i][j] = logsum(target->silent_to_M_e[i][j] ,source->silent_to_M_e[i][j]);//  org->silent_to_M_e[i];

				      
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





/** \fn struct model_bag* init_model_bag(struct parameters* param,double* back)
 
 \brief Initializes whole HMM model based on user input.
 
  
 \param param @ref parameters   Model to hold the sums.
 \param back Backgound nucleotide probabilities. 
 
 */

struct model_bag* init_model_bag(struct parameters* param,double* back)
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
	
	mb->random_scores = malloc(sizeof(double) * param->num_query);
	mb->num_random_scores = param->num_query;
	assert(mb->random_scores);
	
	for(i= 0;i < mb->num_random_scores;i++){
		mb->random_scores[i] = 0.0f;
	}
	
	mb->f_score = prob2scaledprob(0.0f);
	mb->b_score = prob2scaledprob(0.0f);
	mb->num_models = param->read_structure->num_segments;
	// get read length estimate...
	read_length = param->average_read_length;
	//fprintf(stderr,"READlength: %d\n",read_length);
	for(i = 0; i < mb->num_models;i++){
		//mb->model[i] = malloc_model_according_to_read_structure(param->read_structure,i);
		//fprintf(stderr," %d\n",read_length );
		if(param->read_structure->type[i] == 'G'){
			read_length = read_length -2;
		}else if(param->read_structure->type[i] == 'R'){
		}else if(param->read_structure->type[i] == 'P'){
			read_length = read_length - (int)strlen(param->read_structure->sequence_matrix[i][0])/2; // Initial guess - we don't know how much of the linker is present at this stage
		}else{
		//	fprintf(stderr,"%s : %d \n", param->read_structure->sequence_matrix[i][0], (int)strlen(param->read_structure->sequence_matrix[i][0]));
			read_length = read_length - (int)strlen(param->read_structure->sequence_matrix[i][0]);
		}
		
	//	fprintf(stderr,"READlength: %d\n",read_length);
		
	}
	//fprintf(stderr,"READlength: %d\n",read_length);
	
	if(read_length < 20){
		read_length = 20; // the expected read length should never be lower than 20!
	}
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
	
	mb->model_multiplier = 1.0f;
	
	c = 0;
	for(i = 0; i < mb->num_models ;i++){
		mb->model_multiplier  *= mb->model[i]->num_hmms;
		for(j = 0; j < mb->model[i]->num_hmms;j++){
			mb->label[c] = (j << 16) | i ;
			if(mb->model[i]->skip != prob2scaledprob(0.0)){
				mb->label[c]  |= 0x80000000;
			}
			//fprintf(stderr,"%d %d	%d %d\n",c,mb->label[c],mb->label[c] & 0xFFFF, (mb->label[c] >> 16) & 0x7FFF);
			c++;
			
		}
	}
	
	mb->model_multiplier = prob2scaledprob(mb->model_multiplier);
	
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



/** \fn void free_model_bag(struct model_bag* mb)
 
 \brief Frees whole HMM model.
 
 
 \param mb @ref model_bag. 
 
 */

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

	free(mb->random_scores);
	
	free(mb->model);// = malloc(sizeof(struct model* ) * param->read_structure->num_segments);
	
	
	free(mb);// = malloc(sizeof(struct model_bag));
}


















