#if HAVE_CONFIG_H
#include "config.h"
#endif
#include "kslib.h"

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

#define MAX_NUM_ARCH 100

int test_architectures(struct parameters* param, int file_num)
{
	struct read_info** ri = 0;
	int status;
	FILE* inarch = 0;
	int index = 0;
	int i,j,c,architecture_found;
	char line[MAX_LINE];
	char* tmp = 0;
	
	MMALLOC(tmp, sizeof(char) * MAX_LINE);
	struct model_bag* mb = 0;
	
	int best_architecture = -1;
	float best_score = -1.0;
	
	init_logsum();
	
#if DEBUG
	param->num_query = 100000;
#else
	param->num_query = 100000;
#endif
	
	ri = malloc_read_info(ri,param->num_query);
	
		
	struct sequence_stats_info* ssi = 0;
	
	struct arch_bag* ab = 0;
	MMALLOC(ab,sizeof(struct arch_bag) );
	ab->archs = 0;
	ab->arch_posterior = 0;
	ab->command_line = 0;
	MMALLOC(ab->archs,sizeof(struct model_bag*) * MAX_NUM_ARCH);
	MMALLOC(ab->arch_posterior,sizeof(float) * MAX_NUM_ARCH);
	MMALLOC(ab->command_line,sizeof(char*) *MAX_NUM_ARCH);
	
	for(i = 0; i < MAX_NUM_ARCH;i++){
		ab->command_line[i] = 0;
		ab->archs[i] = 0;
	}
	
	ab->num_arch = 0;
	//1) file 0 = architecture file (tagdust commands on  each line)
	//    file 1 = read file to test .
	sprintf(param->buffer,"Looking at file:%s\n", param->infile[file_num]);
	param->messages = append_message(param->messages, param->buffer);
	sprintf(param->buffer,"Searching for best architecture in file '%s'\n", param->arch_file);
	param->messages = append_message(param->messages, param->buffer);
	if((inarch = fopen(param->arch_file, "r")) == NULL) KSLIB_XEXCEPTION_SYS(kslEWRT,"Failed to open file:%s", param->arch_file);
	
	while(fgets(line, MAX_LINE, inarch)){
		
		if( byg_end("tagdust",line)){
			if(param->read_structure){
				free_read_structure(param->read_structure);
				param->read_structure = NULL;
			}
			param->read_structure=malloc_read_structure();
			architecture_found = 1;
			
			for(c = 0;c < 10;c++){
				
				tmp[0] = '-';
				tmp[1] = (char)(c+49);
				tmp[2] = 0;
				
				index = byg_end(tmp, line  );
				if(index){
					j = 0;
					while(isspace((int)line[index])){
						index++;
					}
					for(i = index ;i < MAX_LINE;i++){
						if(isspace((int)line[i])){
							tmp[j] = 0;
							j++;
							break;
						}
						tmp[j] = line[i];
						j++;
						
					}
					if((status = assign_segment_sequences(param, tmp , c  )) != kslOK) KSLIB_XEXCEPTION(kslFAIL,"Some problem with parsing an HMM segment: %s.\n",optarg);
					//param->read_structure = assign_segment_sequences(param->read_structure, tmp , c );
				}else{
					if(!c){
						architecture_found = 0;
					}
				}
			}
			if(architecture_found){
				
				ssi = get_sequence_stats(param, ri, file_num );
				
				if(QC_read_structure(param)){
					free_param(param);
					exit(EXIT_FAILURE);
				}
				if(! ab->num_arch){
					mb = init_model_bag(param, ssi);
				}
				MMALLOC(ab->command_line[ab->num_arch],sizeof(char) * (strlen(line)+2));
				strcpy(ab->command_line[ab->num_arch] , line);
				//fprintf(stderr,"TESTING:\n%s\n",ab->command_line[ab->num_arch]);
				ab->archs[ab->num_arch] = init_model_bag(param, ssi);
				ab->num_arch++;
				if(ab->num_arch == MAX_NUM_ARCH){
					sprintf(param->buffer,"Error - your architechture file has too many architectures. Currently only %d allowed.\n", MAX_NUM_ARCH);
					param->messages = append_message(param->messages, param->buffer);
					free_param(param);
					exit(EXIT_FAILURE);
				}
				MFREE (ssi);
			}
		}
	}
	
	fclose(inarch);
	
	if(!ab->num_arch){
		sprintf(param->buffer,"Error - could not find any architectures in file: %s\n", param->arch_file);
		param->messages = append_message(param->messages, param->buffer);
		free_param(param);
		exit(EXIT_FAILURE);

	}
	
	if(ab->num_arch > 1){
		for(i = 0 ;i < ab->num_arch;i++){
			for(j = i+1;j < ab->num_arch;j++){
				if(!strcmp(ab->command_line[i], ab->command_line[j])){
					sprintf(param->buffer,"ERROR: two architectures in %s are the same:%s\n%s\n", param->arch_file,ab->command_line[i], ab->command_line[j] );
					param->messages = append_message(param->messages, param->buffer);
					free_param(param);
					exit(EXIT_FAILURE);

				}
			}
		}
		// clean up posteriors!!
		
		for(i= 0; i < ab->num_arch;i++){
			ab->arch_posterior[i] = prob2scaledprob(1.0);
		}
		
		
		
		//2) run models and calculate P(M1|x)  = P(x | M1) / sum(P(Mj| x) for all models j
		int (*fp)(struct read_info** ,struct parameters*,FILE*, int* buffer_count) = 0;
		FILE* file = NULL;
		int numseq;
		file =  io_handler(file, file_num,param);
		
		if(param->sam == 0){
			fp = &read_fasta_fastq;
		}else {
			fp = &read_sam_chunk;
		}
		
		
		if((status = fp(ri, param,file,&numseq)) != kslOK) KSLIB_XFAIL(kslFAIL, param->errmsg,"Reading chunk from file %snumseq failed.\n", param->infile[file_num]);
		
		if((status =run_pHMM(ab,mb,ri,param,0,numseq,MODE_ARCH_COMP)) != kslOK) KSLIB_XFAIL(kslFAIL,param->errmsg,"run_pHMM failed\n");

		
		//mb =  run_pHMM(ab,mb,ri,param,0,numseq,MODE_ARCH_COMP);
		
		pclose(file);
		
		float sum = prob2scaledprob(0.0f);
		for(i = 0; i < ab->num_arch;i++){
		//	fprintf(stderr,"%d %s %f	%s", file_num, param->infile[file_num] ab->arch_posterior[i],ab->command_line[i]);
			sum = logsum(sum, ab->arch_posterior[i]);
		}
		best_architecture = -1;
		best_score = -1.0;

		for(i = 0; i < ab->num_arch;i++){
			ab->arch_posterior[i] = scaledprob2prob(ab->arch_posterior[i] - sum);
		//	fprintf(stderr,"%f	%s",  ab->arch_posterior[i],  ab->command_line[i]);
			if(ab->arch_posterior[i]  > best_score){
				best_score =ab->arch_posterior[i] ;
				best_architecture = i;
			}
		}
		
		for(i = 0; i < ab->num_arch;i++){
			if(i == best_architecture){
				param->buffer = pretty_print_selected_architecture(ab->command_line[i],param->buffer);
				param->messages = append_message(param->messages, param->buffer);
				sprintf(param->buffer,"%0.2f Confidence.\n", ab->arch_posterior[i]);
				param->messages = append_message(param->messages, param->buffer);
			}
		}
	}else{
		best_architecture = 0;
		param->buffer = pretty_print_selected_architecture(ab->command_line[0],param->buffer);
		param->messages = append_message(param->messages, param->buffer);
		sprintf(param->buffer,"Confidence: %0.2f\n", 1.0);
		param->messages = append_message(param->messages, param->buffer);
	}
	//4) set param->read_structure to best arch...
	if(param->read_structure){
		free_read_structure(param->read_structure);
		param->read_structure = NULL;
	}
	param->read_structure= malloc_read_structure();
	architecture_found = 1;
	
	for(c = 0;c < 10;c++){
		
		tmp[0] = '-';
		tmp[1] = (char)(c+49);
		tmp[2] = 0;
		
		index = byg_end(tmp, ab->command_line[best_architecture] );
		if(index){
			
			
			
			j = 0;
			while(isspace((int)ab->command_line[best_architecture][index])){
				index++;
			}
			for(i = index ;i < MAX_LINE;i++){
				if(isspace((int)ab->command_line[best_architecture][i])){
					tmp[j] = 0;
					j++;
					break;
				}
				tmp[j] = ab->command_line[best_architecture][i];
				j++;
			}
			
			if((status = assign_segment_sequences(param, tmp , c  )) != kslOK) KSLIB_XEXCEPTION(kslFAIL,"Some problem with parsing an HMM segment: %s.\n",optarg);
			
			//param->read_structure = assign_segment_sequences(param->read_structure, tmp , c );
		}else{
			if(!c){
				architecture_found = 0;
			}
		}
	}
	
	
	if(QC_read_structure(param)){
		free_param(param);
		exit(EXIT_FAILURE);
	}
	
	free_read_info(ri, param->num_query);
	
	
	free_model_bag(mb);

	for(i = 0; i < ab->num_arch;i++){
		free_model_bag(ab->archs[i]);
		MFREE(ab->command_line[i]);
	}
	MFREE(ab->command_line);
	MFREE(ab->arch_posterior);
	MFREE(ab->archs);
	MFREE(ab);
	MFREE(tmp);
	return kslOK;
ERROR:
	return kslFAIL;
}


char* pretty_print_selected_architecture(char* command_line, char* buffer)
{
	char* tmp = 0;
	int status;
	int i,j,c,index;
	
	MMALLOC(tmp, sizeof(char) * MAX_LINE);
	buffer[0] = 'U';
	buffer[1] = 's';
	buffer[2] = 'i';
	buffer[3] = 'n';
	buffer[4] = 'g';
	buffer[5] = ':';
	buffer[6] = ' ';
	c = 7;

	for(i = 0;i < 10;i++){
		
		tmp[0] = '-';
		tmp[1] = (char)(i+49);
		tmp[2] = 0;
		
		index = byg_end(tmp, command_line  );
		
		if(index){
			buffer[c] = tmp[0];
			c++;
			buffer[c] = tmp[1];
			c++;
			buffer[c] = ' ';
			c++;
			while(isspace((int)command_line[index])){
				index++;
			}
			for(j = index ;j < MAX_LINE;j++){
				if(isspace((int)command_line[j])){
					buffer[c] = ' ';
					c++;
					break;
				}
				buffer[c] =command_line[j];
				c++;
			}
		}
	}
	buffer[c] = '\n';
	c++;
	buffer[c] = 0;
	MFREE(tmp);
	return buffer;
ERROR:
	KSLIB_MESSAGE(status,"Something wrong in pretty_print_selected_architecture.\n");
	return NULL;
}


