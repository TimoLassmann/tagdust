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

#ifndef MMALLOC
#include "malloc_macro.h"
#endif


#define MAX_NUM_ARCH 100


struct parameters* test_architectures(struct parameters* param, int file_num)
{
	struct read_info** ri = 0;
	FILE* inarch = 0;
	int index = 0;
	int i,j,c,architecture_found;
	char line[MAX_LINE];
	char* tmp = 0;
	
	MMALLOC(tmp, sizeof(char) * 100);
	struct model_bag* mb = 0;
	
	
	int best_architecture = -1;
	float best_score = -1.0;
	
	init_logsum();
	
#if DEBUG
	param->num_query = 11;
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
	
	sprintf(param->buffer,"Searching for best architecture in file '%s'\n", param->arch_file);
	param->messages = append_message(param->messages, param->buffer);
	if (!(inarch = fopen( param->arch_file, "r" ))){
		sprintf(param->buffer,"Cannot open file '%s'\n", param->arch_file);
		param->messages = append_message(param->messages, param->buffer);
		free_param(param);
		exit(EXIT_FAILURE);
	}
	while(fgets(line, MAX_LINE, inarch)){
		
		if( byg_end("tagdust",line)){
			
			free_read_structure(param->read_structure);
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
					param->read_structure = assign_segment_sequences(param->read_structure, tmp , c );
				}else{
					if(!c){
						architecture_found = 0;
					}
				}
			}
			if(architecture_found){
				
				ssi = get_sequence_stats(param, ri, 0 );
				if(QC_read_structure(param)){
					free_param(param);
					exit(EXIT_FAILURE);
				}
				if(! ab->num_arch){
					mb = init_model_bag(param, ssi);
				}
				MMALLOC(ab->command_line[ab->num_arch],sizeof(char) * (strlen(line)+2));
				strcpy(ab->command_line[ab->num_arch] , line);
				ab->archs[ab->num_arch] = init_model_bag(param, ssi);
				ab->num_arch++;
				if(ab->num_arch == MAX_NUM_ARCH){
					fprintf(stderr,"Error - your architechture file has too many architectures. Currently only %d allowed.\n", MAX_NUM_ARCH);
					free_param(param);
					exit(EXIT_FAILURE);
				}
				MFREE (ssi);
			}
		}
	}
	
	
	fclose(inarch);
	
	if(!ab->num_arch){
		fprintf(stderr,"Error - could not find any architectures in file: %s\n", param->arch_file);
		free_param(param);
		exit(EXIT_FAILURE);

	}
	
	
	if(ab->num_arch > 1){
		
		// clean up posteriors!!
		
		for(i= 0; i < ab->num_arch;i++){
			ab->arch_posterior[i] = prob2scaledprob(0.0);
		}
		
		
		
		//2) run models and calculate P(M1|x)  = P(x | M1) / sum(P(Mj| x) for all models j
		int (*fp)(struct read_info** ,struct parameters*,FILE* ) = 0;
		FILE* file = NULL;
		int numseq;
		file =  io_handler(file, file_num,param);
		
		if(param->sam == 0){
			fp = &read_fasta_fastq;
		}else {
			fp = &read_sam_chunk;
		}
		
		
		numseq = fp(ri, param,file);
		
		mb =  run_pHMM(ab,mb,ri,param,0,numseq,MODE_ARCH_COMP);
		
		pclose(file);
		
		//3) print models and scores...
		
		float sum = prob2scaledprob(0.0f);
		for(i = 0; i < ab->num_arch;i++){
			sum = logsum(sum, ab->arch_posterior[i]);
		}
		best_architecture = -1;
		best_score = -1.0;

		for(i = 0; i < ab->num_arch;i++){
			ab->arch_posterior[i] = scaledprob2prob(ab->arch_posterior[i] - sum);
			if(ab->arch_posterior[i]  > best_score){
				best_score =ab->arch_posterior[i] ;
				best_architecture = i;
			}
		}
		
		for(i = 0; i < ab->num_arch;i++){
			if(i == best_architecture){
				//fprintf(stderr,"BEST%d:	%f	%s", i, ab->arch_posterior[i] ,ab->command_line[i]);
				sprintf(param->buffer,"Using: %s", ab->command_line[i]);
				param->messages = append_message(param->messages, param->buffer);
				sprintf(param->buffer,"%0.2f Confidence.\n", ab->arch_posterior[i]);
				param->messages = append_message(param->messages, param->buffer);
			}//else{
				//fprintf(stderr,"Arch%d:	%f	%s", i, ab->arch_posterior[i] ,ab->command_line[i]);
			//}
		}
	}else{
		best_architecture = 0;
		sprintf(param->buffer,"Using: %s", ab->command_line[0]);
		param->messages = append_message(param->messages, param->buffer);
		sprintf(param->buffer,"Confidence: %0.2f\n", 1.0);
		param->messages = append_message(param->messages, param->buffer);

	}
	//4) set param->read_structure to best arch...
	
	free_read_structure(param->read_structure);
	param->read_structure=malloc_read_structure();
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
			param->read_structure = assign_segment_sequences(param->read_structure, tmp , c );
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
	
	MFREE(ab->arch_posterior);// = malloc(sizeof(float) * MAX_NUM_ARCH);
	
	MFREE(ab->archs);
	MFREE(ab);
	MFREE(tmp);
	return param;
}







