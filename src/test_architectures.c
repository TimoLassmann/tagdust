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


void test_architectures(struct parameters* param, int file_num)
{
	struct read_info** ri = 0;
	FILE* inarch = 0;
	int index = 0;
	int i,j,c,architecture_found;
	char line[MAX_LINE];
	char* tmp = malloc(sizeof(char) * 100);
	struct model_bag* mb = 0;
	
	
	
	init_logsum();
	
	//double* back = 0;
	//int average_length = 0;
	
	//back = malloc(sizeof(double)*5);
	
#if DEBUG
	//printf("Debug\n");
	param->num_query = 11;
#else
	//printf("No Debug\n");
	param->num_query = 101;
#endif
	
	//param->num_query = 500000;
	
	
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
	
	
	struct sequence_stats_info* ssi = 0;
	
	
	
	struct arch_bag* ab = malloc(sizeof(struct arch_bag) );
	
	ab->archs = malloc(sizeof(struct model_bag*) * MAX_NUM_ARCH);
	ab->arch_posterior = malloc(sizeof(float) * MAX_NUM_ARCH);
	ab->num_arch = 0;
	
	//1) file 0 = architecture file (tagdust commands on  each line)
	//    file 1 = read file to test .
	
	
	
	
	
	if (!(inarch = fopen( param->arch_file, "r" ))){
		fprintf(stderr,"Cannot open file '%s'\n", param->arch_file);
		exit(EXIT_FAILURE);
	}
	while(fgets(line, MAX_LINE, inarch)){
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
			if(QC_read_structure(param->read_structure)){
				free_param(param);
				exit(EXIT_FAILURE);
			}
			if(! ab->num_arch){
				mb = init_model_bag(param, ssi);
			}
			
			ab->archs[ab->num_arch] = init_model_bag(param, ssi);
			ab->num_arch++;
			if(ab->num_arch == MAX_NUM_ARCH){
				fprintf(stderr,"Error - your architechture file has too many architectures. Currently only %d allowed.\n", MAX_NUM_ARCH);
				free_param(param);
				exit(EXIT_FAILURE);
			}
			free(ssi);
		}
		
	}
	
	
	fclose(inarch);
	
	
	// clean up posteriors!!
	
	for(i= 0; i < ab->num_arch;i++){
		ab->arch_posterior[i] = prob2scaledprob(0.0);
	}
	
	
	
	//2) run models and calculate P(M1|x)  = P(x | M1) / sum(P(Mj| x) for all models j
	int (*fp)(struct read_info** ,struct parameters*,FILE* ) = 0;
	FILE* file = NULL;
	int numseq;
	file =  io_handler(file, 0,param);
	
	if(param->sam == 0){
		fp = &read_fasta_fastq;
	}else {
		fp = &read_sam_chunk;
	}
	
	
	numseq = fp(ri, param,file);
	
	mb =  run_pHMM(ab,mb,ri,param,0,numseq,MODE_ARCH_COMP);
	
	
	float sum = prob2scaledprob(0.0f);
	for(i = 0; i < ab->num_arch;i++){
		sum = logsum(sum, ab->arch_posterior[i]);
		//fprintf(stderr,"%f\n", ab->arch_posterior[i]);
	}
	
	for(i = 0; i < ab->num_arch;i++){
		//sum = logsum(sum, ab->arch_posterior[i]);
		fprintf(stderr,"Arch%d:	%f	%f	%f\n", i, ab->arch_posterior[i] ,ab->arch_posterior[i]- sum,   scaledprob2prob(  ab->arch_posterior[i]- sum ));
	}
	
	
	
	
	//3) print models and scores...
	for(i = 0; i < param->num_query;i++){
		free(ri[i]->strand);
		free(ri[i]->hits);
		
		if(ri[i]->cigar){
			free(ri[i]->cigar);
		}
		if(ri[i]->labels){
			free(ri[i]->labels);
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
	//free(back);
	free(ri);
	
	free_model_bag(mb);

	for(i = 0; i < ab->num_arch;i++){
		free_model_bag(ab->archs[i]);
	}
	
	free(ab->arch_posterior);// = malloc(sizeof(float) * MAX_NUM_ARCH);
	
	free(ab->archs);
	free(ab);
	free(tmp);
	if(param->sam == 2 || param->sam == 1 || param->gzipped ){
		pclose(file);
	}else{
		fclose(file);
	}
	
}







