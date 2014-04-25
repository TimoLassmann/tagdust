#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <getopt.h>

#include "misc.h"
#include "interface.h"
#include "io.h"

#if HAVE_CONFIG_H
#include "config.h"
#endif

struct libraries{
	int total;
	int* count_in_file;
	int lib_id;
};


int main (int argc,char * argv[])
{
	struct parameters* param = 0;
	
	struct read_info** ri =0 ;
	
	struct libraries** libs = 0;
	
	FILE* file = 0;
	
	int i,j,c;
	int numseq;
	
	param = interface(param,argc,argv);
	
	int (*fp)(struct read_info** ,struct parameters*,FILE* ) = 0;
	fp = &read_fasta_fastq;
	
	libs = malloc(sizeof(struct libraries*) *  100);
	
	for(i = 0; i < 100; i++){
		libs[i] = malloc(sizeof(struct libraries));
		libs[i]->count_in_file = malloc(sizeof(int)* param->infiles);
		for(j = 0; j < param->infiles;j++){
			libs[i]->count_in_file[j] = 0;
		}
		libs[i]->lib_id= i;
		libs[i]->total = 0;
	}
		
	ri = malloc_read_info(ri, 1000000);
	
	for(i = 0; i < param->infiles;i++){
	
		file =  io_handler(file, i,param);
		
		while ((numseq = fp(ri, param,file)) != 0){
			for(j = 0;j < numseq;j++){
				//read all reads into librarues ..
				c = byg_end("BARNUM:", ri[j]->name);
				if(c){
					libs[atoi(ri[j]->name + c )]->count_in_file[i]++;
					libs[atoi(ri[j]->name + c )]->total++;
				}
			}
		}
		
		if(param->sam == 2 || param->sam == 1 || param->gzipped ){
			pclose(file);
		}else{
			fclose(file);
		}
	}
	
	free_read_info(ri, 1000000);
		
	int* lib_file_assignment = 0;
	
	lib_file_assignment = malloc(sizeof(int) * param->infiles);
	
	for(i = 0; i < param->infiles;i++){
		lib_file_assignment[i] = -1;
	}
	
	// look for which file contains the most reads of lib X
	// assign file to that library
	// check if there is a one to many assignment - > quit if so.
	int max;
	for(i = 0; i < param->infiles; i++){
		max = -1;
		fprintf(stderr,"lib:%d\t",i);
		for(j = 0;j < param->infiles;j++){
			fprintf(stderr,"%d\t", libs[i]->count_in_file[j] );
			if(libs[i]->count_in_file[j] > max){
				max =libs[i]->count_in_file[j];
				lib_file_assignment[i]  = j;
			}
		}
		fprintf(stderr,"%d\n",   libs[i]->total );
	}
	
	for(i = 0; i < param->infiles; i++){
		for(j = i+1; j < param->infiles; j++){
			if(lib_file_assignment[i] == lib_file_assignment[j]){
				sprintf(param->buffer,"Cannot determine which file belongs to which library...\n");
				param->messages = append_message(param->messages, param->buffer);
				free_param(param);
				exit(EXIT_FAILURE);
			}
		}
		fprintf(stderr,"lib:%d -> file:%d\n",i,lib_file_assignment[i] );
	}
	
	double TP,FP,FN, TN,sensitivity,specificity,precision,kappa,P_e,P_o,sum;
	
	TP = 0.0;
	FP = 0.0;
	FN = 0.0;
	TN = 0.0;
	sum = 0.0;
	
	for(i = 0; i < param->infiles; i++){
		for(j = 0;j < param->infiles;j++){
			sum += libs[i]->count_in_file[j];
			if(i == 0){
				if(lib_file_assignment[i] == j){
					TN += libs[i]->count_in_file[j];
				}else{
					FN +=libs[i]->count_in_file[j];
				}
			}else{
				if(lib_file_assignment[i] == j){
					TP += libs[i]->count_in_file[j];
				}else{
					FP +=libs[i]->count_in_file[j];
				}
			}
			
		}
	}
	
	precision = TP / (TP + FP);
	sensitivity = TP/( TP + FN );
	specificity =  TN / ( TN + FP);
	
	P_e = ((TP+FN) / (double)sum) * ((TP+FP) / (double)sum) +  ( ((FP+TN) / (double)sum  ) * ((FN+TN) / (double)sum));
	P_o =(TP+TN)/(double)sum ;
	
	kappa = (P_o - P_e) / (1.0 - P_e);
	
	fprintf(stderr,"Sen:%f	Spe:%f	Precision:%f	Kappa:%f	TP:%f	FP:%f	FN:%f	TN:%f\n",sensitivity,specificity,precision,kappa,TP,FP,FN,TN);

	for(i = 0; i < param->infiles; i++){
		free(libs[i]->count_in_file);
		free(libs[i]);
	}
	free(libs);
	return EXIT_SUCCESS;
}






