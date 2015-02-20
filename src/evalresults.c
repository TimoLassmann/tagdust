
#if HAVE_CONFIG_H
#include "config.h"
#endif
#include "kslib.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <getopt.h>
#include <ctype.h>

#include "misc.h"
#include "interface.h"
#include "io.h"
#include "nuc_code.h"

#include <sys/stat.h>



struct libraries{
	int total;
	int* count_in_file;
	int lib_id;
};


int main (int argc,char * argv[])
{
	struct parameters* param = NULL;
	int status;
	struct read_info** ri =NULL ;
	
	struct libraries** libs = NULL;
	
	FILE* file = 0;
	
	struct stat buf;
	
	int i,j,c,g;
	int numseq = 0;
	
	int max_num_lib_detected = -1;
	
	int org_read_len = 0;
	
	float num_extracted = 0;
	float average_error_in_extracted = 0;
	char* orgread = 0;
	
	
	MMALLOC(orgread,sizeof(char)* 300);
	
	init_nuc_code();
	
	param = interface(param,argc,argv);
	
	
	if(!param->format){
		sprintf(param->buffer,"Error: You need to specify the name of program with the -name option. \n");
		param->messages = append_message(param->messages, param->buffer);
		free_param(param);
		exit(EXIT_FAILURE);
	}
	
	param->num_query = 1000000;
	
	int (*fp)(struct read_info** ,struct parameters*,FILE* , int* buffer_count) = 0;
	
	fp = &read_fasta_fastq;
	
	MMALLOC(libs,sizeof(struct libraries*) *  100);
	
	for(i = 0; i < 100; i++){
		libs[i] = 0;
		MMALLOC(libs[i],sizeof(struct libraries));
		libs[i]->count_in_file = 0;
		MMALLOC(libs[i]->count_in_file,sizeof(int)* (param->infiles+5));
		for(j = 0; j < param->infiles+5;j++){
			libs[i]->count_in_file[j] = 0;
		}
		libs[i]->lib_id= i;
		libs[i]->total = 0;
	}
		
	ri = malloc_read_info(ri, param->num_query);
	
	for(i = 0; i < param->infiles;i++){
		file =  io_handler(file, i,param);
		while(1){
			if((status = fp(ri, param,file,&numseq)) != kslOK)  exit(status);
			if(!numseq){
				break;
			}
		//while ((numseq = fp(ri, param,file)) != 0){
			for(j = 0;j < numseq;j++){
				//read all reads into librarues ..
				c = byg_end("BARNUM:", ri[j]->name);
				if(c){
					if(atoi(ri[j]->name + c) > max_num_lib_detected ){
						max_num_lib_detected =atoi(ri[j]->name + c) ;
					}
					libs[atoi(ri[j]->name + c )]->count_in_file[i]++;
					libs[atoi(ri[j]->name + c )]->total++;
				}
			}
		}
		pclose(file);
		
		file = 0;
	}
		
	int* lib_file_assignment = 0;
	
	MMALLOC(lib_file_assignment,sizeof(int) * (max_num_lib_detected+1) );
	
	for(i = 0; i <= max_num_lib_detected ;i++){
		lib_file_assignment[i] = -1;
	}
	
	// look for which file contains the most reads of lib X
	// assign file to that library
	// check if there is a one to many assignment - > quit if so.
	int max;
	for(i = 0; i <= max_num_lib_detected; i++){
		max = 0;
		fprintf(stderr,"lib:%d\t",i);
		for(j = 0;j <= max_num_lib_detected;j++){
			fprintf(stderr,"%d\t", libs[i]->count_in_file[j] );
			if(libs[i]->count_in_file[j] > max){
				max =libs[i]->count_in_file[j];
				lib_file_assignment[i]  = j;
			}
		}
		fprintf(stderr,"%d\n",   libs[i]->total );
	}
	for(i = 0; i <= max_num_lib_detected; i++){
		fprintf(stderr,"lib:%d -> file:%d\n",i,lib_file_assignment[i] );
	}
	
	if(param->sim_numseq ){
		i = (int)((float) param->sim_numseq * (1.0-param->sim_random_frac)) + 1 ;
		j = param->sim_numseq - (int)((float) param->sim_numseq * (1.0-param->sim_random_frac)) - 1;
		if(libs[0]->total  < j){
			lib_file_assignment[0] = -1;
		}
	}
		
	//Sanity check for libraries 1 .... X (unmapped reads are excluded!)
	for(i = 0; i <= max_num_lib_detected; i++){
		for(j = i+1; j <= max_num_lib_detected; j++){
			if(lib_file_assignment[i] == lib_file_assignment[j]){
				sprintf(param->buffer,"Cannot determine which file belongs to which library...\n");
				param->messages = append_message(param->messages, param->buffer);
				free_param(param);
				exit(EXIT_FAILURE);
			}
		}
	}
	
	
	for(i = 0; i < param->infiles;i++){
		c = 0;
		for(j = 1; j <= max_num_lib_detected;j++){
			if(lib_file_assignment[j] == i){
				c = 1;
				break;
			}
		}
		if(c){
			file =  io_handler(file, i,param);
			
			//while ((numseq = fp(ri, param,file)) != 0){
			while(1){
				if((status = fp(ri, param,file,&numseq)) != kslOK)  exit(status);
				if(!numseq){
					break;
				}
				for(j = 0;j < numseq;j++){
					c = byg_end("SEQ:", ri[j]->name  );
					if(c){
						org_read_len = 0;
						for(g = c ;g < (int)strlen(ri[j]->name);g++){
							org_read_len++;
							if(isspace((int) ri[j]->name[g]) ||  ri[j]->name[g]  == ';'){
								break;
							}
							
						}
						//orgread = malloc(sizeof(unsigned char)* g);
						org_read_len = 0;
						for(g = c ;g < (int)strlen(ri[j]->name);g++){
							
							if(isspace((int) ri[j]->name[g]) ||  ri[j]->name[g]  == ';'){
								orgread[org_read_len] = 0;
								break;
							}
							
							orgread[org_read_len] = nuc_code[(int)ri[j]->name[g]];
							org_read_len++;
						}
					}
					c =  byg_count("READ", ri[j]->name);
					
					if(c){
						num_extracted++;
						
						
						if(ri[i]->len < org_read_len){
							c = bpm_check_error_global((unsigned char*)ri[j]->seq,(unsigned char*) orgread, ri[j]->len, org_read_len);
						}else{
							c = bpm_check_error_global((unsigned char*) orgread,(unsigned char*)ri[j]->seq,org_read_len, ri[j]->len);
						}
						/*fprintf(stderr,"%s\t%d\t%d\n", ri[j]->name,j,c);
						for(g = 0;g < ri[j]->len;g++){
							fprintf(stderr,"%d %d\n",orgread[g],ri[j]->seq[g] );
						}*/
						 
						 
						g = (org_read_len > ri[j]->len) ?org_read_len : ri[j]->len ;
						average_error_in_extracted += (double) c / (double)g;
					}
				}
			}
			
			pclose(file);
			
		}
	}

	free_read_info(ri, param->num_query);
	
	double TP,FP,FN, TN,sensitivity,specificity,precision,kappa,P_e,P_o,sum;
	
	TP = 0.0;
	FP = 0.0;
	FN = 0.0;
	TN = 0.0;
	sum = 0.0;
	
	for(i = 0; i <= max_num_lib_detected; i++){
		for(j = 0;j <= max_num_lib_detected;j++){
			sum += libs[i]->count_in_file[j];
			if(i == 0){
				if(lib_file_assignment[i] == j){
					TN += libs[i]->count_in_file[j];
				}else{
					FP += libs[i]->count_in_file[j];
				}
			}else{
				if(lib_file_assignment[i] == j){
					TP += libs[i]->count_in_file[j];
				}else{
					FP += libs[i]->count_in_file[j];
				}
			}
		}
	}
	
	if(param->sim_numseq ){
		i = (int)((float) param->sim_numseq * (1.0-param->sim_random_frac));
		j = param->sim_numseq - (int)((float) param->sim_numseq * (1.0-param->sim_random_frac)) ;
		//fprintf(stderr,"%d	%d	%d	%d\n",i,j,i+j,libs[0]->total);
		
		//TN =
		TN +=  (j - libs[0]->total); // total negatives deteced minus number detccted in files...
		sum += (j - libs[0]->total);
		
		//FN:
		FN +=param->sim_numseq - sum; // everything else...
		sum +=param->sim_numseq - sum;
		
		
	}

	precision = TP / (TP + FP);
	sensitivity = TP/( TP + FN );
	specificity =  TN / ( TN + FP);
	
	P_e = ((TP+FN) / (double)sum) * ((TP+FP) / (double)sum) +  ( ((FP+TN) / (double)sum  ) * ((FN+TN) / (double)sum));
	P_o =(TP+TN)/(double)sum ;
	
	kappa = (P_o - P_e) / (1.0 - P_e);
	
	param->buffer[0] = 0;
	sprintf(param->buffer, "%s_results.txt",param->outfile );
	if(!stat ( param->buffer, &buf )){
		if((file = fopen(param->buffer, "a")) == NULL) KSLIB_XEXCEPTION_SYS(kslEWRT,"Failed to open file:%s",param->buffer);
		//if ((file = fopen(param->buffer, "a")) == NULL){
		//	fprintf(stderr,"can't open output\n");
		//	exit(-1);
		//}
		
	}else {
		if((file = fopen(param->buffer, "w")) == NULL) KSLIB_XEXCEPTION_SYS(kslEWRT,"Failed to open file:%s",param->buffer);
		//if ((file = fopen(param->buffer, "w")) == NULL){
		//	fprintf(stderr,"can't open output\n");
		//	exit(-1);
		//}
		fprintf(file,"Program\tSensitivity\tSpecificity\tPrecision\tKappa\tAvgError\tTP\tFP\tFN\tTN\n");
	}
	
	
	
	fprintf(file,"%s\t%0.4f\t%0.4f\t%0.4f\t%0.4f\t%0.4f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\n",param->format,sensitivity,specificity,precision,kappa,average_error_in_extracted / num_extracted,TP,FP,FN,TN);

	
	//fprintf(stderr,"%f	%f\n",average_error_in_extracted,average_error_in_extracted / num_extracted);
	
	fclose(file);
	
	for(i = 0; i < 100; i++){
		MFREE(libs[i]->count_in_file);
		MFREE(libs[i]);
	}
	MFREE(libs);
	
	MFREE(lib_file_assignment);
	MFREE(orgread);
	free_param(param);
	
	return kslOK;
ERROR:
	return kslFAIL;
}






