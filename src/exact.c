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

/*! \file exact.c
 \brief Functions to extract sequences without using errors.
 
 Appears to only use a single barcode sequence. 
 
 \deprecated For comparison fastx_barcode_splitter can be used.  
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

#include "exact.h"



/** \fn void exact_controller(struct parameters* param,int (*fp)(struct read_info** ,struct parameters*,FILE* ),int file_num)
 \brief Wrapper function for exact extraction of reads. 
 
 \param param @parameters.
 
  \param fp Pointer to function used to read sequences (either SAM or fastq).
 \param filenum Number of input files.
 
 */
void exact_controller(struct parameters* param,int (*fp)(struct read_info** ,struct parameters*,FILE* ),int file_num)
{
	struct read_info** ri = 0;
	
	FILE* outfile;
	int i,c;
	int numseq;
	int total_read = 0;
		
	init_logsum();
	init_nuc_code();
	
	char pattern[1000];
	float* back = 0;
	int average_length = 0;
	int plen = 0;
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
	plen = (int)strlen(param->exact5);
	for(i = 0; i < plen;i++ ){
		pattern[i] = nuc_code[(int)param->exact5[i]];
	}
	pattern[strlen(param->exact5)] = 0;
	
	int success = 0;
	int failure = 0;
	total_read = 0;
	
	if(param->outfile){
		if ((outfile = fopen( param->outfile, "w")) == NULL){
			fprintf(stderr,"can't open output\n");
			exit(-1);
		}
	}else{
		outfile= stdout;
	}
	
	while ((numseq = fp(ri, param,file)) != 0){
		//fprintf(stderr,"rread: %d\n",numseq);
		total_read+= numseq;
		for(i = 0; i < numseq;i++){
			
			c  = byg_end_barcode(pattern, ri[i]->seq ,plen, ri[i]->len);
			
			if((ri[i]->len - c ) < param->minlen){
				c =  -1;
			}
			//print_seq_from_position_x(ri[i],outfile,0);
			
			//fprintf(outfile,"position:%d\n",c);
			switch (c) {
				case -1:
					failure++;
					break;
				
				default:
					print_seq_from_position_x(ri[i],outfile,c);
					success++;
					break;
			}
			
			
		}
	}
	fprintf(stderr,"%d\n", total_read);
	fprintf(stderr,"%d	successfully extracted\n" ,success);
	fprintf(stderr,"%d	low probability\n" , failure);
		
	fprintf(stderr,"%0.1f%% extracted\n",  (float) success / (float) total_read  *100.0f);
	
	for(i = 0; i < param->num_query;i++){
		free(ri[i]->strand);
		free(ri[i]->hits);
		
		if(ri[i]->name){
			free(ri[i]->name);
		}
		if(ri[i]->seq){
			free(ri[i]->seq);
		}
		if(ri[i]->qual){
			free(ri[i]->qual );
		}
		if(ri[i]->labels){
			free(ri[i]->labels);
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


/** \fn void print_seq_from_position_x(struct read_info* ri,FILE* out,int x)
 \brief Prots sequences excluding the Barcode. 
 
 \param x Position of printed sequence in read. 
 \param ri @ref read_info - hols sequence. .
 \param out Output file pointer. 
 
 */
void print_seq_from_position_x(struct read_info* ri,FILE* out,int x)
{
	char alphabet[] = "ACGTN";
	int i;
	fprintf(out,"@%s\n",ri->name);
	for(i =x;i < ri->len;i++){
		fprintf(out,"%c",alphabet[(int) ri->seq[i]]);
	}
	fprintf(out,"\n+\n");
	if(ri->qual){
		for(i =x;i < ri->len;i++){
			fprintf(out,"%c",ri->qual[i]);
		}

		//fprintf(out,"%s",ri->qual);
	}else{
		for(i =x;i < ri->len;i++){
			fprintf(out,".");
		}
	}
	fprintf(out,"\n");
}

/** \fn int byg_end_barcode(const char* pattern,const char*text, int m, int n)
 \brief Finds end of barcode in reads. 
 
 \param pattern Barcode sequence. 
 \param text Read sequence. 
 \param m Barcode Length. 
 \param n Read length
 
 */


int byg_end_barcode(const char* pattern,const char*text, int m, int n)// m is pattern length , n is text length....
{
	//const char* tmp = 0;
	int Tc;
	int i  = 0;
	int j = 0;
	int s = 0;
	int T[4];
	for (i = 0;i < 4;i++){
		T[i] = 0;
	}
	
	if(m > 31){
		m = 31;
	}
	
	//int m = (int)strlen(pattern);
	//int n = (int)strlen(text);
	if (m > n){
		return  -1;
	}
	
	int mb = (1 << (m-1));
	
	for (i= 0;i < m;i++){
		if((int) pattern[i] == 4 ){
			for(j = 0; j < 4;j++){
				T[j] |= (1 << i);
			}
		}else{
			T[(int)pattern[i]] |= (1 << i);
		}
	//	fprintf(stderr,"%d\n", pattern[i]);
	}
	/*fprintf(stderr,"PAttern length = %d\n",m );//
	for(j = 0; j < 4;j++){
		fprintf(stderr,"%x\n",T[j] );
	}
	*/
	
	for (i = 0;i < n;i++){
		s <<= 1;
		s |= 1;
		//if(!text[i]){
		//	return -1;
		//}
		Tc = T[(int)text[i]];
		s &= Tc;
		if(s & mb){
			return i+1;
		}
	}
	//exit(0);
	return -1;
}

