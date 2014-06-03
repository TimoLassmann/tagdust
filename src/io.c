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

/*! \file io.c
 \brief functions for reading sequences.
 
 Initializes nucleotide alphabet needed to parse input. Calls parameter parser. Calls functions to process the data. \author Timo Lassmann \bug No known bugs.
 */

#include <ctype.h>
#include "interface.h"
#include "nuc_code.h"
#include "misc.h"

#include "io.h"

#include "tagdust2.h"

#ifndef MMALLOC
#include "malloc_macro.h"
#endif


struct sequence_stats_info* get_sequence_stats(struct parameters* param, struct read_info** ri, int file_num )
{
	struct sequence_stats_info* ssi = 0;
	
	MMALLOC(ssi, sizeof(struct sequence_stats_info));
	FILE* file = 0;
	
	int (*fp)(struct read_info** ,struct parameters*,FILE* ) = 0;
	
	int i,j,c,numseq,total_read;
	
	int five_len = 0;
	int three_len = 0;
	
	double sum = 0.0;
	
	double five_s0 = 0.0;
	double five_s1 = 0.0;
	double five_s2 = 0.0;
	double three_s0 = 0.0;
	double three_s1 = 0.0;
	double three_s2 = 0.0;
	
	char* five_test_sequence = 0;
	char* three_test_sequence = 0;
	
	ssi->average_length = 0;
	for(i = 0; i < 5;i++){
		ssi->background[i] = 1.0;
	}
	
	ssi->expected_5_len = 0;
	ssi->expected_3_len = 0;
	ssi->mean_5_len = 0.0f;
	ssi->stdev_5_len = 0.0f;
	ssi->mean_3_len = 0.0f;
	ssi->stdev_3_len = 0.0f;
	ssi->average_length = 0.0f;
	ssi->max_seq_len = 0;
	
	
	file =  io_handler(file, file_num,param);
	
	if(param->sam == 0){
		fp = &read_fasta_fastq;
	}else {
		fp = &read_sam_chunk;
	}
	numseq = 0;
	total_read= 0;
	
	//Do I need to test for a 5' partial sequence?
	if(param->read_structure->type[0] == 'P'){
		five_len = (int) strlen(param->read_structure->sequence_matrix[0][0]);
		ssi->expected_5_len = five_len;
		MMALLOC(five_test_sequence, sizeof(char) * (five_len+1));
		
		for(i = 0; i < five_len;i++){
			five_test_sequence[i] = nuc_code[(int) param->read_structure->sequence_matrix[0][0][i]];
		}
		five_test_sequence[five_len] = 0;
	}
	//Do I need to test for a 3' partial sequence?
	if(param->read_structure->type[param->read_structure->num_segments-1] == 'P'){
		three_len = (int) strlen(param->read_structure->sequence_matrix[ param->read_structure->num_segments-1][0]);
		ssi->expected_3_len = three_len;
		MMALLOC(three_test_sequence ,sizeof(char) * (three_len+1));
		for(i = 0; i < three_len;i++){
			three_test_sequence[i] = nuc_code[(int) param->read_structure->sequence_matrix[ param->read_structure->num_segments-1][0][i]];
		}
		three_test_sequence[three_len] = 0;
		
	}
	
	while ((numseq = fp(ri, param,file)) != 0){
		for(i = 0; i < numseq;i++){
			if(ri[i]->len > ssi->max_seq_len){
				ssi->max_seq_len = ri[i]->len;
			}
			ssi->average_length += ri[i]->len;
			for(j = 0;j < ri[i]->len;j++){
				//fprintf(stderr,"%d ",(int)ri[i]->seq[j] );
				ssi->background[(int)ri[i]->seq[j]] += 1.0f;
			}
			//check length of exact 5' matching sequence... 
			if(five_len){
				for(j = 0;j <= five_len ;j++){
					for(c = 0;c < five_len-j;c++){
						if(ri[i]->seq[c] != five_test_sequence[j +c]){
							break;
						}
					}
					if(c == five_len-j && c > 3 ){
						five_s0++;
						five_s1 += five_len -j;
						five_s2 += (five_len-j) * (five_len-j);
						break;
					}
					
				}
			}
			//check length of exact 3' matching sequence...
			if(three_len){
				for(j = 0;j <= three_len ;j++){
					
					for(c = 0;c < three_len-j;c++){
						if(ri[i]->seq[ri[i]->len - (three_len-j -c)] != three_test_sequence[c]){
							break;
						}
					}
					if(c == three_len-j  && c > 3){
						three_s0++;
						three_s1 += three_len -j;
						three_s2 += (three_len-j) * (three_len-j);
						break;
					}
				}
			}
			
			
		}
		
		total_read += numseq;
#if DEBUG
		if(total_read > 10001){
			break;
		}
#else
		if(total_read > 1000000){
			break;
		}
#endif
	}
	
	
	if(five_len){
		if(five_s0 <= 1){
			sprintf(param->buffer,"WARNING: there seems to e not a single read containing the 5' partial sequence.\n");
			param->messages = append_message(param->messages, param->buffer);
			ssi->mean_5_len  = ssi->expected_5_len;
			ssi->stdev_5_len  = 1.0;
			
			
		}else{
		
		ssi->mean_5_len = five_s1 / five_s0;
		ssi->stdev_5_len = sqrt(  (five_s0 * five_s2 - pow(five_s1,2.0))   /  (  five_s0 *(five_s0-1.0) )) ;
		if(!ssi->stdev_5_len){
			ssi->stdev_5_len = 10000.0;
		}
		//fprintf(stderr,"5: %f %f	%f\n", ssi->mean_5_len,  ssi->stdev_5_len,five_s0);
		//if(ssi->stdev_5_len < 1){
		//	ssi->stdev_5_len = 1;
		//}
		
		//fprintf(stderr,"5: %f %f	%f\n", ssi->mean_5_len,  ssi->stdev_5_len,five_s0);
			if(ssi->mean_5_len <= 1){
				sprintf(param->buffer,"WARNING: 5' partial segment seems not to be present in the data (length < 1).\n");
				param->messages = append_message(param->messages, param->buffer);
				//free_param(param);
				//exit(EXIT_FAILURE);
			}
		}
	}else{
		ssi->mean_5_len =  -1.0;
		ssi->stdev_5_len = -1.0;
	}
	
	
	
	if(three_len){
		if(three_s0 <= 1){
			sprintf(param->buffer,"WARNING: 3' partial segment seems not to be present in the data.\n");
			param->messages = append_message(param->messages, param->buffer);
			ssi->mean_3_len  = ssi->expected_3_len;
			ssi->stdev_3_len  = 1.0;

		}else{
		
		ssi->mean_3_len = three_s1 / three_s0;
		ssi->stdev_3_len = sqrt(  (three_s0 * three_s2 - pow(three_s1,2.0))   /  (  three_s0 *(three_s0-1.0) )) ;
		if(!ssi->stdev_3_len){
			ssi->stdev_3_len = 10000.0;
		}
		//fprintf(stderr,"3: %f %f	%f\n", ssi->mean_3_len,  ssi->stdev_3_len,three_s0);
		//if(ssi->stdev_3_len < 1){
		//	ssi->stdev_3_len = 1;
		//}
		//fprintf(stderr,"3: %f %f	%f\n", ssi->mean_3_len,  ssi->stdev_3_len,three_s0);
			if(ssi->mean_3_len <= 1){
				sprintf(param->buffer,"WARNING: 3' partial segment seems not to be present in the data (length < 1).\n");
			//	fprintf(stderr,"%s",param->buffer);
				param->messages = append_message(param->messages, param->buffer);
				//free_param(param);
				//exit(EXIT_FAILURE);
			}
		}
	}else{
		ssi->mean_3_len =  -1.0;
		ssi->stdev_3_len = -1.0;
	}
	
	if(param->matchstart!= -1 || param->matchend !=-1){
		ssi->average_length = (param->matchend - param->matchstart )* total_read;
	}
	ssi->average_length =  (int) floor((double)  ssi->average_length / (double) total_read   + 0.5);
	
	sum = 0.0;
	for(i = 0; i < 5;i++){
		sum += ssi->background[i];
	}
	
	for(i = 0; i < 5;i++){
		ssi->background[i] = prob2scaledprob(ssi->background[i]  / sum);
	}
	
	if(five_test_sequence){
		MFREE(five_test_sequence);
	}
	if(three_test_sequence){
		MFREE(three_test_sequence);
	}
#ifdef DEBUG
	fprintf(stderr,"Backgound:\n");
	fprintf(stderr,"A:%f\n",scaledprob2prob( ssi->background[0]));
	fprintf(stderr,"C:%f\n",scaledprob2prob( ssi->background[1]));
	fprintf(stderr,"G:%f\n",scaledprob2prob( ssi->background[2]));
	fprintf(stderr,"T:%f:\n",scaledprob2prob( ssi->background[3]));
	fprintf(stderr,"N:%f\n",scaledprob2prob( ssi->background[4]));
	
	fprintf(stderr,"\nExpected Length:5':%f	3':%f\n",ssi->expected_5_len,ssi->expected_3_len );
	fprintf(stderr,"Observed Length:5':%f	3':%f\n",ssi->mean_5_len,ssi->mean_3_len );
	fprintf(stderr,"STDEV:5':%f	3':%f\n",ssi->stdev_5_len,ssi->stdev_3_len );
	    
	        
#endif
	
	
	
	
	pclose(file);
	return ssi;
}






/** \fn int qsort_ri_prob_compare(const void *a, const void *b)
 \brief Compares reads based their probability.
 Used to sort arrays of string using qsort.
 \param a void pointer to first @ref read_info.
 \param b void pointer to second @ref read_info.
 */
int qsort_ri_barcode_compare(const void *a, const void *b)
{
	
	//struct mys **a = (struct mys **)i1;
	//struct mys **b = (struct mys **)i2;
	//return (*b)->id - (*a)->id;
	
	const struct read_info **elem1 = (const struct read_info**) a;
	
	const struct read_info **elem2 = (const struct read_info**) b;
	
	if ( (*elem1)->read_type <  (*elem2)->read_type){
		return -1;
	}else if ((*elem1)->read_type > (*elem2)->read_type){
		return 1;
	}else{
		if ( (*elem1)->barcode <  (*elem2)->barcode){
			return -1;
		}else if ((*elem1)->barcode > (*elem2)->barcode){
			return 1;
		}else{
			return 0;
		}
	}
	
	
}

/** \fn int qsort_ri_prob_compare(const void *a, const void *b)
 \brief Compares reads based their probability.
 Used to sort arrays of string using qsort.
 \param a void pointer to first @ref read_info.
 \param b void pointer to second @ref read_info.
 */
int qsort_ri_mapq_compare(const void *a, const void *b)
{
	
	//struct mys **a = (struct mys **)i1;
	//struct mys **b = (struct mys **)i2;
	//return (*b)->id - (*a)->id;
	
	const struct read_info **elem1 = (const struct read_info**) a;
	
	const struct read_info **elem2 = (const struct read_info**) b;
	
	if ( (*elem1)->mapq > (*elem2)->mapq)
		return -1;
	
	else if ((*elem1)->mapq < (*elem2)->mapq)
		return 1;
	
	else
		return 0;
}



/** \fn FILE* io_handler(FILE* file, int file_num,struct parameters* param)
 \brief Opens stream to files. 
 
 Recognizes file types by prefix. 
 
 Used to sort arrays of string using qsort.
 \param file empty file pointer.
\param file_num index of input file.
 \param param @ref parameters.

 */

FILE* io_handler(FILE* file, int file_num,struct parameters* param)
{
	char command[1000];
	char  tmp[1000];
	int i = 0; 
	int gzcat = -1;
	if(access("/usr/bin/gzcat", X_OK) == 0){
		gzcat = 1;
	}else if(access("/bin/gzcat", X_OK) == 0){
		gzcat = 1;
	}else if(access("/usr/bin/zcat", X_OK) == 0){
		gzcat = 0;
	}else if(access("/bin/zcat", X_OK) == 0){
		gzcat = 0;
	}

	param->gzipped = 0;
	param->bzipped = 0;
	param->sam = 0;
	param->fasta = 0;
	
	if(!file_exists(param->infile[file_num])){
		sprintf(param->buffer,"Error: Cannot find input file: %s\n",param->infile[file_num] );
		param->messages = append_message(param->messages, param->buffer);
		free_param(param);
		exit(EXIT_FAILURE);
	}
	
	if(!strcmp(".sam", param->infile[file_num] + (strlen(param->infile[file_num] ) - 4))){
		param->sam = 1;
	}else if (!strcmp(".bam", param->infile[file_num] + (strlen(param->infile[file_num] ) - 4))){
		param->sam = 2;
	}else if (!strcmp(".fa", param->infile[file_num] + (strlen(param->infile[file_num] ) - 3))){
		param->sam = 0;
		param->fasta = 1;
	}else if (!strcmp(".fq", param->infile[file_num] + (strlen(param->infile[file_num] ) - 3))){
		param->sam = 0;
	}else if (!strcmp(".fastq", param->infile[file_num] + (strlen(param->infile[file_num] ) - 6))){
		param->sam = 0;
	}else if (!strcmp(".fastaq", param->infile[file_num] + (strlen(param->infile[file_num] ) - 7))){
		param->sam = 0;
	}else if (!strcmp(".fasta", param->infile[file_num] + (strlen(param->infile[file_num] ) - 6))){
		param->sam = 0;
		param->fasta = 1;
	}else if(!strcmp(".sam.gz", param->infile[file_num] + (strlen(param->infile[file_num] ) - 7))){
		param->sam = 1;
		param->gzipped  = 1;
	}else if (!strcmp(".bam.gz", param->infile[file_num] + (strlen(param->infile[file_num] ) - 7))){
		param->sam = 2;
		param->gzipped  = 1;
	}else if (!strcmp(".fa.gz", param->infile[file_num] + (strlen(param->infile[file_num] ) - 6))){
		param->sam = 0;
		param->fasta = 1;
		param->gzipped  = 1;
	}else if (!strcmp(".fq.gz", param->infile[file_num] + (strlen(param->infile[file_num] ) - 6))){
		param->sam = 0;
		param->gzipped  = 1;
	}else if (!strcmp(".fastq.gz", param->infile[file_num] + (strlen(param->infile[file_num] ) - 9))){
		param->sam = 0;
		param->gzipped  = 1;
	}else if (!strcmp(".fastaq.gz", param->infile[file_num] + (strlen(param->infile[file_num] ) - 10))){
		param->sam = 0;
		param->gzipped  = 1;
	}else if (!strcmp(".fasta.gz", param->infile[file_num] + (strlen(param->infile[file_num] ) - 9))){
		param->sam = 0;
		param->gzipped  = 1;
	}else if (!strcmp(".fastq.bz2", param->infile[file_num] + (strlen(param->infile[file_num] ) - 10))){
		param->sam = 0;
		param->bzipped  = 1;
	}else if (!strcmp(".fq.bz2", param->infile[file_num] + (strlen(param->infile[file_num] ) - 7))){
		param->sam = 0;
		param->bzipped  = 1;
	}else{
		param->sam = -1;
	}
	
	
	if(param->gzipped && gzcat == -1){
		sprintf(param->buffer,"Cannot find gzcat / zcat on your system. Try gzcat <infile> | samstat -f sam/bam/fa/fq\n");
		param->messages = append_message(param->messages, param->buffer);
		free_param(param);
		exit(EXIT_FAILURE);
	}
	
	if(file_num == -1){
		if(param->sam == 2){
			command[0] = 0;
			if(!param->filter){
				strcat ( command, "samtools view -F 768 "); 
			}else{
				strcat ( command, "samtools view -F "); 
				i = sprintf (tmp, "%s ",param->filter);
				strcat ( command, tmp);
			}
			i = sprintf (tmp, "%s ","-");
			strcat ( command, tmp);
			if (!(file = popen(command, "r"))) {
				fprintf(stderr,"Cannot open bam file '%s' with command:%s\n",param->infile[file_num],command);
				exit(-1);
			}
		}else if(param->sam == 1){
			command[0] = 0;
			if(!param->filter){
				strcat ( command, "samtools view -SF 768 "); 
			}else{
				strcat ( command, "samtools view -SF "); 
				i = sprintf (tmp, "%s ",param->filter);
				strcat ( command, tmp);
			}
			i = sprintf (tmp, "%s ", "-");
			strcat ( command, tmp);
			if (!(file = popen(command, "r"))) {
				fprintf(stderr,"Cannot open bam file '%s' with command:%s\n",param->infile[file_num],command);
				exit(-1);
			}
		}else{
			file = stdin;
		}
	}else{
		if(param->sam == 2){
			command[0] = 0;
			
			if(param->bzipped){
				strcat ( command, "bzcat ");
				if(!param->filter){
					i = sprintf (tmp, "%s | samtools view -F 768 - ", param->infile[file_num]);
					strcat ( command, tmp);
				}else{
					i = sprintf (tmp, "%s | samtools view -F  ", param->infile[file_num]);
					strcat ( command, tmp);
					i = sprintf (tmp, "%s - ",param->filter);
					strcat ( command, tmp);
				}
				
			}else if(param->gzipped){
				if(gzcat == 1){
					strcat ( command, "gzcat "); 
				}else{
					strcat ( command, "zcat "); 
				}
				if(!param->filter){
					i = sprintf (tmp, "%s | samtools view -F 768 - ", param->infile[file_num]);
					strcat ( command, tmp);
				}else{
					i = sprintf (tmp, "%s | samtools view -F  ", param->infile[file_num]);
					strcat ( command, tmp);
					i = sprintf (tmp, "%s - ",param->filter);
					strcat ( command, tmp);
				}
			}else{
				if(!param->filter){
					strcat ( command, "samtools view -F 768 "); 
				}else{
					strcat ( command, "samtools view -F "); 
					i = sprintf (tmp, "%s ",param->filter);
					strcat ( command, tmp);
				}
				i = sprintf (tmp, "%s ", param->infile[file_num]);
				strcat ( command, tmp);
			}
			if (!(file = popen(command, "r"))) {
				sprintf(param->buffer,"Cannot open bam file '%s' with command:%s\n",param->infile[file_num],command);
				param->messages = append_message(param->messages, param->buffer);
				free_param(param);
				exit(EXIT_FAILURE);
			}
		}else if(param->sam == 1){
			command[0] = 0;
			if(param->gzipped){
				if(gzcat == 1){
					strcat ( command, "gzcat "); 
				}else{
					strcat ( command, "zcat "); 
				}
				if(!param->filter){
					i = sprintf (tmp, "%s | samtools view -SF 768 - ", param->infile[file_num]);
					strcat ( command, tmp);
				}else{
					i = sprintf (tmp, "%s | samtools view -SF  ", param->infile[file_num]);
					strcat ( command, tmp);
					i = sprintf (tmp, "%s - ",param->filter);
					strcat ( command, tmp);
				}
			}else{
				if(!param->filter){
					strcat ( command, "samtools view -SF 768 "); 
				}else{
					strcat ( command, "samtools view -SF "); 
					i = sprintf (tmp, "%s ",param->filter);
					strcat ( command, tmp);
				}
				i = sprintf (tmp, "%s ", param->infile[file_num]);
				strcat ( command, tmp);
			}
			if (!(file = popen(command, "r"))) {
				sprintf(param->buffer,"Cannot open bam file '%s' with command:%s\n",param->infile[file_num],command);
				param->messages = append_message(param->messages, param->buffer);
				free_param(param);
				exit(EXIT_FAILURE);
			}
		}else{
			command[0] = 0;
			if(param->bzipped){
				strcat ( command, "bzcat ");
				
			}else if(param->gzipped){
				if(gzcat == 1){
					strcat ( command, "gzcat ");
				}else{
					strcat ( command, "zcat ");
				}
			}else{
				strcat ( command, "cat ");
			}
			i = sprintf (tmp, "%s ", param->infile[file_num]);
			strcat ( command, tmp);
			//fprintf(stderr,"%s\n",command);
			if (!(file = popen(command, "r"))) {
				sprintf(param->buffer,"Cannot open bam file '%s' with command:%s\n",param->infile[file_num],command);
				param->messages = append_message(param->messages, param->buffer);
				free_param(param);
				exit(EXIT_FAILURE);
			}
		}
	}
	return file;
}



/** \fn void print_sequence(struct read_info* ri,FILE* out)
 \brief Prints sequence and quality string to file.  

 \param ri @ref read_info containing the sequence.
 \param out Pointer to output file. 
 \deprecated use @ref print_seq instead. 

 */

void print_sequence(struct read_info* ri,FILE* out)
{
	int i;
	char alpha[6] = "ACGTNN";
	fprintf(out,"@%s\n",ri->name);
	for(i = 0; i < ri->len;i++){
		fprintf(out,"%c", alpha[(int) ri->seq[i]]);
	}
	fprintf(out,"\n+\n%s\n" ,ri->qual);
}





void print_split_files(struct parameters* param, struct read_info** ri, int numseq)
{
	struct stat buf;
	FILE* out_read1 = NULL;
	FILE* out_read2 = NULL;
	char alphabet[] = "ACGTNN";
	int i,j;
	int old_type = -100;
	int old_bar = -100;
	int read_file = 0;
	void* tmp = 0;
	int start = 0;
	int stop = -1;
	int last = 0;
	int segment = 1;
	int print_segment_name = 0;
	
	static int check_for_files = 1;
	
	char* buffer =  0;
	MMALLOC(buffer,sizeof(char)* 1000 );
	qsort(ri,numseq, sizeof(struct read_info*), qsort_ri_barcode_compare);
	
	//1 ) has barcode or nor
	//	has: barcode >= 0 && extract_successs
	//	has_nor barcore == -1 && extract_success
	// 2) fail : ! extract_success..
	// ri[i]->read_type,ri[i]->barcode
	int h = 0;
	for(i = 0; i < numseq;i++){
		read_file = 0;
		
		if(ri[i]->read_type != old_type || ri[i]->barcode != old_bar  ){
			if(ri[i]->read_type  == EXTRACT_SUCCESS){
				buffer[0] = 0;
				if(ri[i]->barcode != -1){
					if(param->multiread == 2){
						sprintf (buffer, "%s_BC_%s_READ1.fq",param->outfile,param->read_structure->sequence_matrix[(ri[i]->barcode >> 16) &0XFF][ri[i]->barcode &0XFF]);
					}else{
						sprintf (buffer, "%s_BC_%s.fq",param->outfile,param->read_structure->sequence_matrix[(ri[i]->barcode >> 16) &0XFF][ri[i]->barcode &0XFF]);
					}
				}else{
					if(param->multiread == 2){
						sprintf (buffer, "%s_READ1.fq",param->outfile);
					}else{
						sprintf (buffer, "%s.fq",param->outfile);
					}
				}
			}else{
				buffer[0] = 0;
				if(param->multiread == 2){
					sprintf (buffer, "%s_un_READ1.fq",param->outfile);
				}else{
					sprintf (buffer, "%s_un.fq",param->outfile);
				}
				
			}
				
			if(check_for_files){
				check_for_files = 0;
				if(!stat ( buffer, &buf )){
					//file found.
					sprintf(param->buffer,"ERROR: output file: %s already exists.\n", buffer);
					param->messages = append_message(param->messages, param->buffer);
					free_param(param);
					exit(EXIT_FAILURE);
				}
				if ((out_read1 = fopen(buffer, "w")) == NULL){
					fprintf(stderr,"can't open output\n");
					exit(-1);
				}
				buffer[0] = 0;
				if(param->multiread ==2){
					if(ri[i]->read_type  == EXTRACT_SUCCESS){
						buffer[0] = 0;
						if(ri[i]->barcode != -1){
							sprintf (buffer, "%s_BC_%s_READ2.fq",param->outfile,param->read_structure->sequence_matrix[(ri[i]->barcode >> 16) &0XFF][ri[i]->barcode &0XFF]);
						}else{
							sprintf (buffer, "%s_READ2.fq",param->outfile);
							
						}
					}else{
						buffer[0] = 0;
						sprintf (buffer, "%s_un_READ2.fq",param->outfile);
						
					}
					if ((out_read2 = fopen(buffer, "w")) == NULL){
						fprintf(stderr,"can't open output\n");
						exit(-1);
					}
				}
			}else{
				if(i){
					if(param->multiread ==2){
						fclose(out_read2);
					}
					fclose(out_read1);
				}
				
				if ((out_read1 = fopen(buffer, "a")) == NULL){
					fprintf(stderr,"can't open output\n");
					exit(-1);
				}
				if(param->multiread ==2){
					if(ri[i]->read_type  == EXTRACT_SUCCESS){
						buffer[0] = 0;
						if(ri[i]->barcode != -1){
							sprintf (buffer, "%s_BC_%s_READ2.fq",param->outfile,param->read_structure->sequence_matrix[(ri[i]->barcode >> 16) &0XFF][ri[i]->barcode &0XFF]);
						}else{
							sprintf (buffer, "%s_READ2.fq",param->outfile);
							
						}
					}else{
						buffer[0] = 0;
						sprintf (buffer, "%s_un_READ2.fq",param->outfile);
						
					}
					if ((out_read2 = fopen(buffer, "a")) == NULL){
						fprintf(stderr,"can't open output\n");
						exit(-1);
					}
				}
			}
		//	file open / close
			old_type = ri[i]->read_type;
			old_bar = ri[i]->barcode;
		}
		
		if(ri[i]->read_type == EXTRACT_SUCCESS){
			buffer[0] = 0;
			j = 0;
			if(ri[i]->fingerprint != -1){
				j |= 2;
			}
			if(ri[i]->barcode != -1){
				j |= 1;
			}
			switch (j) {
				case 0:
					sprintf (buffer, "%s",ri[i]->name);
					break;
				case 1:
					sprintf (buffer, "%s;BC:%s",ri[i]->name, param->read_structure->sequence_matrix[(ri[i]->barcode >> 16) &0XFF][ri[i]->barcode &0XFF]);
					break;
				case 2:
					sprintf (buffer, "%s;FP:%d",ri[i]->name,ri[i]->fingerprint);
					break;
				case 3:
					sprintf (buffer, "%s;FP:%d;BC:%s",ri[i]->name,ri[i]->fingerprint ,param->read_structure->sequence_matrix[(ri[i]->barcode >> 16) &0XFF][ri[i]->barcode &0XFF]);
					break;
				default:
					break;
			}
			
			//strcat (buffer, tmp);
			//ri[i]->name = realloc(ri[i]->name, sizeof(char) * (strlen(buffer) + 1) );
			MREALLOC(ri[i]->name,tmp,sizeof(char) * (strlen(buffer) + 1)) ;
			assert(ri[i]->name  != NULL);
			strcpy(ri[i]->name, buffer);
			
			
						
			start = 0;
			stop = -1;
			last = 0;
			segment = 1;
			print_segment_name = 0;
			for(j =0;j < ri[i]->len;j++){
				if( ri[i]->seq[j] == 65){
					print_segment_name = 1;
					break;
				}
			}
			if(print_segment_name){
				//if we have a multiread scroll to the first read if starting with 65....
				while(ri[i]->seq[start] == 65 && start < ri[i]->len){
					start++;
					
				}
			}

			while(stop != ri[i]->len){
				last = 0;
				for(j =start;j < ri[i]->len;j++){
					if( ri[i]->seq[j] == 65){
						stop = j;
						last = 1;
						break;
					}
				}
				if(!last){
					stop = ri[i]->len;
				}
				//print to file segment 'X'
				
				//if(print_segment_name){
				//	fprintf(out,"@%sRS:%d\n",ri->name,segment);
				//}else{
				
				if(segment ==1){
					h++;
					fprintf(out_read1,"@%s;RQ:%0.2f\n",ri[i]->name,ri[i]->mapq);
					
					for(j = start; j < stop;j++){
						fprintf(out_read1,"%c",alphabet[(int) ri[i]->seq[j]]);
					}
					fprintf(out_read1,"\n+\n");
					if(ri[i]->qual){
						for(j =start;j < stop;j++){
							fprintf(out_read1,"%c",ri[i]->qual[j]);
						}
					}else{
						for(j =start;j < stop;j++){
							fprintf(out_read1,".");
						}
					}
					fprintf(out_read1,"\n");

				}else if ( segment == 2){
					fprintf(out_read2,"@%s;RQ:%0.2f\n",ri[i]->name,ri[i]->mapq);
					
					
					for(j = start; j < stop;j++){
						fprintf(out_read2,"%c",alphabet[(int) ri[i]->seq[j]]);
					}
					fprintf(out_read2,"\n+\n");
					if(ri[i]->qual){
						for(j =start;j < stop;j++){
							fprintf(out_read2,"%c",ri[i]->qual[j]);
						}
					}else{
						for(j =start;j < stop;j++){
							fprintf(out_read2,".");
						}
					}
					fprintf(out_read2,"\n");
				}
				segment++;
				
				start = stop;
				while(ri[i]->seq[start] == 65 && start < ri[i]->len){
					start++;
				}
				
				if(segment > 100){
					exit(EXIT_FAILURE);
				}
			}
			
			
		}else{
			start = 0;
			stop = -1;
			last = 0;
			segment = 1;
			print_segment_name = 0;
			for(j =0;j < ri[i]->len;j++){
				if( ri[i]->seq[j] == 65){
					print_segment_name = 1;
					
					break;
				}
			}
			
#ifdef DEBUG
			fprintf(stderr,"READ unextracted!:%d\n" , i);
			for(j = 0; j < ri[i]->len;j++){
				fprintf(stderr,"%d,",ri[i]->seq[j]);
			}
			fprintf(stderr,"\n");
#endif
			
			
			if(print_segment_name){
				//if we have a multiread scroll to the first read if starting with 65....
				while(ri[i]->seq[start] == 65 && start < ri[i]->len){
					start++;
					
				}
			}
			while(stop != ri[i]->len){
				last = 0;
				for(j =start;j < ri[i]->len;j++){
					if( ri[i]->seq[j] == 65){
						stop = j;
						last = 1;
						//	segment++;
						break;
					}
				}
				
				if(!last){
					stop = ri[i]->len;
				}
				
				if(segment ==1){
					h++;
					fprintf(out_read1,"@%s;RQ:%0.2f\n",ri[i]->name,ri[i]->mapq);
					
					for(j = start; j < stop;j++){
						fprintf(out_read1,"%c",alphabet[(int) ri[i]->seq[j]]);
					}
					fprintf(out_read1,"\n+\n");
					if(ri[i]->qual){
						for(j =start;j < stop;j++){
							fprintf(out_read1,"%c",ri[i]->qual[j]);
						}
					}else{
						for(j =start;j < stop;j++){
							fprintf(out_read1,".");
						}
					}
					fprintf(out_read1,"\n");
				}else if ( segment == 2){
					fprintf(out_read2,"@%s;RQ:%0.2f\n",ri[i]->name,ri[i]->mapq);
					
					for(j = start; j < stop;j++){
						fprintf(out_read2,"%c",alphabet[(int) ri[i]->seq[j]]);
					}
					fprintf(out_read2,"\n+\n");
					if(ri[i]->qual){
						for(j =start;j < stop;j++){
							fprintf(out_read2,"%c",ri[i]->qual[j]);
						}
					}else{
						for(j =start;j < stop;j++){
							fprintf(out_read2,".");
						}
					}
					fprintf(out_read2,"\n");
				}
				segment++;
				
				start = stop;
				while(ri[i]->seq[start] == 65 && start < ri[i]->len){
					start++;
				}
				
				if(segment > 100){
					exit(EXIT_FAILURE);
				}
			}
		}
	}
	MFREE(buffer);
	
	if(param->multiread ==2){
		fclose(out_read2);
	}
	fclose(out_read1);
}



/** \fn void print_seq(struct read_info* ri,FILE* out)
 \brief Prints sequence and quality string to file.
 If the quality string is empty it is replaced by '.' characters. 
 \param ri @ref read_info containing the sequence.
 \param out Pointer to output file.
 
 */
void print_seq(struct read_info* ri,FILE* out)
{
	char alphabet[] = "ACGTNN";
	int i;
	int start = 0;
	int stop = -1;
	int last = 0;
	int segment = 1;
	int print_segment_name = 0;
	
	for(i =0;i < ri->len;i++){
		if( ri->seq[i] == 65){
			print_segment_name = 1;
			
			break;
		}
	}
	if(print_segment_name){
		//if we have a multiread scroll to the first read if starting with 65....
		while(ri->seq[start] == 65 && start < ri->len){
			start++;
		
		}
	}
	
	while(stop != ri->len){
		last = 0;
		for(i =start;i < ri->len;i++){
			if( ri->seq[i] == 65){
				stop = i;
				last = 1;
			//	segment++;
				break;
			}
		}
		if(!last){
			stop = ri->len;
		}
		if(print_segment_name){
			fprintf(out,"@%sRS:%d\n",ri->name,segment);
		}else{
		
			fprintf(out,"@%s\n",ri->name);
		}
		for(i = start; i < stop;i++){
			fprintf(out,"%c",alphabet[(int) ri->seq[i]]);
		}
		fprintf(out,"\n+\n");
		if(ri->qual){
			for(i =start;i < stop;i++){
				fprintf(out,"%c",ri->qual[i]);
			}
		}else{
			for(i =start;i < stop;i++){
				fprintf(out,".");
			}
		}
		fprintf(out,"\n");
		segment++;
		
		start = stop;
		while(ri->seq[start] == 65 && start < ri->len){
			start++;
		}
		
		if(segment > 100){
			exit(EXIT_FAILURE);
		}
	}
}


/** \fn int read_sam_chunk(struct read_info** ri,struct parameters* param,FILE* file)
 \brief Reads sequences from SAM / BAM file. 
Sequences are converted to 0-4 strings. 
 
 \param ri Array of @ref read_info to hold the sequences.
 \param param @ref parameters .  
 \param file Pointer to input file.
 
 */
int read_sam_chunk(struct read_info** ri,struct parameters* param,FILE* file)
{
	char line[MAX_LINE];
	int column = 0; 
	int i,j,g,tmp;
	unsigned int pos;
	int read = 0;
	int c = 0;
	int hit = 0;
	int strand = 0;
	
	ri = clear_read_info(ri, param->num_query);
	
	while(fgets(line, MAX_LINE, file)){
		if(line[0] != '@'){
			column = 1; //<QNAME> 
			tmp = 0;
			hit = 0;
			pos = 0xFFFFFFFFu;
			for(j = 0;j < MAX_LINE;j++){
				tmp++;
				if(isspace((int)line[j])){
					break;
				}
			}
			
			MMALLOC(ri[c]->name,sizeof(unsigned char)* tmp);
			for(j = 0;j < MAX_LINE;j++){
				
				if(isspace((int)line[j])){
					ri[c]->name[j] = 0;
					break;
				}
				ri[c]->name[j] = line[j];
			}
			
			for(i = 0; i < MAX_LINE;i++){
				if(line[i] == '\n'){
					break;
				}
				if(isspace((int)line[i])){
					column++;
					switch(column){
						case 2: // <FLAG>
							tmp = atoi(line+i+1);
							strand = (tmp & 0x10);

							//WARNING - read should be reverse complemented if mapped to negative strand before tagdusting...
							
							/*tmp = atoi(line+i+1);
							ri[c]->strand[hit] = (tmp & 0x10);
							if(tmp == 4){
								ri[c]->hits[hit] = 0;
							}else{
								ri[c]->hits[hit] = 1;
							}
							hit++;*/
							
							break;
						case 3: // <RNAME> 
							
							break;
						case 4: // <POS>
							
							break;
						case 5: //  <MAPQ>
							
							ri[c]->mapq =  atof(line +i +1); 
							
							break;
						case 6: //  <CIGAR>
							
							break;
						case 7: //  <MRNM>
							break;
						case 8: //  <MPOS>
							break;
						case 9: //  <ISIZE>
							break;
						case 10: // <SEQ>
							
							tmp = 0;
							for(j = i+1;j < MAX_LINE;j++){
								tmp++;
								if(isspace((int)line[j])){
									break;
								}
							}
							
							MMALLOC(ri[c]->seq,sizeof(unsigned char)* tmp);
							MMALLOC(ri[c]->labels,sizeof(unsigned char)* tmp);
							
							g = 0;
							for(j = i+1;j < MAX_LINE;j++){
								
								if(isspace((int)line[j])){
									ri[c]->seq[g] = 0;
									ri[c]->labels[g] = 0;
									break;
								}
								ri[c]->seq[g] = nuc_code[(int)line[j]];
								ri[c]->labels[g] = 0;

								g++;
							}
							
							ri[c]->len = g;
							break;
						case 11: // <QUAL>
							tmp = 0;
							for(j = i+1;j < MAX_LINE;j++){
								tmp++;
								if(isspace((int)line[j])){
									break;
								}
							}
							g= 0;
							MMALLOC(ri[c]->qual,sizeof(unsigned char)* tmp);
							for(j = i+1;j < MAX_LINE;j++){
								
								if(isspace((int)line[j])){
									ri[c]->qual[g] = 0;
									break;
								}
								ri[c]->qual[g] = line[j];
								g++;
							}
							break;
						default: 
							
									
							i = MAX_LINE;
							break;
					}				}

			}
			tmp = byg_end("NM:i:", line  );
			if(tmp){
				ri[c]->read_type = atoi(line+tmp);
			}else{
				ri[c]->read_type = -1;
			}	
			
			if(strand != 0){
			//	ri[c]->seq = reverse_complement(ri[c]->seq,ri[c]->len);
			}

			
			
			
			//ri[c]->hits[hit] = 0xFFFFFFFFu;
			
			c++;
			read++;
			if(c == param->num_query){
				return c;
			}
		}
	}
	return c;
}

//FASTQ files from CASAVA-1.8 Should have the following READ-ID format:
//@<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos> <read>:<is filtered>:<control number>:<index sequence>



/** \fn int read_fasta_fastq(struct read_info** ri,struct parameters* param,FILE *file)
 \brief Reads sequences from fasta / fastq file.
 Sequences are converted to 0-4 strings.
 
 \param ri Array of @ref read_info to hold the sequences.
 \param param @ref parameters .
 \param file Pointer to input file.
 
 */

int read_fasta_fastq(struct read_info** ri,struct parameters* param,FILE *file) 
{
	int park_pos = -1;
	char line[MAX_LINE];
	int i;//,j;
	int seq_p = 0;
	int set = 0;
	int len = 0;
	int size = 0;
	
	ri = clear_read_info(ri, param->num_query);
	while(fgets(line, MAX_LINE, file)){
		if((line[0] == '@' && !set)|| (line[0] == '>' && !set)){
			//set sequence length of previous read
			
			//check if there is still space....
			//if(param->num_query == size){
			//	fseek (file , -  strlen(line) , SEEK_CUR);
			//	return size;
			//}
			park_pos++;
			len = 0;
			seq_p = 1;
			for(i = 1;i < MAX_LINE;i++){
				len++;
				if(iscntrl((int)line[i])){
					break;
				}
				
			}
			
			//ri[park_pos]->hits[0] = 0;
			//ri[park_pos]->strand[0] = 0;
			MMALLOC(ri[park_pos]->name,sizeof(unsigned char)* (len+1));
			for(i = 1;i < MAX_LINE;i++){
				
				if(iscntrl((int)line[i])){
					ri[park_pos]->name[i-1] = 0;
					break;
				}
				if(isspace((int)line[i])){
					ri[park_pos]->name[i-1] = ';';
				}
				
				ri[park_pos]->name[i-1] = line[i];
			}
			//fprintf(stderr,"LEN:%d	%s\n",len,ri[park_pos]->name);
			
			set = 1;
			size++;
			//get ready to read quality if present  
		}else if(line[0] == '+' && !set){
			seq_p = 0;
			set = 1;
			//reading sequence or quality  
		}else{	
			if(set){
				if(seq_p){
					len = 0;
					for(i = 0;i < MAX_LINE;i++){
						len++;
						if(iscntrl((int)line[i])){
							break;
						}
					}
					//fprintf(stderr,"SEQ LEN:%d	%s\n",len,line);
					MMALLOC(ri[park_pos]->seq,sizeof(unsigned char)* (len+1));
					
					MMALLOC(ri[park_pos]->labels, sizeof(unsigned char)* (len+1));
					
					for(i = 0;i < MAX_LINE;i++){
						if(iscntrl((int)line[i])){
							ri[park_pos]->seq[i] = 0;
							ri[park_pos]->labels[i] = 0;
							break;
						}
						ri[park_pos]->seq[i] = nuc_code[(int)line[i]];
						ri[park_pos]->labels[i] = 0;
					}
					ri[park_pos]->len = len-1;
				}else{
					len = 0;
					for(i = 0;i < MAX_LINE;i++){
						len++;
						if(iscntrl((int)line[i])){
							break;
						}
						
					}
					
					if(len-1 != ri[park_pos]->len ){
						sprintf(param->buffer,"ERROR: Length of sequence and base qualities differ!.\n");
						param->messages = append_message(param->messages, param->buffer);
						free_param(param);
						exit(EXIT_FAILURE);
					}
					
					//fprintf(stderr,"QUAL LEN:%d\n",len);
					MMALLOC(ri[park_pos]->qual,sizeof(unsigned char)* (len+1));
					for(i = 0;i < MAX_LINE;i++){
						if(iscntrl((int)line[i])){
							ri[park_pos]->qual[i] = 0;
							break;
						}
						ri[park_pos]->qual[i] = line[i];
					}
				}
			}
			set = 0;
		}
		if(param->num_query == size ){//here I know I am in the last entry AND filled the quality...
			if(!param->fasta && ri[park_pos]->qual){
				return size;
			}
			if(param->fasta && ri[park_pos]->seq){
			   
				return size;
			}
		}
	}
	return size;
}

/** \fn struct fasta* get_fasta(struct fasta* p,char *infile)
 \brief Old Kalign function to read in fasta sequences.
 All sequences are read into one array indexed by pointers. 
 \param p @ref fasta struct to hold sequences.

 \param infile Pointer to input file.
 
 */

struct fasta* get_fasta(struct fasta* p,char *infile)
{
	MMALLOC(p,sizeof(struct fasta));
	p->string = 0;
	p->mer_hash = 0;
	p->s_index = 0;
	p->sn = 0;
	p->string = 0;
	p->boost = 0;
	p->max_len = 0;
	p->numseq = 0;
	p->string_len = 0;
	
	p->string =  get_input_into_string(p->string,infile);
	if(!p->string){
		fprintf(stderr,"Analysing input %s ... nothing\n",infile);
		exit(-1);
	}
	
	p = read_fasta(p);
	
	return p;
}



/** \fn unsigned char* get_input_into_string(unsigned char* string,char* infile)
 \brief Old Kalign function to copy file content into a string. 
 All sequences are read into one array indexed by pointers.
 \param string target string.
 \param infile Pointer to input file.
 
 */

unsigned char* get_input_into_string(unsigned char* string,char* infile)
{
	long int i = 0;
	
	
	FILE *file = 0;
	if (!(file = fopen( infile, "r" ))){
		return 0;
		fprintf(stderr,"Cannot open file '%s'\n", infile);
		exit(-1);
	}
	if (fseek(file,0,SEEK_END) != 0){
		(void)fprintf(stderr, "ERROR: fseek failed\n");
		(void)exit(EXIT_FAILURE);
	}
	i= ftell (file);
	if (fseek(file,0,SEEK_START) != 0){
		(void)fprintf(stderr, "ERROR: fseek failed\n");
		(void)exit(EXIT_FAILURE);
	}
	if(!string){
		MMALLOC(string,(i+1+18)* sizeof(unsigned char));
	
	}
	fread(string,sizeof(unsigned char), i, file);
	string[i] = 0;
	fclose(file);
	
	return string;
}


/** \fn struct fasta* read_fasta(struct fasta* f)
 \brief Old Kalign function to set pointers to start of sequence names, sequences ...
 
 Also sets the sequence lengths. 
 \param string target string.
 \param f @ref fasta .
 
 */


struct fasta* read_fasta(struct fasta* f)
{
	int c = 0;
	int n = 0;
	int i = 0;
	int j = 0;
	int len = 0;
	int stop = 0;
	int nbytes;
	
	nbytes = (int) strlen((char*) f->string);
	

	//aln->org_seq = 0;
	stop = 0;
	
	//count filenames....
	for (i =0;i < nbytes;i++){
		if (f->string[i] == '>'&& stop == 0){
			f->numseq++;
			stop = 1;
		}else if (f->string[i] == '\n'){
			stop = 0;
		}else if(f->string[i] == '\r'){
			f->string[i] = '\n';
		}
	}
	
	
	MMALLOC(f->sn, sizeof(unsigned char*)*f->numseq);
	//aln->c = 0;
	
	MMALLOC(f->s_index, sizeof(int)*(f->numseq+1));
	
	for(i = 0; i < f->numseq;i++){
		f->sn[i] = 0;
		f->s_index[i] = 0;
	}
	f->s_index[f->numseq] = 0;
	
	
	for (i =0;i < nbytes;i++){
		if (f->string[i] == '>'){
			if(f->max_len < len){
				f->max_len = len;
			}
			
			len = 0;
			j = i+1;
			while(f->string[j] != '\n'){
				//	fprintf(stderr,"%c",aln->string[j]);
				j++;
			}
			//fprintf(stderr,"	%d\n",j-i);
			MMALLOC(f->sn[c],sizeof(char)*(j-i));
			f->s_index[c] = n;
			j = i+1;
			len = 0;
			while(f->string[j] != '\n'){
				if(isspace((int)f->string[j])){
					f->sn[c][len] = '_';
				}else{
					f->sn[c][len] = f->string[j];
				}
				len++;
				j++;
			}
			f->sn[c][len] = 0;
			f->string[n] = 'X';
			n++;
			i+= j-i;
			c++;
		}else if(isalnum((int)f->string[i])){// if (aln->string[i] != '\n' && aln->string[i] != 0 ){
			f->string[n] =   nuc_code[(int)f->string[i]];//  toupper(aln->string[i]);
			n++;
			len++;
		}
	}
	
	f->s_index[c] = n;
	f->string[n] = 'X';
	f->string[n+1] = 0;
	
	f->string_len = n+1;
	return f;
}


/** \fn void free_fasta(struct fasta*f)
 
 \brief frees @ref fasta .
 \param f @ref fasta.
 \warning I used to free f->mer_hash even though it is allocated / used in this project.... 
 */
void free_fasta(struct fasta*f)
{
	int i;
	for (i =0;i < f->numseq;i++){
		MFREE(f->sn[i]);
	}
	if(f->mer_hash){
		MFREE(f->mer_hash);
	}
	if(f->suffix){
		MFREE(f->suffix);
	}
	
	//free(aln->mer_hash);
	MFREE(f->s_index);
	MFREE(f->string);
	MFREE(f->sn);
	
	MFREE(f);
}



struct read_info** malloc_read_info(struct read_info** ri, int numseq)
{
	int i;
	MMALLOC(ri, sizeof(struct read_info*) * numseq);
	
	for(i = 0; i < numseq;i++){
		ri[i] = 0;
		MMALLOC(ri[i], sizeof(struct read_info));
		
		ri[i]->seq = 0;
		ri[i]->name = 0;
		ri[i]->qual = 0;
		ri[i]->labels = 0;
		ri[i]->len = 0;
		ri[i]->bar_prob = 0;
		ri[i]->mapq = -1.0;
		ri[i]->barcode = -1;
		ri[i]->fingerprint = -1;
		ri[i]->read_type = 0;
		//ri[i]->strand = malloc(sizeof(unsigned int)* (LIST_STORE_SIZE+1));
		//ri[i]->hits = malloc(sizeof(unsigned int)* (LIST_STORE_SIZE+1));
	}
	return ri;
}

struct read_info** clear_read_info(struct read_info** ri, int numseq)
{
	int i;
	
	for(i = 0; i < numseq;i++){
		if(ri[i]->seq){
			MFREE(ri[i]->seq);
		}
		if(ri[i]->name){
			MFREE(ri[i]->name);
		}
		if(ri[i]->qual){
			MFREE(ri[i]->qual);
		}
		if(ri[i]->labels){
			MFREE(ri[i]->labels);
		}
		ri[i]->seq = 0;
		ri[i]->name = 0;
		ri[i]->qual = 0;
		ri[i]->labels = 0;
		ri[i]->len = 0;
		ri[i]->bar_prob = 0;
		ri[i]->mapq = -1.0;
		ri[i]->barcode = -1;
		ri[i]->fingerprint = -1;
		ri[i]->read_type = 0;
		//ri[i]->strand = malloc(sizeof(unsigned int)* (LIST_STORE_SIZE+1));
		//ri[i]->hits = malloc(sizeof(unsigned int)* (LIST_STORE_SIZE+1));
	}
	return ri;
}

void free_read_info(struct read_info** ri, int numseq)
{
	int i;
	
	for(i = 0; i < numseq;i++){
		//free(ri[i]->strand);
		//free(ri[i]->hits);
		
		
		if(ri[i]->labels){
			MFREE(ri[i]->labels);
		}
		if(ri[i]->name){
			MFREE(ri[i]->name);
		}
		if(ri[i]->seq){
			MFREE(ri[i]->seq);
		}
		if(ri[i]->qual){
			MFREE(ri[i]->qual );
		}
		
		MFREE(ri[i]);
	}
	MFREE(ri);
}













