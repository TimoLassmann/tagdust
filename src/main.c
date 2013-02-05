/*
 
 Copyright (C) 2010 Timo Lassmann <timolassmann@gmail.com>
 
 This file is part of SAMstat.
 
 Delve is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 Delve is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with Delve.  If not, see <http://www.gnu.org/licenses/>.
 
 */


#include "interface.h"
#include "nuc_code.h"
#include "io.h"
#include "tagdust2.h"
#include <math.h>

int main (int argc,char * argv[]) {
	struct parameters* param = 0;
	//struct seq_stats* seq_stats = 0;
	FILE* outfile =0;
	int i;
	
	init_nuc_code();
	
	param = interface(param,argc,argv);
	if(param->summary){
		if ((outfile = fopen(param->summary, "w")) == NULL){
			fprintf(stderr,"can't open output\n");
			exit(-1);
		}
	}
	
	if(!param->infiles && !isatty(0)){
		if(!param->format){
			fprintf(stderr,"No format specified. Use -f <sam | bam | fa | fq > \n");
			exit(-1);
		}
		if(!strcmp("sam", param->format)){
			param->sam = 1;
			//}else if (byg_end(".bam", param->infile[i])  == strlen(param->infile[i])){
		}else if (!strcmp("bam",  param->format)){
			param->sam = 2;
			//}else if (byg_end(".fa", param->infile[i])  == strlen(param->infile[i])){
		}else if (!strcmp("fa", param->format)){
			param->sam = 0;
			//}else if (byg_end(".fq", param->infile[i])  == strlen(param->infile[i])){
		}else if (!strcmp("fq",  param->format)){
			param->sam = 0;
			//}else if (byg_end(".fastq", param->infile[i])  == strlen(param->infile[i])){
		}else if (!strcmp("fastq",  param->format)){
			param->sam = 0;
			//}else if (byg_end(".fastaq", param->infile[i])  == strlen(param->infile[i])){
		}else if (!strcmp("fastaq",  param->format)){
			param->sam = 0;
			//}else if (byg_end(".fasta", param->infile[i])  == strlen(param->infile[i])){
		}else if (!strcmp("fasta",  param->format)){
			param->sam = 0;
		}else{
			param->sam = -1;
		}
		if(param->sam != -1){
		/*	fprintf(stdout,"Working on: stdin\n");
			seq_stats = init_seq_stats(param->kmer_size);
			seq_stats->sam = param->sam;
			if(param->sam == 0){
				seq_stats = collect_data(seq_stats,param,&read_fasta_fastq,-1);
			}else if(param->sam == 2){
				seq_stats = collect_data(seq_stats,param,&read_sam_chunk,-1);
			}else{
				seq_stats = collect_data(seq_stats,param,&read_sam_chunk,-1);
			}
			
			if(sanity_check(seq_stats)){
				if(param->summary){
					print_summary(seq_stats,param,-1,outfile);
				}else{
					print_html_page(seq_stats,param,-1);
				}
			}
			free_seq_stats(seq_stats);*/
		}
	}
	
	for(i = 0; i < param->infiles;i++){
		param->sam = 0;
		//if(byg_end(".sam", param->infile[i])   == strlen(param->infile[i])){
		if(!strcmp(".sam", param->infile[i] + (strlen(param->infile[i] ) - 4))){
			param->sam = 1;
		//}else if (byg_end(".bam", param->infile[i])  == strlen(param->infile[i])){
		}else if (!strcmp(".bam", param->infile[i] + (strlen(param->infile[i] ) - 4))){
			param->sam = 2;
		//}else if (byg_end(".fa", param->infile[i])  == strlen(param->infile[i])){
		}else if (!strcmp(".fa", param->infile[i] + (strlen(param->infile[i] ) - 3))){
			param->sam = 0;
		//}else if (byg_end(".fq", param->infile[i])  == strlen(param->infile[i])){
		}else if (!strcmp(".fq", param->infile[i] + (strlen(param->infile[i] ) - 3))){
			param->sam = 0;
		//}else if (byg_end(".fastq", param->infile[i])  == strlen(param->infile[i])){
		}else if (!strcmp(".fastq", param->infile[i] + (strlen(param->infile[i] ) - 6))){
			param->sam = 0;
		//}else if (byg_end(".fastaq", param->infile[i])  == strlen(param->infile[i])){
		}else if (!strcmp(".fastaq", param->infile[i] + (strlen(param->infile[i] ) - 7))){
			param->sam = 0;
		//}else if (byg_end(".fasta", param->infile[i])  == strlen(param->infile[i])){
		}else if (!strcmp(".fasta", param->infile[i] + (strlen(param->infile[i] ) - 6))){
			param->sam = 0;
		}else if(!strcmp(".sam.gz", param->infile[i] + (strlen(param->infile[i] ) - 7))){
			param->sam = 1;
			param->gzipped  = 1;
			//}else if (byg_end(".bam", param->infile[i])  == strlen(param->infile[i])){
		}else if (!strcmp(".bam.gz", param->infile[i] + (strlen(param->infile[i] ) - 7))){
			param->sam = 2;
			param->gzipped  = 1;
			//}else if (byg_end(".fa", param->infile[i])  == strlen(param->infile[i])){
		}else if (!strcmp(".fa.gz", param->infile[i] + (strlen(param->infile[i] ) - 6))){
			param->sam = 0;
			param->gzipped  = 1;
			//}else if (byg_end(".fq", param->infile[i])  == strlen(param->infile[i])){
		}else if (!strcmp(".fq.gz", param->infile[i] + (strlen(param->infile[i] ) - 6))){
			param->sam = 0;
			param->gzipped  = 1;
			//}else if (byg_end(".fastq", param->infile[i])  == strlen(param->infile[i])){
		}else if (!strcmp(".fastq.gz", param->infile[i] + (strlen(param->infile[i] ) - 9))){
			param->sam = 0;
			param->gzipped  = 1;
			//}else if (byg_end(".fastaq", param->infile[i])  == strlen(param->infile[i])){
		}else if (!strcmp(".fastaq.gz", param->infile[i] + (strlen(param->infile[i] ) - 10))){
			param->sam = 0;
			param->gzipped  = 1;
			//}else if (byg_end(".fasta", param->infile[i])  == strlen(param->infile[i])){
		}else if (!strcmp(".fasta.gz", param->infile[i] + (strlen(param->infile[i] ) - 9))){
			param->sam = 0;
			param->gzipped  = 1;
		}else{
			param->sam = -1;
		}
		//fprintf(stdout,"Loking at on:%s	%d\n",param->infile[i],sam);
		if(param->sam != -1){
		/*	fprintf(stdout,"Working on:%s\n",param->infile[i]);
			seq_stats = init_seq_stats(param->kmer_size);
			seq_stats->sam = param->sam;
			if(param->sam == 0){
				seq_stats = collect_data(seq_stats,param,&read_fasta_fastq,i);
			}else if(param->sam == 2){
				seq_stats = collect_data(seq_stats,param,&read_sam_chunk,i);
			}else{
				seq_stats = collect_data(seq_stats,param,&read_sam_chunk,i);
			}
			
			if(sanity_check(seq_stats)){
				if(param->summary){
					print_summary(seq_stats,param,i,outfile);
				}else{
					print_html_page(seq_stats,param,i);
				}
			}
			free_seq_stats(seq_stats);*/
		}
	}
	
	
	if(param->summary){
		fclose(outfile);
	}
	free_param(param);
	return 0;
}

struct seq_stats* collect_data(struct seq_stats* seq_stats,struct parameters* param,int (*fp)(struct read_info** ,struct parameters*,FILE* ),int file_num)
{
	struct read_info** ri = 0;
	int i,j,c;
	int test = 1;
	int numseq;
	int qual_key = 0;
	int aln_len = 0;
	//char command[1000];
	//char  tmp[1000];
	FILE* file = 0;
	FILE* unmapped = 0;
	
	ri = malloc(sizeof(struct read_info*) * param->num_query);
	
	for(i = 0; i < param->num_query;i++){
		ri[i] = malloc(sizeof(struct read_info));
		ri[i]->seq = 0;
		ri[i]->name = 0;
		ri[i]->qual = 0;
		ri[i]->len = 0;
		ri[i]->cigar = 0;
		ri[i]->md = 0;
		ri[i]->xp = 0;
		ri[i]->priors = 0;// malloc(sizeof(unsigned int)* (LIST_STORE_SIZE+1));
		ri[i]->strand = malloc(sizeof(unsigned int)* (LIST_STORE_SIZE+1));
		ri[i]->hits = malloc(sizeof(unsigned int)* (LIST_STORE_SIZE+1));
		ri[i]->identity = malloc(sizeof(float)* (LIST_STORE_SIZE+1));
	}
	file =  io_handler(file, file_num,param);
	
	if(param->print_unmapped){
		if (!(unmapped = fopen(param->print_unmapped, "w" ))){
			fprintf(stderr,"Cannot open file '%s'\n",param->print_unmapped);
			exit(-1);
		}
	}
	
	while ((numseq = fp(ri, param,file)) != 0){
		//fprintf(stderr,"rread: %d\n",numseq);
		
	}
	
	if(param->print_unmapped){
		fclose(unmapped);
	}
		
	for(i = 0; i < param->num_query;i++){
		free(ri[i]->strand); 
		free(ri[i]->hits);
		free(ri[i]->identity);
		
		if(ri[i]->cigar){
			free(ri[i]->cigar);
		}
		if(ri[i]->md){
			free(ri[i]->md);
		}
		
		free(ri[i]);
	}
	free(ri);
	//fprintf(stderr,"%p\n",file);
	if(param->sam == 2 || param->sam == 1 || param->gzipped ){
		pclose(file);
	}else{
		//if(file_num != -1){
		fclose(file);
		//}
	}
	return seq_stats;
}



int parse_cigar_md(struct read_info* ri,struct seq_stats* seq_stats,int qual_key)
{
	int* read = malloc(sizeof(int)* MAX_SEQ_LEN);
	int* genome = malloc(sizeof(int) * MAX_SEQ_LEN);
	int reverse_int[5]  ={3,2,1,0,4};
	char tmp_num[8];
	int l,i,j,c,rp,gp,sp,exit_loop,aln_len,add; 
	
	for(i = 0; i < MAX_SEQ_LEN;i++){
		genome[i] = 0;
		read[i] = 0;
	}
	
	l = strlen((char*)ri->cigar);
	exit_loop = 0;
	i =0;
	rp = 0;
	sp = 0;
	while(!exit_loop){
		c = 0;
		if(isdigit((int)ri->cigar[i])){
			j = 0;
			while (isdigit(ri->cigar[i])) {
				tmp_num[j] = ri->cigar[i];
				j++;
				i++;
				if(i == l){
					exit_loop =1;
					break;
				}
			}
			tmp_num[j] = 0;
			
			c = atoi(tmp_num);
		}
		if(isalpha((int)ri->cigar[i])){
			switch (ri->cigar[i]) {
				case 'M':
					for(j = 0; j < c;j++){
						read[rp] = ri->seq[sp];
						rp++;
						sp++;
					}
		//			fprintf(stderr,"M:%d\n",c);
					break;
				case 'I':
					for(j = 0; j < c;j++){
						read[rp] = ri->seq[sp];
						genome[rp] = -1;
						rp++;
						sp++;
					}
					
		//			fprintf(stderr,"I:%d\n",c);
					break;
				case 'D':
					for(j = 0; j < c;j++){
						read[rp] = -1;
						rp++;
						
					}
					
		//			fprintf(stderr,"D:%d\n",c);
					break;
				default:
					break;
			}
			c = 0;
		}
		i++;
		if(i == l){
			exit_loop =1;
			break;
		}
		
	}
	aln_len = rp;
	
	i =0;
	rp = 0;
	
	while(read[rp] == -1){
		rp++;
	}
	gp = 0;
	exit_loop = 0;
	add  = 0;
	l = strlen((char*)ri->md);
	
	//int gg;
	
//#1C20^ATT0C3A0
	while(!exit_loop){
		if(isdigit((int)ri->md[i])){
			j = 0;
			while (isdigit(ri->md[i])) {
				tmp_num[j] = ri->md[i];
				j++;
				i++;
				if(i == l){
					exit_loop = 1;
					break;
				}
			}
			tmp_num[j] = 0;
			
			c = atoi(tmp_num);
			
			//fprintf(stderr,"MD:%d\n",c);
			for(j = 0; j < c;j++){
				while(genome[gp] == -1){
					gp++;
					rp++;
				}
				while(read[rp] == -1){
					rp++;
				}
				genome[gp] = read[rp];
				
				gp++;
				rp++;
				//fprintf(stderr,"%d	%d	%d \n",aln_len,rp,i);
				while(read[rp] == -1){
					rp++;
				}
			}
			add = 0;
		}else if(isalpha((int)ri->md[i])){
			//fprintf(stderr,"MD:%c\n",ri->md[i]);
			while(genome[gp] == -1){
				gp++;
			}
			genome[gp] = nuc_code[ri->md[i]];
			gp++;
			i++;
			if(!add){
				rp++;
			}
			
		}else{
			add = 1;
			i++;
		}
		if(i == l){
			exit_loop = 1;
			break;
		}
	}
	
	if(ri->strand[0] == 0){
		gp = 0;
		for(i =0;i < aln_len;i++){
			
			if(read[i] != -1 && genome[i] != -1){
				if(read[i] != genome[i]){
//					seq_stats->mismatches[qual_key][gp][read[i]] += 1;
				}
				gp++;
				//			fprintf(stderr,"Mismatch %d\n",i);
			}else if(read[i] == -1 && genome[i] != -1){
//				seq_stats->deletions[qual_key][gp] += 1;
				
				//			fprintf(stderr,"Deletion %d\n",i);
			}else if(read[i] != -1 && genome[i] == -1){
//				seq_stats->insertions[qual_key][gp][read[i]] += 1;
				gp++;
				//			fprintf(stderr,"Insertion %d\n",i);
			}
		}
	}else{
		gp = ri->len-1;
		for(i = 0;i < aln_len;i++){
			
			if(read[i] != -1 && genome[i] != -1){
				if(read[i] != genome[i]){
//					seq_stats->mismatches[qual_key][gp][reverse_int[read[i]]] += 1;
					//		fprintf(stderr,"Mismatch %d->%d\n",i,gp);
				}
				gp--;
				
			}else if(read[i] == -1 && genome[i] != -1){
//				seq_stats->deletions[qual_key][gp] += 1;
				
				//	fprintf(stderr,"Deletion %d\n",i);
			}else if(read[i] != -1 && genome[i] == -1){
//				seq_stats->insertions[qual_key][gp][reverse_int[read[i]]] += 1;
				gp--;
				//	fprintf(stderr,"Insertion %d\n",i);
			}
		}
	}
	
	free(read);
	free(genome);
	return aln_len;
}












