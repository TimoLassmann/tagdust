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
 along with Delve.  If not, see <http://www.gnu.org/licenses/>.
 
 */


#include "interface.h"
#include "nuc_code.h"
#include "io.h"
#include "misc.h"
#include "tagdust2.h"
#include "barcode_hmm.h"
#include "pst.h"
#include <math.h>


int main (int argc,char * argv[]) {
	struct parameters* param = 0;
	//struct seq_stats* seq_stats = 0;
	FILE* outfile =0;
	int i,j;
	
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
			fprintf(stdout,"Working on:%s\n",param->infile[i]);
			if(param->sam == 0){
				pst_controller(param,&read_fasta_fastq,i);
			}else{
				pst_controller(param,&read_sam_chunk,i);
			}
		
		}
	}
	
	
	if(param->summary){
		fclose(outfile);
	}
	free_param(param);
	return 0;
}












