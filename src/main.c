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



/**
 @mainpage TagDust 2
 @author Timo Lassmann
 
 
 Raw sequences produced by next generation sequencing (NGS) machines may contain adapter, linker, barcode and fingerprint sequences. TagDust2 is a program to extract and correctly label the sequences to be mapped in downstream pipelines.
 
 TagDust allows users to specify the expected architecture of a read and converts it into a hidden Markov model. The latter can assign sequences to a particular barcode (or index) even in the presence of sequencing errors. Sequences not matching the architecture (primer dimers, contaminants etc.) are automatically discarded).

 \image html figure1.jpg
  \image latex figure1.eps "Overview of the TagDust workflow. Sequences are labelled according to the HMM architecture and relevant information written to the output." width=10cm
 
 
Copyright 2013 Timo Lassmann (timolassmann@gmail.com)
 
 This document is free;  you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This document is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with TagDust.
 
 If not, see (http://www.gnu.org/licenses/).

 
 */



/*! \file main.c 
 \brief Figures out the nature of the input and calls the main functions. 

 Calls parameter parser. Calls functions to process the data. \author Timo Lassmann \bug No known bugs.
 */


#include "interface.h"
#include "nuc_code.h"
#include "io.h"
#include "misc.h"
#include "tagdust2.h"
#include "barcode_hmm.h"
#include "exact.h"
#include "sim.h"
#include <math.h>


/*! \brief  Initializes nucleotide alphabet needed to parse input. Calls parameter parser. Calls functions to process the data. 
 * \param argc number of command line parameters
* \param argv command line parameters
 * \return EXIT_SUCCESS */
int main (int argc,char * argv[]) {
	struct parameters* param = 0;

	int i;
	
	init_nuc_code();

	param = interface(param,argc,argv);
	
	
	if(param->sim){
		simulate(param);
	}
	
	for(i = 0; i < param->infiles;i++){
		param->sam = 0;
		if(!strcmp(".sam", param->infile[i] + (strlen(param->infile[i] ) - 4))){
			param->sam = 1;
		}else if (!strcmp(".bam", param->infile[i] + (strlen(param->infile[i] ) - 4))){
			param->sam = 2;
		}else if (!strcmp(".fa", param->infile[i] + (strlen(param->infile[i] ) - 3))){
			param->sam = 0;
			param->fasta = 1;
		}else if (!strcmp(".fq", param->infile[i] + (strlen(param->infile[i] ) - 3))){
			param->sam = 0;
		}else if (!strcmp(".fastq", param->infile[i] + (strlen(param->infile[i] ) - 6))){
			param->sam = 0;
		}else if (!strcmp(".fastaq", param->infile[i] + (strlen(param->infile[i] ) - 7))){
			param->sam = 0;
		}else if (!strcmp(".fasta", param->infile[i] + (strlen(param->infile[i] ) - 6))){
			param->sam = 0;
			param->fasta = 1;
		}else if(!strcmp(".sam.gz", param->infile[i] + (strlen(param->infile[i] ) - 7))){
			param->sam = 1;
			param->gzipped  = 1;
		}else if (!strcmp(".bam.gz", param->infile[i] + (strlen(param->infile[i] ) - 7))){
			param->sam = 2;
			param->gzipped  = 1;
		}else if (!strcmp(".fa.gz", param->infile[i] + (strlen(param->infile[i] ) - 6))){
			param->sam = 0;
			param->fasta = 1;
			param->gzipped  = 1;
		}else if (!strcmp(".fq.gz", param->infile[i] + (strlen(param->infile[i] ) - 6))){
			param->sam = 0;
			param->gzipped  = 1;
		}else if (!strcmp(".fastq.gz", param->infile[i] + (strlen(param->infile[i] ) - 9))){
			param->sam = 0;
			param->gzipped  = 1;
		}else if (!strcmp(".fastaq.gz", param->infile[i] + (strlen(param->infile[i] ) - 10))){
			param->sam = 0;
			param->gzipped  = 1;
		}else if (!strcmp(".fasta.gz", param->infile[i] + (strlen(param->infile[i] ) - 9))){
			param->sam = 0;
			param->gzipped  = 1;
		}else if (!strcmp(".fastq.bz2", param->infile[i] + (strlen(param->infile[i] ) - 10))){
			param->sam = 0;
			param->bzipped  = 1;
		}else if (!strcmp(".fq.bz2", param->infile[i] + (strlen(param->infile[i] ) - 7))){
			param->sam = 0;
			param->bzipped  = 1;
		}else{
			param->sam = -1;
		}
		
		if(param->sam != -1){
			if(param->sam == 0){
				if(param->exact5){
					exact_controller(param,&read_fasta_fastq,i);
				}else{
					hmm_controller(param,&read_fasta_fastq,i);
				}
			}else{
				if(param->exact5){
					exact_controller(param,&read_sam_chunk,i);
				}else{
					hmm_controller(param,&read_sam_chunk,i);
				}
			}
		
		}
	}
	free_param(param);
	return EXIT_SUCCESS;
}












