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
 
 @latexonly
 \begin{figure}[H]
 \includegraphics[scale = 0.7]{../figures/figure1.pdf}
 \caption{Overview of the TagDust workflow. Sequences are labelled according to the HMM architecture and relevant information written to the output.}
 \end{figure}
 
 @endlatexonly
 
 
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

#ifndef MMALLOC
#include "malloc_macro.h"
#endif

#include "interface.h"
#include "nuc_code.h"
#include "io.h"
#include "misc.h"
#include "tagdust2.h"
#include "barcode_hmm.h"
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
	
	sprintf(param->buffer,"Start Run\n--------------------------------------------------\n");
	param->messages = append_message(param->messages, param->buffer);

	// Paired end or single end ?
	if(param->infiles == 0){
		fprintf(stderr,"Sorry - no input file found.\n");
		free_param(param);
		exit(EXIT_FAILURE);
	}else	if(param->infiles == 1){
		sprintf(param->buffer,"Running in single end mode.\n");
		param->messages = append_message(param->messages, param->buffer);
		if(param->arch_file){
			param = test_architectures(param, 0);
		}
		if (param->read_structure->num_segments == 0){
			filter_controller(param,0);
			
		}else{
			hmm_controller(param,0);
		}
	}else if(param->infiles == 2){
		sprintf(param->buffer,"Running in paired end mode.\n");
		param->messages = append_message(param->messages, param->buffer);
		
		if(param->arch_file){
			param = test_architectures(param, 0);
			param->read_structure_R1 = param->read_structure;
			param->read_structure = 0;
			fprintf(stderr,"ARCH1: %d segments\n", param->read_structure_R1->num_segments);
			for(i = 0; i < param->read_structure_R1->num_segments;i++){
				fprintf(stderr,"%c",  param->read_structure_R1->type[i]);
			}
			fprintf(stderr,"\n");
		}
		
		if(param->arch_file){
			param->read_structure = malloc_read_structure();
			param = test_architectures(param, 1);
			param->read_structure_R2 = param->read_structure;
			param->read_structure = 0;
			fprintf(stderr,"ARCH2: %d segments\n", param->read_structure_R2->num_segments);
			for(i = 0; i < param->read_structure_R2->num_segments;i++){
				fprintf(stderr,"%c",  param->read_structure_R2->type[i]);
			}
			fprintf(stderr,"\n");

			
		}
		
		hmm_controller_pe(param);
		
	}else{
		fprintf(stderr,"Sorry - using more than two input files is presently disabled.\n");
		free_param(param);
		exit(EXIT_FAILURE);
	}
	
	if(param->join){
		concatenate_reads(param,&read_fasta_fastq);
		free_param(param);
		return EXIT_SUCCESS;
	}
	if(param->split){
		split(param,&read_fasta_fastq);
		free_param(param);
		return EXIT_SUCCESS;
	}
	
	if(param->infiles > 1){
		
	}
	
	
	
	
	
	
	free_param(param);
	return EXIT_SUCCESS;
	
}












