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


/*! \file interface.c
 \brief Functions to deal with user inputs.
 */


#include "tagdust2.h"
#include "interface.h"
#include "io.h"
#include "misc.h"


/** \fn struct parameters* interface(struct parameters* param,int argc, char *argv[])
 \brief Read command line options into @ref parameters.
 
 \param param nucleotide sequence.
 \param argc number of command line arguments.
  \param argv command line arguments.
 */
struct parameters* interface(struct parameters* param,int argc, char *argv[])
{
	int i,c,f,g;
	int help = 0;
	int last;
	
	if (argc < 2 && isatty(0)){
		usage();
		exit(0);
	}
		
	param = malloc(sizeof(struct parameters));
	param->infiles = 0;
	param->infile = 0;
	param->outfile = 0;
	param->num_threads = 8;

	param->quiet_flag = 0;
	param->num_query = 1000000;
	param->format = 0;
	param->gzipped = 0;
	param->bzipped = 0;
	param->sam = 0;
	param->train = 0;
	param->fasta = 0;
	param->matchstart = -1;
	param->matchend = -1;
	param->minlen = 16;
	param->exact5 = 0;
	param->sim = 0;
	param->dust = 100;
	
	param->sequencer_error_rate = 0.05f;
	param->indel_frequency = 0.1f;
	param->average_read_length = 50;
	param->numbarcode = 8;
	param->confidence_threshold = 20.0;//ence
	
	param->read_structure = 0;
	param->filter_error = 2;
	param->reference_fasta  = 0;
	param->random_prior = 0;
	
	
	
	param->sim_3seq = 0;
	param->sim_5seq = 0;
	param->sim_barlen = 0;
	param->sim_barnum = 0;
	param->sim_error_rate = 0.0f;
	param->sim_InDel_frac = 0.0f;
	param->sim_numseq = 0;
	param->sim_random_frac = 0.0f;
	param->sim_readlen =0;
	param->sim_readlen_mod = 0;
	param->sim_sequenced_len = 0;
	param->log = 0;
	param->print_artifact = 0;
	
	
	
	param->read_structure = malloc(sizeof(struct read_structure));
	assert(param->read_structure !=0);
	param->read_structure->sequence_matrix = malloc(sizeof(char**) * 10 );
	assert(param->read_structure->sequence_matrix !=0);
	
	param->read_structure->numseq_in_segment  = malloc(sizeof(int) * 10);
	assert(param->read_structure->numseq_in_segment !=0);
	param->read_structure->type = malloc(sizeof(char) * 10 );
	
	assert(param->read_structure->type !=0);
	
	for(i = 0;i  <10;i++){
		param->read_structure->sequence_matrix[i] = 0;
		param->read_structure->numseq_in_segment[i] = 0;
		param->read_structure->type[i] = 0;
		
	}
	
	param->read_structure->num_segments = 0;
		
	while (1){	 
		static struct option long_options[] ={
			{"1",required_argument,0, OPT_SEG1},
			{"2",required_argument,0, OPT_SEG2},
			{"3",required_argument,0, OPT_SEG3},
			{"4",required_argument,0, OPT_SEG4},
			{"5",required_argument,0, OPT_SEG5},
			{"6",required_argument,0, OPT_SEG6},
			{"7",required_argument,0, OPT_SEG7},
			{"8",required_argument,0, OPT_SEG8},
			{"9",required_argument,0, OPT_SEG9},
			{"10",required_argument,0, OPT_SEG10},
			{"train",required_argument,0, OPT_TRAIN},
			{"format",required_argument,0, OPT_FORMAT},
			{"minlen",required_argument,0, OPT_MINLEN},
			{"start",required_argument,0, OPT_START},
			{"exact5",required_argument,0, OPT_EXACT5},
			{"simulation",required_argument,0, OPT_SIM},
			{"numbarcode",required_argument,0, OPT_NUMBARCODE},
			{"end",required_argument,0, OPT_END},
			{"threshold",required_argument,0, 'q'},
			{"fe",required_argument,0,OPT_FILTER_ERROR},
			{"ref",required_argument,0,OPT_FILTER_REFERENCE},
			{"dust",required_argument,0,OPT_DUST},
			{"out",required_argument,0, 'o'},
			{"filter",required_argument,0, 'f'},
			{"sim_barlen",required_argument,0,OPT_sim_barlen},
			{"sim_barnum",required_argument,0,OPT_sim_barnum},
			{"sim_5seq",required_argument,0,OPT_sim_5seq},
			{"sim_3seq",required_argument,0,OPT_sim_3seq},
			{"sim_readlen",required_argument,0,OPT_sim_readlen},
			{"sim_readlen_mod",required_argument,0,OPT_sim_readlen_mod},
			{"sim_error_rate",required_argument,0,OPT_sim_error_rate},
			{"sim_InDel_frac",required_argument,0,OPT_sim_InDel_frac},
			{"sim_numseq",required_argument,0,OPT_sim_numseq},
			{"sim_random_frac",required_argument,0,OPT_sim_random_frac},
			{"sim_sequenced_len",required_argument,0,OPT_sim_sequenced_len},
			{"help",0,0,'h'},
			{"log",0,0,'l'},
			{0, 0, 0, 0}
		};
		
		int option_index = 0;
		c = getopt_long_only (argc, argv,"Q:e:o:p:q:hf:t:i:la:",long_options, &option_index);
		
		if (c == -1){
			break;
		}
		
		switch(c) {
			case 0:
				break;
			case OPT_sim_sequenced_len:
				param->sim_sequenced_len = atoi(optarg);
				break;
			case OPT_sim_barlen:
				param->sim_barlen = atoi(optarg);
				break;
			case OPT_sim_barnum:
				param->sim_barnum = atoi(optarg);
				break;
			case OPT_sim_5seq:
				param->sim_5seq = optarg;
				break;
			case OPT_sim_3seq:
				param->sim_3seq = optarg;
				break;
			case OPT_sim_readlen:
				param->sim_readlen = atoi(optarg);
				break;
			case OPT_sim_readlen_mod:
				param->sim_readlen_mod = atoi(optarg);
				break;
			case OPT_sim_error_rate:
				param->sim_error_rate = atof(optarg);
				break;
			case OPT_sim_InDel_frac:
				param->sim_InDel_frac = atof(optarg);
				break;
			case OPT_sim_numseq:
				param->sim_numseq = atoi(optarg);
				param->sim = 1;
				break;
			case OPT_sim_random_frac:
				param->sim_random_frac = atof(optarg);
				break;
				
			case OPT_SEG1:
				param = assign_segment_sequences(param, optarg , 0 );
				break;
			case OPT_SEG2:
				param = assign_segment_sequences(param, optarg , 1 );
				break;
			case OPT_SEG3:
				param = assign_segment_sequences(param, optarg , 2 );
				break;
			case OPT_SEG4:
				param = assign_segment_sequences(param, optarg , 3 );
				break;
			case OPT_SEG5:
				param = assign_segment_sequences(param, optarg , 4 );
				break;
			case OPT_SEG6:
				param = assign_segment_sequences(param, optarg , 5 );
				break;
			case OPT_SEG7:
				param = assign_segment_sequences(param, optarg , 6 );
				break;
			case OPT_SEG8:
				param = assign_segment_sequences(param, optarg , 7 );
				break;
			case OPT_SEG9:
				param = assign_segment_sequences(param, optarg , 8 );
				break;
			case OPT_SEG10:
				param = assign_segment_sequences(param, optarg , 9 );
				break;
			case OPT_TRAIN:
				param->train = optarg;
				break;
			case OPT_FORMAT:
				param->format = optarg;
				break;
			case OPT_START:
				param->matchstart = atoi(optarg)-1;
				break;
			case OPT_END:
				param->matchend = atoi(optarg);
				break;
			case OPT_THRESHOLD:
				param->confidence_threshold = atof(optarg);
				break;
			case OPT_EXACT5:
				param->exact5 = optarg;
				break;
			case OPT_MINLEN:
				param->minlen = atoi(optarg);
				break;
			case OPT_SIM:
				param->sim = atoi(optarg);
				break;
			case OPT_NUMBARCODE:
				param->numbarcode = atoi(optarg);
				break;
			case OPT_FILTER_ERROR:
				param->filter_error = atoi(optarg);
				break;
			case OPT_FILTER_REFERENCE:
				param->reference_fasta = optarg;
				break;
			case OPT_DUST:
				param->dust = atoi(optarg);
				break;
			case 'l':
			case 'L':
				param->log = 1;
				break;
			case 'f':
				param->filter = optarg;
				break;
			case 'a':
				param->print_artifact = optarg;
				break;
				
		
			case 'o':
				param->outfile = optarg;
				break;
			case 'e':
				param->sequencer_error_rate = atof(optarg); //0.01f;
				break;
			case 'i':
				param->indel_frequency = atof(optarg); //0.01f;
				break;
			case 't':
				param->num_threads = atoi(optarg);
				break;
			case 'q':
				param->confidence_threshold = atof(optarg);
				//param->quiet_flag = 1;
				break;
				
			case 'Q':
				param->confidence_threshold = atof(optarg);
				//param->quiet_flag = 1;
				break;
			case 'h':
				help = 1;
				break;
			case '?':
				exit(1);
				break;
			default:
				fprintf(stderr,"default\n\n\n\n");
				abort ();
		}
	}
	
	//if(param->matchstart)
	//fprintf(stderr,"Viterbi: %d\n",param->viterbi);
	last = -1;
	c = 0;
	for(i = 0; i < param->read_structure->num_segments;i++){
		if(param->read_structure->type[i] == 'R'){
			c = 1;
		}
		
		if(param->read_structure->sequence_matrix[i]){
			if(last +1 != i){
				fprintf(stderr,"ERROR: a hmm building lock was skipped??\n");
				free_param(param);
				exit(-1);
			}
			
			//serious checking...
			for(g = 0;g < param->read_structure->numseq_in_segment[i];g++){
				for(f = g+1;f < param->read_structure->numseq_in_segment[i];f++){
					if(strlen(param->read_structure->sequence_matrix[i][g]) != strlen(param->read_structure->sequence_matrix[i][f])){
						fprintf(stderr,"ERROR: the sequences in the same segment have to have the same length\n");
						fprintf(stderr,"%dseq\n%s	%d\n%s	%d\n",param->read_structure->numseq_in_segment[i],param->read_structure->sequence_matrix[i][g],g,param->read_structure->sequence_matrix[i][f],f );
						free_param(param);
						exit(-1);
					}
				}
			}
			
			fprintf(stderr,"Found %c segment %d with %d sequences\n",param->read_structure->type[i] ,i, param->read_structure->numseq_in_segment[i] );
			for(g = 0;g < param->read_structure->numseq_in_segment[i];g++){
				fprintf(stderr,"%s, ",param->read_structure->sequence_matrix[i][g] );
			}
			fprintf(stderr,"\n" );
			last = i;
			
		}
	}
	
	/*if(!c  && param->sim == 0){
		fprintf(stderr,"ERROR: no read segment specified!%d\n",param->read_structure->num_segments);
		
		free_param(param);
		exit(EXIT_FAILURE);
	}*/
	
	if(help){
		usage();
		free(param);
		exit(EXIT_SUCCESS);
	}
	param->infile = malloc(sizeof(char*)* (argc-optind));
	
	c = 0;
	while (optind < argc){
		param->infile[c] =  argv[optind++];
		c++;
	}
	param->infiles = c;
	return param;
}

/** \fn struct parameters* assign_segment_sequences(struct parameters* param, char* tmp, int segment)
 \brief Assigns sequences to segments.
 
 \param param @ref parameters 
 \param tmp sequences from user input.
 \param segment number of segment.
 */
struct parameters* assign_segment_sequences(struct parameters* param, char* tmp, int segment)
{
	int i,f,g;
	int count;
	int len;
	//tmp = optarg;
	count = byg_count(",", tmp);
	
	if(tmp[0] == 'B'){ // add extra space for all N barcode....
		count++;
	}
	f = 0;
	
	//fprintf(stderr,"Segment %d: %d	sequences\n",segment,count+1);
	param->read_structure->numseq_in_segment[segment] = count+1;
	param->read_structure->sequence_matrix[segment] = malloc(sizeof(char*) * (count+1));
	assert(param->read_structure->sequence_matrix[segment] !=0);
	for(i = 0; i < count+1;i++){
		param->read_structure->sequence_matrix[segment][i] = malloc(sizeof(char)* strlen(tmp));
	}
	param->read_structure->type[segment] = tmp[0];
	
	if(tmp[0] == 'R'){
		param->read_structure->sequence_matrix[segment][0][0]  = 'N';
		param->read_structure->sequence_matrix[segment][0][1]  = 0;
		
	}else{
	
		f = 0;
		g = 0;
		len = (int)strlen(tmp) ;
		//fprintf(stderr,"%d\n",len);
		for(i = 2;i <  len;i++){
			if(tmp[i] != ','){
				param->read_structure->sequence_matrix[segment][f][g] = tmp[i];
				g++;
			}else{
				param->read_structure->sequence_matrix[segment][f][g] = 0;
				f++;
				g = 0;
			}
		}
		param->read_structure->sequence_matrix[segment][f][g] = 0;
	}
	
	if(tmp[0] == 'B'){
		f = f + 1;
		g = 0;
		for(i = 0; i < strlen(param->read_structure->sequence_matrix[segment][0]);i++){
			param->read_structure->sequence_matrix[segment][f][g] = 'N';
			g++;
		}
		param->read_structure->sequence_matrix[segment][f][g] = 0;
	}
	
	
	if(segment+1 >param->read_structure->num_segments  ){
		param->read_structure->num_segments = segment+1;
	}
	return param;
}


/** \fn void usage()
 \brief Prints usage.
 */
void usage()
{
	fprintf(stdout, "\n%s %s, Copyright (C) 2013 Timo Lassmann <%s>\n",PACKAGE_NAME, PACKAGE_VERSION,PACKAGE_BUGREPORT);
	fprintf(stdout, "\n");
	fprintf(stdout, "Usage:   tagdust [options] <file>  .... \n\n");
	fprintf(stdout, "Options:\n");
	
	fprintf(stdout, "         -Q         FLT     confidence threshold [20].\n");
	fprintf(stdout, "         -l         NA      write a log file.\n");
	fprintf(stdout, "         -start     INT     start of search area [0].\n");
	fprintf(stdout, "         -end       INT     end of search area [length of sequence].\n");
	fprintf(stdout, "         -format    STR     format of input sequence file.\n");
	fprintf(stdout, "         -minlen    INT     minimal accepted read length [16].\n");
	fprintf(stdout, "         -ref       STR     reference fasta file to be compared against[].\n");
	fprintf(stdout, "         -fe        INT     number of errors allowed when comparing to reference[2].\n");
	fprintf(stdout, "         -dust      INT     remove low complexity sequences. [100].\n");
	fprintf(stdout, "         -e         FLT     expected sequencer error rate [0.05].\n");
	fprintf(stdout, "         -o         STR     output file name.\n");
	fprintf(stdout, "         -a         STR     output file for artifacts [NA].\n");
	fprintf(stdout, "         -t         INT     number of threads [8].\n");
	fprintf(stdout, "         -1         STR     type of the first HMM building block.\n");
	fprintf(stdout, "         -2         STR     type of the second HMM building block.\n");
	fprintf(stdout, "         -...       STR     type of the . . . HMM building block.\n");
	fprintf(stdout, "\n");

}


/** \fn void free_param(struct parameters* param)
 \brief free @ref parameters.
 */
void free_param(struct parameters* param)
{
	int i,j;
	for(i = 0; i < 10;i++){
		if(param->read_structure->sequence_matrix[i]){
			for(j = 0; j < param->read_structure->numseq_in_segment[i];j++){
				free(param->read_structure->sequence_matrix[i][j]);
			}
			
			free(param->read_structure->sequence_matrix[i]);
		}
	}
	free(param->read_structure->sequence_matrix);
	
	free(param->read_structure->numseq_in_segment );
	free(param->read_structure->type);
	free(param->read_structure);
	free(param->infile);
	free(param);
}

