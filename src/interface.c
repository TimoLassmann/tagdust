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


#include "tagdust2.h"
#include "interface.h"
#include "io.h"
#include "misc.h"



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
	param->sam = 0;
	param->train = 0;
	param->fasta = 0;
	param->matchstart = -1;
	param->matchend = -1;
	param->minlen = 16;
	
	param->sequencer_error_rate = 0.01f;
	param->indel_frequency = 0.0f;
	param->average_read_length = 50;
	
	param->confidence_threshold = 0.99;//ence
	
	param->read_structure = 0;
	
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
			{"end",required_argument,0, OPT_END},
			{"threshold",required_argument,0, OPT_THRESHOLD},
			{"out",required_argument,0, 'o'},
			//{"format",required_argument,0, OPT_FORMAT},
			
			{"filter",required_argument,0, 'f'},
			
			{"quiet",0,0,'q'},
			{"help",0,0,'h'},
			{0, 0, 0, 0}
		};
		
		int option_index = 0;
		c = getopt_long_only (argc, argv,"e:o:p:qhf:t:",long_options, &option_index);
		
		if (c == -1){
			break;
		}
		
		switch(c) {
			case 0:
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
				
			case 'f':
				param->filter = optarg;
				break;
			case 'o':
				param->outfile = optarg;
				break;
			case 'e':
				param->sequencer_error_rate = atof(optarg); //0.01f;
				break;
			case 't':
				param->num_threads = atoi(optarg);
				break;
			case 'q':
				param->quiet_flag = 1;
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
	for(i = 0; i < 10;i++){
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
						free_param(param);
						exit(-1);
					}
					//assert(strlen(param->read_structure->sequence_matrix[i][g]) == strlen(param->read_structure->sequence_matrix[i][f]));
					//fprintf(stderr,"\t%s\n",param->read_structure->sequence_matrix[i][g] );
				}
				//fprintf(stderr,"\t%s\n",param->read_structure->sequence_matrix[i][g] );
			}
			
			
			
			fprintf(stderr,"Found %c segment %d with %d sequences\n",param->read_structure->type[i] ,i, param->read_structure->numseq_in_segment[i] );
			for(g = 0;g < param->read_structure->numseq_in_segment[i];g++){
				fprintf(stderr,"\t%s\n",param->read_structure->sequence_matrix[i][g] );
			}
			last = i;
			
		}
	}
	
	//exit(0);
	
	if(help){
		usage();
		free(param);
		exit(0);
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

struct parameters* assign_segment_sequences(struct parameters* param, char* tmp, int segment)
{
	int i,f,g;
	int count;
	//tmp = optarg;
	count = byg_count(",", tmp);
	fprintf(stderr,"Segment %d: %d	sequences\n",segment,count+1);
	param->read_structure->numseq_in_segment[segment] = count+1;
	param->read_structure->sequence_matrix[segment] = malloc(sizeof(char*) * (count+1));
	assert(param->read_structure->sequence_matrix[segment] !=0);
	for(i = 0; i < count+1;i++){
		param->read_structure->sequence_matrix[segment][i] = malloc(sizeof(char)* 32);
	}
	param->read_structure->type[segment] = tmp[0];
	
	if(tmp[0] == 'R'){
		param->read_structure->sequence_matrix[segment][0][0]  = 'N';
		param->read_structure->sequence_matrix[segment][0][1]  = 0;
		
	}else{
	
	f = 0;
	g = 0;
	for(i = 2;i < strlen(tmp);i++){
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
	if(segment+1 >param->read_structure->num_segments  ){
		param->read_structure->num_segments = segment+1;
	}
	return param;
}


void usage()
{
	fprintf(stdout, "\nTagDust %0.2f, Copyright (C) 2010, 2011 Timo Lassmann <timolassmann@gmail.com>\n", VERSION);
	fprintf(stdout, "\n");
	fprintf(stdout, "Usage:   tagdust <file.sam> <file.bam> <file.fa> <file.fq> .... \n\n");
	fprintf(stdout, "Options:\n");
	fprintf(stdout, "         -s             STR    prints summary statistics of multiple libraries into one FILE [NA].\n");
	fprintf(stdout, "         -k             INT    treat all reads mapping more than k errors as unmapped [off].\n");
	fprintf(stdout, "         -n             STR    name of library and output file when reading from stdin\n");
	fprintf(stdout, "         -f             STR    specifies the format when reading from stdin.\n");
	fprintf(stdout, "         -F             INT    filter used in samtools [768].\n");
	//fprintf(stdout, "Options: -unmapped      STR    print unmapped reads to Length of simulated reads [30]\n");
	fprintf(stdout, "\n");
	fprintf(stdout, "Example: Reading from standard input:\n\n");
	fprintf(stdout, "	cat library.sam | samstat -f sam -n library_one\n");
	fprintf(stdout, "	or:\n");
	fprintf(stdout, "	samtools view -ub  ~/tmp/small.bam   | ./samstat -f bam -n library_one\n\n");
	fprintf(stdout, "	In both cases SAMStat will list the stats in the file \"library_one.html\"\n");
}


void free_param(struct parameters* param)
{
	int i,j;
	for(i = 0; i < 5;i++){
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

