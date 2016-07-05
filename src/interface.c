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



#if HAVE_CONFIG_H
#include "config.h"
#endif

#include "kslib.h"
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
	int i,j,c;
	int help = 0;
	int version = 0;
	
	int status;
	
	//int last;
	
	if (argc < 2 && isatty(0)){
		usage();
		return NULL;
	//	return kslOK;
	}
		
	MMALLOC(param,sizeof(struct parameters));
	param->infiles = 0;
	param->infile = NULL;
	param->outfile = NULL;
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
	param->confidence_threshold = 0.0;//ence
	param->confidence_thresholds = NULL;//ence
	param->confidence_threshold_R1 = 0.0;
	param->confidence_threshold_R2 = 0.0;
	
	param->read_structure = 0;
	param->read_structure_R1 = 0;
	param->read_structure_R2 = 0;
	param->filter_error = 2;
	param->reference_fasta  = 0;
	param->random_prior = 0;
	
	param->print_seq_finger = 0;
	
	
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
	param->sim_end_loss = 0;
	param->log = 0;
	param->print_artifact = 0;
	param->multiread = 0;
	param->join = 0;
	param->split = 0;
	param->messages = 0;
	param->buffer = 0;
	param->seed = 0;
	
	param->arch_file = NULL;
	
	param->read_structures = NULL;
	if((param->read_structure = malloc_read_structure()) == NULL) KSLIB_XEXCEPTION_SYS(kslEMEM,"Malloc of readstructure failed.\n");
	
		
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
			{"name",required_argument,0, OPT_FORMAT},
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
			{"sim_endloss",required_argument,0,OPT_sim_endloss},
			{"arch",required_argument,0,OPT_archfile},
			{"seed",required_argument,0, OPT_seed},
			{"show_finger_seq",0,0,OPT_show_finger_seq},
			{"join",0,0,OPT_join_paired},
			{"split",0,0,OPT_split},
			{"help",0,0,'h'},
			{"version",0,0,'v'},
			{"log",required_argument,0,'l'},
			{0, 0, 0, 0}
		};
		
		int option_index = 0;
		c = getopt_long_only (argc, argv,"Q:e:o:p:q:hvf:t:i:l:L:a:",long_options, &option_index);
		
		if (c == -1){
			break;
		}
		
		switch(c) {
			case 0:
				break;
			case OPT_seed:
				param->seed = atoi(optarg);
				break;
			case OPT_archfile:
				param->arch_file = optarg;
				break;
			case	OPT_split:
				param->split = 1;
				break;
			case OPT_join_paired:
				param->join = 1;
				break;
			case OPT_sim_endloss:
				param->sim_end_loss = atoi(optarg);
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
				if((status = assign_segment_sequences(param, optarg , 0 )) != kslOK) KSLIB_XEXCEPTION(kslFAIL,"Some problem with parsing an HMM segment: %s.\n",optarg);
				//param->read_structure = assign_segment_sequences(param->read_structure, optarg , 0 );
				break;
			case OPT_SEG2:
				if((status = assign_segment_sequences(param, optarg , 1 )) != kslOK) KSLIB_XEXCEPTION(kslFAIL,"Some problem with parsing an HMM segment: %s.\n",optarg);
				//param->read_structure = assign_segment_sequences(param->read_structure, optarg , 1 );
				break;
			case OPT_SEG3:
				if((status = assign_segment_sequences(param, optarg , 2 )) != kslOK) KSLIB_XEXCEPTION(kslFAIL,"Some problem with parsing an HMM segment: %s.\n",optarg);
				//param->read_structure = assign_segment_sequences(param->read_structure, optarg , 2 );
				break;
			case OPT_SEG4:
				if((status = assign_segment_sequences(param, optarg , 3 )) != kslOK) KSLIB_XEXCEPTION(kslFAIL,"Some problem with parsing an HMM segment: %s.\n",optarg);
				//param->read_structure = assign_segment_sequences(param->read_structure, optarg , 3 );
				break;
			case OPT_SEG5:
				if((status = assign_segment_sequences(param, optarg , 4 )) != kslOK) KSLIB_XEXCEPTION(kslFAIL,"Some problem with parsing an HMM segment: %s.\n",optarg);
				//param->read_structure = assign_segment_sequences(param->read_structure, optarg , 4 );
				break;
			case OPT_SEG6:
				if((status = assign_segment_sequences(param, optarg , 5 )) != kslOK) KSLIB_XEXCEPTION(kslFAIL,"Some problem with parsing an HMM segment: %s.\n",optarg);
				//param->read_structure = assign_segment_sequences(param->read_structure, optarg , 5 );
				break;
			case OPT_SEG7:
				if((status = assign_segment_sequences(param, optarg , 6 )) != kslOK) KSLIB_XEXCEPTION(kslFAIL,"Some problem with parsing an HMM segment: %s.\n",optarg);
				//param->read_structure = assign_segment_sequences(param->read_structure, optarg , 6 );
				break;
			case OPT_SEG8:
				if((status = assign_segment_sequences(param, optarg , 7 )) != kslOK) KSLIB_XEXCEPTION(kslFAIL,"Some problem with parsing an HMM segment: %s.\n",optarg);
				//param->read_structure = assign_segment_sequences(param->read_structure, optarg , 7 );
				break;
			case OPT_SEG9:
				if((status = assign_segment_sequences(param, optarg , 8 )) != kslOK) KSLIB_XEXCEPTION(kslFAIL,"Some problem with parsing an HMM segment: %s.\n",optarg);
				//param->read_structure = assign_segment_sequences(param->read_structure, optarg , 8 );
				break;
			case OPT_SEG10:
				if((status = assign_segment_sequences(param, optarg , 9 )) != kslOK) KSLIB_XEXCEPTION(kslFAIL,"Some problem with parsing an HMM segment: %s.\n",optarg);
				//param->read_structure = assign_segment_sequences(param->read_structure, optarg , 9 );
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
				
			case OPT_show_finger_seq:
				param->print_seq_finger = 1;
				break;
			case 'l':
			case 'L':
				param->log = optarg;
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
			case 'v':
				version = 1;
				break;
			case '?':
				exit(1);
				break;
			default:
				fprintf(stderr,"default\n\n\n\n");
				abort ();
		}
	}
	if(help){
		usage();
		free_param(param);
		return NULL;
		//
		//return kslOK;
	}
	
	if(version){
		fprintf(stdout,"%s %s\n",PACKAGE_NAME,PACKAGE_VERSION);
		
		free_param(param);
		return NULL;
	}
	
	
	MMALLOC(param->buffer,sizeof(char) * kslibMSGBUFSIZE);
	
	snprintf(param->buffer , kslibMSGBUFSIZE,"%s %s, Copyright (C) 2013 Timo Lassmann <%s>\n",PACKAGE_NAME, PACKAGE_VERSION,PACKAGE_BUGREPORT);
	param->messages = append_message(param->messages, param->buffer  );
	//command_line[0] = 0;
	//c = 0;
	param->buffer[0] = 'c';
	param->buffer[1] = 'm';
	param->buffer[2] = 'd';
	param->buffer[3] = ':';
	param->buffer[4] = ' ';
	c = 5;
	for(i =0 ; i < argc;i++){
		for(j = 0; j < strlen(argv[i]);j++){
			if(c == kslibMSGBUFSIZE-2){
				break;
			}
			param->buffer[c] = argv[i][j];
			c++;
		}
		if(c == kslibMSGBUFSIZE-2){
			break;
		}
		param->buffer[c] = ' ';
		c++;
		
	}
	param->buffer[c] = '\n';
	param->buffer[c+1] = 0;
	
	param->messages = append_message(param->messages, param->buffer );
	
	
	//if(param->matchstart)
	//fprintf(stderr,"Viterbi: %d\n",param->viterbi);
	
	
	c = 0;
	for(i = 0; i < param->read_structure->num_segments;i++){
		if(param->read_structure->type[i] == 'R'){
			c++;
		}
	}
	
	
	if(c >= 2){
		param->multiread = c;
	}
	/*if(!c  && param->sim == 0){
		fprintf(stderr,"ERROR: no read segment specified!%d\n",param->read_structure->num_segments);
		
		free_param(param);
		exit(EXIT_FAILURE);
	}*/
	
	
	
	
	
	if(param->reference_fasta || param->dust){
		if(param->multiread){
			sprintf(param->buffer,"WARNING: cannot dust or filter sequences by comparison to a known sequence if multiple reads are present in one input seqeunce.\n");
			param->messages = append_message(param->messages, param->buffer);
			param->dust = 0;
			param->reference_fasta = 0;
			
		}
	//param->messages = append_message(param->messages, param->buffer );
	}
	
	MMALLOC(param->infile,sizeof(char*)* (argc-optind));
	
	c = 0;
	while (optind < argc){
		param->infile[c] =  argv[optind++];
		c++;
	}
	param->infiles = c;
	
	


	
	//fprintf(stderr,"%d %p\n",param->read_structure->num_segments, param->arch_file  );
	
	
	
	
	
	//free(command_line);
	return param;
ERROR:
	
	if(param){
		fprintf(stdout,"%s",param->errmsg);
		free_param(param);
	}
	return NULL;
}

/** \fn struct parameters* assign_segment_sequences(struct parameters* param, char* tmp, int segment)
 \brief Assigns sequences to segments.
 
 \param param @ref parameters 
 \param tmp sequences from user input.
 \param segment number of segment.
 */
int assign_segment_sequences(struct parameters* param, char* tmp, int segment)
{
	int i,f,g;
	int count;
	int len;
	int status;
	struct read_structure* read_structure = param->read_structure;
	
	
	
	
	//tmp = optarg;
	
	
	switch (tmp[0]) {
		case 'R':
		case 'G':
		case 'O':
		case 'P':
		case 'S':
		case 'F':
		case 'B':
			break;
			
		default:
			KSLIB_FAIL(kslFAIL,param->errmsg,"Segment type :%c not recognized.\n",tmp[0]);
			break;
	}

	
	count = byg_count(",", tmp);
	
	if(tmp[0] == 'B'){ // add extra space for all N barcode....
		count++;
	}
	
	if(tmp[0] == 'S'){ // add extra space for all N barcode....
		count++;
	}
	f = 0;
	
	//fprintf(stderr,"Segment %d: %d	X%sX sequences\n",segment,count+1,tmp);
	read_structure->numseq_in_segment[segment] = count+1;
	
	MMALLOC(read_structure->sequence_matrix[segment],sizeof(char*) * (count+1));
	for(i = 0; i < count+1;i++){
		read_structure->sequence_matrix[segment][i] = 0;
		MMALLOC(read_structure->sequence_matrix[segment][i],sizeof(char)* strlen(tmp));
	}
	read_structure->type[segment] = tmp[0];
	
	if(tmp[0] == 'R'){
		read_structure->sequence_matrix[segment][0][0]  = 'N';
		read_structure->sequence_matrix[segment][0][1]  = 0;
		
	}else{
	
		f = 0;
		g = 0;
		len = (int)strlen(tmp) ;
		//fprintf(stderr,"%d\n",len);
		for(i = 2;i <  len;i++){
			if(tmp[i] != ','){
				read_structure->sequence_matrix[segment][f][g] = tmp[i];
				g++;
			}else{
				read_structure->sequence_matrix[segment][f][g] = 0;
				f++;
				g = 0;
			}
		}
		read_structure->sequence_matrix[segment][f][g] = 0;
	}
	
	if(tmp[0] == 'B'){
		f = f + 1;
		g = 0;
		for(i = 0; i < strlen(read_structure->sequence_matrix[segment][0]);i++){
			read_structure->sequence_matrix[segment][f][g] = 'N';
			g++;
		}
		read_structure->sequence_matrix[segment][f][g] = 0;
	}
	
	if(tmp[0] == 'S'){
		f = f + 1;
		g = 0;
		for(i = 0; i < strlen(read_structure->sequence_matrix[segment][0]);i++){
			read_structure->sequence_matrix[segment][f][g] = 'N';
			g++;
		}
		read_structure->sequence_matrix[segment][f][g] = 0;
	}
	
	
	if(segment+1 >read_structure->num_segments  ){
		read_structure->num_segments = segment+1;
	}
	return kslOK;
ERROR:
	if(read_structure->sequence_matrix[segment]){
	
		for(i = 0; i < count+1;i++){
			MFREE(read_structure->sequence_matrix[segment][i]);//,sizeof(char)* strlen(tmp));
		}
		MFREE(read_structure->sequence_matrix[segment]);//,sizeof(char*) * (count+1));
	}
	
	return status;
}


/** \fn void usage()
 \brief Prints usage.
 */

#ifdef TAGDUST
void usage()
{
	fprintf(stdout, "\n%s %s, Copyright (C) 2015 Timo Lassmann <%s>\n",PACKAGE_NAME, PACKAGE_VERSION,PACKAGE_BUGREPORT);
	fprintf(stdout, "\n");
	fprintf(stdout, "Usage:   tagdust [options] <file>  -o <output prefix> \n\n");
	fprintf(stdout, "Options:\n");
	
	
	fprintf (stdout,"\t%-17s%10s%7s%-30s\n","-Q","FLT","", "confidence threshold [20].");
	fprintf (stdout,"\t%-17s%10s%7s%-30s\n","-l","STR","", "log file directory name.");
	fprintf (stdout,"\t%-17s%10s%7s%-30s\n","-start","INT","", "start of search area [0].");
	fprintf (stdout,"\t%-17s%10s%7s%-30s\n","-end","INT","", "end of search area [length of sequence].");
	fprintf (stdout,"\t%-17s%10s%7s%-30s\n","-format","STR","", "format of input sequence file.");
	fprintf (stdout,"\t%-17s%10s%7s%-30s\n"," -minlen","INT","", "minimal accepted read length [16].");
	fprintf (stdout,"\t%-17s%10s%7s%-30s\n","-ref","STR","", "reference fasta file to be compared against[].");
	fprintf (stdout,"\t%-17s%10s%7s%-30s\n","-fe","INT","", "number of errors allowed when comparing to reference[2].");
	fprintf (stdout,"\t%-17s%10s%7s%-30s\n","-dust","INT","", "remove low complexity sequences. [100].");
	fprintf (stdout,"\t%-17s%10s%7s%-30s\n","-e","FLT","", "expected sequencer error rate [0.05].");
	fprintf (stdout,"\t%-17s%10s%7s%-30s\n","-o","STR","", "output file name.");
	fprintf (stdout,"\t%-17s%10s%7s%-30s\n","-a","STR","", "output file for artifacts [NA].");
	fprintf (stdout,"\t%-17s%10s%7s%-30s\n","-t","INT","", "number of threads [8].");
	fprintf (stdout,"\t%-17s%10s%7s%-30s\n","-show_finger_seq","NA","", "print fingerprint as sequence (default is as base 4 number).");
	
	fprintf (stdout,"\t%-17s%10s%7s%-30s\n","-h/help","NA","", "print help.");
	fprintf (stdout,"\t%-17s%10s%7s%-30s\n","-v/version","NA","", "print version number.");
	fprintf (stdout,"\t%-17s%10s%7s%-30s\n","-1","STR","", "type of the first HMM building block.");
	fprintf (stdout,"\t%-17s%10s%7s%-30s\n","-2","STR","", "type of the second HMM building block.");
	fprintf (stdout,"\t%-17s%10s%7s%-30s\n","-...","STR","", "type of the . . . HMM building block.");
	fprintf(stdout, "\n");

}
#endif

#ifdef SIMREADS
void usage()
{
	fprintf(stdout, "\n%s %s, Copyright (C) 2015 Timo Lassmann <%s>\n",PACKAGE_NAME, PACKAGE_VERSION,PACKAGE_BUGREPORT);
	fprintf(stdout, "\n");
	fprintf(stdout, "Usage:   simreads  [options] <barcodefile from EDITTAG>-o <file>  .... \n\n");
	fprintf(stdout, "Options:\n");
	fprintf (stdout,"\t%-17s%10s%7s%-30s\n","-sim_barlen","INT","", "Barcode length.");
	fprintf (stdout,"\t%-17s%10s%7s%-30s\n","-sim_barnum","INT" ,"","Number of samples.");
	fprintf (stdout,"\t%-17s%10s%7s%-30s\n","-sim_5seq","STR" ,"", "Sequence of 5' linker.");
	fprintf (stdout,"\t%-17s%10s%7s%-30s\n","-sim_3seq","STR" ,"", "Sequence of 3' linker.");
	fprintf (stdout,"\t%-17s%10s%7s%-30s\n","-sim_readlen","INT" ,"", "Length of read.");
	fprintf (stdout,"\t%-17s%10s%7s%-30s\n","-sim_readlen_mod","INT" ,"", "+/- mod of read length.");
	fprintf (stdout,"\t%-17s%10s%7s%-30s\n","-sim_error_rate","FLT" ,"", "Simulated error rate.");
	fprintf (stdout,"\t%-17s%10s%7s%-30s\n","-sim_InDel_frac","FLT" ,"", "INDEL fraction.");
	fprintf (stdout,"\t%-17s%10s%7s%-30s\n","-sim_numseq","INT" ,"", "Number of simulated sequences.");
	fprintf (stdout,"\t%-17s%10s%7s%-30s\n","-sim_random_frac","FLT" ,"", "Fraction of totally random sequences.");
	fprintf (stdout,"\t%-17s%10s%7s%-30s\n","-sim_endloss","INT" ,"", "mean number of nucleotides lost on either end of the read.");

	
	fprintf(stdout, "\n");
	
}

#endif

#ifdef MERGE
void usage()
{
	fprintf(stdout, "\n%s %s, Copyright (C) 2015 Timo Lassmann <%s>\n",PACKAGE_NAME, PACKAGE_VERSION,PACKAGE_BUGREPORT);
	fprintf(stdout, "\n");
	fprintf(stdout, "Usage:   merge <file>  .... \n\n");
	fprintf(stdout, "Options:\n");
	
	
	fprintf(stdout, "\n");
	
}
#endif

#ifdef RENAME
void usage()
{
	fprintf(stdout, "\n%s %s, Copyright (C) 2015 Timo Lassmann <%s>\n",PACKAGE_NAME, PACKAGE_VERSION,PACKAGE_BUGREPORT);
	fprintf(stdout, "\n");
	fprintf(stdout, "Usage:   rename_qiime <map file> <file>  .... \n\n");
	fprintf(stdout, "Options:\n");
	fprintf(stdout, "\n");
}
#endif

#ifdef EVALRES
void usage()
{
	fprintf(stdout, "\n%s %s, Copyright (C) 2015 Timo Lassmann <%s>\n",PACKAGE_NAME, PACKAGE_VERSION,PACKAGE_BUGREPORT);
	fprintf(stdout, "\n");
	fprintf(stdout, "Usage:   evalres <file>  .... \n\n");
	fprintf(stdout, "Options:\n");
	
	
	fprintf(stdout, "\n");
	
}
#endif



/** \fn void free_param(struct parameters* param)
 \brief free @ref parameters.
 */
int free_param(struct parameters* param)
{
	char logfile[200];
	FILE* outfile = 0;
	int i,status;
	//if(param->log){
	if(param->outfile){
		sprintf (logfile, "%s_logfile.txt",param->outfile);
		if((outfile = fopen(logfile, "w")) == NULL) KSLIB_XEXCEPTION_SYS(kslEWRT,"Failed to open file:%s",logfile);
		//if ((outfile = fopen( logfile, "w")) == NULL){
		//	fprintf(stderr,"can't open logfile\n");
		//	exit(-1);
		//}
		fprintf(outfile,"%s\n",param->messages);
	
		fclose(outfile);
	
	}
	if(param->read_structures){
		for(i = 0; i < param->infiles;i++){
			if(param->read_structures[i]){
				free_read_structure(param->read_structures[i]);
			}
		}
		MFREE(param->read_structures);
	}
	
	if(param->read_structure){
		free_read_structure(param->read_structure);
	}
	if(param->read_structure_R1){
		free_read_structure(param->read_structure_R1);
	}
	if(param->read_structure_R2){
		free_read_structure(param->read_structure_R2);
	}
	if(param->confidence_thresholds){
		MFREE(param->confidence_thresholds);
	}
	
	MFREE (param->infile);
	MFREE(param->messages);
	MFREE(param->buffer);
	MFREE(param);
	return kslOK;
ERROR:
	return status;
}


int QC_read_structure(struct parameters* param )
{
	struct read_structure* read_structure = param->read_structure;
	int i,g,f,min_error,errors;
	int last = -1;
	int num_pairs = 0;
	for(i = 0; i < read_structure->num_segments;i++){
		
		
		if(read_structure->sequence_matrix[i]){
			if(last +1 != i){
				sprintf(param->buffer,"ERROR: a hmm building lock was skipped??\n");
				param->messages = append_message(param->messages, param->buffer);

				return 1;
			}
			
			//serious checking...
			for(g = 0;g < read_structure->numseq_in_segment[i];g++){
				for(f = g+1;f < read_structure->numseq_in_segment[i];f++){
					if(strlen(read_structure->sequence_matrix[i][g]) != strlen(read_structure->sequence_matrix[i][f])){
						
						sprintf(param->buffer,"ERROR: the sequences in the same segment have to have the same length.\n");
						param->messages = append_message(param->messages, param->buffer);
						
						sprintf(param->buffer,"Segment %d\n%s	%d\n%s	%d\n",read_structure->numseq_in_segment[i], read_structure->sequence_matrix[i][g],g, read_structure->sequence_matrix[i][f],f );
						param->messages = append_message(param->messages, param->buffer);
						
						
						//fprintf(stderr,"ERROR: the sequences in the same segment have to have the same length\n");
						//fprintf(stderr,"%dseq\n%s	%d\n%s	%d\n",read_structure->numseq_in_segment[i], read_structure->sequence_matrix[i][g],g, read_structure->sequence_matrix[i][f],f );
						return 1;
					}
				}
			}
			
			/*fprintf(stderr,"Found %c segment %d with %d sequences\n", read_structure->type[i] ,i, read_structure->numseq_in_segment[i] );
			for(g = 0;g < read_structure->numseq_in_segment[i];g++){
				fprintf(stderr,"%s, ",read_structure->sequence_matrix[i][g] );
			}
			fprintf(stderr,"\n" );*/
			last = i;
			
		}
		
		
		
		if(read_structure->type[i] == 'B'){
			min_error = 1000;
			for(g = 0;g <  read_structure->numseq_in_segment[i];g++){
				for(f = g+1;f < read_structure->numseq_in_segment[i];f++){
					errors = bpm(read_structure->sequence_matrix[i][g], read_structure->sequence_matrix[i][f], (int)strlen(read_structure->sequence_matrix[i][0]),(int)strlen(read_structure->sequence_matrix[i][0]));
					
					if(errors < min_error){
						min_error = errors;
						num_pairs = 1;
						//	fprintf(stderr,"%s\n%s\t%d\n",param->read_structure->sequence_matrix[i][j],param->read_structure->sequence_matrix[i][c] ,numseq);
					}else if(errors == min_error ){
						num_pairs++;
					}
				}
			}
			//if(min_error != 1000){
				//sprintf(param->buffer,"Minumum edit distance among barcodes: %d, %d pairs\n", min_error,num_pairs);
				
				//param->messages = append_message(param->messages, param->buffer);

				
			//}
			
		}
	}
	return kslOK;
}

struct read_structure* malloc_read_structure(void)
{
	struct read_structure* read_structure = 0;
	int i;
	int status;
	MMALLOC(read_structure, sizeof(struct read_structure));
	read_structure->sequence_matrix = 0;
	read_structure->numseq_in_segment = 0;
	read_structure->type = 0;
	MMALLOC(read_structure->sequence_matrix ,sizeof(char**) * 10 );
	MMALLOC(read_structure->numseq_in_segment, sizeof(int) * 10);
	MMALLOC(read_structure->type ,sizeof(char) * 10 );
	
	
	for(i = 0;i  <10;i++){
		read_structure->sequence_matrix[i] = NULL;
		read_structure->numseq_in_segment[i] = 0;
		read_structure->type[i] = 0;
		
	}
	
	read_structure->num_segments = 0;
	return read_structure;
ERROR:
	KSLIB_MESSAGE(status,"Something went wrong in malloc_read_structure.\n");
	return NULL;
}

void free_read_structure(struct read_structure* read_structure)
{
	int i,j;
	for(i = 0; i < 10;i++){
		if(read_structure->sequence_matrix[i]){
			for(j = 0; j < read_structure->numseq_in_segment[i];j++){
				MFREE(read_structure->sequence_matrix[i][j]);
			}
			
			MFREE(read_structure->sequence_matrix[i]);
		}
	}
	MFREE(read_structure->sequence_matrix);
	
	MFREE(read_structure->numseq_in_segment );
	MFREE(read_structure->type);
	MFREE(read_structure);
}




