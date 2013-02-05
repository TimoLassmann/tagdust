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



struct parameters* interface(struct parameters* param,int argc, char *argv[])
{
	int c;
	int help = 0;
	
	if (argc < 2 && isatty(0)){
		usage();
		exit(0);
	}
		
	param = malloc(sizeof(struct parameters));
	param->infiles = 0;
	param->infile = 0;
	param->outfile = 0;

	param->quiet_flag = 0;
	param->num_query = 100000;
	param->kmer_size = 2;
	param->print_unmapped = 0;
	param->solexa = 0;
	param->print_qual = 0;
	param->print_posteriors = 0;
	param->summary = 0;
	param->k_errors_allowed = -1;
	param->format = 0;
	param->filter = 0;
	param->alt_lib_name = 0;
	param->gzipped = 0;
	param->sam = 0;
		
	while (1){	 
		static struct option long_options[] ={
			{"unmapped",required_argument,0, OPT_UNMAPPED},
			{"quiet",0,0,'q'},
			{"help",0,0,'h'},
			{0, 0, 0, 0}
		};
		
		int option_index = 0;
		c = getopt_long_only (argc, argv,"qhk:s:f:n:F:",long_options, &option_index);
		
		if (c == -1){
			break;
		}
		
		switch(c) {
			case 0:
				
				break;
			case OPT_UNMAPPED:
				param->print_unmapped = optarg;
				break;
			case 'q':
				param->quiet_flag = 1;
				break;
			case 'f':
				param->format = optarg;
				break;
			case 'F':
				param->filter = optarg;
				break;
			case 's':
				param->summary = optarg;
				break;
			case 'k':
				param->k_errors_allowed = atoi(optarg);
				break;
			case 'h':
				help = 1;
				break;
			case 'n':
				param->alt_lib_name = optarg;
				break;
			case '?':
				exit(1);
				break;
			default:
				fprintf(stderr,"default\n\n\n\n");
				abort ();
		}
	}
	//fprintf(stderr,"Viterbi: %d\n",param->viterbi);
	
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


void usage()
{
	fprintf(stdout, "\nSAMStat %0.2f, Copyright (C) 2010, 2011 Timo Lassmann <timolassmann@gmail.com>\n", VERSION);
	fprintf(stdout, "\n");
	fprintf(stdout, "Usage:   samstat <file.sam> <file.bam> <file.fa> <file.fq> .... \n\n");
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
	free(param->infile);
	free(param);
}

