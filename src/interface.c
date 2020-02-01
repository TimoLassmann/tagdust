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

#include <unistd.h>
#include "interface.h"

#include "tlmisc.h"


#define OPT_SHOWW 1
#define OPT_SEED 2
#define OPT_RECIPE 3
#define OPT_BOOK 4
#define OPT_FILTER 5
#define OPT_ERROR 6

static int print_tagdust_warranty(void);

static int print_tagdust_header(void);

int print_tagdust_header(void)
{
        fprintf(stdout,"\n");
        fprintf(stdout,"Tagdust (%s)\n", PACKAGE_VERSION);
        fprintf(stdout,"\n");
        fprintf(stdout,"Copyright (C) 2009,2020 Timo Lassmann\n");
        fprintf(stdout,"\n");
        fprintf(stdout,"This program comes with ABSOLUTELY NO WARRANTY; for details type:\n");
        fprintf(stdout,"`tagdust -showw'.\n");
        fprintf(stdout,"This is free software, and you are welcome to redistribute it\n");
        fprintf(stdout,"under certain conditions; consult the COPYING file for details.\n");
        fprintf(stdout,"\n");
        fprintf(stdout,"Please cite:\n");


        fprintf(stdout,"  Lassmann, Timo.\n");
        fprintf(stdout,"  \"TagDust2: a generic method to extract reads from sequencing data.\"\n");
        fprintf(stdout,"  BMC bioinformatics (2015) \n");
        fprintf(stdout,"  https://doi.org/10.1186/s12859-015-0454-y\n");
        fprintf(stdout,"\n");

        /*fprintf(stdout,"  Lassmann, Timo, Oliver Frings, and Erik LL Sonnhammer.\n");
        fprintf(stdout,"  \"Kalign2: high-performance multiple alignment of protein and\n");
        fprintf(stdout,"  nucleotide sequences allowing external features.\"\n");
        fprintf(stdout,"  Nucleic acids research 37.3 (2008): 858-865.\n");
        fprintf(stdout,"\n");
        fprintf(stdout,"  Lassmann, Timo, and Erik LL Sonnhammer. \"Kalign–an accurate and\n");
        fprintf(stdout,"  fast multiple sequence alignment algorithm.\"\n  BMC bioinformatics 6.1 (2005): 298.\n");
        fprintf(stdout,"\n");*/

        return OK;
}

int print_tagdust_warranty(void)
{
        fprintf(stdout,"Here is the Disclaimer of Warranty section of the GNU General Public License (GPL):\n");
        fprintf(stdout,"\n");
        fprintf(stdout,"15. Disclaimer of Warranty.\n");
        fprintf(stdout,"THERE IS NO WARRANTY FOR THE PROGRAM, TO THE EXTENT PERMITTED BY\n");
        fprintf(stdout,"APPLICABLE LAW.  EXCEPT WHEN OTHERWISE STATED IN WRITING THE COPYRIGHT\n");
        fprintf(stdout,"HOLDERS AND/OR OTHER PARTIES PROVIDE THE PROGRAM \"AS IS\" WITHOUT WARRANTY\n");
        fprintf(stdout,"OF ANY KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO,\n");
        fprintf(stdout,"THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR\n");
        fprintf(stdout,"PURPOSE.  THE ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE PROGRAM\n");
        fprintf(stdout,"IS WITH YOU.  SHOULD THE PROGRAM PROVE DEFECTIVE, YOU ASSUME THE COST OF\n");
        fprintf(stdout,"ALL NECESSARY SERVICING, REPAIR OR CORRECTION.\n");
        fprintf(stdout,"\n");
        fprintf(stdout,"A complete copy of the GPL can be found in the \"COPYING\" file.\n");
        return OK;
}

#ifndef HAVE_AVX2
static int print_AVX_warning(void);

int print_AVX_warning(void)
{
        fprintf(stdout,"\n");
        fprintf(stdout,"WARNING: AVX2 instruction set not found!\n");
        fprintf(stdout,"         Tagdust will not run optimally.\n");
        fprintf(stdout,"\n");

        return OK;
}
#endif

static void usage(void);

int interface(struct parameters** p,int argc, char *argv[])
{
        struct parameters* param = NULL;
        int i,c;
        int help = 0;
        int version = 0;
        int showw = 0;

        print_tagdust_header();
        if (argc < 2 && isatty(0)){
                usage();
                return OK;
        }



#ifndef HAVE_AVX2
        RUN(print_AVX_warning());
#endif

        //int showw = 0;
        MMALLOC(param,sizeof(struct parameters));
        //param->segments = NULL;
        //param->num_segments = 0;
        //MMALLOC(param->segments, sizeof(char*) * 10);
        //for(i = 0;i < 10;i++){
        //param->segments[i] = NULL;
        //}
        param->num_infiles = 0;
        param->infile = NULL;
        param->outfile = NULL;
        param->filter_fasta = NULL;
        param->filter_error = -1;
        param->num_threads = 8;
        param->bam = 1;
        param->seed = 42;


        param->recipe = NULL;
        param->book_file = NULL;
        //param->read_structures = NULL;

        //RUN(malloc_read_structure(&param->read_structure));


        while (1){
                static struct option long_options[] ={
                        {"showw", 0,0,OPT_SHOWW },
                        {"out",required_argument,0, 'o'},
                        {"seed",required_argument,0, OPT_SEED},
                        {"book",required_argument,0, OPT_BOOK },
                        {"filter", required_argument,0,OPT_FILTER},
                        {"error", required_argument, 0, OPT_ERROR},
                        {"help",0,0,'h'},
                        {"version",0,0,'v'},
                        {"recipe",required_argument,0,'r'},
                        {0, 0, 0, 0}
                };

                int option_index = 0;
                c = getopt_long_only (argc, argv,"o:t:r:qhv",long_options, &option_index);

                if (c == -1){
                        break;
                }

                switch(c) {
                case 0:
                        break;
                case OPT_ERROR:
                        param->filter_error = atoi(optarg);
                        break;
                case OPT_FILTER:
                        param->filter_fasta = optarg;
                        break;
                case OPT_BOOK:
                        param->book_file = optarg;
                        break;
                case OPT_SEED:
                        param->seed = atoi(optarg);
                        break;
                case 'r':
                        param->recipe = optarg;
                        break;
                case 'o':
                        param->outfile = optarg;
                        break;
                case 't':
                        param->num_threads = atoi(optarg);
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
                        ERROR_MSG("Option not recognised");
                }
        }
        if(help){
                usage();
                free_param(param);
                return OK;
        }

        if(version){
                fprintf(stdout,"%s %s\n",PACKAGE_NAME,PACKAGE_VERSION);
                free_param(param);
                return OK;
        }

        if(showw){
                print_tagdust_warranty();
                free_param(param);
                return OK;

        }

        //MMALLOC(param->buffer,sizeof(char) * kslibMSGBUFSIZE);

        //snprintf(param->buffer , kslibMSGBUFSIZE,"%s %s, Copyright (C) 2013-2019 Timo Lassmann <%s>\n",PACKAGE_NAME, PACKAGE_VERSION,PACKAGE_BUGREPORT);
        //param->messages = append_message(param->messages, param->buffer  );
        //command_line[0] = 0;
        //c = 0;
        /*
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
        */
        //param->messages = append_message(param->messages, param->buffer );


        //if(param->matchstart)
        //fprintf(stderr,"Viterbi: %d\n",param->viterbi);


        /*c = 0;
        for(i = 0; i < param->read_structure->num_segments;i++){
                if(param->read_structure->type[i] == 'R'){
                        c++;
                }
        }


        if(c >= 2){
                param->multiread = c;
        }
        */
        if(argc-optind){
                MMALLOC(param->infile,sizeof(char*)* (argc-optind));

                c = 0;
                while (optind < argc){
                        param->infile[c] =  argv[optind++];
                        c++;
                }
                param->num_infiles = c;
        }

        if(param->num_infiles){
                for(i = 0; i < param->num_infiles;i++){
                        ASSERT(my_file_exists(param->infile[i]),"Could not find file %s.",param->infile[i]);
                }
        }
        *p = param;
        return OK;
ERROR:
        if(param){
                free_param(param);
        }
        return FAIL;
}

void usage()
{
        fprintf(stdout, "Usage:   tagdust [options] <file>  -o <output prefix> \n\n");
        fprintf(stdout, "Options:\n");


        fprintf (stdout,"\t%-17s%10s%7s%-30s\n","-Q","FLT","", "confidence threshold [20].");
        fprintf (stdout,"\t%-17s%10s%7s%-30s\n","-l","STR","", "log file directory name.");
        fprintf (stdout,"\t%-17s%10s%7s%-30s\n","-start","INT","", "start of search area [0].");
        fprintf (stdout,"\t%-17s%10s%7s%-30s\n","-end","INT","", "end of search area [length of sequence].");
        fprintf (stdout,"\t%-17s%10s%7s%-30s\n","-format","STR","", "format of input sequence file.");
        fprintf (stdout,"\t%-17s%10s%7s%-30s\n","-minlen","INT","", "minimal accepted read length [16].");
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

#ifdef SIMREADS
void usage()
{
        fprintf(stdout, "\n%s %s, Copyright (C) 2015-2019 Timo Lassmann <%s>\n",PACKAGE_NAME, PACKAGE_VERSION,PACKAGE_BUGREPORT);
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
        fprintf(stdout, "\n%s %s, Copyright (C) 2015-2019 Timo Lassmann <%s>\n",PACKAGE_NAME, PACKAGE_VERSION,PACKAGE_BUGREPORT);
        fprintf(stdout, "\n");
        fprintf(stdout, "Usage:   merge <file>  .... \n\n");
        fprintf(stdout, "Options:\n");


        fprintf(stdout, "\n");

}
#endif

#ifdef RENAME
void usage()
{
        fprintf(stdout, "\n%s %s, Copyright (C) 2015-2019 Timo Lassmann <%s>\n",PACKAGE_NAME, PACKAGE_VERSION,PACKAGE_BUGREPORT);
        fprintf(stdout, "\n");
        fprintf(stdout, "Usage:   rename_qiime <map file> <file>  .... \n\n");
        fprintf(stdout, "Options:\n");
        fprintf(stdout, "\n");
}
#endif

#ifdef EVALRES
void usage()
{
        fprintf(stdout, "\n%s %s, Copyright (C) 2015-2019 Timo Lassmann <%s>\n",PACKAGE_NAME, PACKAGE_VERSION,PACKAGE_BUGREPORT);
        fprintf(stdout, "\n");
        fprintf(stdout, "Usage:   evalres <file>  .... \n\n");
        fprintf(stdout, "Options:\n");


        fprintf(stdout, "\n");

}
#endif

int free_param(struct parameters* param)
{
        //char logfile[200];
        //FILE* outfile = 0;
        if(param){
                //int i;
                //if(param->log){
                //if(param->outfile){
                        //sprintf (logfile, "%s_logfile.txt",param->outfile);
                        //RUNP(outfile = fopen(logfile, "w"));
                        //if((outfile = fopen(logfile, "w")) == NULL) KSLIB_XEXCEPTION_SYS(kslEWRT,"Failed to open file:%s",logfile);
                        //if ((outfile = fopen( logfile, "w")) == NULL){
                        //	fprintf(stderr,"can't open logfile\n");
                        //	exit(-1);
                        //}
                        //fprintf(outfile,"%s\n",param->messages);

                        //fclose(outfile);
                //}
                /*if(param->read_structures){
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
                  }*/

                MFREE(param->infile);
                MFREE(param);
        }


        return OK;
//ERROR:
        //return FAIL;
}





#ifdef INTERFACE_TEST
/* valgrind ./interface_test  -1 P:agggaggacgatgcgg -2 B:TGCT,AAAA,AACC,AAGG,AATT,ACAC,ACCA,ACGT -3 R:N -4 P:gtgtcagtcacttccagcgg in.fq -o out.fq */
int main(int argc, char *argv[])
{
        struct parameters* param = NULL;
        RUN(interface(&param, argc, argv));
        free_param(param);
        return EXIT_SUCCESS;
ERROR:
        return EXIT_FAILURE;
}

#endif
