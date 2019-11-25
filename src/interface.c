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

#include "interface.h"

#include "misc.h"



int interface(struct parameters** p,int argc, char *argv[])
{
        struct parameters* param = NULL;
        int c;
        int help = 0;
        int version = 0;

        if (argc < 2 && isatty(0)){
                usage();
                return OK;
        }

        MMALLOC(param,sizeof(struct parameters));
        param->segments = NULL;
        MMALLOC(param->segments, sizeof(char*) * 10);
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
        param->confidence_threshold = 0.0;
        param->confidence_thresholds = NULL;
        param->confidence_threshold_R1 = 0.0;
        param->confidence_threshold_R2 = 0.0;

        //param->read_structure = 0;
        //param->read_structure_R1 = 0;
        //param->read_structure_R2 = 0;
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

        //param->read_structures = NULL;

        //RUN(malloc_read_structure(&param->read_structure));


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
                        param->segments[0] = optarg;
                        break;
                case OPT_SEG2:
                        param->segments[1] = optarg;
                        break;
                case OPT_SEG3:
                        param->segments[2] = optarg;
                        break;
                case OPT_SEG4:
                        param->segments[3] = optarg;
                        break;
                case OPT_SEG5:
                        param->segments[4] = optarg;
                        break;
                case OPT_SEG6:
                        param->segments[5] = optarg;
                        break;
                case OPT_SEG7:
                        param->segments[6] = optarg;
                        break;
                case OPT_SEG8:
                        param->segments[7] = optarg;
                        break;
                case OPT_SEG9:
                        param->segments[8] = optarg;
                        break;
                case OPT_SEG10:
                        param->segments[9] = optarg;
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
                return OK;
        }

        if(version){
                fprintf(stdout,"%s %s\n",PACKAGE_NAME,PACKAGE_VERSION);
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
        if(param->reference_fasta || param->dust){
                if(param->multiread){
                        WARNING_MSG("WARNING: cannot dust or filter sequences by comparison to a known sequence if multiple reads are present in one input seqeunce.\n");
                        param->dust = 0;
                        param->reference_fasta = 0;

                }
        }

        MMALLOC(param->infile,sizeof(char*)* (argc-optind));

        c = 0;
        while (optind < argc){
                param->infile[c] =  argv[optind++];
                c++;
        }
        param->infiles = c;


        *p = param;
        return OK;
ERROR:

        if(param){
                //fprintf(stdout,"%s",param->errmsg);
                free_param(param);
        }
        return FAIL;
}

void usage()
{
        fprintf(stdout, "\n%s %s, Copyright (C) 2015-2019 Timo Lassmann <%s>\n",PACKAGE_NAME, PACKAGE_VERSION,PACKAGE_BUGREPORT);
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
        char logfile[200];
        FILE* outfile = 0;
        if(param){
                //int i;
                //if(param->log){
                if(param->outfile){
                        sprintf (logfile, "%s_logfile.txt",param->outfile);
                        RUNP(outfile = fopen(logfile, "w"));
                        //if((outfile = fopen(logfile, "w")) == NULL) KSLIB_XEXCEPTION_SYS(kslEWRT,"Failed to open file:%s",logfile);
                        //if ((outfile = fopen( logfile, "w")) == NULL){
                        //	fprintf(stderr,"can't open logfile\n");
                        //	exit(-1);
                        //}
                        fprintf(outfile,"%s\n",param->messages);

                        fclose(outfile);

                }
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
                if(param->confidence_thresholds){
                        MFREE(param->confidence_thresholds);
                }
                MFREE(param->segments);

                MFREE(param->infile);
                if(param->messages){
                        MFREE(param->messages);
                }
                if(param->buffer){
                        MFREE(param->buffer);
                }
                MFREE(param);
        }
        return OK;
ERROR:
        return FAIL;
}



#ifdef INTERFACE_TEST

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
