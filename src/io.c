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

/*! \file io.c
  \brief functions for reading sequences.
n
n  Initializes nucleotide alphabet needed to parse input. Calls parameter parser. Calls functions to process the data. \author Timo Lassmann \bug No known bugs.
*/


#include <ctype.h>
#include <unistd.h>
#include <string.h>

#if HAVE_CONFIG_H
#include "config.h"
#endif

#include "tldevel.h"

//#include "interface.h"
#include "nuc_code.h"
//#include "misc.h"

#include "io.h"




#define MAX_LINE 10000


int write_all(struct assign_struct* as,char* prefix)
{
        char filename[256];
        char alphabet[] = "ACGTNN";
        int i,j,gg;


        /* should check before running if files exist */
        for(i = 0; i < as->total;i++){
                fprintf(stdout,"READ %d %s (PASS: %d)  ",i, as->bits[i]->name, as->bits[i]->pass);
                for(j = 0; j < as->num_files;j++){

                        fprintf(stdout,"%f ", as->bits[i]->Q[j]);
                }
                fprintf(stdout,"\n");

                fprintf(stdout,"%s\n%s\n", as->bits[i]->name,as->bits[i]->bc);
                for(j = 0; j < as->bits[i]->num_bit;j++){
                        switch(as->bits[i]->bits[j]->type){
                        case READ_TYPE:
                                fprintf(stdout,"READ (file: %d): ", as->bits[i]->bits[j]->file);
                                for(gg = 0; gg < as->bits[i]->bits[j]->len;gg++){
                                        fprintf(stdout,"%c", alphabet[(int)as->bits[i]->bits[j]->p[gg]]);
                                }
                                fprintf(stdout,"\n");
                                fprintf(stdout,"QUAL (file: %d): ", as->bits[i]->bits[j]->file);
                                for(gg = 0; gg < as->bits[i]->bits[j]->len;gg++){
                                        fprintf(stdout,"%c",as->bits[i]->bits[j]->q[gg]);
                                }
                                fprintf(stdout,"\n");
                                break;
                        case BAR_TYPE:
                                fprintf(stdout,"BAR (file: %d): ", as->bits[i]->bits[j]->file);
                                fprintf(stdout,"%s  and then...:", as->bits[i]->bits[j]->p);

                                fprintf(stdout,"\n");
                                break;
                        case UMI_TYPE:
                                fprintf(stdout,"UMI (file: %d): ", as->bits[i]->bits[j]->file);

                                for(gg = 0; gg < as->bits[i]->bits[j]->len;gg++){
                                        fprintf(stdout,"%c", alphabet[(int)as->bits[i]->bits[j]->p[gg]]);
                                }
                                fprintf(stdout,"\n");
                                break;
                        }
                }
        }
        return OK;
}


int io_handler(struct file_handler** fh, char*filename)
{
        struct file_handler* f_handle = NULL;
        FILE* file;

        char command[1000];
        char tmp[1000];

        int gzcat = -1;

        if(!*fh){
                MMALLOC(f_handle, sizeof(struct file_handler));
                *fh = f_handle;
        }

        f_handle = *fh;
        f_handle->bzipped = 0;
        f_handle->fasta = 0;
        f_handle->gzipped = 0;
        f_handle->sam = 0;
        if(access("/usr/bin/gzcat", X_OK) == 0){
                gzcat = 1;
        }else if(access("/bin/gzcat", X_OK) == 0){
                gzcat = 1;
        }else if(access("/usr/bin/zcat", X_OK) == 0){
                gzcat = 0;
        }else if(access("/bin/zcat", X_OK) == 0){
                gzcat = 0;
        }


        /*if(!file_exists(filename)){
                ERROR_MSG("Error: Cannot find input file: %s\n",filename );
                }*/

        if(!strcmp(".sam", filename + (strlen(filename ) - 4))){
                f_handle->sam = 1;
        }else if (!strcmp(".bam", filename + (strlen(filename ) - 4))){
                f_handle->sam = 2;
        }else if (!strcmp(".fa", filename + (strlen(filename ) - 3))){
                f_handle->sam = 0;
                f_handle->fasta = 1;
        }else if (!strcmp(".fq", filename + (strlen(filename ) - 3))){
                f_handle->sam = 0;
        }else if (!strcmp(".fastq", filename + (strlen(filename ) - 6))){
                f_handle->sam = 0;
        }else if (!strcmp(".fastaq", filename + (strlen(filename ) - 7))){
                f_handle->sam = 0;
        }else if (!strcmp(".fasta", filename + (strlen(filename ) - 6))){
                f_handle->sam = 0;
                f_handle->fasta = 1;
        }else if(!strcmp(".sam.gz", filename + (strlen(filename ) - 7))){
                f_handle->sam = 1;
                f_handle->gzipped  = 1;
        }else if (!strcmp(".bam.gz", filename + (strlen(filename ) - 7))){
                f_handle->sam = 2;
                f_handle->gzipped  = 1;
        }else if (!strcmp(".fa.gz", filename + (strlen(filename ) - 6))){
                f_handle->sam = 0;
                f_handle->fasta  = 1;
                f_handle->gzipped  = 1;
        }else if (!strcmp(".fq.gz", filename + (strlen(filename ) - 6))){
                f_handle->sam = 0;
                f_handle->gzipped  = 1;
        }else if (!strcmp(".fastq.gz", filename + (strlen(filename ) - 9))){
                f_handle->sam = 0;
                f_handle->gzipped  = 1;
        }else if (!strcmp(".fastaq.gz", filename + (strlen(filename ) - 10))){
                f_handle->sam = 0;
                f_handle->gzipped  = 1;
        }else if (!strcmp(".fasta.gz", filename + (strlen(filename ) - 9))){
                f_handle->sam = 0;
                f_handle->gzipped  = 1;
        }else if (!strcmp(".fastq.bz2", filename + (strlen(filename ) - 10))){
                f_handle->sam = 0;
                f_handle->bzipped  = 1;
        }else if (!strcmp(".fq.bz2", filename + (strlen(filename ) - 7))){
                f_handle->sam = 0;
                f_handle->bzipped  = 1;
        }else{
                f_handle->sam = -1;
        }


        if(f_handle->gzipped && gzcat == -1){
                ERROR_MSG("Cannot find gzcat / zcat on your system. Try gzcat <infile> | samstat -f sam/bam/fa/fq\n");}

        if(filename == NULL){
                if(f_handle->sam == 2){
                        command[0] = 0;
                        //if(!param->filter){
                        strcat ( command, "samtools view -F 768 ");
                        //}else{
                        //strcat ( command, "samtools view -F ");
                        //sprintf (tmp, "%s ",param->filter);
                        //strcat ( command, tmp);
                        //}
                        sprintf (tmp, "%s ","-");
                        strcat ( command, tmp);
                        if (!(file = popen(command, "r"))) {
                                fprintf(stderr,"Cannot open bam file '%s' with command:%s\n",filename,command);
                                exit(-1);
                        }
                }else if(f_handle->sam == 1){
                        command[0] = 0;
                        //if(!param->filter){
                                strcat ( command, "samtools view -SF 768 ");
                                //}else{
                                //strcat ( command, "samtools view -SF ");
                                //sprintf (tmp, "%s ",param->filter);
                                //strcat ( command, tmp);
                                //}
                        sprintf (tmp, "%s ", "-");
                        strcat ( command, tmp);
                        if (!(file = popen(command, "r"))) {
                                fprintf(stderr,"Cannot open bam file '%s' with command:%s\n",filename,command);
                                exit(-1);
                        }
                }else{
                        file = stdin;
                }
        }else{
                if(f_handle->sam == 2){
                        command[0] = 0;

                        if(f_handle->bzipped){
                                strcat ( command, "bzcat ");
                                //if(!param->filter){
                                        sprintf (tmp, "%s | samtools view -F 768 - ", filename);
                                        strcat ( command, tmp);
                                        //}else{
                                        //sprintf (tmp, "%s | samtools view -F  ", filename);
                                        //strcat ( command, tmp);
                                        //sprintf (tmp, "%s - ",param->filter);
                                        //strcat ( command, tmp);
                                        //}

                        }else if(f_handle->gzipped){
                                if(gzcat == 1){
                                        strcat ( command, "gzcat ");
                                }else{
                                        strcat ( command, "zcat ");
                                }
                                //if(!param->filter){
                                        sprintf (tmp, "%s | samtools view -F 768 - ", filename);
                                        strcat ( command, tmp);
                                        //}else{
                                        //sprintf (tmp, "%s | samtools view -F  ", filename);
                                        //strcat ( command, tmp);
                                        //sprintf (tmp, "%s - ",param->filter);
                                        //strcat ( command, tmp);
                                        //}
                        }else{
                                //if(!param->filter){
                                        strcat ( command, "samtools view -F 768 ");
                                        //}else{
                                        //strcat ( command, "samtools view -F ");
                                        //sprintf (tmp, "%s ",param->filter);
                                        //strcat ( command, tmp);
                                        //}
                                sprintf (tmp, "%s ", filename);
                                strcat ( command, tmp);
                        }

                        RUNP(file = popen(command, "r"));
                }else if(f_handle->sam == 1){
                        command[0] = 0;
                        if(f_handle->gzipped){
                                if(gzcat == 1){
                                        strcat ( command, "gzcat ");
                                }else{
                                        strcat ( command, "zcat ");
                                }
                                //if(!param->filter){
                                        sprintf (tmp, "%s | samtools view -SF 768 - ", filename);
                                        strcat ( command, tmp);
                                        //}else{
                                        //sprintf (tmp, "%s | samtools view -SF  ", filename);
                                        //strcat ( command, tmp);
                                        //sprintf (tmp, "%s - ",param->filter);
                                        //strcat ( command, tmp);
                                        //}
                        }else{
                                //if(!param->filter){
                                        strcat ( command, "samtools view -SF 768 ");
                                        //}else{
                                        //strcat ( command, "samtools view -SF ");
                                        //sprintf (tmp, "%s ",param->filter);
                                        //strcat ( command, tmp);
                                        //}
                                sprintf (tmp, "%s ", filename);
                                strcat ( command, tmp);
                        }
                        RUNP(file = popen(command, "r"));

                }else{
                        command[0] = 0;
                        if(f_handle->bzipped){
                                strcat ( command, "bzcat ");

                        }else if(f_handle->gzipped){
                                if(gzcat == 1){
                                        strcat ( command, "gzcat ");
                                }else{
                                        strcat ( command, "zcat ");
                                }
                        }else{
                                strcat ( command, "cat ");
                        }
                        sprintf (tmp, "%s ", filename);
                        strcat ( command, tmp);
                        RUNP(file = popen(command, "r"));
                }
        }
        //return file;

        f_handle->f_ptr = file;

        *fh = f_handle;
        return OK;
ERROR:
        return FAIL;
}

void print_sequence(struct read_info* ri,FILE* out)
{
        int i;
        char alpha[6] = "ACGTNN";
        fprintf(out,"@%s\n",ri->name);
        for(i = 0; i < ri->len;i++){
                fprintf(out,"%c", alpha[(int) ri->seq[i]]);
        }
        fprintf(out,"\n+\n%s\n" ,ri->qual);
}


FILE* open_file(struct parameters* param, char* buffer, char* mode)
{
        FILE* file = NULL;
        //int status;
        //sprintf (buffer, "%s_READ%d.fq",param->outfile,i+1);
        //if ((file = fopen(buffer, mode )) == NULL){
        RUNP(file = fopen(buffer, mode));
        //if((file = fopen(buffer, mode)) == NULL) KSLIB_XEXCEPTION_SYS(kslEWRT,"Failed to open file:%s",buffer);
        //	sprintf(param->buffer,"ERROR: cannot open file %s for writing.\n", buffer);
        //	param->messages = append_message(param->messages, param->buffer);
        //	free_param(param);
        //	exit(EXIT_FAILURE);
        //}
#ifdef DEBUG
        fprintf(stderr,"Opening:%s in %s mode\n",buffer,mode );
#endif

        return file;
ERROR:
        //KSLIB_MESSAGE(status,"Something wrong in open file.\n");
        return NULL;
}


/*
  void print_split_files(struct parameters* param, struct read_info** ri, int numseq)
  {
	struct stat buf;
	FILE* out_read1 = NULL;
	FILE* out_read2 = NULL;
	char alphabet[] = "ACGTNN";
	int i,j;
	int old_type = -100;
	int old_bar = -100;
	//int read_file = 0;
	void* tmp = 0;
	int start = 0;
	int stop = -1;
	int last = 0;
	int segment = 1;
	int print_segment_name = 0;

	static int check_for_files = 1;

	char* buffer =  0;
	MMALLOC(buffer,sizeof(char)* 1000 );
	qsort(ri,numseq, sizeof(struct read_info*), qsort_ri_barcode_compare);

	//1 ) has barcode or nor
	//	has: barcode >= 0 && extract_successs
	//	has_nor barcore == -1 && extract_success
	// 2) fail : ! extract_success..
	// ri[i]->read_type,ri[i]->barcode
	int h = 0;
	for(i = 0; i < numseq;i++){
  //fprintf(stderr,"BARCODE: %d	%d	%d	%p\n", ri[i]->barcode,(ri[i]->barcode >> 16) &0XFF,ri[i]->barcode &0XFF,param->read_structure);
  if(ri[i]->read_type != old_type || ri[i]->barcode != old_bar  ){
  if(ri[i]->read_type  == EXTRACT_SUCCESS){
  buffer[0] = 0;
  if(ri[i]->barcode != -1){
  if(param->multiread == 2){
  sprintf (buffer, "%s_BC_%s_READ1.fq",param->outfile,param->read_structure->sequence_matrix[(ri[i]->barcode >> 16) &0XFF][ri[i]->barcode &0XFF]);
  }else{
  sprintf (buffer, "%s_BC_%s.fq",param->outfile,param->read_structure->sequence_matrix[(ri[i]->barcode >> 16) &0XFF][ri[i]->barcode &0XFF]);
  }
  }else{
  if(param->multiread == 2){
  sprintf (buffer, "%s_READ1.fq",param->outfile);
  }else{
  sprintf (buffer, "%s.fq",param->outfile);
  }
  }
  }else{
  buffer[0] = 0;
  if(param->multiread == 2){
  sprintf (buffer, "%s_un_READ1.fq",param->outfile);
  }else{
  sprintf (buffer, "%s_un.fq",param->outfile);
  }

  }

  if(check_for_files){
  check_for_files = 0;
  if(!stat ( buffer, &buf )){
  //file found.
  sprintf(param->buffer,"ERROR: output file: %s already exists.\n", buffer);
  param->messages = append_message(param->messages, param->buffer);
  free_param(param);
  exit(EXIT_FAILURE);
  }
  if ((out_read1 = fopen(buffer, "w")) == NULL){
  fprintf(stderr,"can't open output\n");
  exit(-1);
  }
  buffer[0] = 0;
  if(param->multiread ==2){
  if(ri[i]->read_type  == EXTRACT_SUCCESS){
  buffer[0] = 0;
  if(ri[i]->barcode != -1){
  sprintf (buffer, "%s_BC_%s_READ2.fq",param->outfile,param->read_structure->sequence_matrix[(ri[i]->barcode >> 16) &0XFF][ri[i]->barcode &0XFF]);
  }else{
  sprintf (buffer, "%s_READ2.fq",param->outfile);

  }
  }else{
  buffer[0] = 0;
  sprintf (buffer, "%s_un_READ2.fq",param->outfile);

  }
  if ((out_read2 = fopen(buffer, "w")) == NULL){
  fprintf(stderr,"can't open output\n");
  exit(-1);
  }
  }
  }else{
  if(i){
  if(param->multiread ==2){
  fclose(out_read2);
  }
  fclose(out_read1);
  }

  if ((out_read1 = fopen(buffer, "a")) == NULL){
  fprintf(stderr,"can't open output\n");
  exit(-1);
  }
  if(param->multiread ==2){
  if(ri[i]->read_type  == EXTRACT_SUCCESS){
  buffer[0] = 0;
  if(ri[i]->barcode != -1){
  sprintf (buffer, "%s_BC_%s_READ2.fq",param->outfile,param->read_structure->sequence_matrix[(ri[i]->barcode >> 16) &0XFF][ri[i]->barcode &0XFF]);
  }else{
  sprintf (buffer, "%s_READ2.fq",param->outfile);

  }
  }else{
  buffer[0] = 0;
  sprintf (buffer, "%s_un_READ2.fq",param->outfile);

  }
  if ((out_read2 = fopen(buffer, "a")) == NULL){
  fprintf(stderr,"can't open output\n");
  exit(-1);
  }
  }
  }
  //	file open / close
  old_type = ri[i]->read_type;
  old_bar = ri[i]->barcode;
  }

  if(ri[i]->read_type == EXTRACT_SUCCESS){
  buffer[0] = 0;
  j = 0;
  if(ri[i]->fingerprint != -1){
  j |= 2;
  }
  if(ri[i]->barcode != -1){
  j |= 1;
  }
  switch (j) {
  case 0:
  sprintf (buffer, "%s",ri[i]->name);
  break;
  case 1:
  sprintf (buffer, "%s;BC:%s",ri[i]->name, param->read_structure->sequence_matrix[(ri[i]->barcode >> 16) &0XFF][ri[i]->barcode &0XFF]);
  break;
  case 2:
  sprintf (buffer, "%s;FP:%d",ri[i]->name,ri[i]->fingerprint);
  break;
  case 3:
  sprintf (buffer, "%s;FP:%d;BC:%s",ri[i]->name,ri[i]->fingerprint ,param->read_structure->sequence_matrix[(ri[i]->barcode >> 16) &0XFF][ri[i]->barcode &0XFF]);
  break;
  default:
  break;
  }

  //strcat (buffer, tmp);
  //ri[i]->name = realloc(ri[i]->name, sizeof(char) * (strlen(buffer) + 1) );
  MREALLOC(ri[i]->name,tmp,sizeof(char) * (strlen(buffer) + 1)) ;
  assert(ri[i]->name  != NULL);
  strcpy(ri[i]->name, buffer);



  start = 0;
  stop = -1;
  last = 0;
  segment = 1;
  print_segment_name = 0;
  for(j =0;j < ri[i]->len;j++){
  if( ri[i]->seq[j] == 65){
  print_segment_name = 1;
  break;
  }
  }
  if(print_segment_name){
  //if we have a multiread scroll to the first read if starting with 65....
  while(ri[i]->seq[start] == 65 && start < ri[i]->len){
  start++;

  }
  }

  while(stop != ri[i]->len){
  last = 0;
  for(j =start;j < ri[i]->len;j++){
  if( ri[i]->seq[j] == 65){
  stop = j;
  last = 1;
  break;
  }
  }
  if(!last){
  stop = ri[i]->len;
  }
  //print to file segment 'X'

  //if(print_segment_name){
  //	fprintf(out,"@%sRS:%d\n",ri->name,segment);
  //}else{

  if(segment ==1){
  h++;
  fprintf(out_read1,"@%s;RQ:%0.2f\n",ri[i]->name,ri[i]->mapq);

  for(j = start; j < stop;j++){
  fprintf(out_read1,"%c",alphabet[(int) ri[i]->seq[j]]);
  }
  fprintf(out_read1,"\n+\n");
  if(ri[i]->qual){
  for(j =start;j < stop;j++){
  fprintf(out_read1,"%c",ri[i]->qual[j]);
  }
  }else{
  for(j =start;j < stop;j++){
  fprintf(out_read1,".");
  }
  }
  fprintf(out_read1,"\n");

  }else if ( segment == 2){
  fprintf(out_read2,"@%s;RQ:%0.2f\n",ri[i]->name,ri[i]->mapq);


  for(j = start; j < stop;j++){
  fprintf(out_read2,"%c",alphabet[(int) ri[i]->seq[j]]);
  }
  fprintf(out_read2,"\n+\n");
  if(ri[i]->qual){
  for(j =start;j < stop;j++){
  fprintf(out_read2,"%c",ri[i]->qual[j]);
  }
  }else{
  for(j =start;j < stop;j++){
  fprintf(out_read2,".");
  }
  }
  fprintf(out_read2,"\n");
  }
  segment++;

  start = stop;
  while(ri[i]->seq[start] == 65 && start < ri[i]->len){
  start++;
  }

  if(segment > 100){
  exit(EXIT_FAILURE);
  }
  }


  }else{
  start = 0;
  stop = -1;
  last = 0;
  segment = 1;
  print_segment_name = 0;
  for(j =0;j < ri[i]->len;j++){
  if( ri[i]->seq[j] == 65){
  print_segment_name = 1;

  break;
  }
  }

  #ifdef DEBUG
  fprintf(stderr,"READ unextracted!:%d\n" , i);
  for(j = 0; j < ri[i]->len;j++){
  fprintf(stderr,"%d,",ri[i]->seq[j]);
  }
  fprintf(stderr,"\n");
  #endif


  if(print_segment_name){
  //if we have a multiread scroll to the first read if starting with 65....
  while(ri[i]->seq[start] == 65 && start < ri[i]->len){
  start++;

  }
  }
  while(stop != ri[i]->len){
  last = 0;
  for(j =start;j < ri[i]->len;j++){
  if( ri[i]->seq[j] == 65){
  stop = j;
  last = 1;
  //	segment++;
  break;
  }
  }

  if(!last){
  stop = ri[i]->len;
  }

  if(segment ==1){
  h++;
  fprintf(out_read1,"@%s;RQ:%0.2f\n",ri[i]->name,ri[i]->mapq);

  for(j = start; j < stop;j++){
  fprintf(out_read1,"%c",alphabet[(int) ri[i]->seq[j]]);
  }
  fprintf(out_read1,"\n+\n");
  if(ri[i]->qual){
  for(j =start;j < stop;j++){
  fprintf(out_read1,"%c",ri[i]->qual[j]);
  }
  }else{
  for(j =start;j < stop;j++){
  fprintf(out_read1,".");
  }
  }
  fprintf(out_read1,"\n");
  }else if ( segment == 2){
  fprintf(out_read2,"@%s;RQ:%0.2f\n",ri[i]->name,ri[i]->mapq);

  for(j = start; j < stop;j++){
  fprintf(out_read2,"%c",alphabet[(int) ri[i]->seq[j]]);
  }
  fprintf(out_read2,"\n+\n");
  if(ri[i]->qual){
  for(j =start;j < stop;j++){
  fprintf(out_read2,"%c",ri[i]->qual[j]);
  }
  }else{
  for(j =start;j < stop;j++){
  fprintf(out_read2,".");
  }
  }
  fprintf(out_read2,"\n");
  }
  segment++;

  start = stop;
  while(ri[i]->seq[start] == 65 && start < ri[i]->len){
  start++;
  }

  if(segment > 100){
  exit(EXIT_FAILURE);
  }
  }
  }
	}
	MFREE(buffer);

	if(param->multiread ==2){
  fclose(out_read2);
	}
	fclose(out_read1);
  }

*/

/** \fn void print_seq(struct read_info* ri,FILE* out)
    \brief Prints sequence and quality string to file.
    If the quality string is empty it is replaced by '.' characters.
    \param ri @ref read_info containing the sequence.
    \param out Pointer to output file.

*/
void print_seq(struct read_info* ri,FILE* out)
{
        char alphabet[] = "ACGTNN";
        int i;
        int start = 0;
        int stop = -1;
        int last = 0;
        int segment = 1;
        int print_segment_name = 0;

        for(i =0;i < ri->len;i++){
                if( ri->seq[i] == 65){
                        print_segment_name = 1;

                        break;
                }
        }
        if(print_segment_name){
                //if we have a multiread scroll to the first read if starting with 65....
                while(ri->seq[start] == 65 && start < ri->len){
                        start++;

                }
        }

        while(stop != ri->len){
                last = 0;
                for(i =start;i < ri->len;i++){
                        if( ri->seq[i] == 65){
                                stop = i;
                                last = 1;
                                //	segment++;
                                break;
                        }
                }
                if(!last){
                        stop = ri->len;
                }
                if(print_segment_name){
                        fprintf(out,"@%sRS:%d\n",ri->name,segment);
                }else{

                        fprintf(out,"@%s\n",ri->name);
                }
                for(i = start; i < stop;i++){
                        fprintf(out,"%c",alphabet[(int) ri->seq[i]]);
                }
                fprintf(out,"\n+\n");
                if(ri->qual){
                        for(i =start;i < stop;i++){
                                fprintf(out,"%c",ri->qual[i]);
                        }
                }else{
                        for(i =start;i < stop;i++){
                                fprintf(out,".");
                        }
                }
                fprintf(out,"\n");
                segment++;

                start = stop;
                while(ri->seq[start] == 65 && start < ri->len){
                        start++;
                }

                if(segment > 100){
                        exit(EXIT_FAILURE);
                }
        }
}


/** \fn int read_sam_chunk(struct read_info** ri,struct parameters* param,FILE* file)
    \brief Reads sequences from SAM / BAM file.
    Sequences are converted to 0-4 strings.

    \param ri Array of @ref read_info to hold the sequences.
    \param param @ref parameters .
    \param file Pointer to input file.

*/
int read_sam_chunk(struct read_info_buffer* rb, struct file_handler* f_handle)//  struct read_info** ri,struct parameters* param,FILE* file,int* buffer_count)
{
        FILE* file = NULL;
        struct read_info** ri = NULL;
        //int status;
        char* line = NULL;
        size_t line_len = 0;
        ssize_t nread;

        //char line[MAX_LINE];
        int column = 0;
        int i,j,g,tmp;
        //unsigned int pos;
        int read = 0;
        int c = 0;
        //int hit = 0;
        //int strand = 0;

        ri = rb->ri;
        rb->num_seq = 0;

        //*buffer_count = 0;

        ri = clear_read_info(ri,   rb->num_alloc);

        file = f_handle->f_ptr;

        //while(fgets(line, MAX_LINE, file)){
        while ((nread = getline(&line, &line_len, file)) != -1) {
                if(line[0] != '@'){
                        column = 1; //<QNAME>
                        tmp = 0;
                        //		hit = 0;
                        //		pos = 0xFFFFFFFFu;
                        for(j = 0;j < nread;j++){
                                tmp++;
                                if(isspace((int)line[j])){
                                        break;
                                }
                        }

                        MMALLOC(ri[c]->name,sizeof(unsigned char)* tmp);
                        for(j = 0;j < nread;j++){

                                if(isspace((int)line[j])){
                                        ri[c]->name[j] = 0;
                                        break;
                                }
                                ri[c]->name[j] = line[j];
                        }

                        for(i = 0; i < nread;i++){
                                if(line[i] == '\n'){
                                        break;
                                }
                                if(isspace((int)line[i])){
                                        column++;
                                        switch(column){
                                        case 2: // <FLAG>
                                                tmp = atoi(line+i+1);
                                                //strand = (tmp & 0x10);

                                                //WARNING - read should be reverse complemented if mapped to negative strand before tagdusting...

                                                /*tmp = atoi(line+i+1);
                                                  ri[c]->strand[hit] = (tmp & 0x10);
                                                  if(tmp == 4){
                                                  ri[c]->hits[hit] = 0;
                                                  }else{
                                                  ri[c]->hits[hit] = 1;
                                                  }
                                                  hit++;*/

                                                break;
                                        case 3: // <RNAME>

                                                break;
                                        case 4: // <POS>

                                                break;
                                        case 5: //  <MAPQ>

                                                //ri[c]->mapq =  atof(line +i +1);

                                                break;
                                        case 6: //  <CIGAR>

                                                break;
                                        case 7: //  <MRNM>
                                                break;
                                        case 8: //  <MPOS>
                                                break;
                                        case 9: //  <ISIZE>
                                                break;
                                        case 10: // <SEQ>

                                                tmp = 0;
                                                for(j = i+1;j < nread;j++){
                                                        tmp++;
                                                        if(isspace((int)line[j])){
                                                                break;
                                                        }
                                                }

                                                MMALLOC(ri[c]->seq,sizeof(unsigned char)* tmp);
                                                MMALLOC(ri[c]->labels,sizeof(unsigned char)* tmp);

                                                g = 0;
                                                for(j = i+1;j < nread;j++){

                                                        if(isspace((int)line[j])){
                                                                ri[c]->seq[g] = 0;
                                                                ri[c]->labels[g] = 0;
                                                                break;
                                                        }
                                                        ri[c]->seq[g] = nuc_code[(int)line[j]];
                                                        ri[c]->labels[g] = 0;

                                                        g++;
                                                }

                                                ri[c]->len = g;
                                                break;
                                        case 11: // <QUAL>
                                                tmp = 0;
                                                for(j = i+1;j < nread;j++){
                                                        tmp++;
                                                        if(isspace((int)line[j])){
                                                                break;
                                                        }
                                                }
                                                g= 0;
                                                MMALLOC(ri[c]->qual,sizeof(unsigned char)* tmp);
                                                for(j = i+1;j < nread;j++){

                                                        if(isspace((int)line[j])){
                                                                ri[c]->qual[g] = 0;
                                                                break;
                                                        }
                                                        ri[c]->qual[g] = line[j];
                                                        g++;
                                                }
                                                break;
                                        default:


                                                i = nread;
                                                break;
                                        }				}

                        }
                        /*tmp = byg_end("NM:i:", line  );
                        if(tmp){
                                ri[c]->read_type = atoi(line+tmp);
                        }else{
                                ri[c]->read_type = -1;
                                }*/

                        ///if(strand != 0){
                                //	ri[c]->seq = reverse_complement(ri[c]->seq,ri[c]->len);
                        //}




                        //ri[c]->hits[hit] = 0xFFFFFFFFu;

                        c++;
                        read++;
                        if(c ==  rb->num_alloc){
                                rb->num_seq = c;
                                MFREE(line);
                                return OK;
                                //*buffer_count = c;
                                //return c;
                        }
                }
        }
        rb->num_seq = c;
        MFREE(line);
        //*buffer_count = c;
        return OK;
        //return c;
ERROR:
        return FAIL;
        //return status;
}

//FASTQ files from CASAVA-1.8 Should have the following READ-ID format:
//@<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos> <read>:<is filtered>:<control number>:<index sequence>



/** \fn int read_fasta_fastq(struct read_info** ri,struct parameters* param,FILE *file)
    \brief Reads sequences from fasta / fastq file.
    Sequences are converted to 0-4 strings.

    \param ri Array of @ref read_info to hold the sequences.
    \param param @ref parameters .
    \param file Pointer to input file.

*/
//int read_sam_chunk(struct read_info_buffer* rb, struct file_handler* f_handle)//  struct read_info** ri,struct parameters* param,FILE* file,int* buffer_count)
int read_fasta_fastq(struct read_info_buffer* rb, struct file_handler* f_handle)//struct read_info** ri,struct parameters* param,FILE *file,int* buffer_count)
{
        FILE* file;
        struct read_info** ri = NULL;
        int park_pos = -1;
        char* line = NULL;
        size_t line_len = 0;
        ssize_t nread;
        //char line[MAX_LINE];
        int i;//,j;
        int seq_p = 0;
        int set = 0;
        int len = 0;
        int size = 0;
        //int status;

        ri = rb->ri;
        rb->num_seq = 0;
        //*buffer_count = 0;

        ri = clear_read_info(ri, rb->num_alloc);

        file = f_handle->f_ptr;

        while ((nread = getline(&line, &line_len, file)) != -1) {
                //while(fgets(line, MAX_LINE, file)){
                if((line[0] == '@' && !set)|| (line[0] == '>' && !set)){
                        //set sequence length of previous read

                        //check if there is still space....
                        //if(param->num_query == size){
                        //	fseek (file , -  strlen(line) , SEEK_CUR);
                        //	return size;
                        //}
                        park_pos++;
                        len = 0;
                        seq_p = 1;
                        for(i = 1;i < nread;i++){
                                len++;
                                if(iscntrl((int)line[i])){
                                        break;
                                }

                        }

                        //ri[park_pos]->hits[0] = 0;
                        //ri[park_pos]->strand[0] = 0;
                        MMALLOC(ri[park_pos]->name,sizeof(unsigned char)* (len+1));
                        for(i = 1;i < nread;i++){

                                if(iscntrl((int)line[i])){
                                        ri[park_pos]->name[i-1] = 0;
                                        break;
                                }
                                if(isspace((int)line[i])){
                                        ri[park_pos]->name[i-1] = ';';
                                }

                                ri[park_pos]->name[i-1] = line[i];
                        }
                        //fprintf(stderr,"LEN:%d	%s\n",len,ri[park_pos]->name);

                        set = 1;
                        size++;
                        //get ready to read quality if present
                }else if(line[0] == '+' && !set){
                        seq_p = 0;
                        set = 1;
                        //reading sequence or quality
                }else{
                        if(set){
                                if(seq_p){
                                        len = 0;
                                        for(i = 0;i < nread;i++){
                                                len++;
                                                if(iscntrl((int)line[i])){
                                                        break;
                                                }
                                        }
                                        //fprintf(stderr,"SEQ LEN:%d	%s\n",len,line);
                                        MMALLOC(ri[park_pos]->seq,sizeof(unsigned char)* (len+1));

                                        MMALLOC(ri[park_pos]->labels, sizeof(unsigned char)* (len+1));

                                        for(i = 0;i < nread;i++){
                                                if(iscntrl((int)line[i])){
                                                        ri[park_pos]->seq[i] = 0;
                                                        ri[park_pos]->labels[i] = 0;
                                                        break;
                                                }
                                                ri[park_pos]->seq[i] = nuc_code[(int)line[i]];
                                                ri[park_pos]->labels[i] = 0;
                                        }
                                        ri[park_pos]->len = len-1;
                                }else{
                                        len = 0;
                                        for(i = 0;i < nread;i++){
                                                len++;
                                                if(iscntrl((int)line[i])){
                                                        break;
                                                }

                                        }

                                        if(len-1 != ri[park_pos]->len ){
                                                MFREE(line);
                                                ERROR_MSG("Length of sequence and base qualities differ!.\n");

                                        }

                                        //fprintf(stderr,"QUAL LEN:%d\n",len);
                                        MMALLOC(ri[park_pos]->qual,sizeof(unsigned char)* (len+1));
                                        for(i = 0;i < nread;i++){
                                                if(iscntrl((int)line[i])){
                                                        ri[park_pos]->qual[i] = 0;
                                                        break;
                                                }
                                                ri[park_pos]->qual[i] = line[i];
                                        }
                                }
                        }
                        set = 0;
                }
                if( rb->num_alloc == size ){//here I know I am in the last entry AND filled the quality...
                        if(! f_handle->fasta && ri[park_pos]->qual){
                                //return size;
                                rb->num_seq = size;
                                //*buffer_count = size;
                                MFREE(line);
                                return OK;
                        }
                        if(f_handle->fasta && ri[park_pos]->seq){
                                //*buffer_count = size;
                                rb->num_seq = size;
                                MFREE(line);
                                return OK;

                        }
                }
        }
        rb->num_seq = size;
        MFREE(line);
        return OK;
ERROR:
        return FAIL;

}

/** \fn struct fasta* get_fasta(struct fasta* p,char *infile)
    \brief Old Kalign function to read in fasta sequences.
    All sequences are read into one array indexed by pointers.
    \param p @ref fasta struct to hold sequences.

    \param infile Pointer to input file.

*/

struct fasta* get_fasta(struct fasta* p,char *infile)
{
        //int status;
        MMALLOC(p,sizeof(struct fasta));
        p->string = 0;
        p->mer_hash = 0;
        p->s_index = 0;
        p->sn = 0;
        p->string = 0;
        p->boost = 0;
        p->max_len = 0;
        p->numseq = 0;
        p->string_len = 0;
        p->suffix = 0;

        p->string =  get_input_into_string(p->string,infile);
        if(!p->string){
                fprintf(stderr,"Analysing input %s ... nothing\n",infile);
                exit(-1);
        }

        p = read_fasta(p);

        return p;
ERROR:
        //KSLIB_MESSAGE(status,"Something went wrong in get_fasta.\n ");
        return NULL;
}



/** \fn unsigned char* get_input_into_string(unsigned char* string,char* infile)
    \brief Old Kalign function to copy file content into a string.
    All sequences are read into one array indexed by pointers.
    \param string target string.
    \param infile Pointer to input file.

*/

unsigned char* get_input_into_string(unsigned char* string,char* infile)
{
        long int i = 0;
        //int status;
        size_t bytes_read;

        FILE *file = NULL;
        RUNP(file = fopen(infile, "r"));
        //if((file = fopen(infile, "r")) == NULL) KSLIB_XEXCEPTION_SYS(kslEWRT,"Failed to open file:%s",infile);

        if (fseek(file,0,SEEK_END) != 0){
                (void)fprintf(stderr, "ERROR: fseek failed\n");
                (void)exit(EXIT_FAILURE);
        }
        i= ftell (file);
        if (fseek(file,0,SEEK_START) != 0){
                (void)fprintf(stderr, "ERROR: fseek failed\n");
                (void)exit(EXIT_FAILURE);
        }
        if(!string){
                MMALLOC(string,(i+1+18)* sizeof(unsigned char));

        }
        bytes_read = fread(string,sizeof(unsigned char), i, file);
        if(!bytes_read){
                fprintf(stderr,"Reading from file:%s failed \n", infile );
                exit(EXIT_FAILURE);
        }
        string[i] = 0;
        fclose(file);

        return string;
ERROR:
        //KSLIB_MESSAGE(status,"Something went wrong in get_input_into_string");
        return NULL;
}


/** \fn struct fasta* read_fasta(struct fasta* f)
    \brief Old Kalign function to set pointers to start of sequence names, sequences ...

    Also sets the sequence lengths.
    \param string target string.
    \param f @ref fasta .

*/


struct fasta* read_fasta(struct fasta* f)
{
        //int status;
        int c = 0;
        int n = 0;
        int i = 0;
        int j = 0;
        int len = 0;
        int stop = 0;
        int nbytes;

        nbytes = (int) strlen((char*) f->string);


        //aln->org_seq = 0;
        stop = 0;

        //count filenames....
        for (i =0;i < nbytes;i++){
                if (f->string[i] == '>'&& stop == 0){
                        f->numseq++;
                        stop = 1;
                }else if (f->string[i] == '\n'){
                        stop = 0;
                }else if(f->string[i] == '\r'){
                        f->string[i] = '\n';
                }
        }


        MMALLOC(f->sn, sizeof(unsigned char*)*f->numseq);
        //aln->c = 0;

        MMALLOC(f->s_index, sizeof(int)*(f->numseq+1));

        for(i = 0; i < f->numseq;i++){
                f->sn[i] = 0;
                f->s_index[i] = 0;
        }
        f->s_index[f->numseq] = 0;


        for (i =0;i < nbytes;i++){
                if (f->string[i] == '>'){
                        if(f->max_len < len){
                                f->max_len = len;
                        }

                        len = 0;
                        j = i+1;
                        while(f->string[j] != '\n'){
                                //	fprintf(stderr,"%c",aln->string[j]);
                                j++;
                        }
                        //fprintf(stderr,"	%d\n",j-i);
                        MMALLOC(f->sn[c],sizeof(char)*(j-i));
                        f->s_index[c] = n;
                        j = i+1;
                        len = 0;
                        while(f->string[j] != '\n'){
                                if(isspace((int)f->string[j])){
                                        f->sn[c][len] = '_';
                                }else{
                                        f->sn[c][len] = f->string[j];
                                }
                                len++;
                                j++;
                        }
                        f->sn[c][len] = 0;
                        f->string[n] = 'X';
                        n++;
                        i+= j-i;
                        c++;
                }else if(isalnum((int)f->string[i])){// if (aln->string[i] != '\n' && aln->string[i] != 0 ){
                        f->string[n] =   nuc_code[(int)f->string[i]];//  toupper(aln->string[i]);
                        n++;
                        len++;
                }
        }

        f->s_index[c] = n;
        f->string[n] = 'X';
        f->string[n+1] = 0;

        f->string_len = n+1;
        return f;
ERROR:
        //KSLIB_MESSAGE(status,"Something wrong in read_fasta.\n");
        return NULL;
}


/** \fn void free_fasta(struct fasta*f)

    \brief frees @ref fasta .
    \param f @ref fasta.
    \warning I used to free f->mer_hash even though it is allocated / used in this project....
*/
void free_fasta(struct fasta*f)
{
        int i;
        for (i =0;i < f->numseq;i++){
                MFREE(f->sn[i]);
        }
        if(f->mer_hash){
                MFREE(f->mer_hash);
        }
        if(f->suffix){
                MFREE(f->suffix);
        }

        //free(aln->mer_hash);
        MFREE(f->s_index);
        MFREE(f->string);
        MFREE(f->sn);

        MFREE(f);
}


int alloc_read_info_buffer(struct read_info_buffer** rb, int size)
{
        struct read_info_buffer* read_buffer = NULL;

        MMALLOC(read_buffer, sizeof(struct read_info_buffer));
        read_buffer->num_alloc = size;
        read_buffer->num_seq = 0;
        read_buffer->offset = 0;
        read_buffer->ri = NULL;

        RUNP(read_buffer->ri = malloc_read_info(read_buffer->ri, read_buffer->num_alloc));

        *rb = read_buffer;
        return OK;

ERROR:
        return FAIL;

}

void free_read_info_buffer(struct read_info_buffer* rb)
{
        if(rb){
                free_read_info(rb->ri, rb->num_alloc);
                MFREE(rb);
        }
}

struct read_info** malloc_read_info(struct read_info** ri, int numseq)
{
        int i;
        //int status;
        MMALLOC(ri, sizeof(struct read_info*) * numseq);

        for(i = 0; i < numseq;i++){
                ri[i] = 0;
                MMALLOC(ri[i], sizeof(struct read_info));

                ri[i]->seq = 0;
                ri[i]->name = 0;
                ri[i]->qual = 0;
                ri[i]->labels = 0;
                ri[i]->len = 0;
                ri[i]->bar_prob = 0;
                //ri[i]->mapq = -1.0;
                //ri[i]->barcode = -1;
                //ri[i]->fingerprint = -1;
                //ri[i]->read_type = 0;
                //ri[i]->barcode_string = NULL;
                //ri[i]->strand = malloc(sizeof(unsigned int)* (LIST_STORE_SIZE+1));
                //ri[i]->hits = malloc(sizeof(unsigned int)* (LIST_STORE_SIZE+1));
        }
        return ri;
ERROR:
        //KSLIB_MESSAGE(status,"Somehting wrong in malloc_read_info.\n");
        return NULL;
}

struct read_info** clear_read_info(struct read_info** ri, int numseq)
{
        int i;

        for(i = 0; i < numseq;i++){
                if(ri[i]->seq){
                        MFREE(ri[i]->seq);
                }
                if(ri[i]->name){
                        MFREE(ri[i]->name);
                }
                if(ri[i]->qual){
                        MFREE(ri[i]->qual);
                }
                if(ri[i]->labels){
                        MFREE(ri[i]->labels);
                }
                ri[i]->seq = 0;
                ri[i]->name = 0;
                ri[i]->qual = 0;
                ri[i]->labels = 0;
                ri[i]->len = 0;
                ri[i]->bar_prob = 0;
                //ri[i]->mapq = -1.0;
                //ri[i]->barcode = -1;
                //ri[i]->fingerprint = -1;
                //ri[i]->read_type = 0;
                //ri[i]->barcode_string = NULL;
                //ri[i]->strand = malloc(sizeof(unsigned int)* (LIST_STORE_SIZE+1));
                //ri[i]->hits = malloc(sizeof(unsigned int)* (LIST_STORE_SIZE+1));
        }
        return ri;
}

void free_read_info(struct read_info** ri, int numseq)
{
        int i;

        for(i = 0; i < numseq;i++){
                //free(ri[i]->strand);
                //free(ri[i]->hits);


                if(ri[i]->labels){
                        MFREE(ri[i]->labels);
                }
                if(ri[i]->name){
                        MFREE(ri[i]->name);
                }
                if(ri[i]->seq){
                        MFREE(ri[i]->seq);
                }
                if(ri[i]->qual){
                        MFREE(ri[i]->qual );
                }

                MFREE(ri[i]);
        }
        MFREE(ri);
}





int compare_read_names(char* name1, char* name2)
{

#ifdef UTEST
        int detected = -1;
#else
        static int detected = -1;
#endif
        char instrument_R1[100];
        int run_id_R1 = 0;
        char flowcell_R1[100];
        int flowcell_lane_R1= 0;
        int tile_number_R1= 0;
        int x_coordinate_R1= 0;
        int y_coordinate_R1= 0;

        char instrument_R2[100];
        int run_id_R2= 0;
        char flowcell_R2[100];
        int flowcell_lane_R2= 0;
        int tile_number_R2= 0;
        int x_coordinate_R2= 0;
        int y_coordinate_R2= 0;

        int i;

        int number_of_values_found = 0;

        instrument_R1[0] = 0;
        instrument_R2[0] = 0;

        flowcell_R1[0] = 0;
        flowcell_R2[0] = 0;

        if(detected == -1){
                //option 1: casava 1.8
                // name should look like this:@EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG
                number_of_values_found =sscanf(name1,"%[^:]:%d:%[^:]:%d:%d:%d:%d ", instrument_R1,&run_id_R1,flowcell_R1,&flowcell_lane_R1,&tile_number_R1,&x_coordinate_R1,&y_coordinate_R1 );
                if(number_of_values_found == 7){
                        detected = 1;
                        LOG_MSG("Detected casava 1.8 format.\n");
                        //param->messages = append_message(param->messages, param->buffer);
                }
                //fprintf(stderr,"casava 1.8?:%d %s\n",number_of_values_found, name1);
                //fprintf(stderr,"%s\n%d\n%s\n%d\n%d\n%d\n%d\n", instrument_R1,run_id_R1,flowcell_R1,flowcell_lane_R1,tile_number_R1,x_coordinate_R1,y_coordinate_R1);

        }

        if(detected == -1){
                //option 2: casava 1.7
                // name should look like this:@HWUSI-EAS100R:6:73:941:1973#0/1
                //HWUSI-EAS747_0040_FC64GRTAAXX:8:1:3268:1065#0/1
                number_of_values_found =sscanf(name1,"%[^:]:%d:%d:%d:%d", instrument_R1,&flowcell_lane_R1,&tile_number_R1,&x_coordinate_R1,&y_coordinate_R1);
                //fprintf(stderr,"casava 1.7?:%d %s\n", number_of_values_found,name1);
                //fprintf(stderr,"%s\n%d\n%d\n%d\n%d\n", instrument_R1,flowcell_lane_R1,tile_number_R1,x_coordinate_R1,y_coordinate_R1);

                if(number_of_values_found == 5){
                        detected = 2;
                        LOG_MSG("Detected casava <1.7 format.\n");
                        //param->messages = append_message(param->messages, param->buffer);
                }
        }

        if(detected == -1){
                detected = 1000;
        }

        if(detected == 1){
                number_of_values_found =sscanf(name1,"%[^:]:%d:%[^:]:%d:%d:%d:%d ", instrument_R1,&run_id_R1,flowcell_R1,&flowcell_lane_R1,&tile_number_R1,&x_coordinate_R1,&y_coordinate_R1 );
                if(number_of_values_found != 7){
                        ERROR_MSG("File name %s\n does not match detected casava 1.8 format.\n",name1);
                        //param->messages = append_message(param->messages, param->buffer);
                        return 1;

                }

                number_of_values_found =sscanf(name2,"%[^:]:%d:%[^:]:%d:%d:%d:%d ", instrument_R2,&run_id_R2,flowcell_R2,&flowcell_lane_R2,&tile_number_R2,&x_coordinate_R2,&y_coordinate_R2 );
                if(number_of_values_found != 7){
                        ERROR_MSG("File name %s\n does not match detected casava 1.8 format.\n",name2);
#ifdef UTEST
                        fprintf(stderr,"@%s:%d:%d:%d:%d\n", instrument_R1,flowcell_lane_R1,tile_number_R1,x_coordinate_R1,y_coordinate_R1);
                        fprintf(stderr,"@%s:%d:%d:%d:%d\n", instrument_R2,flowcell_lane_R2,tile_number_R2,x_coordinate_R2,y_coordinate_R2);
#endif
                        //param->messages = append_message(param->messages, param->buffer);
                        return 1;
                }

                if(y_coordinate_R1 != y_coordinate_R2){
                        ERROR_MSG("Files seem to contain reads in different order:\n%s\n%s\n",name1, name2);
#ifdef UTEST
                        fprintf(stderr,"@%s:%d:%d:%d:%d\n", instrument_R1,flowcell_lane_R1,tile_number_R1,x_coordinate_R1,y_coordinate_R1);
                        fprintf(stderr,"@%s:%d:%d:%d:%d\n", instrument_R2,flowcell_lane_R2,tile_number_R2,x_coordinate_R2,y_coordinate_R2);
#endif
                        //param->messages = append_message(param->messages, param->buffer);
                        return 1;
                }

                if(x_coordinate_R1 != x_coordinate_R2){
                        ERROR_MSG("Files seem to contain reads in different order:\n%s\n%s\n",name1, name2);
#ifdef UTEST
                        fprintf(stderr,"@%s:%d:%d:%d:%d\n", instrument_R1,flowcell_lane_R1,tile_number_R1,x_coordinate_R1,y_coordinate_R1);
                        fprintf(stderr,"@%s:%d:%d:%d:%d\n", instrument_R2,flowcell_lane_R2,tile_number_R2,x_coordinate_R2,y_coordinate_R2);
#endif
                        //param->messages = append_message(param->messages, param->buffer);
                        return 1;
                }

                if(tile_number_R1 !=  tile_number_R2){
                        ERROR_MSG("Files seem to contain reads in different order:\n%s\n%s\n",name1, name2);
#ifdef UTEST
                        fprintf(stderr,"@%s:%d:%d:%d:%d\n", instrument_R1,flowcell_lane_R1,tile_number_R1,x_coordinate_R1,y_coordinate_R1);
                        fprintf(stderr,"@%s:%d:%d:%d:%d\n", instrument_R2,flowcell_lane_R2,tile_number_R2,x_coordinate_R2,y_coordinate_R2);
#endif
                        //param->messages = append_message(param->messages, param->buffer);
                        return 1;
                }

                if(flowcell_lane_R1 !=  flowcell_lane_R2){
                        ERROR_MSG("Files seem to contain reads in different order:\n%s\n%s\n",name1, name2);
#ifdef UTEST
                        fprintf(stderr,"@%s:%d:%d:%d:%d\n", instrument_R1,flowcell_lane_R1,tile_number_R1,x_coordinate_R1,y_coordinate_R1);
                        fprintf(stderr,"@%s:%d:%d:%d:%d\n", instrument_R2,flowcell_lane_R2,tile_number_R2,x_coordinate_R2,y_coordinate_R2);
#endif
                        //param->messages = append_message(param->messages, param->buffer);
                        return 1;
                }

                if(strcmp(flowcell_R1,flowcell_R2)){
                        ERROR_MSG("Files seem to contain reads in different order:\n%s\n%s\n",name1, name2);
#ifdef UTEST
                        fprintf(stderr,"@%s:%d:%d:%d:%d\n", instrument_R1,flowcell_lane_R1,tile_number_R1,x_coordinate_R1,y_coordinate_R1);
                        fprintf(stderr,"@%s:%d:%d:%d:%d\n", instrument_R2,flowcell_lane_R2,tile_number_R2,x_coordinate_R2,y_coordinate_R2);
#endif
                        //param->messages = append_message(param->messages, param->buffer);
                        return 1;
                }
                if(run_id_R1 !=  run_id_R2){
                        ERROR_MSG("Files seem to contain reads in different order:\n%s\n%s\n",name1, name2);
#ifdef UTEST
                        fprintf(stderr,"@%s:%d:%d:%d:%d\n", instrument_R1,flowcell_lane_R1,tile_number_R1,x_coordinate_R1,y_coordinate_R1);
                        fprintf(stderr,"@%s:%d:%d:%d:%d\n", instrument_R2,flowcell_lane_R2,tile_number_R2,x_coordinate_R2,y_coordinate_R2);
#endif
                        //param->messages = append_message(param->messages, param->buffer);
                        return 1;
                }
                if(strcmp(instrument_R1,instrument_R2)){
                        ERROR_MSG("Files seem to contain reads in different order:\n%s\n%s\n",name1, name2);
#ifdef UTEST
                        fprintf(stderr,"@%s:%d:%d:%d:%d\n", instrument_R1,flowcell_lane_R1,tile_number_R1,x_coordinate_R1,y_coordinate_R1);
                        fprintf(stderr,"@%s:%d:%d:%d:%d\n", instrument_R2,flowcell_lane_R2,tile_number_R2,x_coordinate_R2,y_coordinate_R2);
#endif
                        //param->messages = append_message(param->messages, param->buffer);
                        return 1;
                }



        }

        if(detected == 2){
                number_of_values_found =sscanf(name1,"%[^:]:%d:%d:%d:%d", instrument_R1,&flowcell_lane_R1,&tile_number_R1,&x_coordinate_R1,&y_coordinate_R1);

                if(number_of_values_found != 5){
                        ERROR_MSG("File name %s\n does not match detected casava <1.8 format.\n",name1);
#ifdef UTEST
                        fprintf(stderr,"@%s:%d:%d:%d:%d\n", instrument_R1,flowcell_lane_R1,tile_number_R1,x_coordinate_R1,y_coordinate_R1);
                        fprintf(stderr,"@%s:%d:%d:%d:%d\n", instrument_R2,flowcell_lane_R2,tile_number_R2,x_coordinate_R2,y_coordinate_R2);
#endif
                        //param->messages = append_message(param->messages, param->buffer);
                        return 1;
                }
                number_of_values_found =sscanf(name2,"%[^:]:%d:%d:%d:%d", instrument_R2,&flowcell_lane_R2,&tile_number_R2,&x_coordinate_R2,&y_coordinate_R2);

                if(number_of_values_found != 5){
                        ERROR_MSG("File name %s\n does not match detected casava <1.8 format.\n",name2);
#ifdef UTEST
                        fprintf(stderr,"@%s:%d:%d:%d:%d\n", instrument_R1,flowcell_lane_R1,tile_number_R1,x_coordinate_R1,y_coordinate_R1);
                        fprintf(stderr,"@%s:%d:%d:%d:%d\n", instrument_R2,flowcell_lane_R2,tile_number_R2,x_coordinate_R2,y_coordinate_R2);
#endif
                        //param->messages = append_message(param->messages, param->buffer);
                        return 1;
                }

                if(y_coordinate_R1 != y_coordinate_R2){
                        ERROR_MSG("Files seem to contain reads in different order:\n%s\n%s\n",name1, name2);
#ifdef UTEST
                        fprintf(stderr,"@%s:%d:%d:%d:%d\n", instrument_R1,flowcell_lane_R1,tile_number_R1,x_coordinate_R1,y_coordinate_R1);
                        fprintf(stderr,"@%s:%d:%d:%d:%d\n", instrument_R2,flowcell_lane_R2,tile_number_R2,x_coordinate_R2,y_coordinate_R2);
#endif
                        //param->messages = append_message(param->messages, param->buffer);
                        return 1;
                }

                if(x_coordinate_R1 != x_coordinate_R2){
                        ERROR_MSG("Files seem to contain reads in different order:\n%s\n%s\n",name1, name2);
#ifdef UTEST
                        fprintf(stderr,"@%s:%d:%d:%d:%d\n", instrument_R1,flowcell_lane_R1,tile_number_R1,x_coordinate_R1,y_coordinate_R1);
                        fprintf(stderr,"@%s:%d:%d:%d:%d\n", instrument_R2,flowcell_lane_R2,tile_number_R2,x_coordinate_R2,y_coordinate_R2);
#endif
                        //param->messages = append_message(param->messages, param->buffer);
                        return 1;
                }

                if(tile_number_R1 !=  tile_number_R2){
                        ERROR_MSG("Files seem to contain reads in different order:\n%s\n%s\n",name1, name2);
#ifdef UTEST
                        fprintf(stderr,"@%s:%d:%d:%d:%d\n", instrument_R1,flowcell_lane_R1,tile_number_R1,x_coordinate_R1,y_coordinate_R1);
                        fprintf(stderr,"@%s:%d:%d:%d:%d\n", instrument_R2,flowcell_lane_R2,tile_number_R2,x_coordinate_R2,y_coordinate_R2);
#endif
                        //param->messages = append_message(param->messages, param->buffer);
                        return 1;
                }

                if(flowcell_lane_R1 !=  flowcell_lane_R2){
                        ERROR_MSG("Files seem to contain reads in different order:\n%s\n%s\n",name1, name2);
#ifdef UTEST
                        fprintf(stderr,"@%s:%d:%d:%d:%d\n", instrument_R1,flowcell_lane_R1,tile_number_R1,x_coordinate_R1,y_coordinate_R1);
                        fprintf(stderr,"@%s:%d:%d:%d:%d\n", instrument_R2,flowcell_lane_R2,tile_number_R2,x_coordinate_R2,y_coordinate_R2);
#endif
                        //param->messages = append_message(param->messages, param->buffer);
                        //return 1;
                }

                if(strcmp(instrument_R1,instrument_R2)){
                        ERROR_MSG("Files seem to contain reads in different order:\n%s\n%s\n",name1, name2);
#ifdef UTEST
                        fprintf(stderr,"@%s:%d:%d:%d:%d\n", instrument_R1,flowcell_lane_R1,tile_number_R1,x_coordinate_R1,y_coordinate_R1);
                        fprintf(stderr,"@%s:%d:%d:%d:%d\n", instrument_R2,flowcell_lane_R2,tile_number_R2,x_coordinate_R2,y_coordinate_R2);
#endif
                        //param->messages = append_message(param->messages, param->buffer);
                        //return 1;
                }
        }
        if(detected == 1000){
                number_of_values_found = 0;
                for(i = 0; i < strlen(name1);i++){
                        if(isspace(name1[i]) || name1[i] == ';'){
                                break;
                        }
                        if(name1[i] != name2[i]){
                                number_of_values_found = 1;
                                break;
                        }
                        /*if(isspace(name1[i])){
                          name1[i] = 0;
                          name2[i] = 0;
                          }

                          if(name1[i] == ';'){
                          name1[i] = 0;
                          name2[i] = 0;
                          }*/
                }
                if(number_of_values_found){
                        ERROR_MSG("Files seem to contain reads in different order:\n%s\n%s\n",name1,name2);
#ifdef UTEST
                        fprintf(stderr,"@%s:%d:%d:%d:%d\n", instrument_R1,flowcell_lane_R1,tile_number_R1,x_coordinate_R1,y_coordinate_R1);
                        fprintf(stderr,"@%s:%d:%d:%d:%d\n", instrument_R2,flowcell_lane_R2,tile_number_R2,x_coordinate_R2,y_coordinate_R2);
#endif
                        //param->messages = append_message(param->messages, param->buffer);
                        return FAIL;
                }
        }


        return OK;
ERROR:
        return FAIL;
}



#ifdef UTEST
int main (int argc,char * argv[]) {
        struct parameters* param = 0;
        int status;
        char* name1 = 0;

        char* name2 = 0;
        MMALLOC(param,sizeof(struct parameters));

        param->read_structure = 0;
        param->read_structure_R1 = 0;
        param->read_structure_R2 = 0;
        param->read_structures = NULL;
        param->confidence_thresholds = NULL;//ence

        param->outfile = 0;
        param->infile =0 ;
        param->buffer =0;
        param->messages = 0;

        MMALLOC(param->buffer,sizeof(char) * MAX_LINE);

        MMALLOC(name1, sizeof(char)* 1000);
        MMALLOC(name2, sizeof(char)* 1000);

        sprintf (name1, "EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG");
        sprintf (name2, "EAS139:136:FC706VJ:2:2104:15343:197393 2:Y:18:ATCACG");

        int i = 0;

        fprintf(stdout,"Running I/O Unit tests.\n");
        i = compare_read_names(param, name1, name2);
        fprintf(stdout,"	cmp1:\n	%s\n	%s\n	Result:%d\n",name1,name2,i);
        if(i == 1){
                sprintf(param->buffer , "ERROR - names should match\n");
                param->messages = append_message(param->messages, param->buffer  );
                MFREE(name1);
                MFREE(name2);
                free_param(param);
                exit(EXIT_FAILURE);
        }
        name1[0] = 0;
        name2[0] = 0;

        sprintf (name1, "EAS139:136:FC706VJ:2:2104:15343:197393 2:Y:18:ATCACG");
        sprintf (name2, "EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG");

        i = compare_read_names(param, name1, name2);
        fprintf(stdout,"	cmp2:\n	%s\n	%s\n	Result:%d\n",name1,name2,i);
        if(i == 1){
                sprintf(param->buffer , "ERROR - names should match\n");
                param->messages = append_message(param->messages, param->buffer  );
                MFREE(name1);
                MFREE(name2);
                free_param(param);
                exit(EXIT_FAILURE);
        }


        name1[0] = 0;
        name2[0] = 0;

        sprintf (name1, "HWUSI-EAS100R:6:73:941:1973#0/1");
        sprintf (name2, "HWUSI-EAS100R:6:73:941:1973#0/2");

        i = compare_read_names(param, name1, name2);
        fprintf(stdout,"	cmp3:\n	%s\n	%s\n	Result:%d\n",name1,name2,i);
        if(i == 1){
                sprintf(param->buffer , "ERROR - names should match\n");
                param->messages = append_message(param->messages, param->buffer  );
                MFREE(name1);
                MFREE(name2);
                free_param(param);
                exit(EXIT_FAILURE);
        }

        name1[0] = 0;
        name2[0] = 0;

        sprintf (name1, "HWUSI-EAS100R:6:73:941:1973#0/2");
        sprintf (name2, "HWUSI-EAS100R:6:73:941:1973#0/1");

        i = compare_read_names(param, name1, name2);
        fprintf(stdout,"	cmp4:\n	%s\n	%s\n	Result:%d\n",name1,name2,i);
        if(i == 1){
                sprintf(param->buffer , "ERROR - names should match\n");
                param->messages = append_message(param->messages, param->buffer  );
                MFREE(name1);
                MFREE(name2);
                free_param(param);
                exit(EXIT_FAILURE);
        }

        name1[0] = 0;
        name2[0] = 0;

        sprintf (name1, "EAS139:136:FC706VJ:2:2104:15343:197393 2:Y:18:ATCACG");
        sprintf (name2, "HWUSI-EAS100R:6:73:941:1973#0/1");

        i = compare_read_names(param, name1, name2);
        fprintf(stdout,"	cmp5:\n	%s\n	%s\n	Result:%d\n",name1,name2,i);
        if(i == 0){
                sprintf(param->buffer , "ERROR - names should not match\n");
                param->messages = append_message(param->messages, param->buffer  );
                MFREE(name1);
                MFREE(name2);
                free_param(param);
                exit(EXIT_FAILURE);
        }

        name1[0] = 0;
        name2[0] = 0;

        sprintf (name1, "EAS139:136:FC706VJ:2:2104:15343:197393 2:Y:18:ATCACG");
        sprintf (name2, "EAS139:136:FC706VJ:2:2104:15343:197393 1:N:18:GGGACG");

        i = compare_read_names(param, name1, name2);
        fprintf(stdout,"	cmp6:\n	%s\n	%s\n	Result:%d\n",name1,name2,i);
        if(i == 1){
                sprintf(param->buffer , "ERROR - names should match\n");
                param->messages = append_message(param->messages, param->buffer  );
                MFREE(name1);
                MFREE(name2);
                free_param(param);
                exit(EXIT_FAILURE);
        }

        name1[0] = 0;
        name2[0] = 0;

        sprintf (name1, "EAS139:136:FC706VJ:2:2104:15343:197393 2:Y:18:ATCACG");
        sprintf (name2, "EAS139:136:FC706VJ:2:2104:15344:197393 1:N:18:GGGACG");

        i = compare_read_names(param, name1, name2);
        fprintf(stdout,"	cmp7:\n	%s\n	%s\n	Result:%d\n",name1,name2,i);
        if(i == 0){
                sprintf(param->buffer , "ERROR - names should not match\n");
                param->messages = append_message(param->messages, param->buffer  );
                MFREE(name1);
                MFREE(name2);
                free_param(param);
                exit(EXIT_FAILURE);
        }

        name1[0] = 0;
        name2[0] = 0;


        sprintf (name1, "HWUSI-EAS747_0040_FC64GRTAAXX:8:1:3268:1065#0/1");
        sprintf (name2, "HWUSI-EAS747_0040_FC64GRTAAXX:8:1:3268:1065#0/2");
        i = compare_read_names(param, name1, name2);
        fprintf(stdout,"	cmp8:\n	%s\n	%s\n	Result:%d\n",name1,name2,i);
        if(i == 1){
                sprintf(param->buffer , "ERROR - names should  match\n");
                param->messages = append_message(param->messages, param->buffer  );
                MFREE(name1);
                MFREE(name2);
                free_param(param);
                exit(EXIT_FAILURE);
        }



        MFREE(name1);
        MFREE(name2);

        free_param(param);
        return kslOK;
ERROR:
        KSLIB_MESSAGE(status,"Something wrong in main of io itest.\n");
        return kslFAIL;
}


#endif
