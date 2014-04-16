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

/*! \file sim.c
 \brief Functions to simulate sequences.
 
 \author Timo Lassmann
 \bug No known bugs.
 */

#include <stdio.h>
#include <time.h>
#include "interface.h"
#include <math.h>
#include <assert.h>
#include <ctype.h>
#include "misc.h"
#include "io.h"
#include "sim.h"
#include "nuc_code.h"
#include <errno.h>
#include <unistd.h>

/** \fn void simulation_for_benchmark(struct parameters* param)
 \brief Prints out simulates sequences for comparison to Btrim,cutadapt and fastx..
 
 \param param @parameters.
 
 
 */
void simulation_for_benchmark(struct parameters* param)
{
	struct read_info** ri = 0;
	//char template[35], dirpath[22], *dp;
	//int fd = -1;
	//./tagdust -sim_numseq 10 -sim_5seq AAAAAA -sim_3seq TTTT -sim_readlen 10 -sim_readlen_mod 5 -sim_random_frac 0.001 -o fafa -sim_sequenced_len 100  -sim_barlen 6 -sim_barnum 4
	//./tagdust -sim_numseq 10  -sim_5seq agggaggacgatgcgg   -sim_3seq gtgtcagtcacttccagcgg  -sim_readlen 25 -sim_readlen_mod 5 -sim_random_frac 0.5 -o ~/tmp  -sim_sequenced_len 100   -sim_error_rate 0.02 -sim_InDel_frac 0.1 -sim_sequenced_len 1098 -sim_barlen 10 -sim_barnum 24
		
	int i,j,c,n,n_dash,barcode_used,numseq;//,no5,no3;
	float r = 0.0;
	FILE* file = 0;
	FILE* file2 = 0;
	int errors_allowed = 6;

	char* runid;
	char* basedir;
	char* dp;
	char* read = 0;
	char* sequenced_read = 0;
	char* sequenced_read_mutated = 0;
	
	char* outfile = 0;
	
	char* command = 0;
	char alpha[6] = "ACGTNN";
	
	
	//char btrimfiller[] = "ZZZZZZZ";
	
	unsigned int seed = (unsigned int) (time(NULL) * ( 42));
	
	if(!param->outfile){
		fprintf(stderr,"No output file name given.\n");
		
		free_param(param);
		exit(EXIT_FAILURE);
	}
	if(param->sim_readlen > 64){
		fprintf(stderr,"Warning: setting sim_readlen to 64.\n");
	}
	if(param->sim_sequenced_len < 10 ){
		fprintf(stderr,"Warning: sim_sequenced_len to short.\n");
		free_param(param);
		exit(EXIT_FAILURE);

	}
	
	if(param->sim_sequenced_len > 1000){
		param->sim_sequenced_len  = param->sim_sequenced_len  - 1000;
		//
		// this adds totally random sequences..
		//
		c = 0; // random sequence length
		if(param->sim_5seq){
			c += (int) strlen (param->sim_5seq);
		}
		if(param->sim_3seq){
			c += (int) strlen (param->sim_3seq);
		}
		
		c +=param->sim_barlen ;
		
		c+= param->sim_readlen;
		
		param->sim_sequenced_len  = (int)((float)c * (float) param->sim_sequenced_len / 100.0);
		
		
	}
	
	outfile = malloc(sizeof(char)*(200));
	assert(outfile != 0);
	read = malloc(sizeof(char)* 200);
	sequenced_read = malloc(sizeof(char)* 200);
	sequenced_read_mutated = malloc(sizeof(char)* 220);
	
	runid  = malloc(sizeof(char)* 220);
	basedir =  malloc(sizeof(char)* 220);
	int code = rand_r(&seed);
	if(param->sim_5seq || param->sim_3seq){
		c = sprintf(runid, "%d_%f_%s_%s_%d_%d_%d_%d_%f_%f_%d",param->sim_numseq,
			   param->sim_random_frac,
			   param->sim_5seq,
			   param->sim_3seq,
			   param->sim_barnum,
			   param->sim_barlen,
			   param->sim_readlen,
			   param->sim_readlen_mod,
			   param->sim_error_rate,
			   param->sim_InDel_frac,
			   param->sim_sequenced_len);
	}else {
		c = sprintf(runid, "%d_%f_%d_%d_%d_%d_%f_%f_%d",param->sim_numseq,
			   param->sim_random_frac,
			   param->sim_barnum,
			   param->sim_barlen,
			   param->sim_readlen,
			   param->sim_readlen_mod,
			   param->sim_error_rate,
			   param->sim_InDel_frac,
			   param->sim_sequenced_len);
	}

	fprintf(stderr,"%s\n",runid);
	
	for(i = 0; i < strlen(runid);i++){
		code = code ^ ((int)runid[i] + i);
	}
	
	fprintf(stderr,"%d	%X\n",code,code);
	c = sprintf(runid,"%X",code);
	
	sprintf (basedir, "%s%s%s",param->outfile,runid,"tmpXXXXXX");
	if ((dp = mkdtemp(basedir)) != NULL) {
		fprintf(stderr,"NEW TMP DIR:%s\n",dp);
	}else{
		fprintf(stderr, "Failed to create temp directory.\n%s	%s\n",strerror(errno),basedir);
	}
	//exit(0);
	//if((param->outfile = )
	/*param->outfile
	strcpy(dirpath, "/tmp/tempdir.XXXXXXXX");
	if ((dp = mkdtemp(dirpath)) != NULL) {
		sprintf(template, "%s/%s", dp, "tempXXXXXXXX");
		if ((fd = mkstemp(template)) != -1) {
			unlink(template);
			close(fd);
		} else
			fprintf(stderr, "Failed to create temp file.\n%s\n",
			        strerror(errno));
		rmdir(dirpath);
	} else{
		fprintf(stderr, "Failed to create temp directory.\n%s\n",
		        strerror(errno));
	}*/
	

	
	
	char** barcode = 0;
	if(param->sim_barnum){
		int num_barcode = param->sim_barnum;
		if(param->sim_barlen < 16){
			if( (int) powf(4.0, (float)param->sim_barlen) <= num_barcode){
				num_barcode =  (int) powf(4.0, (float)param->sim_barlen) ;
			}
			//fprintf(stderr,"%f\n",  powf(4.0, (float)param->sim));
		}
		//int num_barcode = (int) powf(4.0, (float)param->sim) ;
		
		
		//fprintf(stderr,"Number of Barcodes:%d\n",num_barcode );
		
		barcode = malloc(sizeof(char*)* num_barcode);
		assert(barcode!=0);
		
		for(i = 0; i < num_barcode;i++){
			barcode[i] = malloc(sizeof(char) *(param->sim_barlen+1) );
			assert(barcode[i]!=0);
		}
		
		n_dash = 100;
		c = 0;
		
		while(n_dash){
			n_dash = 0;
			//fprintf(stderr,"Errors allowed: %d\n");
			while(c != num_barcode){
				for(i = 0; i < param->sim_barlen;i++){
					r = (float)rand_r(&seed)/(float)RAND_MAX;
					if(r < 0.25){
						barcode[c][i] = 'A';
					}else if(r < 0.5){
						barcode[c][i] = 'C';
					}else if(r < 0.75){
						barcode[c][i] = 'G';
					}else{
						barcode[c][i] = 'T';
					}
				}
				barcode[c][param->sim_barlen] = 0;
				j = 1;
				for(i = 0; i < c;i++){
					if(bpm(barcode[c], barcode[i], param->sim_barlen, param->sim_barlen)  <= errors_allowed){
						j = 0;
						break;
					}
				}
				if(j){
					c++;
				}
				//
				n_dash++;
				if(n_dash == 1000000){
					break;
				}
			}
			
			if(n_dash == 1000000){
				errors_allowed--;
			}
			
		}
		
		if(errors_allowed == 0){// two identical barcodes were found.... I will simply pick the fist numbarcode sequences...
			for(i = 0; i < num_barcode ;i++){
				for(j = 0; j < param->sim_barlen;j++){
					barcode[i][j] = alpha[0x3 & (i >> (j*2))];
				}
			}
			
			errors_allowed = 100;
			for(i = 0; i < num_barcode-1 ;i++){
				
				for(j= i+1; j < num_barcode ;j++){
					c = bpm(barcode[i], barcode[j], param->sim_barlen, param->sim_barlen);
					if(c < errors_allowed){
						errors_allowed = c;
					}
				}
			}
			
			
		}
		fprintf(stderr,"Edit Distance between barcodes: %d\n", errors_allowed);
		for(i = 0; i < num_barcode;i++){
			fprintf(stderr,"%d: %s\n",i,barcode[i]);
		}
		
		sprintf (outfile, "%s/%sfastxbarcodefile.txt",dp,runid);
		
		
		if ((file = fopen(outfile, "w")) == NULL){
			fprintf(stderr,"can't open output\n");
			exit(-1);
		}
		
		for(i = 0 ;i < param->sim_barnum;i++){
			fprintf(file,"BC%d %s\n",i,barcode[i]);
			
		}
		fclose(file);

	}
	
	sprintf (outfile, "%s/%sbtrim_pattern.txt",dp,runid);
	if ((file = fopen(outfile, "w")) == NULL){
		fprintf(stderr,"can't open output\n");
		exit(-1);
	}
	/*no5 = 0;
	no3 = 0;
	if(!param->sim_5seq){
		param->sim_5seq = btrimfiller;
		no5 = 1;
	}
	
	if(!param->sim_3seq){
		param->sim_3seq = btrimfiller;
		no3 = 1;
	}*/
	
	if(param->sim_barnum){
		for(i = 0 ;i < param->sim_barnum;i++){
			fprintf(file,"%s%s %s\n",param->sim_5seq,barcode[i], param->sim_3seq);
		
		}
	}else{
		fprintf(file,"%s %s\n",param->sim_5seq, param->sim_3seq);

	}
	/*if(no5){
		param->sim_5seq = 0;
	}
	if(no3){
		param->sim_3seq = 0;
	}*/
	
	fclose(file);
	
	sprintf (outfile, "%s/%sread.fq",dp,runid);
	if ((file = fopen(outfile, "w")) == NULL){
		fprintf(stderr,"can't open output\n");
		exit(-1);
	}

	
	//Here I simulate sequences containing the read....
	for(i = 0; i < (int)((float) param->sim_numseq * (1.0-param->sim_random_frac));i++){
		for(j = 0; j < 200;j++){
			sequenced_read[j] = 0;
		}
		
		if(param->sim_5seq){
			strcat ( sequenced_read,  param->sim_5seq);
		}
		barcode_used = 0;
		if(param->sim_barnum){
			barcode_used = (int) (rand_r(&seed) % (int) (param->sim_barnum)) ;
		
			strcat ( sequenced_read, barcode[barcode_used]);
		}
		
		c = param->sim_readlen - param->sim_readlen_mod +  (int) (rand_r(&seed) % (int) (param->sim_readlen_mod*2)) ;
		//fprintf(stderr,"%d\n",c);
		for(j = 0; j < c;j++){
			r = (float)rand_r(&seed)/(float)RAND_MAX;
			if(r < 0.25){
				n = 'A';
			}else if(r < 0.5){
				n = 'C';
			}else if(r < 0.75){
				n = 'G';
			}else{
				n = 'T';
			}
			read[j] = n;
		}
		read[c] = 0;
		//fprintf(stderr,"%s\n", read);

		strcat ( sequenced_read, read);
		
		if(param->sim_3seq){
			strcat ( sequenced_read,  param->sim_3seq);
		}
		//fprintf(stderr,"%s\n", sequenced_read);
		c = 0;
		for(j = 0;j <  strlen(sequenced_read) ;j++){
			r = (float)rand_r(&seed)/(float)RAND_MAX;
			if(r <= param->sim_error_rate){
				//we have an error
				
				
				r = (float)rand_r(&seed)/(float)RAND_MAX;
				if(r <= param->sim_InDel_frac){
					//indel++;
					// we have an indel (only considering single nucleotide.....
					r = (float)rand_r(&seed)/(float)RAND_MAX;
					if(r <= 0.5){
						//insertion
						//n_dash = read[j];
						r = (float)rand_r(&seed)/(float)RAND_MAX;
						if(r < 0.25){
							n = 'A';
						}else if(r < 0.5){
							n = 'C';
						}else if(r < 0.75){
							n = 'G';
						}else{
							n = 'T';
						}
						sequenced_read_mutated[c] = sequenced_read[j];
						c++;
						sequenced_read_mutated[c] = n;
						c++;
						
						
						
						
					}else{
						//deletion
					}
					
				}else{
					//mismatches++;
					n = sequenced_read[j];
					
					while(n == sequenced_read[j]){
						r = (float)rand_r(&seed)/(float)RAND_MAX;
						if(r < 0.25){
							n = 'A';
						}else if(r < 0.5){
							n = 'C';
						}else if(r < 0.75){
							n = 'G';
						}else{
							n = 'T';
						}
					}
					sequenced_read_mutated[c] = n;
					c++;
				}
			}else{
				sequenced_read_mutated[c] = sequenced_read[j];
				c++;
			}
		}
		if(param->sim_sequenced_len < c){
			c = param->sim_sequenced_len;
		}
		sequenced_read_mutated[c]  =0;
		
		
		//sequenced_read_mutated[param->sim_sequenced_len] = 0;
		
		//fprintf(stderr,"%d	%s\n",i, sequenced_read_mutated);
		fprintf(file,"@READ%d;SEQ:%s;",i, read);
		if(param->sim_barnum){
			fprintf(file,"RBC:%s;",barcode[barcode_used]);
		
		}
		fprintf(file,"\n");
		fprintf(file,"%s\n+\n",sequenced_read_mutated);
		for(j = 0; j < c;j++){
			fprintf(file,"I");
		}
		fprintf(file,"\n");
		
	}
	//
	// this adds totally random sequences.. 
	//
	c = 0; // random sequence length
	if(param->sim_5seq){
		c += (int) strlen (param->sim_5seq);
	}
	if(param->sim_3seq){
		c += (int) strlen (param->sim_3seq);
	}
	
	c +=param->sim_barlen ;
	
	c+= param->sim_readlen;

	
	
	for(i = (int)((float) param->sim_numseq * (1.0-param->sim_random_frac)); i < param->sim_numseq   ;i++){
		if(param->sim_5seq){
			strcat ( sequenced_read,  param->sim_5seq);
		}
		for(j = 0; j < c;j++){
			r = (float)rand_r(&seed)/(float)RAND_MAX;
			if(r < 0.25){
				n = 'A';
			}else if(r < 0.5){
				n = 'C';
			}else if(r < 0.75){
				n = 'G';
			}else{
				n = 'T';
			}
			sequenced_read[j] = n;
		}
		if(param->sim_sequenced_len < c){
			c = param->sim_sequenced_len;
		}
		sequenced_read[c] = 0;
		//sequenced_read[param->sim_sequenced_len] = 0;
		//fprintf(stderr,"%s\n", sequenced_read);
		fprintf(file,"@RAND%d;SEQ:NONE;",i);
		if(param->sim_barnum){
			fprintf(file,"RBC:NONE;");
			
		}
		fprintf(file,"\n");
		fprintf(file,"%s\n+\n",sequenced_read);
		for(j = 0; j < c;j++){
			fprintf(file,"I");
		}
		fprintf(file,"\n");
		
		
	}
	fclose(file);

	
	
	//run tests:
	command = malloc(sizeof(char) * 10000);
	ri = malloc(sizeof(struct read_info*) * param->num_query);
	
	assert(ri !=0);
	
	for(i = 0; i < param->num_query;i++){
		ri[i] = malloc(sizeof(struct read_info));
		ri[i]->seq = 0;
		ri[i]->name = 0;
		ri[i]->qual = 0;
		ri[i]->labels = 0;
		ri[i]->len = 0;
		ri[i]->cigar = 0;
		ri[i]->bar_prob = 0;
		ri[i]->md = 0;
		ri[i]->strand = malloc(sizeof(unsigned int)* (LIST_STORE_SIZE+1));
		ri[i]->hits = malloc(sizeof(unsigned int)* (LIST_STORE_SIZE+1));
	}
	
	// btrim
	
	if(param->sim_3seq){
		sprintf (outfile, "%s/%sbtrimout.fq",dp,runid);
		if ((file2 = fopen(outfile, "w")) == NULL){
			fprintf(stderr,"can't open output\n");
			exit(-1);
		}
		
		
		if(param->sim_barnum){
			c = sprintf(command, "btrim -p %s/%sbtrim_pattern.txt -t %s/%sread.fq -o %s/%sbtrimout.fq -l %d -B", dp,runid,dp,runid,dp,runid, param->sim_readlen - param->sim_readlen_mod - 2);
			fprintf(stderr,"%s\n",command);
//btrim -p ~/tmp/testtag_btrim_pattern.txt  -t ~/tmp/testtag_read.fq -o  ~/tmp/btrimtest -s ~/tmp/btrimtestlog -l %d  -B
		}else{
			c = sprintf(command, "btrim -p %s/%sbtrim_pattern.txt -t %s/%sread.fq -o %s/%sbtrimout.fq -l %d ", dp,runid,dp,runid,dp,runid, param->sim_readlen - param->sim_readlen_mod - 2);
			fprintf(stderr,"%s\n",command);
		}
		//run btrim!!!
		if (!(file = popen(command, "r"))) {
			fprintf(stderr,"Cannot execute btrim with command:%s\n",command);
			exit(-1);
		}
		while ((c = fgetc (file)) != EOF){
			putchar (c);
		}
		pclose(file);
		
		//collect output in one file!
		if(param->sim_barnum){
			for(i = 0; i < param->sim_barnum;i++){
				c = sprintf(command,"cat %s/%sbtrimout.fq.%d  ", dp,runid,i);
				if (!(file = popen(command, "r"))) {
					fprintf(stderr,"Cannot execute command:%s\n",command);
					exit(-1);
				}
				while ((numseq = read_fasta_fastq(ri, param,file)) != 0){
					for(j = 0; j < numseq;j++){
						//fprintf(stderr,"%sBC:%s\n",ri[j]->name,   barcode[i]  );
						fprintf(file2,"@%s;BC:%s;\n",ri[j]->name,   barcode[i] );
						//for(i = 0; i < ri->len;i++){
						for(c = 0; c < ri[j]->len;c++){
							fprintf(file2,"%c", alpha[(int) ri[j]->seq[c]]);
						}
						fprintf(file2,"\n+\n%s\n" ,ri[j]->qual);
						
					}
				}
				
				pclose(file);
								
			}
		}
		fclose(file2);
	}
		
	//AdapterRemoval
	
	if(!param->sim_5seq && !param->sim_barnum && param->sim_3seq){
		c = sprintf(command, "cat %s/%sread.fq | AdapterRemoval --pcr1 %s  > %s/%sAdapterRemovalout.fq ",dp,runid,param->sim_3seq,dp,runid);
		fprintf(stderr,"%s\n",command);
		if (!(file = popen(command, "r"))) {
			fprintf(stderr,"Cannot execute cutadapt with command:%s\n",command);
			exit(-1);
		}
		while ((c = fgetc (file)) != EOF){
			putchar (c);
		}
		pclose(file);
	}
	
	
	
	//cutadapt!!!
	
	
	
	if(param->sim_5seq && param->sim_3seq){
		c = sprintf(command, "cutadapt --discard-untrimmed -n 2 --front %s    --adapter=%s  %s/%sread.fq   -o %s/%scutadaptout.fq ", param->sim_5seq,param->sim_3seq,dp,runid,dp,runid);
		fprintf(stderr,"%s\n",command);
		if (!(file = popen(command, "r"))) {
			fprintf(stderr,"Cannot execute cutadapt with command:%s\n",command);
			exit(-1);
		}
		while ((c = fgetc (file)) != EOF){
			putchar (c);
		}
		pclose(file);
		
		
	 //cutadapt -n 2 --front test=agggaggacgatgcgg   --adapter=gtgtcagtcacttccagcgg ~/tmp/testtag_read.fq    > ~/tmp/testtag_cutadapt.fq --discard-untrimmed
	}
	
	
	//fastx_barcode...
	
	//split reads which were cut by cutadapt.
	
	
	if(param->sim_5seq && param->sim_3seq && param->sim_barnum){
		sprintf (outfile, "%s/%scutadaptfastxout.fq",dp,runid);
		if ((file2 = fopen(outfile, "w")) == NULL){
			fprintf(stderr,"can't open output\n");
			exit(-1);
		}
		c = sprintf(command, "cat %s/%scutadaptout.fq | fastx_barcode_splitter.pl --bcfile %s/%sfastxbarcodefile.txt --prefix %s/%sfastx -bol", dp,runid,dp,runid,dp,runid);
		fprintf(stderr,"%s\n",command);
		if (!(file = popen(command, "r"))) {
			fprintf(stderr,"Cannot execute cutadapt with command:%s\n",command);
			exit(-1);
		}
		while ((c = fgetc (file)) != EOF){
			putchar (c);
		}
		pclose(file);
		for(i = 0; i < param->sim_barnum;i++){
			c = sprintf(command,"cat %s/%sfastxBC%d  ", dp,runid,i);
			if (!(file = popen(command, "r"))) {
				fprintf(stderr,"Cannot execute command:%s\n",command);
				exit(-1);
			}
			while ((numseq = read_fasta_fastq(ri, param,file)) != 0){
				for(j = 0; j < numseq;j++){
					//fprintf(stderr,"%sBC:%s\n",ri[j]->name,   barcode[i]  );
					fprintf(file2,"@%s;BC:%s;\n",ri[j]->name,   barcode[i] );
					
					//need to start from the position after barcode .... because fastX does not trim the barcode!
					
					for(c = param->sim_barlen ; c < ri[j]->len;c++){
						fprintf(file2,"%c", alpha[(int) ri[j]->seq[c]]);
					}
					fprintf(file2,"\n+\n%s\n" ,ri[j]->qual);
					
				}
			}
			
			
			pclose(file);
			
			
		}
		fclose(file2);
	}
	
		
	//split reads without adapters.
	
	
	if(!param->sim_5seq && !param->sim_3seq && param->sim_barnum){
		sprintf (outfile, "%s/%sfastxout.fq",dp,runid);
		if ((file2 = fopen(outfile, "w")) == NULL){
			fprintf(stderr,"can't open output\n");
			exit(-1);
		}
		c = sprintf(command, "cat %s/%sread.fq | fastx_barcode_splitter.pl --bcfile %s/%sfastxbarcodefile.txt --prefix %s/%sfastx -bol", dp,runid,dp,runid,dp,runid);
		fprintf(stderr,"%s\n",command);
		if (!(file = popen(command, "r"))) {
			fprintf(stderr,"Cannot execute cutadapt with command:%s\n",command);
			exit(-1);
		}
		while ((c = fgetc (file)) != EOF){
			putchar (c);
		}
		pclose(file);
		for(i = 0; i < param->sim_barnum;i++){
			c = sprintf(command,"cat %s/%sfastxBC%d  ", dp,runid,i);
			if (!(file = popen(command, "r"))) {
				fprintf(stderr,"Cannot execute command:%s\n",command);
				exit(-1);
			}
			while ((numseq = read_fasta_fastq(ri, param,file)) != 0){
				for(j = 0; j < numseq;j++){
					//fprintf(stderr,"%sBC:%s\n",ri[j]->name,   barcode[i]  );
					fprintf(file2,"@%s;BC:%s;\n",ri[j]->name,   barcode[i] );
					
					//need to start from the position after barcode .... because fastX does not trim the barcode!
					
					for(c = param->sim_barlen ; c < ri[j]->len;c++){
						fprintf(file2,"%c", alpha[(int) ri[j]->seq[c]]);
					}
					fprintf(file2,"\n+\n%s\n" ,ri[j]->qual);
					
				}
			}
			
			pclose(file);
		}
		fclose(file2);
	}
	
	
	
	
	//TAGDUST !!!!!!
	if(param->sim_5seq && param->sim_barnum && param->sim_3seq){
		c = sprintf(command,"tagdust -threshold %f -t 80 -1 P:%s -2 B:", param->confidence_threshold,param->sim_5seq);
		for(i = 0; i < param->sim_barnum-1;i++){
			c = sprintf (read, "%s,",barcode[i]);
			strcat ( command, read);
		}
		c = sprintf (read, "%s ",barcode[param->sim_barnum-1]);
		strcat ( command, read);
		c = sprintf (read, "-3 R:N ");
		strcat ( command, read);
		c = sprintf (read, "-4 P:%s ", param->sim_3seq);
		strcat ( command, read);
		c = sprintf (read, "%s/%sread.fq  -o %s/%stagdustout.fq ",dp,runid,dp,runid);
		strcat ( command, read);
		fprintf(stderr,"%s\n",command);
		if (!(file = popen(command, "r"))) {
			fprintf(stderr,"Cannot execute tagdust with command:%s\n",command);
			exit(-1);
		}
		while ((c = fgetc (file)) != EOF){
			putchar (c);
		}
		pclose(file);

	}
	
	if(!param->sim_5seq && param->sim_barnum && param->sim_3seq){
		c = sprintf(command,"tagdust -threshold %f -t 80  -1 B:", param->confidence_threshold);
		for(i = 0; i < param->sim_barnum-1;i++){
			c = sprintf (read, "%s,",barcode[i]);
			strcat ( command, read);
		}
		c = sprintf (read, "%s ",barcode[param->sim_barnum-1]);
		strcat ( command, read);
		c = sprintf (read, "-2 R:N ");
		strcat ( command, read);
		c = sprintf (read, "-3 P:%s ", param->sim_3seq);
		strcat ( command, read);
		c = sprintf (read, "%s/%sread.fq  -o %s/%stagdustout.fq ",dp,runid,dp,runid);
		strcat ( command, read);
		fprintf(stderr,"%s\n",command);
		if (!(file = popen(command, "r"))) {
			fprintf(stderr,"Cannot execute tagdust with command:%s\n",command);
			exit(-1);
		}
		while ((c = fgetc (file)) != EOF){
			putchar (c);
		}
		pclose(file);
		
	}
	
	if(param->sim_5seq && !param->sim_barnum && param->sim_3seq){
		c = sprintf(command,"tagdust -threshold %f -t 80 -1 P:%s ", param->confidence_threshold, param->sim_5seq);
		/*for(i = 0; i < param->sim_barnum-1;i++){
			c = sprintf (read, "%s,",barcode[i]);
			strcat ( command, read);
		}
		c = sprintf (read, "%s ",barcode[param->sim_barnum-1]);
		strcat ( command, read);
		*/
		c = sprintf (read, "-2 R:N ");
		strcat ( command, read);
		c = sprintf (read, "-3 P:%s ", param->sim_3seq);
		strcat ( command, read);
		c = sprintf (read, "%s/%sread.fq  -o %s/%stagdustout.fq ", dp,runid,dp,runid);
		strcat ( command, read);
		fprintf(stderr,"%s\n",command);
		if (!(file = popen(command, "r"))) {
			fprintf(stderr,"Cannot execute tagdust with command:%s\n",command);
			exit(-1);
		}
		while ((c = fgetc (file)) != EOF){
			putchar (c);
		}
		pclose(file);
		
	}
	
	if(param->sim_5seq && param->sim_barnum && !param->sim_3seq){
		c = sprintf(command,"tagdust -threshold %f -t 80 -1 P:%s -2 B:", param->confidence_threshold, param->sim_5seq);
		for(i = 0; i < param->sim_barnum-1;i++){
			c = sprintf (read, "%s,",barcode[i]);
			strcat ( command, read);
		}
		c = sprintf (read, "%s ",barcode[param->sim_barnum-1]);
		strcat ( command, read);
		c = sprintf (read, "-3 R:N ");
		strcat ( command, read);
		c = sprintf (read, "%s/%sread.fq  -o %s/%stagdustout.fq ", dp,runid,dp,runid);
		strcat ( command, read);
		fprintf(stderr,"%s\n",command);
		if (!(file = popen(command, "r"))) {
			fprintf(stderr,"Cannot execute tagdust with command:%s\n",command);
			exit(-1);
		}
		while ((c = fgetc (file)) != EOF){
			putchar (c);
		}
		pclose(file);
		
	}
	
	if(!param->sim_5seq && param->sim_barnum && !param->sim_3seq){
		c = sprintf(command,"tagdust -threshold %f -t 80  -1 B:", param->confidence_threshold);
		for(i = 0; i < param->sim_barnum-1;i++){
			c = sprintf (read, "%s,",barcode[i]);
			strcat ( command, read);
		}
		c = sprintf (read, "%s ",barcode[param->sim_barnum-1]);
		strcat ( command, read);
		c = sprintf (read, "-2 R:N ");
		strcat ( command, read);
		c = sprintf (read, "%s/%sread.fq  -o %s/%stagdustout.fq ", dp,runid,dp,runid);
		strcat ( command, read);
		fprintf(stderr,"%s\n",command);
	
		if (!(file = popen(command, "r"))) {
			fprintf(stderr,"Cannot execute tagdust with command:%s\n",command);
			exit(-1);
		}
		while ((c = fgetc (file)) != EOF){
			putchar (c);
		}
		pclose(file);
		
	}
	
	
	if(!param->sim_5seq && !param->sim_barnum && param->sim_3seq){
		c = sprintf(command,"tagdust -threshold %f -t 80 ", param->confidence_threshold);
		/*for(i = 0; i < param->sim_barnum-1;i++){
		 c = sprintf (read, "%s,",barcode[i]);
		 strcat ( command, read);
		 }
		 c = sprintf (read, "%s ",barcode[param->sim_barnum-1]);
		 strcat ( command, read);
		 */
		c = sprintf (read, "-1 R:N ");
		strcat ( command, read);
		c = sprintf (read, "-2 P:%s ", param->sim_3seq);
		strcat ( command, read);
		c = sprintf (read, "%s/%sread.fq  -o %s/%stagdustout.fq ", dp,runid,dp,runid);
		strcat ( command, read);
		fprintf(stderr,"%s\n",command);
		if (!(file = popen(command, "r"))) {
			fprintf(stderr,"Cannot execute tagdust with command:%s\n",command);
			exit(-1);
		}
		while ((c = fgetc (file)) != EOF){
			putchar (c);
		}
		pclose(file);
		
	}


	//EVALUATION!!!
	
	struct eval_results* evalres = 0;
	
	evalres = malloc(sizeof(struct eval_results));
	
	//sprintf (outfile, "%s/read.fq",param->outfile);
	sprintf (outfile, "%s/%sbtrimout.fq",dp,runid );
	if( access( outfile, F_OK ) != -1 ) {
		evalres=  get_results(evalres, param, outfile, "BTRIM");
		unlink(outfile);
	}
	//exit(0);
	sprintf (outfile, "%s/%sAdapterRemovalout.fq",dp,runid );
	if( access( outfile, F_OK ) != -1 ) {
		evalres=  get_results(evalres, param, outfile, "AdapterRemoval");
		unlink(outfile);
	}
	
		
	
	sprintf (outfile, "%s/%stagdustout.fq",dp,runid);
	if( access( outfile, F_OK ) != -1 ) {
		evalres=  get_results(evalres, param, outfile, "TagDust");
		unlink(outfile);
	}
	
	sprintf (outfile, "%s/%sfastxout.fq",dp,runid);
	if( access( outfile, F_OK ) != -1 ) {
		evalres=  get_results(evalres, param, outfile, "FASTX");
		unlink(outfile);
	}
	
	sprintf (outfile, "%s/%scutadaptfastxout.fq",dp,runid);
	if( access( outfile, F_OK ) != -1 ) {
		evalres=  get_results(evalres, param, outfile, "CUTADAPTFASTX");
		unlink(outfile);
	}
	
	sprintf (outfile, "%s/%scutadaptout.fq",dp,runid);
	if( access( outfile, F_OK ) != -1 ) {
		evalres=  get_results(evalres, param, outfile, "CUTADAPT");
		unlink(outfile);
	}
	//exit(0);

	
	for(i = 0; i < param->sim_barnum;i++){
		sprintf (outfile, "%s/%sbtrimout.fq.%d", dp,runid,i);
		
		//fprintf(stderr,"DELETING: %s\n",outfile);
		if( access( outfile, F_OK ) != -1 ) {
			fprintf(stderr,"DELETING: %s\n",outfile);
			c = unlink(outfile);
		//	fprintf(stderr,"%d\n",c);
		}
		sprintf (outfile, "%s/%sfastxBC%d",dp,runid,i);
		
		//fprintf(stderr,"DELETING: %s\n",outfile);
		if( access( outfile, F_OK ) != -1 ) {
			fprintf(stderr,"DELETING: %s\n",outfile);
			c = unlink(outfile);
		//	fprintf(stderr,"%d\n",c);
		}
		
		
	
	}
	sprintf (outfile, "%s/%sfastxunmatched",dp,runid);
	if( access( outfile, F_OK ) != -1 ) {
		fprintf(stderr,"DELETING: %s\n",outfile);
		c = unlink(outfile);
		//fprintf(stderr,"%d\n",c);
	}
	//54F01FF5fastxunmatched
	
	sprintf (outfile, "%s/%sread.fq",dp,runid);
	if( access( outfile, F_OK ) != -1 ) {
		c = unlink(outfile);
	}
	
	sprintf (outfile, "%s/%sfastxbarcodefile.txt",dp,runid);
	if( access( outfile, F_OK ) != -1 ) {
		c = unlink(outfile);
	}
	
	sprintf (outfile, "%s/%sbtrim_pattern.txt",dp,runid);
	if( access( outfile, F_OK ) != -1 ) {
		c = unlink(outfile);
	}

	rmdir(basedir);

	free(sequenced_read_mutated);
	free(sequenced_read);
	free(read);
	free(runid);
	
	free(command);
	for(i = 0; i < param->num_query;i++){
		free(ri[i]->strand);
		free(ri[i]->hits);
		
		if(ri[i]->cigar){
			free(ri[i]->cigar);
		}
		if(ri[i]->md){
			free(ri[i]->md);
		}
		if(ri[i]->name){
			free(ri[i]->name);
		}
		if(ri[i]->seq){
			free(ri[i]->seq);
		}
		if(ri[i]->qual){
			free(ri[i]->qual );
		}
		if(ri[i]->labels){
			free(ri[i]->labels);
		}
		
		free(ri[i]);
	}
	
	free(ri);
	free(basedir);
	if(param->sim_barnum){
		for(i = 0 ;i < param->sim_barnum;i++){
			
			free(barcode[i]);
		}
		free(barcode);
	}
	
	free_param(param);
	
	exit(EXIT_SUCCESS);
}

struct eval_results* get_results(struct eval_results* eval,struct parameters* param, char* filename ,char* program)
{
	FILE* file;
	struct read_info** ri  = 0;
	int numseq;
	int i,j,c,tmp,org_read_len;
	
	//char alpha[] = "ACTGN";
	char* orgread = 0;
	char* bc = 0;
	char* orgbc = 0;
	org_read_len = 0;
	
	double TP,FP,FN, TN,sensitivity,specificity,precision,kappa,P_e,P_o;
	
	
	orgread = malloc(sizeof(char)* 300);
	orgbc = malloc(sizeof(char)* 300);
	bc = malloc(sizeof(char)* 300);
	
	eval->num_extracted = 0;
	eval->average_read_similarity= 0.0;
	eval->num_rand_extracted = 0;
	eval->num_wrong_bc = 0;
	if ((file = fopen(filename, "r")) == NULL){
		sprintf(param->buffer,"can't open output\n");
		fprintf(stderr,"%s",param->buffer);
		param->messages = append_message(param->messages, param->buffer);
		free_param(param);
		exit(EXIT_FAILURE);
	}
	
	ri = malloc(sizeof(struct read_info*) * param->num_query);
	
	assert(ri !=0);
	
	for(i = 0; i < param->num_query;i++){
		ri[i] = malloc(sizeof(struct read_info));
		ri[i]->seq = 0;
		ri[i]->name = 0;
		ri[i]->qual = 0;
		ri[i]->labels = 0;
		ri[i]->len = 0;
		ri[i]->cigar = 0;
		ri[i]->bar_prob = 0;
		ri[i]->md = 0;
		ri[i]->strand = malloc(sizeof(unsigned int)* (LIST_STORE_SIZE+1));
		ri[i]->hits = malloc(sizeof(unsigned int)* (LIST_STORE_SIZE+1));
	}

	
	
	while ((numseq = read_fasta_fastq(ri, param,file)) != 0){
		for(i = 0; i < numseq;i++){
			//fprintf(stderr,"@%s\n",ri[i]->name);
			tmp = byg_end("SEQ:", ri[i]->name  );
			if(tmp){
				org_read_len = 0;
				for(j = tmp ;j < (int)strlen(ri[i]->name);j++){
					org_read_len++;
					if(isspace((int) ri[i]->name[j]) ||  ri[i]->name[j]  == ';'){
						break;
					}
					
				}
				//orgread = malloc(sizeof(unsigned char)* g);
				org_read_len = 0;
				for(j = tmp ;j < (int)strlen(ri[i]->name);j++){
					
					if(isspace((int) ri[i]->name[j]) ||  ri[i]->name[j]  == ';'){
						orgread[org_read_len] = 0;
						break;
					}
					
					orgread[org_read_len] = nuc_code[(int)ri[i]->name[j]];
					org_read_len++;
				}
			}
			
			tmp = byg_end("RBC:", ri[i]->name  );
			if(tmp){
				c = 0;
				for(j = tmp ;j < (int)strlen(ri[i]->name);j++){
					if(ri[i]->name[j]  == ';' ){
						orgbc[c] = 0;
						break;
					}
					orgbc[c] = ri[i]->name[j];
					c++;
				}
			}else{
				orgbc[0] = 0;
			}

			
			
			tmp = byg_end(";BC:", ri[i]->name  );
			
			if(tmp){
				
				c = 0;
				for(j = tmp ;j < (int)strlen(ri[i]->name);j++){
					if(ri[i]->name[j]  == ';' ){
						bc[c] = 0;
						break;
					}
					bc[c] = ri[i]->name[j];
					c++;
				}
			}else{
				bc[0] = 0;
			}
			
			c =  byg_count("READ", ri[i]->name);
			
			if(c){
				eval->num_extracted++;
				

				if(strcmp(orgbc, bc)){
					eval->num_wrong_bc++;
				}
				if(ri[i]->len < org_read_len){
					c = bpm_check_error_global((unsigned char*)ri[i]->seq,(unsigned char*) orgread, ri[i]->len, org_read_len);
				}else{
					c = bpm_check_error_global((unsigned char*) orgread,(unsigned char*)ri[i]->seq,org_read_len, ri[i]->len);
				}
				/*fprintf(stderr,"%s\t%d\n", program,c);
				for(j = 0;  j < org_read_len;j++){
					fprintf(stderr,"%c",  alpha[orgread[j]]);
				}
				fprintf(stderr,"\n");
				
				for(j = 0;  j < ri[i]->len;j++){
					fprintf(stderr,"%c",  alpha[ri[i]->seq[j]]);
				}
				fprintf(stderr,"\n");
				
				*/
				j = (org_read_len > ri[i]->len) ?org_read_len : ri[i]->len ;
				eval->average_read_similarity += (double) c / (double)j ;
			}
			
			c =  byg_count("RAND", ri[i]->name);
			
			if(c){
				//eval->num_extracted++;
				eval->num_rand_extracted++;
			}
			
			
			
			
			//fprintf(stderr,"%sBC:%s\n",ri[j]->name,   barcode[i]  );
			//if(
						
			
			
			//need to start from the position after barcode .... because fastX does not trim the barcode!
			
			/**/
			
			//fprintf(stderr,"\n+\n%s\n" ,ri[j]->qual);
			
		}
	}
	
	TP = eval->num_extracted - eval->num_wrong_bc;
	FP = eval->num_wrong_bc + eval->num_rand_extracted;
	FN = (int)((float) param->sim_numseq * (1.0-param->sim_random_frac)) - eval->num_extracted;
	TN =  param->sim_numseq  - (TP + FP + FN);

	
	precision = TP / (TP + FP);
	sensitivity = TP/( TP + FN );
	specificity =  TN / ( TN + FP);
	
	P_e = ((TP+FN) / (double)param->sim_numseq) * ((TP+FP) / (double)param->sim_numseq) +  ( ((FP+TN) / (double)param->sim_numseq  ) * ((FN+TN) / (double)param->sim_numseq));
	P_o =(TP+TN)/(double)param->sim_numseq ;
	
	kappa = (P_o - P_e) / (1.0 - P_e);
	
	fprintf(stderr,"%s	Sen:%f	Spe:%f	Precision:%f	Kappa:%f	TP:%f	FP:%f	FN:%f	TN:%f\n",program,sensitivity,specificity,precision,kappa,TP,FP,FN,TN);
	sprintf (orgread, "%s/summary.csv",param->outfile);
	
	if ((file = fopen(orgread, "a")) == NULL){
		sprintf(param->buffer,"can't open output\n");
		fprintf(stderr,"%s",param->buffer);
		param->messages = append_message(param->messages, param->buffer);
		free_param(param);
		exit(EXIT_FAILURE);
	}
	struct tm *ptr;
	int hour;
	char am_or_pm;
	//char logfile[100];
	
	time_t current = time(NULL);
	ptr = localtime(&current);
	hour = ptr->tm_hour;
	if (hour <= 11)
		am_or_pm = 'a';
	else {
		hour -= 12;
		am_or_pm = 'p';
	}
	if (hour == 0){
		hour = 12;
	}
	
	if(!param->sim_5seq){
		param->sim_5seq = "NA";
	}
	
	if(!param->sim_3seq){
		param->sim_3seq = "NA";
	}
	if(!(eval->num_extracted +  eval->num_rand_extracted)){
		fprintf(file,"%.2d-%.2d-%d;%d:%d%cm\t%s\t%f\t%f\t%f\t%f\t%0.0f\t%0.0f\t%0.0f\t%0.0f\t%f\t%d\t%f\t%s\t%s\t%d\t%d\t%d\t%d\t%f\t%f\t%d\t%f\n",
		        ptr->tm_mon + 1,ptr->tm_mday, ptr->tm_year + 1900,hour,ptr->tm_min, am_or_pm,
		        program,sensitivity,specificity,precision,kappa,TP,FP,FN,TN,
		        -1.0,
		        param->sim_numseq,
		        param->sim_random_frac,
		        param->sim_5seq,
		        param->sim_3seq,
		        param->sim_barnum,
		        param->sim_barlen,
		        param->sim_readlen,
		        param->sim_readlen_mod,
		        param->sim_error_rate,
		        param->sim_InDel_frac,
		        param->sim_sequenced_len,
		        -1.0
		        );
	}else{
	
	fprintf(file,"%.2d-%.2d-%d;%d:%d%cm\t%s\t%f\t%f\t%f\t%f\t%0.0f\t%0.0f\t%0.0f\t%0.0f\t%f\t%d\t%f\t%s\t%s\t%d\t%d\t%d\t%d\t%f\t%f\t%d\t%f\n",
	        ptr->tm_mon + 1,ptr->tm_mday, ptr->tm_year + 1900,hour,ptr->tm_min, am_or_pm,
	        program,sensitivity,specificity,precision,kappa,TP,FP,FN,TN,
	        eval->average_read_similarity / (double)(eval->num_extracted +  eval->num_rand_extracted),
	        param->sim_numseq,
	        param->sim_random_frac,
	        param->sim_5seq,
	        param->sim_3seq,
	        param->sim_barnum,
	        param->sim_barlen,
	        param->sim_readlen,
	        param->sim_readlen_mod,
	        param->sim_error_rate,
	        param->sim_InDel_frac,
	        param->sim_sequenced_len,
	        (double)(eval->num_rand_extracted+eval->num_wrong_bc) / (double)(eval->num_extracted+ eval->num_rand_extracted)
	        
	        );
	}
	
	
	
	//fprintf(stderr,"%d extracted\n",eval->num_extracted);
	//fprintf(stderr,"%d rand extracted\n",eval->num_rand_extracted);
	//fprintf(stderr,"%d wrong bc\n",eval->num_wrong_bc);
	//fprintf(stderr,"%f(%f) error\n",eval->average_read_similarity / (double)(eval->num_extracted +  eval->num_rand_extracted), eval->average_read_similarity);
	fclose(file);
	free(orgread);
	free(orgbc);
	free(bc);
	for(i = 0; i < param->num_query;i++){
		free(ri[i]->strand);
		free(ri[i]->hits);
		
		if(ri[i]->cigar){
			free(ri[i]->cigar);
		}
		if(ri[i]->md){
			free(ri[i]->md);
		}
		if(ri[i]->name){
			free(ri[i]->name);
		}
		if(ri[i]->seq){
			free(ri[i]->seq);
		}
		if(ri[i]->qual){
			free(ri[i]->qual );
		}
		if(ri[i]->labels){
			free(ri[i]->labels);
		}
		
		free(ri[i]);
	}
	
	free(ri);

	return eval;
}


/** \fn void simulate(struct parameters* param)
 \brief Prints out simulates sequences.
 
 Prints out simulated sequences containing a barcode at the 5' end. Parameters are the length of the barcode, number of barcodes and error rates. The function attempts to find a set of bacodes with a maximum pairwise edit distance. All parameters are passed to the function via the @ref parameters struct. In total 900000 sequence are simulated and 100000 random sequences are added. 
 
 In addition this function prints out a shell script to run Tagdust2 and fastx_barcode_splitter on the simulated data. 
 
 \param param @parameters.
 
 
 */
void simulate(struct parameters* param)
{
	int i,j,c,n,n_dash,read_pos, mismatches, indel;
	float r = 0.0;
	FILE* file = 0;
	int errors_allowed = 6;
	char alpha[6] = "ACGTNN";
	
	char* outfile = 0;
	outfile = malloc(sizeof(char)*(200));
	assert(outfile != 0);
	

	//sprintf (outfile, "%s_tagdust_command.sh",param->outfile);
	
	
	//char seq_a[1000];
	
	char** barcode = 0;
	int num_barcode = param->numbarcode;
	if(param->sim < 16){
		if( (int) powf(4.0, (float)param->sim) <= num_barcode){
			num_barcode =  (int) powf(4.0, (float)param->sim) ;
		}
		//fprintf(stderr,"%f\n",  powf(4.0, (float)param->sim));
	}
	//int num_barcode = (int) powf(4.0, (float)param->sim) ;
	
	
	fprintf(stderr,"Number of Barcodes:%d\n",num_barcode );
	barcode = malloc(sizeof(char*)* num_barcode);
	assert(barcode!=0);
	for(i = 0; i < num_barcode;i++){
		barcode[i] = malloc(sizeof(char) *(param->sim+1) );
		assert(barcode[i]!=0);
	}
	
	unsigned int seed = (unsigned int) (time(NULL) * ( 42));
	
	
	
	
	n_dash = 100;
	c = 0;
	
	while(n_dash){
		n_dash = 0;
		//fprintf(stderr,"Errors allowed: %d\n");
		while(c != num_barcode){
			for(i = 0; i < param->sim;i++){
				r = (float)rand_r(&seed)/(float)RAND_MAX;
				if(r < 0.25){
					barcode[c][i] = 0;
				}else if(r < 0.5){
					barcode[c][i] = 1;
				}else if(r < 0.75){
					barcode[c][i] = 2;
				}else{
					barcode[c][i] = 3;
				}
			}
			j = 1;
			for(i = 0; i < c;i++){
				if(bpm(barcode[c], barcode[i], param->sim, param->sim)  <= errors_allowed){
					j = 0;
					break;
				}
			}
			if(j){
				c++;
			}
			//
			n_dash++;
			if(n_dash == 1000000){
				break;
			}
		}
		
		if(n_dash == 1000000){
			errors_allowed--;
		}
		
	}
	//fprintf(stderr,"Edit Distance between barcodes: %d\n", errors_allowed);

	if(errors_allowed == 0){// two identical barcodes were found.... I will simply pick the fist numbarcode sequences... 
		for(i = 0; i < num_barcode ;i++){
			for(j = 0; j < param->sim;j++){
				barcode[i][j] = 0x3 & (i >> (j*2)) ;
			}
		}
		
		errors_allowed = 100;
		for(i = 0; i < num_barcode-1 ;i++){

			for(j= i+1; j < num_barcode ;j++){
				c = bpm(barcode[i], barcode[j], param->sim, param->sim);
				if(c < errors_allowed){
					errors_allowed = c;
				}
			}
		}

		
	}
	fprintf(stderr,"Edit Distance between barcodes: %d\n", errors_allowed);

	
	
	sprintf (outfile, "%s_tagdust_command.sh",param->outfile);
	
	
	if ((file = fopen(outfile, "w")) == NULL){
		fprintf(stderr,"can't open output\n");
		exit(-1);
	}
	
	sprintf (outfile, "%s_read.fq",param->outfile);
	
	fprintf(file,"tagdust %s -t 80 -1 B:",outfile);
	
	
	for(i = 0 ;i < num_barcode-1;i++){
		//fprintf(stderr,"%d\t",i);
		for(j = 0; j < param->sim;j++){
			fprintf(file,"%c",alpha[(int)barcode[i][j]]);
		}
		fprintf(file,",");
		
	}
	for(j = 0; j < param->sim;j++){
		fprintf(file,"%c",alpha[(int)barcode[num_barcode-1][j]]);
	}
	sprintf (outfile, "%s_read_extracted20.fq",param->outfile);
	
	fprintf(file," -2 R:N -q 20 -o %s  \n", outfile);
	
	fprintf(file,"grep READ  %s  |  awk -v numbarcode=%d -v sim=%d -v errorrate=%f  -v indelrate=%f 'BEGIN{p=0;n=0}{x = split($0,a,\"[;,:]\");if(a[x] == a[3]){p++}else{n++}}END{printf \"%%d\\t%%d\\t%%f\\t%%f\\t%%d\\t%%d\\t%%d\\n\",numbarcode,sim,errorrate,indelrate,p,1000000 - (p + n),n}'  >> tagdust_benchmark20.csv &\n", outfile,param->numbarcode, param->sim,param->sequencer_error_rate, param->indel_frequency );

	
	
	sprintf (outfile, "%s_read.fq",param->outfile);
	
	fprintf(file,"tagdust %s -t 80 -1 B:",outfile);
	
	
	for(i = 0 ;i < num_barcode-1;i++){
		//fprintf(stderr,"%d\t",i);
		for(j = 0; j < param->sim;j++){
			fprintf(file,"%c",alpha[(int)barcode[i][j]]);
		}
		fprintf(file,",");
		
	}
	for(j = 0; j < param->sim;j++){
		fprintf(file,"%c",alpha[(int)barcode[num_barcode-1][j]]);
	}
	sprintf (outfile, "%s_read_extracted10.fq",param->outfile);
	
	fprintf(file," -2 R:N -q 10 -o %s  \n", outfile);
	

	
	
	
	//fprintf(file,"echo -n %d\t%d\t%f\t;\t",param->numbarcode, param->sim,param->sequencer_error_rate);
	fprintf(file,"grep READ  %s  |  awk -v numbarcode=%d -v sim=%d -v errorrate=%f  -v indelrate=%f 'BEGIN{p=0;n=0}{x = split($0,a,\"[;,:]\");if(a[x] == a[3]){p++}else{n++}}END{printf \"%%d\\t%%d\\t%%f\\t%%f\\t%%d\\t%%d\\t%%d\\n\",numbarcode,sim,errorrate,indelrate,p,1000000 - (p + n),n}'  >> tagdust_benchmark10.csv &\n", outfile,param->numbarcode, param->sim,param->sequencer_error_rate, param->indel_frequency );
	
	
	
	//sprintf (outfile, "%s_read.fq",param->outfile);
	//fprintf(file,"cat %s | ", outfile);
	
	
	//sprintf (outfile, "%s_barcodefile.txt",param->outfile);
	
	//fprintf(file,"./fastx_barcode_splitter.pl  --bcfile %s --bol --mismatches 2 --prefix ~/tmp/bla2_ --suffix \".txt\"  | awk -v numbarcode=%d -v sim=%d -v errorrate=%f -v indelrate=%f  '{printf\"%%d\\t%%d\\t%%f\\t%%f\\t%%d\\t%%d\\t%%d\\n\" ,numbarcode,sim,errorrate,indelrate,$1,$2,$3 }' >> fastx_benchmark_2.csv &\n",outfile,param->numbarcode, param->sim,param->sequencer_error_rate, param->indel_frequency  );
	
	sprintf (outfile, "%s_read.fq",param->outfile);
	fprintf(file,"cat %s | ", outfile);
	
	
	sprintf (outfile, "%s_barcodefile.txt",param->outfile);

	
	fprintf(file,"./fastx_barcode_splitter.pl  --bcfile %s --bol --mismatches 1 --prefix ~/tmp/bla1_ --suffix \".txt\"  | awk -v numbarcode=%d -v sim=%d -v errorrate=%f -v indelrate=%f  '{printf\"%%d\\t%%d\\t%%f\\t%%f\\t%%d\\t%%d\\t%%d\\n\" ,numbarcode,sim,errorrate,indelrate,$1,$2,$3 }' >> fastx_benchmark_1.csv &\n",outfile,param->numbarcode, param->sim,param->sequencer_error_rate,param->indel_frequency  );
	
	sprintf (outfile, "%s_read.fq",param->outfile);
	fprintf(file,"cat %s | ", outfile);
	
	
	sprintf (outfile, "%s_barcodefile.txt",param->outfile);

	
	fprintf(file,"./fastx_barcode_splitter.pl  --bcfile %s --bol --mismatches 0 --prefix ~/tmp/bla0_ --suffix \".txt\"  | awk -v numbarcode=%d -v sim=%d -v errorrate=%f -v indelrate=%f  '{printf\"%%d\\t%%d\\t%%f\\t%%f\\t%%d\\t%%d\\t%%d\\n\" ,numbarcode,sim,errorrate,indelrate,$1,$2,$3 }' >> fastx_benchmark_0.csv \n",outfile,param->numbarcode, param->sim,param->sequencer_error_rate,  param->indel_frequency );
	
	
	fclose(file);
	
	sprintf (outfile, "%s_barcodefile.txt",param->outfile);
	if ((file = fopen(outfile, "w")) == NULL){
		fprintf(stderr,"can't open output\n");
		exit(-1);
	}
	
	for(i = 0 ;i < num_barcode-1;i++){
		fprintf(file,"BC%d ",i);
		for(j = 0; j < param->sim;j++){
			fprintf(file,"%c",alpha[(int)barcode[i][j]]);
		}
		fprintf(file,"\n");
		
	}
	fprintf(file,"BC%d ",num_barcode-1);

	for(j = 0; j < param->sim;j++){
		fprintf(file,"%c",alpha[(int)barcode[num_barcode-1][j]]);
	}
	fprintf(file,"\n");


	//free(barcode);
	fclose(file);
	
	char* read = malloc(sizeof(char)* (param->sim *2 + 20));
	char* read_mutated = malloc(sizeof(char)* (param->sim *2 + 20));

	assert(read!=0);
	
	int num_bar_correct = 0;
	int num_bar_misassigned = 0;
	int num_bar_missed = 0;
	int no_bar = 0;
	
	//int num_finger_correct = 0;
	//int num_finger_misassigned = 0;
	//int num_finger_missed = 0;
	
	sprintf (outfile, "%s_read.fq",param->outfile);
	if ((file = fopen(outfile, "w")) == NULL){
		fprintf(stderr,"can't open output\n");
		exit(-1);
	}
	
	
	for(i = 0; i < 1000000*0.1;i++){
		read_pos = 0;
		for(j = 0; j < 20+param->sim;j++){
			r = (float)rand_r(&seed)/(float)RAND_MAX;
			if(r < 0.25){
				n = 0;
			}else if(r < 0.5){
				n = 1;
			}else if(r < 0.75){
				n = 2;
			}else{
				n = 3;
			}
			read_mutated[read_pos] = n;
			read_pos++;
			//fprintf(file,"%c",alpha[(int) n ]);// = barcode[c][j];
		}
		fprintf(file,"@READ%d;RBC:NOBAR\n",i);
		for(j = 0;j < read_pos;j++){
			fprintf(file,"%c",alpha[(int) read_mutated[j]]);// = barcode[c][j];
		}
		fprintf(file,"\n");
		
		fprintf(file,"+\n");
		for(j = 0; j < read_pos;j++){
			fprintf(file,"I");
		}
		fprintf(file,"\n");
		no_bar++;
		
	}
	
	for(i = 1000000*0.1; i < 1000000;i++){
		//construct fake read
		c = (int) (rand_r(&seed) % (int) (num_barcode)) ;
		
		for(j = 0;j < param->sim;j++){
			read[j] = barcode[c][j];
		}
		
		fprintf(file,"@READ%d;RBC:",i);
		for(j = 0;j < param->sim;j++){
			fprintf(file,"%c",alpha[(int) read[j]]);//
		}
		
		
		
		// apply errors
		
		read_pos = 0;
		mismatches =0 ;
		indel = 0;
		
		for(j = 0;j < param->sim;j++){
			r = (float)rand_r(&seed)/(float)RAND_MAX;
			if(r <= param->sequencer_error_rate){
				//we have an error
				
				
				r = (float)rand_r(&seed)/(float)RAND_MAX;
				if(r <= param->indel_frequency){
					indel++;
					// we have an indel (only considering single nucleotide.....
					r = (float)rand_r(&seed)/(float)RAND_MAX;
					if(r <= 0.5){
						//insertion
						//n_dash = read[j];
						r = (float)rand_r(&seed)/(float)RAND_MAX;
						if(r < 0.25){
							n = 0;
						}else if(r < 0.5){
							n = 1;
						}else if(r < 0.75){
							n = 2;
						}else{
							n = 3;
						}
						read_mutated[read_pos] = read[j];
						read_pos++;
						read_mutated[read_pos] = n;
						read_pos++;
						
						
						
						
					}else{
						//deletion
					}
					
				}else{
					mismatches++;
					n = read[j];
				
					while(n == read[j]){
						r = (float)rand_r(&seed)/(float)RAND_MAX;
						if(r < 0.25){
							n = 0;
						}else if(r < 0.5){
							n = 1;
						}else if(r < 0.75){
							n = 2;
						}else{
							n = 3;
						}
					}
					read_mutated[read_pos] = n;
					read_pos++;
				}
			}else{
				read_mutated[read_pos] = read[j];
				read_pos++;
			}
		}
		
		/*if(bpm(barcode[c], read_mutated, param->sim, read_pos)  != 0){
		fprintf(stderr,"%d ",i);
			for(j = 0;j < param->sim;j++){
				fprintf(stderr,"%d", (int) barcode[c][j]);
			}
			fprintf(stderr,"\t");
			
			for(j = 0;j < read_pos;j++){
				fprintf(stderr,"%d", (int) read_mutated[j]);
			}
			fprintf(stderr," Mis:%d Indel:%d\n", mismatches,indel  );
			
			
		}*/
		for(j = 0; j < 20;j++){
			r = (float)rand_r(&seed)/(float)RAND_MAX;
			if(r < 0.25){
				n = 0;
			}else if(r < 0.5){
				n = 1;
			}else if(r < 0.75){
				n = 2;
			}else{
				n = 3;
			}
			read_mutated[read_pos] = n;
			read_pos++;
			//fprintf(file,"%c",alpha[(int) n ]);// = barcode[c][j];
		}
		//fprintf(file,"\n");
		
		
		// test if read is recognizeable...
		n = -1;
		for(j = 0;j < num_barcode;j++){
			if(bpm(barcode[j], read_mutated, param->sim,  param->sim)  == 0){
				if(j == c){
					n = 1;
				}else{
					n = 0;
				}
				break;
			}
		}
		switch(n){
			case -1:
				fprintf(file,";missedBC\n");
				num_bar_missed++;
				break;
			case 0:
				fprintf(file,";wrongBC\n");
				num_bar_misassigned++;
				break;
			case 1:
				fprintf(file,";correctBC\n");
				num_bar_correct++;
				break;
		}
		
		
		for(j = 0;j < read_pos;j++){
			fprintf(file,"%c",alpha[(int) read_mutated[j]]);// = barcode[c][j];
		}
		/*for(j = 0; j < 20;j++){
			r = (float)rand_r(&seed)/(float)RAND_MAX;
			if(r < 0.25){
				n = 0;
			}else if(r < 0.5){
				n = 1;
			}else if(r < 0.75){
				n = 2;
			}else{
				n = 3;
			}
			fprintf(file,"%c",alpha[(int) n ]);// = barcode[c][j];
		}*/
		fprintf(file,"\n");
		
		fprintf(file,"+\n");
		for(j = 0; j < read_pos;j++){
			fprintf(file,"I");
		}
		fprintf(file,"\n");
		
		
		
		//c = (int) (rand_r(&seed) % (int) (num_barcode)) ;
		/*c = 0;
		for(j = 0;j < param->sim;j++){
			r = (float)rand_r(&seed)/(float)RAND_MAX;
			if(r < 0.25){
				n = 0;
			}else if(r < 0.5){
				n = 1;
			}else if(r < 0.75){
				n = 2;
			}else{
				n = 3;
			}
			read[j] = n;
			c = (c << 2) | n;
		}
		
		// apply errors
		for(j = 0;j < param->sim;j++){
			r = (float)rand_r(&seed)/(float)RAND_MAX;
			if(r <= param->sequencer_error_rate){
				n = read[j];
				while(n == read[j]){
					r = (float)rand_r(&seed)/(float)RAND_MAX;
					if(r < 0.25){
						n = 0;
					}else if(r < 0.5){
						n = 1;
					}else if(r < 0.75){
						n = 2;
					}else{
						n = 3;
					}
				}
				read[j] = n;
			}
		}
		
		n_dash = 0;
		
		for(j = 0;j < param->sim;j++){
			n_dash = (n_dash << 2) | read[j];
		}

		if(c != n_dash){
			num_finger_misassigned++;
		}else{
			num_finger_correct++;
		}*/
	}
	fclose(file);
	
	
	sprintf (outfile, "simulation.log");
	if ((file = fopen(outfile, "a")) == NULL){
		fprintf(stderr,"can't open output\n");
		exit(-1);
	}
	fprintf(file,"%d\t%d\t%f\t%d\t%d\t%d\t%d\n",param->numbarcode, param->sim,param->sequencer_error_rate,num_bar_correct,num_bar_missed,num_bar_misassigned,no_bar);
	
	fclose(file);
	fprintf(stderr,"BARCODE:\n");
	fprintf(stderr,"%d(%0.1f) correct\n", num_bar_correct,(float)num_bar_correct / 1000000.0* 100.0);
	fprintf(stderr,"%d(%0.1f) num_bar_misassigned\n", num_bar_misassigned,(float)num_bar_misassigned / 1000000.0* 100.0 );
	fprintf(stderr,"%d(%0.1f) num_bar_missed\n", num_bar_missed,(float)num_bar_missed / 1000000.0 * 100.0);
	fprintf(stderr,"%d(%0.1f) no_bar\n", no_bar,(float)no_bar / 1000000.0 * 100.0);
	
	//fprintf(stderr,"Fingerprint:\n");
	//fprintf(stderr,"%d(%0.1f) correct\n", num_finger_correct,(float)num_finger_correct / 1000000.0* 100.0);
	//fprintf(stderr,"%d(%0.1f) num_bar_misassigned\n", num_finger_misassigned,(float)num_finger_misassigned / 1000000.0* 100.0 );
	
	free(read_mutated);
	free(read);
	
	for(i = 0 ;i < num_barcode;i++){
		
		free(barcode[i]);
	}
	free(barcode);
	
	free(outfile);

}


