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


#include <ctype.h>

#include "interface.h"
#include "nuc_code.h"
#include "misc.h"
#include "io.h"
#include "tagdust2.h"

FILE* io_handler(FILE* file, int file_num,struct parameters* param)
{
	char command[1000];
	char  tmp[1000];
	int i = 0; 
	int gzcat = -1;
	if(access("/usr/bin/gzcat", X_OK) == 0){
		gzcat = 1;
	}else if(access("/bin/gzcat", X_OK) == 0){
		gzcat = 1;
	}else if(access("/usr/bin/zcat", X_OK) == 0){
		gzcat = 0;
	}else if(access("/bin/zcat", X_OK) == 0){
		gzcat = 0;
	}
	
	if(param->gzipped && gzcat == -1){
		fprintf(stderr,"Cannot find gzcat / zcat on your system. Try gzcat <infile> | samstat -f sam/bam/fa/fq\n");
		exit(-1);
	}
	//fprintf(stderr,"%d	%d	%s\n",sam,file_num,param->infile[file_num]);
	if(file_num == -1){
		if(param->sam == 2){
			command[0] = 0;
			if(!param->filter){
				strcat ( command, "samtools view -F 768 "); 
			}else{
				strcat ( command, "samtools view -F "); 
				i = sprintf (tmp, "%s ",param->filter);
				strcat ( command, tmp);
			}
			i = sprintf (tmp, "%s ","-");
			strcat ( command, tmp);
			if (!(file = popen(command, "r"))) {
				fprintf(stderr,"Cannot open bam file '%s' with command:%s\n",param->infile[file_num],command);
				exit(-1);
			}
		}else if(param->sam == 1){
			command[0] = 0;
			if(!param->filter){
				strcat ( command, "samtools view -SF 768 "); 
			}else{
				strcat ( command, "samtools view -SF "); 
				i = sprintf (tmp, "%s ",param->filter);
				strcat ( command, tmp);
			}
			
			//strcat ( command, "samtools view -SF 768 ");     // Copy name into full name
			i = sprintf (tmp, "%s ", "-");
			strcat ( command, tmp);
			if (!(file = popen(command, "r"))) {
				fprintf(stderr,"Cannot open bam file '%s' with command:%s\n",param->infile[file_num],command);
				exit(-1);
			}
		}else{
			file = stdin;
		}
	}else{
		if(param->sam == 2){
			command[0] = 0;
			if(param->gzipped){
				if(gzcat == 1){
					strcat ( command, "gzcat "); 
				}else{
					strcat ( command, "zcat "); 
				}
				if(!param->filter){
					i = sprintf (tmp, "%s | samtools view -F 768 - ", param->infile[file_num]);
					strcat ( command, tmp);
				}else{
					i = sprintf (tmp, "%s | samtools view -F  ", param->infile[file_num]);
					strcat ( command, tmp);
					i = sprintf (tmp, "%s - ",param->filter);
					//strcat ( command, "samtools view -F "); 
					//i = sprintf (tmp, "%s ",param->filter);
					strcat ( command, tmp);
				}
				//i = sprintf (tmp, "%s | samtools view -F 768 - ", param->infile[file_num]);
				//strcat ( command, tmp);
			}else{
				if(!param->filter){
					strcat ( command, "samtools view -F 768 "); 
				}else{
					strcat ( command, "samtools view -F "); 
					i = sprintf (tmp, "%s ",param->filter);
					strcat ( command, tmp);
				}
				//strcat ( command, "samtools view -F 768 ");     // Copy name into full name
				i = sprintf (tmp, "%s ", param->infile[file_num]);
				strcat ( command, tmp);
			}
			if (!(file = popen(command, "r"))) {
				fprintf(stderr,"Cannot open bam file '%s' with command:%s\n",param->infile[file_num],command);
				exit(-1);
			}
		}else if(param->sam == 1){
			command[0] = 0;
			if(param->gzipped){
				if(gzcat == 1){
					strcat ( command, "gzcat "); 
				}else{
					strcat ( command, "zcat "); 
				}
				if(!param->filter){
					i = sprintf (tmp, "%s | samtools view -SF 768 - ", param->infile[file_num]);
					strcat ( command, tmp);
				}else{
					i = sprintf (tmp, "%s | samtools view -SF  ", param->infile[file_num]);
					strcat ( command, tmp);
					i = sprintf (tmp, "%s - ",param->filter);
					//strcat ( command, "samtools view -F "); 
					//i = sprintf (tmp, "%s ",param->filter);
					strcat ( command, tmp);
				}
			}else{
				if(!param->filter){
					strcat ( command, "samtools view -SF 768 "); 
				}else{
					strcat ( command, "samtools view -SF "); 
					i = sprintf (tmp, "%s ",param->filter);
					strcat ( command, tmp);
				}
				//strcat ( command, "samtools view -SF 768 ");     // Copy name into full name
				i = sprintf (tmp, "%s ", param->infile[file_num]);
				strcat ( command, tmp);
			}
			if (!(file = popen(command, "r"))) {
				fprintf(stderr,"Cannot open bam file '%s' with command:%s\n",param->infile[file_num],command);
				exit(-1);
			}
		}else{
			if (!(file = fopen(param->infile[file_num] , "r" ))){
				fprintf(stderr,"Cannot open sam file '%s'\n",param->infile[file_num]);
				exit(-1);
			}
		}
	}
	return file;
}


void print_just_seq(struct read_info* ri,FILE* out)
{
	char alphabet[] = "ACGTN";
	int i;
	fprintf(out,"@%s\n",ri->name);
	for(i =0;i < ri->len;i++){
		fprintf(out,"%c",alphabet[(int) ri->seq[i]]);
	}
	fprintf(out,"\n+\n");
	if(ri->qual){
		fprintf(out,"%s",ri->qual);
	}else{
		for(i =0;i < ri->len;i++){
			fprintf(out,".");
		}
	}
	fprintf(out,"\n");
}


void print_seq(struct read_info* ri,FILE* out)
{
	char alphabet[] = "ACGTN";
	int i;
	fprintf(out,"@%s\n",ri->name);
	for(i =0;i < ri->len;i++){
		fprintf(out,"%c",alphabet[(int) ri->seq[i]]);
	}
	fprintf(out,"\n+\n");
	if(ri->qual){
		fprintf(out,"%s",ri->qual);
	}else{
		for(i =0;i < ri->len;i++){
			fprintf(out,".");
		}
	}
	fprintf(out,"\n");
}


int read_sam_chunk(struct read_info** ri,struct parameters* param,FILE* file)
{
	char line[MAX_LINE];
	int column = 0; 
	int i,j,g,tmp;
	unsigned int pos;
	int read = 0;
	int c = 0;
	int hit = 0;
	
	for(i = 0; i < param->num_query;i++){
		free(ri[i]->seq);
		free(ri[i]->name);
		free(ri[i]->qual);
		if(ri[i]->priors){
			free(ri[i]->priors);
		}
		
		if(ri[i]->cigar){
			free(ri[i]->cigar);
		}
		if(ri[i]->md){
			free(ri[i]->md);
		}
		
		if(ri[i]->xp){
			free(ri[i]->xp);
		}
		
		ri[i]->priors = 0;
		ri[i]->seq = 0;
		ri[i]->name = 0;
		ri[i]->qual = 0;
		ri[i]->len = 0;
		ri[i]->errors = -1;
		ri[i]->cigar = 0;
		ri[i]->md = 0;
		ri[i]->xp = 0;
		ri[i]->read_start = -1;
		ri[i]->read_end = -1;
		
	}
	
	while(fgets(line, MAX_LINE, file)){
		if(line[0] != '@'){
			column = 1; //<QNAME> 
			tmp = 0;
			hit = 0;
			pos = 0xFFFFFFFFu;
			for(j = 0;j < MAX_LINE;j++){
				tmp++;
				if(isspace((int)line[j])){
					break;
				}
			}
			
			ri[c]->name = malloc(sizeof(unsigned char)* tmp);
			for(j = 0;j < MAX_LINE;j++){
				
				if(isspace((int)line[j])){
					ri[c]->name[j] = 0;
					break;
				}
				ri[c]->name[j] = line[j];
			}
			
			for(i = 0; i < MAX_LINE;i++){
				if(line[i] == '\n'){
					break;
				}
				if(isspace((int)line[i])){
					column++;
					switch(column){
						case 2: // <FLAG>
							tmp = atoi(line+i+1);
							ri[c]->strand[hit] = (tmp & 0x10);
							if(tmp == 4){
								ri[c]->hits[hit] = 0;
							}else{
								ri[c]->hits[hit] = 1;
							}
							hit++;
							break;
						case 3: // <RNAME> 
							
							break;
						case 4: // <POS>
							
							break;
						case 5: //  <MAPQ>
							
							ri[c]->mapq =  atof(line +i +1); 
							
							break;
						case 6: //  <CIGAR>
							tmp = 0;
							for(j = i+1;j < MAX_LINE;j++){
								tmp++;
								if(isspace((int)line[j])){
									break;
								}
							}
							
							ri[c]->cigar = malloc(sizeof(unsigned char)* tmp);
							g = 0;
							for(j = i+1;j < MAX_LINE;j++){
								if(isspace((int)line[j])){
									ri[c]->cigar[g] = 0;
									break;
								}
								ri[c]->cigar[g] = line[j];
								g++;
							}
							break;
						case 7: //  <MRNM>
							break;
						case 8: //  <MPOS>
							break;
						case 9: //  <ISIZE>
							break;
						case 10: // <SEQ>
							
							tmp = 0;
							for(j = i+1;j < MAX_LINE;j++){
								tmp++;
								if(isspace((int)line[j])){
									break;
								}
							}
							
							ri[c]->seq = malloc(sizeof(unsigned char)* tmp);
							g = 0;
							for(j = i+1;j < MAX_LINE;j++){
								
								if(isspace((int)line[j])){
									ri[c]->seq[g] = 0;
									break;
								}
								ri[c]->seq[g] = nuc_code5[(int)line[j]];
								g++;
							}
							
							ri[c]->len = g;
							break;
						case 11: // <QUAL>
							tmp = 0;
							for(j = i+1;j < MAX_LINE;j++){
								tmp++;
								if(isspace((int)line[j])){
									break;
								}
							}
							g= 0;
							ri[c]->qual = malloc(sizeof(unsigned char)* tmp);
							for(j = i+1;j < MAX_LINE;j++){
								
								if(isspace((int)line[j])){
									ri[c]->qual[g] = 0;
									break;
								}
								ri[c]->qual[g] = line[j];
								g++;
							}
							break;
						default: 
							
									
							i = MAX_LINE;
							break;
					}				}

			}
			tmp = byg_end("NM:i:", line  );
			if(tmp){
				ri[c]->errors = atoi(line+tmp);
				//if(ri[c]->errors > 20){
				///fprintf(stderr,"%s\n,%c,%c,%c,%d\n",line, *(line +tmp), *(line +tmp+1),*(line +tmp+2), ri[c]->errors); 
				//}
				
			}else{
				ri[c]->errors = -1;								
			}	
			
			tmp = byg_end("MD:Z:", line  );
			if(tmp){
				g = 0;
				for(j = tmp ;j < MAX_LINE;j++){
					g++;
					if(isspace((int)line[j])){
						break;
					}
					
				}
				ri[c]->md = malloc(sizeof(unsigned char)* g);
				g = 0;
				for(j = tmp ;j < MAX_LINE;j++){
					
					if(isspace((int)line[j])){
						ri[c]->md[g] = 0;
						break;
					}
					ri[c]->md[g] = line[j];
					g++;
				}
			}
			
			tmp = byg_end("XP:Z:", line  );
			if(tmp){
				g = 0;
				for(j = tmp ;j < MAX_LINE;j++){
					g++;
					if(isspace((int)line[j])){
						break;
					}
					
				}
				ri[c]->xp = malloc(sizeof(unsigned char)* g);
				g = 0;
				for(j = tmp ;j < MAX_LINE;j++){
					
					if(isspace((int)line[j])){
						ri[c]->xp[g] = 0;
						break;
					}
					ri[c]->xp[g] = line[j];
					g++;
				}
			}
			
			ri[c]->hits[hit] = 0xFFFFFFFFu;
			ri[c]->identity[hit] = -2.0f;
			
			c++;
			read++;
			if(c == param->num_query){
				return c;
			}
		}
	}
	return c;
}

int read_fasta_fastq(struct read_info** ri,struct parameters* param,FILE *file) 
{
	int park_pos = -1;
	char line[MAX_LINE];
	int i;//,j;
	int seq_p = 0;
	int set = 0;
	int len = 0;
	int size = 0;
	
	for(i = 0; i < param->num_query;i++){
		free(ri[i]->seq);
		free(ri[i]->name);
		free(ri[i]->qual);
		if(ri[i]->priors){
			free(ri[i]->priors);
		}
		//free(ri[i]->strand);
		ri[i]->priors = 0;
		ri[i]->seq = 0;
		ri[i]->name = 0;
		ri[i]->qual = 0;
		ri[i]->len = 0;
		ri[i]->md = 0;
		ri[i]->xp = 0;
		ri[i]->cigar = 0;
		ri[i]->errors = 0;
		ri[i]->read_start = -1;
		ri[i]->read_end = -1;
		//ri[i]->strand = 0;
	}
	
	while(fgets(line, MAX_LINE, file)){
		if((line[0] == '@' && !set)|| (line[0] == '>' && !set)){
			//set sequence length of previous read
			
			//check if there is still space.... 
			if(param->num_query == size){
				fseek (file , -  strlen(line) , SEEK_CUR);
				
				return size;
			}
			park_pos++;
			len = 0;
			seq_p = 1;
			for(i = 1;i < MAX_LINE;i++){
				len++;
				if(isspace((int)line[i])){
					break;
				}
				
			}
			
			ri[park_pos]->hits[0] = 0;
			ri[park_pos]->strand[0] = 0;
			ri[park_pos]->name = malloc(sizeof(unsigned char)* (len+1));
			for(i = 1;i < MAX_LINE;i++){
				
				if(isspace((int)line[i])){
					ri[park_pos]->name[i-1] = 0;
					break;
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
					for(i = 0;i < MAX_LINE;i++){
						len++;
						if(isspace((int)line[i])){
							break;
						}
					}
					//fprintf(stderr,"SEQ LEN:%d	%s\n",len,line);
					ri[park_pos]->seq = malloc(sizeof(unsigned char)* (len+1));
					for(i = 0;i < MAX_LINE;i++){
						if(isspace((int)line[i])){
							ri[park_pos]->seq[i] = 0;
							break;
						}
						ri[park_pos]->seq[i] = nuc_code5[(int)line[i]];
					}
					ri[park_pos]->len = len-1;
					
				}else{
					len = 0;
					for(i = 0;i < MAX_LINE;i++){
						len++;
						if(isspace((int)line[i])){
							break;
						}
						
					}
					//fprintf(stderr,"QUAL LEN:%d\n",len);
					ri[park_pos]->qual = malloc(sizeof(unsigned char)* (len+1));
					for(i = 0;i < MAX_LINE;i++){
						if(isspace((int)line[i])){
							ri[park_pos]->qual[i] = 0;
							break;
						}
						ri[park_pos]->qual[i] = line[i];
					}
				}
			}
			set = 0;
		}
	}
	return size;
}

