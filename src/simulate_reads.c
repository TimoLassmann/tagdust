#if HAVE_CONFIG_H
#include "config.h"
#endif
#include "kslib.h"

#include "interface.h"
#include "nuc_code.h"
#include "misc.h"
#include <time.h>

#include <sys/stat.h>




char* mutate(struct parameters* param, char* seq,char* seq_mutated,int len,unsigned int my_rand_max );

int main (int argc,char * argv[]) {
	struct parameters* param = NULL;
	int status;
	unsigned int seed = 0;
	
	init_nuc_code();
	
	param = interface(param,argc,argv);
	
	if(!param){
		return kslOK;
	}
	if(param->seed){
		seed = param->seed;
	}else{
		seed = (unsigned int) (time(NULL) * ( 42));
	}
	
	srand(seed);
	
#ifdef RTEST
	unsigned int my_rand_max = 32768;
#else
	unsigned int my_rand_max = RAND_MAX;
#endif
	
	
	FILE* file = stdout;
	int i,j,c,n,start,barcode_used;
	int read_barcodes = 0;
	float r = 0.0;
	char* read = 0;
	char* sequenced_read = 0;
	char* sequenced_read_mutated = 0;
	char* tmp_str = 0;
	
	char** barcode = 0;
	
	struct stat buf;
		
	if(param->infiles == 0){
		exit(EXIT_FAILURE);
	}
	
	MMALLOC(barcode, sizeof(char*) * 1000);
	for(i = 0; i < 1000;i++){
		barcode[i] = 0;
	}
	
	if ((file = fopen(param->infile[0], "r")) == NULL){
		fprintf(stderr,"can't open barcode file. \n");
		exit(-1);
	}
	char line[10000];
	n =0;
	read_barcodes = 0;
	while(fgets(line, 10000, file)){
		c = byg_end(":", line  );
		if(c){
			MMALLOC(barcode[read_barcodes],sizeof(char) * 100);
			j = 0;
			for(i = c; i< strlen(line);i++){
				
				switch(line[i]){
					case 'A':
					case 'C':
					case 'G':
					case 'T':
					case 'a':
					case 'c':
					case 'g':
					case 't':
						barcode[read_barcodes][j] = line[i];
						j++;
						n =1;
						break;
					default:
						barcode[read_barcodes][j] = 0;
						j++;
						n = 0;
						break;
						
				}
				if(!n){
					break;
				}
			}
			read_barcodes++;
		}
	}
	fclose(file);
	
	if(read_barcodes < param->sim_barnum){
		sprintf(param->buffer,"File contains too few barcodes.\n");
		param->messages = append_message(param->messages, param->buffer);
		free_param(param);
		exit(EXIT_FAILURE);
	}
	
	//for(i = 0; i < read_barcodes;i++){
	//	fprintf(stderr,"%d	%s\n",i,barcode[i]);
	//}
	//exit(EXIT_SUCCESS);
	
	
	MMALLOC(read,sizeof(char)* 1000);
	
	
	MMALLOC(tmp_str,sizeof(char)* 1000);
	MMALLOC(sequenced_read,sizeof(char)* 1000);
	MMALLOC(sequenced_read_mutated,sizeof(char)* 1020);
	
	
	
	if(param->outfile){
		if ((file = fopen(param->outfile, "w")) == NULL){
			fprintf(stderr,"can't open barcode file. \n");
			exit(-1);
		}

	}else{
		file = stdout;
	}
	//Here I simulate sequences containing the read....
	for(i = 0; i < (int)((float) param->sim_numseq * (1.0-param->sim_random_frac));i++){
		for(j = 0; j < 1000;j++){
			sequenced_read[j] = 0;
			sequenced_read_mutated[j]  = 0;
			tmp_str[j] = 0;
			read[j]  =0;
		}
		
		if(param->sim_5seq){
			tmp_str[0] = 0;
			strcat ( tmp_str,  param->sim_5seq);
		}
		barcode_used = 0;
		if(param->sim_barnum){
			barcode_used = (int) (rand() % (int) (param->sim_barnum)) ;
			
			strcat ( tmp_str, barcode[barcode_used]);
		}
		
		sequenced_read_mutated = mutate(param, tmp_str, sequenced_read_mutated, (int)strlen(tmp_str), my_rand_max);
		//fprintf(stderr,"%s\n%s\n",sequenced_read,sequenced_read_mutated);
		
		
		strcat ( sequenced_read, sequenced_read_mutated);

		sequenced_read_mutated[0] = 0;
		
		
		if(param->sim_readlen_mod){
			c = param->sim_readlen - param->sim_readlen_mod +  (int) (rand() % (int) (param->sim_readlen_mod*2)) ;
		}else{
			c = param->sim_readlen;
		}
		//fprintf(stderr,"%d\n",c);
		for(j = 0; j < c;j++){
			r = (float)rand()/(float)my_rand_max;
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
		
		strcat ( sequenced_read, read);
		
		
		
		if(param->sim_3seq){
			tmp_str[0] = 0;
			strcat ( tmp_str,  param->sim_3seq);
			sequenced_read_mutated = mutate(param, tmp_str, sequenced_read_mutated, (int)strlen(tmp_str), my_rand_max);
			//fprintf(stderr,"%s\n%s\n",read,sequenced_read_mutated);
		}
		strcat ( sequenced_read, sequenced_read_mutated);
		
		
		if(param->sim_end_loss){
			start = (int) (rand() % (int) (param->sim_end_loss*2));
			n = 0;
			c = (int)strlen(sequenced_read);
			for(j = start; j < c;j++ ){
				sequenced_read[n]  = sequenced_read[j];
				n++;
			}
			sequenced_read[n]  = 0;
			
			start = (int) (rand() % (int) (param->sim_end_loss*2));
			
			for(j = 0; j < start;j++){
				sequenced_read[n-1-j]  = 0;
			
			}
			
			
		}
		
		        
		        
		        
		
		
		//sequenced_read_mutated[param->sim_sequenced_len] = 0;
		
		//fprintf(stderr,"%d	%s\n",i, sequenced_read_mutated);
		//fprintf(file,"@READ%d;SEQ:%s;",i, read);
		if(param->sim_barnum){
			fprintf(file,"@READ%d;SEQ:%s;RBC:%s;BARNUM:%d",i, read,barcode[barcode_used],barcode_used+1);
			//fprintf(file,"RBC:%s;BARNUM:%d;",barcode[barcode_used],barcode_used+1);
		}else{
			fprintf(file,"@READ%d;SEQ:%s;BARNUM:%d",i, read,1);
		}
		fprintf(file,"\n");
		fprintf(file,"%s\n+\n",sequenced_read);
		c = (int)strlen(sequenced_read);
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
	
	
	
	for(i =  (int)((float) param->sim_numseq * (1.0-param->sim_random_frac)) ; i < param->sim_numseq   ;i++){
		if(param->sim_5seq){
			strcat ( sequenced_read,  param->sim_5seq);
		}
		for(j = 0; j < c;j++){
			r = (float)rand()/(float)my_rand_max;
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
	
		sequenced_read[c] = 0;
		if(param->sim_end_loss){
			start = (int) (rand() % (int) (param->sim_end_loss*2));
			n = 0;
			for(j = start; j < c;j++ ){
				sequenced_read_mutated[n]  = sequenced_read_mutated[j];
				n++;
			}
			sequenced_read_mutated[n]  = 0;
			
			start = (int) (rand() % (int) (param->sim_end_loss*2));
			
			for(j = 0; j < start;j++){
				sequenced_read_mutated[n-1-j]  = 0;
				
			}
			
			
		}
		
		//sequenced_read[param->sim_sequenced_len] = 0;
		//fprintf(stderr,"%s\n", sequenced_read);
		//fprintf(file,"@RAND%d;SEQ:NONE;",i);
		if(param->sim_barnum){
			fprintf(file,"@RAND%d;SEQ:NONE;RBC:NONE;BARNUM:0",i);
		}else{
			fprintf(file,"@RAND%d;SEQ:NONE;BARNUM:0",i);

		}
		fprintf(file,"\n");
		fprintf(file,"%s\n+\n",sequenced_read);
		for(j = 0; j < c;j++){
			fprintf(file,"I");
		}
		fprintf(file,"\n");
		
		
	}
	if(param->outfile){
		fclose(file);
	}
	
	
	param->buffer[0] = 0;
	sprintf(param->buffer, "%s_tagdust_arch.txt",param->outfile );
	
	
	
	
	if(!stat ( param->buffer, &buf )){
		if ((file = fopen(param->buffer, "a")) == NULL){
			fprintf(stderr,"can't open output\n");
			exit(-1);
		}

	}else {
		if ((file = fopen(param->buffer, "w")) == NULL){
			fprintf(stderr,"can't open output\n");
			exit(-1);
		}
	}
	
	
	
	c = 1;
	
	fprintf(file,"tagdust ");
	
	if(param->sim_5seq){
		fprintf(file,"-%d ",c);
		c++;
		fprintf(file,"P:%s ",param->sim_5seq);
	}
	
	if(param->sim_barnum){
		fprintf(file,"-%d ",c);
		c++;
		
		fprintf(file,"B:");
		for(i = 0; i < param->sim_barnum-1;i++){
			fprintf(file,"%s,",barcode[i]);
		}
		fprintf(file,"%s ",barcode[param->sim_barnum-1]);
	}
	
	fprintf(file,"-%d ",c);
	c++;
	
	fprintf(file,"R:N ");
	
	if(param->sim_3seq){
		fprintf(file,"-%d ",c);
		c++;
		fprintf(file,"P:%s ",param->sim_3seq);
	}
	
	fprintf(file,"%s ", "in.fq");
	
	fprintf(file,"-o out.fq");

	
	
	fprintf(file,"\n");
	
	
	
	
	if(param->outfile){
		fclose(file);
	}
	
	
	if(param->outfile){
		param->buffer[0] = 0;
		sprintf(param->buffer, "%s_btrim_pattern.txt",param->outfile );
		if ((file = fopen(param->buffer, "w")) == NULL){
			fprintf(stderr,"can't open barcode file. \n");
			exit(-1);
		}
		
	}else{
		file = stdout;
	}
	
	
	if(param->sim_barnum){
		for(i = 0 ;i < param->sim_barnum;i++){
			if(param->sim_5seq){
				if(param->sim_3seq){
					fprintf(file,"%s%s %s\n",param->sim_5seq,barcode[i], param->sim_3seq);
				}else{
					fprintf(file,"%s%s\n",param->sim_5seq,barcode[i]);
				}
			}else{
				if(param->sim_3seq){
					fprintf(file,"%s %s\n",barcode[i], param->sim_3seq);
				}else{
					fprintf(file,"%s\n",barcode[i]);
				}
			}
		}
	}else{
		if(param->sim_5seq){
			if(param->sim_3seq){
				fprintf(file,"%s %s\n",param->sim_5seq, param->sim_3seq);
			}else{
				fprintf(file,"%s\n",param->sim_5seq);
			}
		}else{
			if(param->sim_3seq){
				fprintf(file,"%s\n",param->sim_3seq);
			}else{
				fprintf(file,"\n");
			}
		}
		
		
	}
	
	if(param->outfile){
		fclose(file);
	}
	if(param->sim_barnum){
		
		if(param->outfile){
			param->buffer[0] = 0;
			sprintf(param->buffer, "%s_fastxbarcodefile.txt",param->outfile );
			if ((file = fopen(param->buffer, "w")) == NULL){
				fprintf(stderr,"can't open barcode file. \n");
				exit(-1);
			}
			
		}else{
			file = stdout;
		}
		
		
		
		for(i = 0 ;i < param->sim_barnum;i++){
			fprintf(file,"BC%d %s\n",i,barcode[i]);
			
		}
		if(param->outfile){
			fclose(file);
		}
		
		
	}
	
	for(i = 0; i <read_barcodes;i++){
		MFREE(barcode[i]);
	}
	MFREE(barcode);
	MFREE(tmp_str);
	MFREE(read);// = malloc(sizeof(char)* 200);
	MFREE(sequenced_read);// = malloc(sizeof(char)* 200);
	MFREE(sequenced_read_mutated);// = malloc(sizeof(char)* 220);
	free_param(param);
	return kslOK;
ERROR:
	KSLIB_MESSAGE(status,"Something wrong in main...\n");
	return kslFAIL;
}


char* mutate(struct parameters* param, char* seq,char* seq_mutated,int len,unsigned int my_rand_max )
{
	
	int j;
	float r = 0.0f;
	char n = 0;
	int newlen = 0;
	
	float cutoff = 0.5;
	
	for(j = 0;j <  strlen(seq) ;j++){
		r = (float)rand()/(float)my_rand_max;
		if(r <= param->sim_error_rate){
			//we have an error
			
			
			r = (float)rand()/(float)my_rand_max;
			if(r <= param->sim_InDel_frac){
				//indel++;
				// we have an indel (only considering single nucleotide.....
				r = (float)rand()/(float)my_rand_max;
				if(j ==  strlen(seq)  -1){
					cutoff = 0.0;
				}else{
					cutoff = 0.5;

				}
				if(r <= cutoff){
					//insertion
					//n_dash = read[j];
					r = (float)rand()/(float)my_rand_max;
					if(r < 0.25){
						n = 'A';
					}else if(r < 0.5){
						n = 'C';
					}else if(r < 0.75){
						n = 'G';
					}else{
						n = 'T';
					}
					seq_mutated[newlen] = seq[j];
					newlen++;
					seq_mutated[newlen] = n;
					newlen++;
					
					
					
					
				}else{
					//deletion
				}
				
			}else{
				//mismatches++;
				n = seq[j];
				
				while(n == seq[j]){
					r = (float)rand()/(float)my_rand_max;
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
				seq_mutated[newlen] = n;
				newlen++;
			}
		}else{
			seq_mutated[newlen] = seq[j];
			newlen++;
		}
	}
	seq_mutated[newlen]  =0;
	
	return seq_mutated;
}
