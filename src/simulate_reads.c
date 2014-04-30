
#include "interface.h"
#include "nuc_code.h"
#include "misc.h"
#include <time.h>

#include <sys/stat.h>


#ifndef MMALLOC
#include "malloc_macro.h"
#endif

/*! \brief  Initializes nucleotide alphabet needed to parse input. Calls parameter parser. Calls functions to process the data.
 * \param argc number of command line parameters
 * \param argv command line parameters
 * \return EXIT_SUCCESS */

/*char* datadir;
char datafile[1000];
int i;

init_nuc_code();

datadir = getenv("testdatafiledir");

fprintf(stdout,"%s\n",datadir);
TOMEDB* db = 0;

if(db_exists ("testdb")){
	fprintf(stderr,"A database with the same name already exists.\n");
	abort();
}

i = sprintf(datafile, "%s/%s",datadir,"hg19_chr22_M.fa");
*/

int main (int argc,char * argv[]) {
	struct parameters* param = 0;
	
	unsigned int seed = 0;
	
	init_nuc_code();
	
	param = interface(param,argc,argv);
	
	
	if(param->seed){
		seed = param->seed;
	}else{
		seed = (unsigned int) (time(NULL) * ( 42));
	}
	
	FILE* file = stdout;
	int i,j,c,n,start,barcode_used;
	int read_barcodes = 0;
	float r = 0.0;
	char* read = 0;
	char* sequenced_read = 0;
	char* sequenced_read_mutated = 0;
	
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
	
	
	MMALLOC(read,sizeof(char)* 200);
	MMALLOC(sequenced_read,sizeof(char)* 200);
	MMALLOC(sequenced_read_mutated,sizeof(char)* 220);
	
	
	
	if(param->outfile){
		if ((file = fopen(param->outfile, "w")) == NULL){
			fprintf(stderr,"can't open barcode file. \n");
			exit(-1);
		}

	}else{
		file = stdout;
	}
	//Here I simulate sequences containing the read....
	for(i = 0; i <= (int)((float) param->sim_numseq * (1.0-param->sim_random_frac));i++){
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
		if(param->sim_readlen_mod){
			c = param->sim_readlen - param->sim_readlen_mod +  (int) (rand_r(&seed) % (int) (param->sim_readlen_mod*2)) ;
		}else{
			c = param->sim_readlen;
		}
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
		sequenced_read_mutated[c]  =0;
		
		
		
		if(param->sim_end_loss){
			start = (int) (rand_r(&seed) % (int) (param->sim_end_loss*2));
			n = 0;
			for(j = start; j < c;j++ ){
				sequenced_read_mutated[n]  = sequenced_read_mutated[j];
				n++;
			}
			sequenced_read_mutated[n]  = 0;
			
			start = (int) (rand_r(&seed) % (int) (param->sim_end_loss*2));
			
			for(j = 0; j < start;j++){
				sequenced_read_mutated[n-1-j]  = 0;
			
			}
			
			
		}
		
		        
		        
		        
		
		
		//sequenced_read_mutated[param->sim_sequenced_len] = 0;
		
		//fprintf(stderr,"%d	%s\n",i, sequenced_read_mutated);
		//fprintf(file,"@READ%d;SEQ:%s;",i, read);
		if(param->sim_barnum){
			fprintf(file,"@READ%d;SEQ:%s;RBC:%s;BARNUM:%d;",i, read,barcode[barcode_used],barcode_used+1);
			//fprintf(file,"RBC:%s;BARNUM:%d;",barcode[barcode_used],barcode_used+1);
		}else{
			fprintf(file,"@READ%d;SEQ:%s;BARNUM:%d;",i, read,1);
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
	
	
	
	for(i = (int)((float) param->sim_numseq * (1.0-param->sim_random_frac)) + 1; i < param->sim_numseq   ;i++){
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
	
		sequenced_read[c] = 0;
		if(param->sim_end_loss){
			start = (int) (rand_r(&seed) % (int) (param->sim_end_loss*2));
			n = 0;
			for(j = start; j < c;j++ ){
				sequenced_read_mutated[n]  = sequenced_read_mutated[j];
				n++;
			}
			sequenced_read_mutated[n]  = 0;
			
			start = (int) (rand_r(&seed) % (int) (param->sim_end_loss*2));
			
			for(j = 0; j < start;j++){
				sequenced_read_mutated[n-1-j]  = 0;
				
			}
			
			
		}
		
		//sequenced_read[param->sim_sequenced_len] = 0;
		//fprintf(stderr,"%s\n", sequenced_read);
		//fprintf(file,"@RAND%d;SEQ:NONE;",i);
		if(param->sim_barnum){
			fprintf(file,"@RAND%d;SEQ:NONE;RBC:NONE;BARNUM:0;",i);
		}else{
			fprintf(file,"@RAND%d;SEQ:NONE;BARNUM:0;",i);

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
	
	fprintf(file,"%s ", param->outfile);
	
	fprintf(file,"-o bananas");

	
	
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
			fprintf(file,"%s%s %s\n",param->sim_5seq,barcode[i], param->sim_3seq);
			
		}
	}else{
		fprintf(file,"%s %s\n",param->sim_5seq, param->sim_3seq);
		
	}
	
	if(param->outfile){
		fclose(file);
	}
	
	
	
	for(i = 0; i <read_barcodes;i++){
		MFREE(barcode[i]);
	}
	MFREE(barcode);
	
	MFREE(read);// = malloc(sizeof(char)* 200);
	MFREE(sequenced_read);// = malloc(sizeof(char)* 200);
	MFREE(sequenced_read_mutated);// = malloc(sizeof(char)* 220);
	free_param(param);
	return EXIT_SUCCESS;
}
