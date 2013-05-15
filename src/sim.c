//
//  sim.c
//  tagdust2
//
//  Created by lassmann on 5/14/13.
//  Copyright (c) 2013 lassmann. All rights reserved.
//

#include <stdio.h>
#include <time.h>
#include "interface.h"
#include <math.h>
#include <assert.h>
#include "misc.h"
//unsigned int seed = (unsigned int) (time(NULL) * ( data->threadID + 1));

//function = (int) (rand_r(split->seed) % (int) (stp->feather)) ;


void simulate(struct parameters* param)
{
	int i,j,c,n,n_dash;
	float r = 0.0;
	FILE* file = 0;
	int errors_allowed = 6;
	char alpha[5] = "ACGTN";
	
	char* outfile = 0;
	outfile = malloc(sizeof(char)*(200));
	assert(outfile != 0);
	

	sprintf (outfile, "%s_tagdust_command.sh",param->outfile);
	
	
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
	sprintf (outfile, "%s_read_extracted.fq",param->outfile);
	
	fprintf(file," -2 R:N -o %s -threshold 0.5 -e 0.05 \n", outfile);
	
	//fprintf(file,"echo -n %d\t%d\t%f\t;\t",param->numbarcode, param->sim,param->sequencer_error_rate);
	fprintf(file,"grep READ  %s  |  awk -v numbarcode=%d -v sim=%d -v errorrate=%f 'BEGIN{p=0;n=0}{x = split($0,a,\"[;,:]\");if(a[x] == a[3]){p++}else{n++}}END{printf \"%%d\\t%%d\\t%%f\\t%%d\\t%%d\\n\",numbarcode,sim,errorrate,p,n}'  >> tagdust_benchmark.csv &\n", outfile,param->numbarcode, param->sim,param->sequencer_error_rate );
	
	sprintf (outfile, "%s_read.fq",param->outfile);
	fprintf(file,"cat %s | ", outfile);
	
	
	sprintf (outfile, "%s_barcodefile.txt",param->outfile);
	
	fprintf(file,"./fastx_barcode_splitter.pl  --bcfile %s --bol --mismatches 2 --prefix ~/tmp/bla2_ --suffix \".txt\"  | awk -v numbarcode=%d -v sim=%d -v errorrate=%f  '{printf\"%%d\\t%%d\\t%%f\\t%%d\\t%%d\\t%%d\\n\" ,numbarcode,sim,errorrate,$1,$2,$3 }' >> fastx_benchmark_2.csv &\n",outfile,param->numbarcode, param->sim,param->sequencer_error_rate );
	
	sprintf (outfile, "%s_read.fq",param->outfile);
	fprintf(file,"cat %s | ", outfile);
	
	
	sprintf (outfile, "%s_barcodefile.txt",param->outfile);

	
	fprintf(file,"./fastx_barcode_splitter.pl  --bcfile %s --bol --mismatches 1 --prefix ~/tmp/bla1_ --suffix \".txt\"  | awk -v numbarcode=%d -v sim=%d -v errorrate=%f  '{printf\"%%d\\t%%d\\t%%f\\t%%d\\t%%d\\t%%d\\n\" ,numbarcode,sim,errorrate,$1,$2,$3 }' >> fastx_benchmark_1.csv &\n",outfile,param->numbarcode, param->sim,param->sequencer_error_rate );
	
	sprintf (outfile, "%s_read.fq",param->outfile);
	fprintf(file,"cat %s | ", outfile);
	
	
	sprintf (outfile, "%s_barcodefile.txt",param->outfile);

	
	fprintf(file,"./fastx_barcode_splitter.pl  --bcfile %s --bol --mismatches 0 --prefix ~/tmp/bla0_ --suffix \".txt\"  | awk -v numbarcode=%d -v sim=%d -v errorrate=%f  '{printf\"%%d\\t%%d\\t%%f\\t%%d\\t%%d\\t%%d\\n\" ,numbarcode,sim,errorrate,$1,$2,$3 }' >> fastx_benchmark_0.csv \n",outfile,param->numbarcode, param->sim,param->sequencer_error_rate );
	
	
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
	assert(read!=0);
	
	int num_bar_correct = 0;
	int num_bar_misassigned = 0;
	int num_bar_missed = 0;
	
	int num_finger_correct = 0;
	int num_finger_misassigned = 0;
	//int num_finger_missed = 0;
	
	sprintf (outfile, "%s_read.fq",param->outfile);
	if ((file = fopen(outfile, "w")) == NULL){
		fprintf(stderr,"can't open output\n");
		exit(-1);
	}
	
	for(i = 0; i < 1000000;i++){
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
		
		
		// test if read is recognizeable...
		n = -1;
		for(j = 0;j < num_barcode;j++){
			if(bpm(barcode[j], read, param->sim, param->sim)  == 0){
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
		
		
		for(j = 0;j < param->sim;j++){
			fprintf(file,"%c",alpha[(int) read[j]]);// = barcode[c][j];
		}
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
			fprintf(file,"%c",alpha[(int) n ]);// = barcode[c][j];
		}
		fprintf(file,"\n");
		
		fprintf(file,"+\n");
		for(j = 0; j < 20+param->sim;j++){
			fprintf(file,".");
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
	fprintf(file,"%d\t%d\t%f\t%d\t%d\t%d\n",param->numbarcode, param->sim,param->sequencer_error_rate,num_bar_correct,num_bar_missed,num_bar_misassigned);
	
	fclose(file);
	fprintf(stderr,"BARCODE:\n");
	fprintf(stderr,"%d(%0.1f) correct\n", num_bar_correct,(float)num_bar_correct / 1000000.0* 100.0);
	fprintf(stderr,"%d(%0.1f) num_bar_misassigned\n", num_bar_misassigned,(float)num_bar_misassigned / 1000000.0* 100.0 );
	fprintf(stderr,"%d(%0.1f) num_bar_missed\n", num_bar_missed,(float)num_bar_missed / 1000000.0 * 100.0);
	
	//fprintf(stderr,"Fingerprint:\n");
	//fprintf(stderr,"%d(%0.1f) correct\n", num_finger_correct,(float)num_finger_correct / 1000000.0* 100.0);
	//fprintf(stderr,"%d(%0.1f) num_bar_misassigned\n", num_finger_misassigned,(float)num_finger_misassigned / 1000000.0* 100.0 );
	
	
	for(i = 0 ;i < num_barcode;i++){
		
		free(barcode[i]);
	}
	free(barcode);
	
	free(outfile);

}


