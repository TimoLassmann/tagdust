//
//  barcode_hmm.c
//  tagdust2
//
//  Created by lassmann on 2/5/13.
//  Copyright (c) 2013 lassmann. All rights reserved.
//

#include <stdio.h>
#include "tagdust2.h"
#include "interface.h"
#include "io.h"
#include "misc.h"


#include "barcode_hmm.h"


void hmm_controller(struct parameters* param,int (*fp)(struct read_info** ,struct parameters*,FILE* ),int file_num)
{
	struct read_info** ri = 0;
	int i;
	int numseq;
	
	
	init_logsum();
	
	
	float* back = 0;
	
	back = malloc(sizeof(float)*5);
	for(i = 0; i < 5;i++){
		back[i]= prob2scaledprob( 0.2);
	}
	
	struct model* model =  malloc_model(16, 10, 64);
		
	model = init_model(model , back,24);
	
	struct model* model2 = copy_and_malloc_model(model);
	
	free_model(model);
	free_model(model2);
	free(back);
	
	free_param(param);
	exit(0);
	//char command[1000];
	//char  tmp[1000];
	FILE* file = 0;
	FILE* unmapped = 0;
	
	ri = malloc(sizeof(struct read_info*) * param->num_query);
	
	for(i = 0; i < param->num_query;i++){
		ri[i] = malloc(sizeof(struct read_info));
		ri[i]->seq = 0;
		ri[i]->name = 0;
		ri[i]->qual = 0;
		ri[i]->len = 0;
		ri[i]->cigar = 0;
		ri[i]->md = 0;
		ri[i]->xp = 0;
		ri[i]->priors = 0;// malloc(sizeof(unsigned int)* (LIST_STORE_SIZE+1));
		ri[i]->strand = malloc(sizeof(unsigned int)* (LIST_STORE_SIZE+1));
		ri[i]->hits = malloc(sizeof(unsigned int)* (LIST_STORE_SIZE+1));
		ri[i]->identity = malloc(sizeof(float)* (LIST_STORE_SIZE+1));
		ri[i]->read_start = -1;
		ri[i]->read_end = -1;
	}
	file =  io_handler(file, file_num,param);
	
	if(param->print_unmapped){
		if (!(unmapped = fopen(param->print_unmapped, "w" ))){
			fprintf(stderr,"Cannot open file '%s'\n",param->print_unmapped);
			exit(-1);
		}
	}
	
	while ((numseq = fp(ri, param,file)) != 0){
		//fprintf(stderr,"rread: %d\n",numseq);
		for(i = 0; i < numseq;i++){
			
		}
	}
	
	if(param->print_unmapped){
		fclose(unmapped);
	}
	
	for(i = 0; i < param->num_query;i++){
		free(ri[i]->strand);
		free(ri[i]->hits);
		free(ri[i]->identity);
		
		if(ri[i]->cigar){
			free(ri[i]->cigar);
		}
		if(ri[i]->md){
			free(ri[i]->md);
		}
		
		free(ri[i]);
	}
	free(ri);
	//fprintf(stderr,"%p\n",file);
	if(param->sam == 2 || param->sam == 1 || param->gzipped ){
		pclose(file);
	}else{
		//if(file_num != -1){
		fclose(file);
		//}
	}

}

struct model* malloc_model(int main_length, int sub_length, int number_sub_models)
{
	struct model* model = NULL;
	int i = 0;
	int j = 0;
	
	
	assert(number_sub_models  <=MAX_NUM_SUB_MODELS );
	
	model = malloc(sizeof(struct model));
	assert(model != 0);
	
	model->num_hmms =  (1+ number_sub_models);
	model->hmms = malloc(sizeof(struct hmm*) * (1+ number_sub_models));
	assert(model->hmms !=0);
	for(i = 0; i < model->num_hmms;i++){
		model->hmms[i] = malloc(sizeof(struct hmm) );
		assert(model->hmms[i]  != 0);
	}
	
	
		
	
	
	model->hmms[0]->num_columns = main_length;
	model->hmms[0]->hmm_column = malloc(sizeof(struct hmm_column*) * main_length);
	assert(model->hmms[0]->hmm_column !=0);
	
	for(j = 0; j < main_length;j++){
		model->hmms[0]->hmm_column[j] = malloc(sizeof(struct hmm_column));
		assert(model->hmms[0]->hmm_column[j] !=0);
	}

	for(i = 1; i < model->num_hmms;i++){
		model->hmms[i]->num_columns = sub_length;
		model->hmms[i]->hmm_column = malloc(sizeof(struct hmm_column*) * sub_length);
		assert(model->hmms[i]->hmm_column !=0);
		for(j = 0; j < sub_length;j++){
			model->hmms[i]->hmm_column[j] = malloc(sizeof(struct hmm_column));
			assert(model->hmms[i]->hmm_column[j] !=0);
			//model->hmms[i]->hmm_column[j]->identifier = -1;
		}
	}
	
	
	
	return model;
}



void free_model(struct model* model)
{
	int i = 0;
	int j = 0;

		
	for(i = 1; i < model->num_hmms;i++){
		
		for(j = 0; j < model->hmms[i]->num_columns;j++){
			free(model->hmms[i]->hmm_column[j]);// = malloc(sizeof(struct hmm_column));
			//assert(model->hmms[i]->hmm_column[j] !=0);
			//model->hmms[i]->hmm_column[j]->identifier = -1;
		}
		//model->hmms[i]->num_columns = sub_length;
		free(model->hmms[i]->hmm_column);// = malloc(sizeof(struct hmm_column*) * sub_length);
		//assert(model->hmms[i]->hmm_column !=0);
	}
	for(j = 0; j < 	model->hmms[0]->num_columns;j++){
		free(model->hmms[0]->hmm_column[j]);// = malloc(sizeof(struct hmm_column));
		//assert(model->hmms[0]->hmm_column[j] !=0);
	}

	//model->hmms[0]->num_columns = main_length;
	free(model->hmms[0]->hmm_column);// = malloc(sizeof(struct hmm_column*) * main_length);
	//assert(model->hmms[0]->hmm_column !=0);
	
	///assert(model->hmms !=0);
	

	for(i = 0; i < model->num_hmms;i++){
		free(model->hmms[i]);// = malloc(sizeof(struct hmm) );
		//assert(model->hmms[i]  != 0);
	}
	free(model->hmms);// = malloc(sizeof(struct hmm*) * (1+ number_sub_models));
	
	free(model);// = malloc(sizeof(struct model));
	//assert(model != 0);

	
	//return model;
}





struct model* copy_and_malloc_model(struct model* org)
{
	struct model* model = NULL;
	int i = 0;
	int j = 0;
	int len = 0;
	int c = 0;

	struct hmm_column* col =0;

	
	
	model = malloc(sizeof(struct model));
	assert(model != 0);
	
	model->num_hmms = org->num_hmms;
	model->hmms = malloc(sizeof(struct hmm*) * org->num_hmms);
	assert(model->hmms !=0);
	for(i = 0; i < model->num_hmms;i++){
		model->hmms[i] = malloc(sizeof(struct hmm) );
		assert(model->hmms[i]  != 0);
	}
	
	
	
	
	
	model->hmms[0]->num_columns = org->hmms[0]->num_columns;
	model->hmms[0]->hmm_column = malloc(sizeof(struct hmm_column*) * model->hmms[0]->num_columns );
	assert(model->hmms[0]->hmm_column !=0);
	
	for(j = 0; j < model->hmms[0]->num_columns ;j++){
		model->hmms[0]->hmm_column[j] = malloc(sizeof(struct hmm_column));
		assert(model->hmms[0]->hmm_column[j] !=0);
	}
	
	for(i = 1; i < model->num_hmms;i++){
		model->hmms[i]->num_columns = org->hmms[i]->num_columns;
		model->hmms[i]->hmm_column = malloc(sizeof(struct hmm_column*) * model->hmms[i]->num_columns);
		assert(model->hmms[i]->hmm_column !=0);
		for(j = 0; j < model->hmms[i]->num_columns;j++){
			model->hmms[i]->hmm_column[j] = malloc(sizeof(struct hmm_column));
			assert(model->hmms[i]->hmm_column[j] !=0);
			//model->hmms[i]->hmm_column[j]->identifier = -1;
		}
	}
	
	len = model->hmms[0]->num_columns;
	
	for(i = 0; i < len;i++){
		col = model->hmms[0]->hmm_column[i];
		
		//col->identifier = -1;
		for(j = 0; j < 5;j++){
			col->i_emit[j] =  org->hmms[0]->hmm_column[i]->i_emit[j];// background[j];
			col->m_emit[j] =  org->hmms[0]->hmm_column[i]->m_emit[j];
			col->i_emit_e[j] = org->hmms[0]->hmm_column[i]->i_emit_e[j];
			col->m_emit_e[j] = org->hmms[0]->hmm_column[i]->m_emit_e[j];
		}
		col->short_transition[NEXT] = org->hmms[0]->hmm_column[i]->short_transition[NEXT]; //  prob2scaledprob( max);
		col->short_transition[SELF] =  org->hmms[0]->hmm_column[i]->short_transition[SELF];// prob2scaledprob(1.0 - max) +  prob2scaledprob(0.5);
		
		col->short_transition_e[NEXT] = org->hmms[0]->hmm_column[i]->short_transition_e[NEXT]; //  prob2scaledprob( max);
		col->short_transition_e[SELF] =  org->hmms[0]->hmm_column[i]->short_transition_e[SELF];// prob2scaledprob(1.0 - max) +  prob2scaledprob(0.5);
		
		
		for(j = 0; j< model->num_hmms-1;j++){
			col->long_transition[j] = org->hmms[0]->hmm_column[i]->long_transition[j];// prob2scaledprob(1.0 - max) +  prob2scaledprob(0.5  /  (float) (model->num_hmms-1));
			col->long_transition_e[j] = org->hmms[0]->hmm_column[i]->long_transition_e[j];
		}
		
		
	}
	
	/*
	 
	 sub models model....
	 
	 */
	len = model->hmms[1]->num_columns;
	for(i = 1; i< model->num_hmms;i++){
		for(j = 0; j < len;j++){
			col = model->hmms[i]->hmm_column[j];
			///col->identifier = -1;
			for(c = 0; c < 5;c++){
				col->i_emit[c] = org->hmms[i]->hmm_column[j]->i_emit[c];
				col->m_emit[c] = org->hmms[i]->hmm_column[j]->m_emit[c];
				
				col->i_emit_e[c] = org->hmms[i]->hmm_column[j]->i_emit_e[c];
				col->m_emit_e[c] =  org->hmms[i]->hmm_column[j]->m_emit_e[c];
			}
			
			
			
			
			col->short_transition[MM] =   org->hmms[i]->hmm_column[j]->short_transition[MM];//  prob2scaledprob( 0.999);
			col->short_transition[MI] = org->hmms[i]->hmm_column[j]->short_transition[MI];// prob2scaledprob(1.0 - 0.999) +  prob2scaledprob(0.5);
			col->short_transition[MD] = org->hmms[i]->hmm_column[j]->short_transition[MD];// prob2scaledprob(1.0 - 0.999) +  prob2scaledprob(0.5);
			
			col->short_transition[II] = org->hmms[i]->hmm_column[j]->short_transition[II];// prob2scaledprob(1.0 - 0.999);
			col->short_transition[IM] = org->hmms[i]->hmm_column[j]->short_transition[IM];// prob2scaledprob(0.999);
			
			col->short_transition[DD] = org->hmms[i]->hmm_column[j]->short_transition[DD];// prob2scaledprob(1.0 - 0.999);
			col->short_transition[DM] = org->hmms[i]->hmm_column[j]->short_transition[DM];// prob2scaledprob(0.999);
			
			col->short_transition_e[MM] =   org->hmms[i]->hmm_column[j]->short_transition_e[MM];//  prob2scaledprob( 0.999);
			col->short_transition_e[MI] = org->hmms[i]->hmm_column[j]->short_transition_e[MI];// prob2scaledprob(1.0 - 0.999) +  prob2scaledprob(0.5);
			col->short_transition_e[MD] = org->hmms[i]->hmm_column[j]->short_transition_e[MD];// prob2scaledprob(1.0 - 0.999) +  prob2scaledprob(0.5);
			
			col->short_transition_e[II] = org->hmms[i]->hmm_column[j]->short_transition_e[II];// prob2scaledprob(1.0 - 0.999);
			col->short_transition_e[IM] = org->hmms[i]->hmm_column[j]->short_transition_e[IM];// prob2scaledprob(0.999);
			
			col->short_transition_e[DD] = org->hmms[i]->hmm_column[j]->short_transition_e[DD];// prob2scaledprob(1.0 - 0.999);
			col->short_transition_e[DM] = org->hmms[i]->hmm_column[j]->short_transition_e[DM];// prob2scaledprob(0.999);
			
			
			
			for(c = 0;c <  (model->num_hmms-1);c++){
				col->long_transition[c] = org->hmms[i]->hmm_column[j]->long_transition[c];
				col->long_transition_e[c] = org->hmms[i]->hmm_column[j]->long_transition_e[c];
			}
			
			
			
		}
		
	}
	
	return model;
}


struct model* add_estimates_to_model(struct model* target, struct model* source)
{
	
	int i = 0;
	int j = 0;
	int len = 0;
	int c = 0;
	
	struct hmm_column* col_target =0;
	struct hmm_column* col_source =0;
	len = target->hmms[0]->num_columns;
	
	for(i = 0; i < len;i++){
		col_target = target->hmms[0]->hmm_column[i];
		col_source = source->hmms[0]->hmm_column[i];
		
		//col->identifier = -1;
		for(j = 0; j < 5;j++){
			col_target->i_emit_e[j] = logsum(col_target->m_emit_e[j] ,col_source->m_emit[j] ); /// org->hmms[0]->hmm_column[i]->i_emit_e[j];
			col_target->i_emit_e[j] = logsum(col_target->i_emit_e[j] ,col_source->i_emit[j] ); ///  org->hmms[0]->hmm_column[i]->m_emit_e[j];
		}
		col_target->short_transition_e[NEXT] = logsum(col_target->short_transition_e[NEXT] ,col_source->short_transition_e[NEXT] );
		col_target->short_transition_e[SELF] = logsum(col_target->short_transition_e[SELF] ,col_source->short_transition_e[SELF] );

		for(j = 0; j < target->num_hmms-1;j++){
			//col_target->long_transition[j] = org->hmms[0]->hmm_column[i]->long_transition[j];// prob2scaledprob(1.0 - max) +  prob2scaledprob(0.5  /  (float) (model->num_hmms-1));
			col_target->long_transition_e[j] = logsum(col_target->long_transition_e[j], col_source->long_transition_e[j]);// org->hmms[0]->hmm_column[i]->long_transition_e[j];
		}
		
		
	}
	
	/*
	 
	 sub models model....
	 
	 */
	len = target->hmms[1]->num_columns;
	for(i = 1; i< target->num_hmms;i++){
		for(j = 0; j < len;j++){
			col_target = target->hmms[i]->hmm_column[i];
			col_source = source->hmms[i]->hmm_column[i];
			//col = model->hmms[i]->hmm_column[j];
			//col->identifier = -1;
			for(c = 0; c < 5;c++){
				//col->i_emit[c] = org->hmms[i]->hmm_column[j]->i_emit[c];
				//col->m_emit[c] = org->hmms[i]->hmm_column[j]->m_emit[c];
				
				col_target->i_emit_e[c] = logsum(col_target->i_emit_e[c] , col_source->i_emit_e[c] ); //  org->hmms[i]->hmm_column[j]->i_emit_e[c];
				col_target->m_emit_e[c] =  logsum(col_target->m_emit_e[c] , col_source->m_emit_e[c] );/// org->hmms[i]->hmm_column[j]->m_emit_e[c];
			}
			
			
			
			
						
			col_target->short_transition_e[MM] = logsum(col_target->short_transition_e[MM] ,col_source->short_transition_e[MM] );//
			col_target->short_transition_e[MI] = logsum(col_target->short_transition_e[MI] ,col_source->short_transition_e[MI] );//
			col_target->short_transition_e[MD] = logsum(col_target->short_transition_e[MD] ,col_source->short_transition_e[MD] );//
			
			col_target->short_transition_e[II] = logsum(col_target->short_transition_e[II] ,col_source->short_transition_e[II] );//
			col_target->short_transition_e[IM] = logsum(col_target->short_transition_e[IM] ,col_source->short_transition_e[IM] );//
			
			col_target->short_transition_e[DD] = logsum(col_target->short_transition_e[DD] ,col_source->short_transition_e[DD] );//
			col_target->short_transition_e[DM] = logsum(col_target->short_transition_e[DM] ,col_source->short_transition_e[DM] );//
		
			
			
			
			for(c = 0;c <  (target->num_hmms-1);c++){

				col_target->long_transition_e[c] = logsum(col_target->long_transition_e[c],col_source->long_transition_e[c] ); //  org->hmms[i]->hmm_column[j]->long_transition_e[c];
			}
			
			
			
		}
		
	}
	return target;
}



struct model* init_model(struct model* model ,float* background,int average_length)
{
	
	int i,j,c;
	double max;
	int len = 0;
	
	int kmer_len = log(model->num_hmms-1) / log(4);
	int x;
	
	struct hmm_column* col =0;
	assert(model != 0);
	len = model->hmms[0]->num_columns;
	c = 0;
	max = -1;
	for(i = 0; i< 1000;i++){
		if(binomial_distribution((double)i / 1000.0 ,average_length ,len)  > max){
			max = binomial_distribution((double)i / 1000.0 ,average_length,len );
			c = i;
		}
	}
	max = (double)c / 1000.0;
/*
 
 Init main model....
 
 */
	
	
	for(i = 0; i < len;i++){
		col = model->hmms[0]->hmm_column[i];
		
		col->identifier = -1;
		for(j = 0; j < 5;j++){
			col->i_emit[j] = background[j];
			col->m_emit[j] = background[j];
			col->i_emit_e[j] =  prob2scaledprob(0.0f);
			col->m_emit_e[j] =  prob2scaledprob(0.0f);
		}
		col->short_transition[NEXT] = prob2scaledprob( max);
		col->short_transition[SELF] = prob2scaledprob(1.0 - max) +  prob2scaledprob(0.5);
		
		col->short_transition_e[NEXT] = prob2scaledprob( max);
		col->short_transition_e[SELF] = prob2scaledprob(1.0 - max) +  prob2scaledprob(0.5);
		
		
		for(j = 0; j< model->num_hmms-1;j++){
			col->long_transition[j] = prob2scaledprob(1.0 - max) +  prob2scaledprob(0.5  /  (float) (model->num_hmms-1));
			col->long_transition_e[j] = prob2scaledprob(1.0 - max) +  prob2scaledprob(0.5  /  (float) (model->num_hmms-1));
		}
		
		
	}
	
	/*
	 
	 sub models model....
	 
	 */
	len = model->hmms[1]->num_columns;
	for(i = 1; i< model->num_hmms;i++){
		for(j = 0; j < len;j++){
			col = model->hmms[i]->hmm_column[j];
			col->identifier = -1;
			for(c = 0; c < 5;c++){
				col->i_emit[c] = background[c];
				
				col->i_emit_e[c] =  prob2scaledprob(0.0f);
				col->m_emit_e[c] =  prob2scaledprob(0.0f);
			}
			
			/*
			
			 This horrible piece of code checks if we are at the start of the sub - model - if so it assigns a high probability to nucleotides of the #kmer given by i-1 -> the current sub-hmm (i.e. first hmm is AAAXXXXX, second AACXXXXX and so forth....
			 
			*/
			if(j < kmer_len){
				for(c = 0; c < 5;c++){
					col->m_emit[c] =prob2scaledprob(0.0001);
				}
				x = (kmer_len - 1 - j)  <<1  ;
				
				col->m_emit[  ((i-1) >> x) & 0x3] = prob2scaledprob(1.0 - (0.0001*4.0) );
			}else{
				for(c = 0; c < 5;c++){
					col->m_emit[c] = background[c];
				}
			}
			
			//for(c = 0; c < 5;c++){
			//	fprintf(stderr,"%d	%d	%f	%f	%f	%f\n",i,j, scaledprob2prob(col->m_emit[0]),  scaledprob2prob(col->m_emit[1]),  scaledprob2prob(col->m_emit[2]),  scaledprob2prob(col->m_emit[3]));
			//}
			
			
			col->short_transition[MM] = prob2scaledprob( 0.999);
			col->short_transition[MI] = prob2scaledprob(1.0 - 0.999) +  prob2scaledprob(0.5);
			col->short_transition[MD] = prob2scaledprob(1.0 - 0.999) +  prob2scaledprob(0.5);
			
			col->short_transition[II] = prob2scaledprob(1.0 - 0.999);
			col->short_transition[IM] = prob2scaledprob(0.999);
			
			col->short_transition[DD] = prob2scaledprob(1.0 - 0.999);
			col->short_transition[DM] = prob2scaledprob(0.999);
			
			
			col->short_transition_e[MM] =  prob2scaledprob(0.0);
			col->short_transition_e[MI] =  prob2scaledprob(0.0);
			col->short_transition_e[MD] =  prob2scaledprob(0.0);
			
			col->short_transition_e[II] =  prob2scaledprob(0.0);
			col->short_transition_e[IM] =  prob2scaledprob(0.0);
			
			col->short_transition_e[DD] =  prob2scaledprob(0.0);
			col->short_transition_e[DM] =  prob2scaledprob(0.0);

			
			
			
			for(c = 0;c <  (model->num_hmms-1);c++){
				col->long_transition[c] = prob2scaledprob(0.0);
				col->long_transition_e[c] = prob2scaledprob(0.0);
			}
			
			
			
		}
		
	}
	
	
	return model;
}



















