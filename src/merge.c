#if HAVE_CONFIG_H
#include "config.h"
#endif

#include "kslib.h"

#include "interface.h"
#include "nuc_code.h"
#include "io.h"
#include "misc.h"
#include "tagdust2.h"
#include <math.h>
#include <pthread.h>


struct profiles{
	float **a;
	float **b;
	float* a_mem;
	float* b_mem;
	char* output;
	char* qual;
	int out_len;
	int malloced;
};

struct thread_data{
	struct profiles* profile;
	struct parameters* param;
	struct read_info*** ri;
	int numseq;
	int start;
	int end;
};


int merge(struct parameters* param);
int resize_profiles(struct profiles* p, int len);
int overlap_reads(struct profiles* p, struct read_info* f, struct read_info* r,int min_overlap, float threshold);
void free_profiles(struct profiles* p);

void* do_merging(void *threadarg);
int run_merge(struct parameters* param, struct read_info*** read_info_container, int numseq);

int main (int argc,char * argv[]) {
	struct parameters* param = NULL;
	init_nuc_code();
	
	param = interface(param,argc,argv);
	if(param){
		merge(param);
	}
	free_param(param);
	return EXIT_SUCCESS;
}



int merge(struct parameters* param)
{
	struct read_info*** read_info_container = NULL;
	
	int* numseqs = NULL;
	
	int status;
	

	//void* tmp = 0;
	
	int i,j,c;//,f;
	
	
	char* read_present = NULL;
	
	MMALLOC(read_present,sizeof(char)* param->infiles);
	
	for(i = 0; i < param->infiles;i++){
		read_present[i] = 0;
	}
	
	if(!param->minlen){
		param->minlen = 0;
	}
	//long long int read_present = 0;
	
	FILE** file_container = NULL;
	int (*fp)(struct read_info** ,struct parameters*,FILE* ,int* buffer_count) = NULL;
	
	
	
	MMALLOC(read_info_container, sizeof(struct read_info**) * param->infiles);
	MMALLOC(file_container,sizeof(FILE*) * param->infiles);
	
	MMALLOC(numseqs, sizeof(int) * param->infiles);
	for(i= 0; i < param->infiles;i++){
		read_info_container[i] = NULL;
		file_container[i] = NULL;
		//param->read_structures[i] = NULL;
		//param->confidence_thresholds[i] = 0.0f;
	}

	// Get ready for log space...
	init_logsum();
	
#if DEBUG
	param->num_query = 1001;
	
#else
#if  RTEST
	param->num_query = 1000;
#else
	param->num_query = 1000000;
#endif
	
#endif
	
	//malloc read_info structure (buffer for reads) AND collect basic sequence statistics from each input file...
	
	for(i =0 ; i < param->infiles;i++){
		read_info_container[i] = malloc_read_info(read_info_container[i],param->num_query);
	}
	param->read_structure = 0;
	
	for(i = 0;i < param->infiles;i++){
		file_container[i] =  io_handler(file_container[i] , i,param);
		numseqs[i] = 0;
	}
	
	/*struct profiles* profiles = NULL;
	MMALLOC(profiles, sizeof(struct profiles));
	profiles->output = NULL;
	profiles->out_len = 0;
	profiles->a_mem = NULL;
	profiles->b_mem = NULL;
	profiles->a = NULL;
	profiles->b = NULL;
	profiles->malloced = 1;*/
	
	
	//if((status = resize_profiles(profiles, 512)) != kslOK) KSLIB_XEXCEPTION_SYS(kslEWRT,"Failed to allocate profiles...");
	
	if(param->sam == 0){
		fp = &read_fasta_fastq;
	}else {
		fp = &read_sam_chunk;
	}
	
	int total_read = 0;
	while(1){
		// read batches of sequences from input file(s)
		c = 0;
		for(i = 0; i < param->infiles;i++){
			if((status = fp( read_info_container[i], param,file_container[i],&numseqs[i])) != kslOK) KSLIB_XEXCEPTION_SYS(kslEWRT,"Failed to read data chunk from file: %s", param->infile[i]);
			//numseqs[i] = fp( read_info_container[i], param,file_container[i]);
			c+=numseqs[i];
		}
		if(!c){
			break;
		}
		
		// stop if files are not sorted...
		
		//check if number of reads read from each file are identical.
		for(i = 0; i < param->infiles-1;i++){
			for(j = i +1; j < param->infiles;j++){
				if(numseqs[i] != numseqs[j]){
					sprintf(param->buffer,"Input File:%s and %s differ in number of entries.\n", param->infile[i],param->infile[j]);
					param->messages = append_message(param->messages, param->buffer);
					free_param(param);
					exit(EXIT_FAILURE);
					
				}
			}
		}
		//fprintf(stderr,"%d	%d\n",numseqs[0], li->total_read);
		
		// On first iteration compare top 1000 read names in all files...
		if(!total_read){
			for(i = 0; i < param->infiles-1;i++){
				for(j = i +1; j < param->infiles;j++){
					for(c = 0;c < HMMER3_MIN(1000,numseqs[0] );c++){
						//fprintf(stderr,"%s\n%s\n", read_info_container[i][c]->name,read_info_container[j][c]->name);
						if(compare_read_names(param,read_info_container[i][c]->name,read_info_container[j][c]->name) ){
							sprintf(param->buffer,"Files seem to contain reads in different order:\n%s\n%s\n", read_info_container[i][c]->name,read_info_container[j][c]->name);
							param->messages = append_message(param->messages, param->buffer);
							free_param(param);
							exit(EXIT_FAILURE);
						}
						//fprintf(stderr,"%s\n%s\n", read_info_container[i][c]->name,read_info_container[j][c]->name);
					}
				}
			}
		}
		
		
		run_merge(param, read_info_container, numseqs[0]);
		
		
		total_read += numseqs[0];
		fprintf(stderr,"%d\n", total_read);
		
		
	}
	
	for(i = 0; i < param->infiles;i++){
		free_read_info(read_info_container[i], param->num_query);
	}
	//free_profiles(profiles);
	MFREE(read_info_container);
	MFREE(file_container);
	MFREE(numseqs);
	return kslOK;
ERROR:
	fprintf(stderr,"%s\n",param->errmsg);
	return kslFAIL;
}

int run_merge(struct parameters* param, struct read_info*** read_info_container, int numseq)
{
	struct thread_data* thread_data = 0;
	int status;
	
	
	pthread_t threads[param->num_threads];
	pthread_attr_t attr;
	int t;
	
	int interval = 0;
	int rc;
	
	
	MMALLOC(thread_data,sizeof(struct thread_data)* param->num_threads);
	
	interval =  (int)((double)numseq /(double)param->num_threads);
	
	for(t = 0;t < param->num_threads ;t++) {
		thread_data[t].ri = read_info_container;
		thread_data[t].profile =  NULL;
		thread_data[t].param = param;
		
		MMALLOC(thread_data[t].profile, sizeof(struct profiles));
		thread_data[t].profile->output = NULL;
		thread_data[t].profile->qual = NULL;
		thread_data[t].profile->out_len = 0;
		thread_data[t].profile->a_mem = NULL;
		thread_data[t].profile->b_mem = NULL;
		thread_data[t].profile->a = NULL;
		thread_data[t].profile->b = NULL;
		thread_data[t].profile->malloced = 1;
		
		
		if((status = resize_profiles(thread_data[t].profile, 512)) != kslOK) KSLIB_XEXCEPTION_SYS(kslEWRT,"Failed to allocate profiles...");
		
		
		thread_data[t].start = t*interval;
		thread_data[t].end = t*interval + interval;
	}
	thread_data[param->num_threads-1].end = numseq;
	
	rc = pthread_attr_init(&attr);
	if(rc){
		sprintf(param->buffer,"ERROR; return code from pthread_attr_init() is %d\n", rc);
		param->messages = append_message(param->messages, param->buffer);
		
		free_param(param);
		exit(EXIT_FAILURE);
	}
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	
	for(t = 0;t < param->num_threads;t++) {
		rc = pthread_create(&threads[t], &attr, do_merging, (void *) &thread_data[t]);
	}
	
	pthread_attr_destroy(&attr);
	
	for (t = 0;t < param->num_threads;t++){
		rc = pthread_join(threads[t], NULL);
		if (rc){
			sprintf(param->buffer,"ERROR; return code from pthread_join()is %d\n", rc);
			param->messages = append_message(param->messages, param->buffer);
			free_param(param);
			exit(EXIT_FAILURE );
		}
	}
	
	for (t = 0;t < param->num_threads;t++){
		free_profiles(thread_data[t].profile);
	}
	MFREE(thread_data);
	

	
	return kslOK;
ERROR:
	return kslFAIL;
}

void* do_merging(void *threadarg)
{
	struct thread_data *data;
	data = (struct thread_data *) threadarg;
	
	struct read_info*** ri  = data->ri;
	
	struct profiles* profiles = data->profile;
	struct parameters* param = data->param;

	int start = data->start;
	int end = data->end;
	int i;
	int status;
	for(i = start; i < end;i++){
		
		ri[1][i]->seq = (char* )reverse_complement((unsigned char* ) ri[1][i]->seq, ri[1][i]->len);
		reverse_sequence(ri[1][i]->qual,ri[1][i]->len);
		if(ri[0][i]->len >= profiles->malloced){
			if((status = resize_profiles(profiles, ri[0][i]->len)) != kslOK) KSLIB_XEXCEPTION_SYS(kslEWRT,"Failed to allocate profiles...");
		}
		
		if(ri[1][i]->len >= profiles->malloced){
			if((status = resize_profiles(profiles, ri[1][i]->len)) != kslOK) KSLIB_XEXCEPTION_SYS(kslEWRT,"Failed to allocate profiles...");
		}
		
		
		
		
		if((status = overlap_reads(profiles, ri[0][i],ri[1][i],param->minlen , param->confidence_threshold)) != kslOK) KSLIB_XEXCEPTION_SYS(kslEWRT,"Failed to merge reads");
		
		if(profiles->out_len){
			fprintf(stdout,"@%s\n%s\n+\n%s\n", ri[0][i]->name,profiles->output,profiles->qual);
		}

		

	}
ERROR:
	
	
	pthread_exit((void *) 0);
}


void free_profiles(struct profiles* p)
{
	MFREE(p->a_mem);//, sizeof(float) * p->malloced * 4);
	MFREE(p->b_mem);//, sizeof(float) * p->malloced * 4);
	MFREE(p->a);//, sizeof(float*) * p->malloced);
	MFREE(p->b);//, sizeof(float*) * p->malloced);
	MFREE(p->output);
	MFREE(p);
}

int resize_profiles(struct profiles* p, int len)
{
	int i,j,status;
	
	if (len > p->malloced){
		while(len > p->malloced){
			p->malloced = p->malloced << 1;
		}
		
		if(p->a_mem == NULL){
			MMALLOC(p->a_mem, sizeof(float) * p->malloced * 4);
			MMALLOC(p->b_mem, sizeof(float) * p->malloced * 4);
			MMALLOC(p->a, sizeof(float*) * p->malloced);
			MMALLOC(p->b, sizeof(float*) * p->malloced);
			
		}else{
			MREALLOC(p->a_mem, sizeof(float) * p->malloced * 4);
			MREALLOC(p->b_mem, sizeof(float) * p->malloced * 4);
			MREALLOC(p->a, sizeof(float*) * p->malloced);
			MREALLOC(p->b, sizeof(float*) * p->malloced);
			
		}
		
		if(p->output){
			MMALLOC(p->output, sizeof(char) * 2 * p->malloced );
			MMALLOC(p->qual, sizeof(char) * 2 * p->malloced );
		}else{
			MREALLOC(p->output, sizeof(char) * 2 * p->malloced );
			MREALLOC(p->qual, sizeof(char) * 2 * p->malloced );
		}
		
		j = 0;
		for(i =0 ; i < p->malloced;i++){
			p->a[i] = p->a_mem + j;
			p->b[i] = p->b_mem + j;
			j +=  4;
		}
		
	}
	return kslOK;
ERROR:
	return kslFAIL;
}


int overlap_reads(struct profiles* p, struct read_info* f, struct read_info* r,int min_overlap, float threshold)
{
	int i,j,c;
	float** f_accuracy = p->a;
	float** r_accuracy = p->b;
	float random = 0.0;
	int len_f = f->len;
	int len_r = r->len;
	
	
	int d = 0;
	int local_i,local_j;
	int best_d = -1;
	int nuc = 0;
	float score = 0.0;
	float max_score = 0.0;
	float sum = 0.0;
	
	//float tgpe = 1.0- (1.0 /( (  ((double) len_r + len_f )/ 2.0) *(1.0 - overlap)));
	
	//fprintf(stderr,"TGPE: %f    %f\n", tgpe, (  ((double) len_r + len_f )/ 2.0) *(1.0 - overlap));
	
	//random = prob2scaledprob(1.0);
	//for(i = 0; i < len_r + len_f;i++){
	//	random = random +prob2scaledprob(tgpe);
	//}
	
	for(i = 0; i < f->len;i++){

		score =1.0 - pow(10.0, - ((int)f->qual[i] -33)   /10.0);
	//	fprintf(stderr,"%d %c ",(int)f->seq[i],f->qual[i]);
		//sum = prob2scaledprob(0.0);
		
		if(f->seq[i] > 3){
			for(j = 0;j <4;j++){
				p->a[i][j] = (0.25);
	//			fprintf(stderr," %f",scaledprob2prob(f_accuracy[i][j]) );
				//sum = logsum(sum, f_accuracy[i][j]);
			}
		}else{
		
			for(j = 0;j <4;j++){
				if(j == f->seq[i]){
					f_accuracy[i][j] = (score);
				}else{
					f_accuracy[i][j]  = ((1.0 -score) / 3.0);
				}
				//sum = logsum(sum, f_accuracy[i][j]);
	//			fprintf(stderr," %f",scaledprob2prob(f_accuracy[i][j]) );
			}
		}
	//	fprintf(stderr," %f\n",scaledprob2prob(sum));
	}
	//fprintf(stderr,"\n\n");
	for(i = 0; i < r->len;i++){
		score =1.0 - pow(10.0, - ((int)r->qual[i] -33)   /10.0);
	//	fprintf(stderr,"%d %c ",(int)r->seq[i],r->qual[i]);
		//sum = prob2scaledprob(0.0);
		if(r->seq[i] > 3){
			for(j = 0;j <4;j++){
				r_accuracy[i][j] = (0.25);
	//			fprintf(stderr," %f",scaledprob2prob(r_accuracy[i][j]) );
			//	sum = logsum(sum, r_accuracy[i][j]);
			}
		}else{
			for(j = 0;j <4;j++){
				if(j == r->seq[i]){
					r_accuracy[i][j] = (score);
				}else{
					r_accuracy[i][j]  = ((1.0 -score) / 3.0);
				}
			//	sum = logsum(sum, r_accuracy[i][j]);
	//			fprintf(stderr," %f",scaledprob2prob(r_accuracy[i][j]) );
			}
		}
	//	fprintf(stderr," %f\n",scaledprob2prob(sum));
	}

	
	max_score = prob2scaledprob(0.0);
	
	for(i = 0; i <  len_f;i++){
		local_j = 0;
		local_i = i;
		score = prob2scaledprob(1.0);
		if(len_f - local_i > min_overlap  &&len_r -local_j > min_overlap){
			/*c = 0;
			while(c != local_i){
				score = score + prob2scaledprob(tgpe);
				c++;
			}*/
			//score = score + prob2scaledprob(1.0 - tgpe);
			while(local_i != len_f && local_j != len_r){
				sum = 0.0;
				for(c = 0; c<4;c++){
					sum +=f_accuracy[local_i][c] * r_accuracy[local_j][c];
					//sum = logsum(sum,f_accuracy[local_i][c] + r_accuracy[local_j][c]);
				}
				score = score + prob2scaledprob(sum);
				local_j++;
				local_i++;
			}
			/*while(local_j < len_r ){
				score = score + prob2scaledprob(tgpe);
				local_j++;
			}
			while(local_i < len_f ){
				score = score + prob2scaledprob(tgpe);
				local_i++;
			}*/
		//	fprintf(stderr,"d:%d i:%d j:%d %f\n",d,i,0,score-random );
			if(score > max_score){
				max_score = score-random;
				best_d = d;
			}
		}
		d++;
	}
	for(j = 0;j < len_r;j++){
		
		local_j = j;
		local_i = 0;
		score = prob2scaledprob(1.0);
		if(len_f - local_i > min_overlap  &&len_r -local_j > min_overlap){
			
			/*c = 0;
			while(c != local_j){
				score = score + prob2scaledprob(tgpe);
				c++;
			}*/
			
			while(local_i != len_f && local_j != len_r){
				sum = 0.0;
				for(c = 0; c<4;c++){
					sum += f_accuracy[local_i][c] * r_accuracy[local_j][c];
					//sum = logsum(sum,f_accuracy[local_i][c] + r_accuracy[local_j][c]);
				}
				score = score + prob2scaledprob(sum);
				local_j++;
				local_i++;
			}
			/*while(local_j < len_r ){
				score = score + prob2scaledprob(tgpe);
				local_j++;
			}
			while(local_i < len_f ){
				score = score + prob2scaledprob(tgpe);
				local_i++;
			}*/

			
			
		//	fprintf(stderr,"d:%d i:%d j:%d %f\n",d,0,j,score-random );
			if(score > max_score){
				max_score = score-random;
				best_d = d;
			}
		}
		d++;
	}
	
	
	d = 0;
	
	float id;
	
	float aligned;
	
	if(best_d < len_f){
		local_j = 0;
		local_i = best_d;
		for(i = 0; i < local_i;i++) {
			
			p->output[d] = "ACGTC"[(int)f->seq[i]];
			p->qual[d] = f->qual[i];
			d++;
			//fprintf(stderr,"%c-\n", f->seq[i]+65);
		}
		id = 0.0;
		aligned = 0.0;
		while(local_i != len_f && local_j != len_r){
			if( f->seq[local_i] == r->seq[local_j]){
				p->output[d] = "ACGTC"[(int)f->seq[local_i]];
				
				id += 1.0;
			}else{
				best_d = 0;
				max_score = prob2scaledprob(0.0);
				for(c = 0;c < 4;c++){
					if(p->a[local_i][c] > max_score){
						max_score =p->a[local_i][c];
						nuc  = c;
					}
					if(p->b[local_j][c] > max_score){
						max_score = p->b[local_j][c];
						nuc  = c;
					}
					
				}
				p->output[d] = "ACGTC"[(int)nuc];
			}
			p->qual[d] =  (f->qual[local_i] > r->qual[local_j]) ? f->qual[local_i] : r->qual[local_j];
			d++;
			aligned+= 1.0;
			//fprintf(stderr,"%c%c\n", f->seq[local_i]+65, r->seq[local_j]+65);
			local_j++;
			local_i++;
		}
		for(i = local_i; i < len_f;i++) {
			//fprintf(stderr,"%c-\n", f->seq[i]+65);
			p->output[d] = "ACGTC"[(int)f->seq[i]];
			p->qual[d] = f->qual[i];
			d++;
		}
		for(j = local_j; j < len_r;j++) {
			//fprintf(stderr,"-%c\n", r->seq[j]+65);
			p->output[d] = "ACGTC"[(int)r->seq[j]];
			p->qual[d] = r->qual[j];
			d++;
		}
		
		
		
		
	}else{
		local_j = best_d - len_f;
		local_i = 0;
		
		for(j = 0; j < local_j;j++) {
			p->output[d] = "ACGTC"[(int)r->seq[j]];
			p->qual[d] = r->qual[j];
			d++;
	//		fprintf(stderr,"-%c\n", r->seq[j]+65);
		}
		id = 0.0;
		aligned = 0.0;
		while(local_i != len_f && local_j != len_r){
			if( f->seq[local_i] == r->seq[local_j]){
				p->output[d] = "ACGTC"[(int)f->seq[local_i]];
				id += 1.0;
			}else{
				best_d = 0;
				max_score = prob2scaledprob(0.0);
				for(c = 0;c < 4;c++){
					if(p->a[local_i][c] > max_score){
						max_score =p->a[local_i][c];
						nuc  = c;
					}
					if(p->b[local_j][c] > max_score){
						max_score = p->b[local_j][c];
						nuc  = c;
					}
					
				}
				p->output[d] = "ACGTC"[(int)nuc];

			}
			aligned+= 1.0;
			p->qual[d] =  (f->qual[local_i] > r->qual[local_j]) ? f->qual[local_i] : r->qual[local_j];
			d++;
	//		fprintf(stderr,"%c%c\n", f->seq[local_i]+65, r->seq[local_j]+65);
			local_j++;
			local_i++;
		}
		for(i = local_i; i < len_f;i++) {
	//		fprintf(stderr,"%c-\n", f->seq[i]+65);
			p->output[d] = "ACGTC"[(int)f->seq[i]];
			p->qual[d] = f->qual[i];
			d++;
		}
		for(j = local_j; j < len_r;j++) {
	//		fprintf(stderr,"-%c\n", r->seq[j]+65);
			p->output[d] = "ACGTC"[(int)r->seq[j]];
			p->qual[d] = r->qual[j];
			d++;
		}
		
	}
	p->output[d]  = 0;
	p->qual[d]  = 0;
	//fprintf(stderr,"ID: %f\n", id / aligned *100);
	
	if(id / aligned >= threshold){
		p->out_len = d;
	}else{
		p->out_len = 0;
	}
	return kslOK;

}


