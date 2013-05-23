//
//  fly.c
//  tagdust2
//
//  Created by lassmann on 5/21/13.
//  Copyright (c) 2013 lassmann. All rights reserved.
//

#include <stdio.h>

#include "misc.h"
#include <time.h>
#include <pthread.h>
#include "fly.h"

void flytest()
{
	float* data = 0;
	float* min = 0;
	float* max = 0;
	float* out = 0;
	//float r;
	int i;
	unsigned int seed = (unsigned int) (time(NULL) * ( 42));
	
	//r = (float)rand_r(&seed)/(float)RAND_MAX;
	
	data = malloc(sizeof(float) * 5000);
	
	min = malloc(sizeof(float) * 2);
	max = malloc(sizeof(float) * 2);
	
	min[0] = -20;
	min[1] = 0;
	
	max[0] = 20;
	max[1] = 20;
	
	for(i = 0; i < 5000;i++){
		data[i] = 5+ (0.5 - (float)rand_r(&seed)/(float)RAND_MAX)* 0.1;
	}
	
	//data = firefly(data,5000,2,max,min, &extreme_value_distribution_eval);
	out = run_firefly_thread(data,5000,2,max,min, &extreme_value_distribution_eval,80);
	
	for(i = 0; i < 2;i++){
		fprintf(stderr,"%f ",out[i]);
	}
	fprintf(stderr,"\n");
	for(i = 0; i < 5000;i++){
		fprintf(stdout," %f	%f	%f\n", data[i] , exp(-1 * exp(-1 * out[EVD_lambda] *(data[i] - out[ EVD_mu]))),   out[EVD_lambda]  * exp(-1.0f *out[EVD_lambda] * (data[i] - out[ EVD_mu]) - exp(-1.0f * out[EVD_lambda]  *  (data[i] - out[ EVD_mu]))  )); // - exp(
	}
	
	
	exit(0);
}



float*  run_firefly_thread(float* data,int n,int num_var,float* max, float* min, float (*fp)(float*,int ,  float* ,int ) ,int num_threads)
{
	struct fly_thread_data* thread_data = 0;
	thread_data = malloc(sizeof(struct fly_thread_data)* num_threads);
	pthread_t threads[num_threads];
	pthread_attr_t attr;
	int t,i,j;
	//int interval = 0;
	int rc;
	int num_fly = 20;
	
	float* result  = 0;
	float best_light = 0;
	
	
	result = malloc(sizeof(float) * num_var);
	
	//unsigned int seed = (unsigned int) (time(NULL) * ( 42));
	
	//interval =  (int)((double)numseq /(double)param->num_threads);
	
	for(t = 0;t < num_threads ;t++) {
		thread_data[t].swarm = init_swarm(num_fly,num_var );
		thread_data[t].data = data; //model_bag* mb;
		thread_data[t].n = n;
		thread_data[t].num_var = num_var;
		thread_data[t].fp = fp;
		thread_data[t].thread_num = t;
		thread_data[t].seed = (unsigned int) (time(NULL) * ( 42 * t));

		//thread_data[t].best = malloc(sizeof(float) * num_var);
		
		
		//init flies....
		
		for(i = 0; i < thread_data[t].swarm->num_fly;i++){
			//flies[i] = malloc(sizeof(float) * num_var);
			for (j = 0; j < num_var;j++){
				thread_data[t].swarm->fly[i][j] =  ( (float)rand_r(&thread_data[t].seed)/(float)RAND_MAX * (max[j] - min[j])) + min[j];
			}
		}
		
	}
	
	
	
	
	//fprintf(stderr,"got here...\n");
	//exit(0);
	rc = pthread_attr_init(&attr);
	if(rc){
		fprintf(stderr,"ERROR; return code from pthread_attr_init() is %d\n", rc);
		exit(-1);
	}
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	
	for(t = 0;t < num_threads;t++) {
		
		rc = pthread_create(&threads[t], &attr, do_fly , (void *) &thread_data[t]);
		
		if (rc) {
			fprintf(stderr,"ERROR; return code from pthread_create() is %d\n", rc);
			exit(-1);
		}
	}
	
	pthread_attr_destroy(&attr);
	
	for (t = 0;t < num_threads;t++){
		rc = pthread_join(threads[t], NULL);
		if (rc){
			fprintf(stderr,"ERROR; return code from pthread_join()is %d\n", rc);
			exit(-1);
		}
	}
	
	
	best_light = -FLT_MAX;
	for (t = 0;t < num_threads;t++){
		
		fprintf(stderr,"%d %f ",t,thread_data[t].swarm->best_light );
		for(i = 0; i < num_var;i++){
			fprintf(stderr,"%f ",thread_data[t].swarm->best[i] );
		}
		fprintf(stderr,"\n");
		
		if(thread_data[t].swarm->best_light  > best_light){
			best_light = thread_data[t].swarm->best_light ;
			for(i = 0; i < num_var;i++){
				result[i]  =thread_data[t].swarm->best[i];
			}
			
			
		}
		
	}
	
	fprintf(stderr,"got here...\n");
	
	for(t = 0;t < num_threads;t++) {
		free_swarm(thread_data[t].swarm );
	}
	
	free(thread_data);
	
	return result;
	//return (&-1.0);
}



void* do_fly(void *threadarg)
{
	struct fly_thread_data *thread_data;
	thread_data = (struct fly_thread_data *) threadarg;
	float* data = thread_data->data;
	
	//float* best = thread_data->best;
	//float* max = thread_data->max;
	//float* min = thread_data->min;
	int num_var = thread_data->num_var;
	int num_data_points = thread_data->n;
	
	struct swarm* swarm = thread_data->swarm;
	//float*
	
	
	unsigned int seed = thread_data->seed;
	

	int i,j,c,f;
	
	int max_gen = 1000;
	float d;
	float lambda = 1.0; // absorption coefficient.
	float alpha = 0.2;
	float beta0 = 1.0;

	for(i = 0; i < swarm->num_fly;i++){
		//light[i] =  fp(data,n,flies[i],num_var );
		swarm->light[i] =  thread_data->fp(data,num_data_points,swarm->fly[i],num_var );
		
		if(swarm->light[i] > swarm->best_light){
			swarm->best_light = swarm->light[i];
			
			for(j = 0; j < num_var;j++){
				swarm->best[j] = swarm->fly[i][j];
			}
			
		}
		
	}
	
	for(i = 0; i < max_gen;i++){
		/*fprintf(stderr,"GENERATION:%d\n",i);
		for(j = 0; j < swarm->num_fly;j++){
			//light[i] =  fp(data,n,flies[i],num_var );
			fprintf(stderr,"%d %f	",j, swarm->light[j]);
			for (f = 0; f < num_var;f++){
				fprintf(stderr,"%f ",swarm->fly[j][f]);
			}
			fprintf(stderr,"\n");
		}*/
		for(j = 0; j < swarm->num_fly;j++){
			swarm->moved[j] = 0;
		}
		
		for(j = 0; j < swarm->num_fly;j++){
			for(c = 0; c < swarm->num_fly;c++){
				if(swarm->light[j] >= swarm->light[c]){
					swarm->moved[c] = 1;
					d = 0;
					for(f = 0;f < num_var;f++){
						d+= (swarm->fly[c][f] - swarm->fly[j][f]) * (swarm->fly[c][f] - swarm->fly[j][f]) ;
					}
					//d = sqrtf(d);
					for(f = 0;f < num_var;f++){
						swarm->fly[c][f] = swarm->fly[c][f] + beta0* exp(-1.0 * lambda * d ) * (swarm->fly[j][f] - swarm->fly[c][f]) + alpha * (0.5 -  (float)rand_r(&seed)/(float)RAND_MAX);
					}
					
					swarm->light[c] =  thread_data->fp(data,num_data_points,swarm->fly[c],num_var );
					if(swarm->light[c] > swarm->best_light){
						swarm->best_light = swarm->light[c];
						
						for(f = 0; f < num_var;f++){
							swarm->best[f] = swarm->fly[c][f];
						}
						
					}
				}
			}
		}
		// move not moved flies randomly...
		for(c= 0; c < swarm->num_fly;c++){
			if(!swarm->moved[c]){
				for(f = 0;f < num_var;f++){
					swarm->fly[c][f] = swarm->fly[c][f] + alpha * (0.5 -  (float)rand_r(&seed)/(float)RAND_MAX);
				}
				
				swarm->light[c] =  thread_data->fp(data,num_data_points,swarm->fly[c],num_var );
				if(swarm->light[c] > swarm->best_light){
					swarm->best_light = swarm->light[c];
					
					for(f = 0; f < num_var;f++){
						swarm->best[f] = swarm->fly[c][f];
					}
					
				}

			}
		}
		
		
		//if(i % 10 == 0){
		//	fprintf(stderr,"Thread %d Best: %f\n", thread_data->thread_num,swarm->best_light );
			//fprintf(stderr,"%")
		//}
	}
	

	pthread_exit((void *) 0);
}


float* firefly(float* data, int n,int num_var,float* max,float* min, float (*fp)(float*,int ,  float* ,int ))
{
	float* var = 0;
	
	float** flies = 0;
	float* light = 0;

	float r = 0;
	int i,j,c,f;
	
	int max_gen = 1000;
	int num_fly = 20;
		
	float d;
	float lambda = 1.0; // absorption coefficient.
	float alpha = 0.2;
	float beta0 = 1.0;
	unsigned int seed = (unsigned int) (time(NULL) * ( 42));
	
	r = (float)rand_r(&seed)/(float)RAND_MAX;
	
	flies = malloc(sizeof(float*) * num_fly);
	light = malloc(sizeof(float) * num_fly);
	for(i = 0; i < num_fly;i++){
		flies[i] = malloc(sizeof(float) * num_var);
		for (j = 0; j < num_var;j++){
			flies[i][j] =  ( (float)rand_r(&seed)/(float)RAND_MAX * (max[j] - min[j])) + min[j];
		}
	}
	
	for(i = 0; i < num_fly;i++){
		light[i] =  fp(data,n,flies[i],num_var );
		fprintf(stderr,"%d %f	",i, light[i]);
		for (j = 0; j < num_var;j++){
			fprintf(stderr,"%f ",flies[i][j]);
		}
		fprintf(stderr,"\n");
	}
	
	for(i = 0; i < max_gen;i++){
		fprintf(stderr,"GENERATION:%d\n",i);
		for(j = 0; j < num_fly;j++){
			//light[i] =  fp(data,n,flies[i],num_var );
			fprintf(stderr,"%d %f	",j, light[j]);
			for (f = 0; f < num_var;f++){
				fprintf(stderr,"%f ",flies[j][f]);
			}
			fprintf(stderr,"\n");
		}
		for(j = 0; j < num_fly;j++){
			for(c = 0; c < num_fly;c++){
				if(light[j] >= light[c]){
					d = 0;
					for(f = 0;f < num_var;f++){
						d+= (flies[c][f] - flies[j][f]) * (flies[c][f] - flies[j][f]) ;
					}
					//d = sqrtf(d);
					for(f = 0;f < num_var;f++){
						flies[c][f] = flies[c][f] + beta0* exp(-1.0 * lambda * d ) * (flies[j][f] - flies[c][f]) + alpha * (0.5 -  (float)rand_r(&seed)/(float)RAND_MAX);
					}
					
					light[c] =  fp(data,n,flies[c],num_var );
				}
			}
		}
	}
	
	
	
	
	return var;
}



float extreme_value_distribution_eval(float*data,int n, float* var, int num_var)
{
	float sum = prob2scaledprob(1.0);
	int i;
	float mu = var[EVD_mu];
	float lambda = var[EVD_lambda];
	float px = 0;
	for(i = 0; i < n;i++){
		px = lambda * exp(-1.0f *lambda* (data[i] - mu) - exp(-1.0f * lambda *  (data[i] - mu))  );
		if(px <0){
			px = 0.0;
		}
		//if(px > 1){
		//	fprintf(stderr," Huh?? %f %f %f %f\n", px, lambda,mu,data[i] );
		//	px = 0.0;
		//}
		
		//fprintf(stderr,"%f\n", px);
		sum += prob2scaledprob(px);
	}
	
	return sum;
	
}


struct swarm* init_swarm(int num_fly, int num_var)
{
	struct swarm* swarm = 0;
	int i,j;
	swarm = malloc(sizeof(struct swarm));
	
	assert(swarm != 0);
	
	swarm->fly = malloc(sizeof(float*) * num_fly);
	assert(swarm->fly  != 0);
	swarm->light = malloc(sizeof(float) * num_fly);
	assert(swarm->light  != 0);
	
	swarm->moved = malloc(sizeof(int) * num_fly);
	assert(swarm->light  != 0);
	
	swarm->num_fly = num_fly;
	swarm->best_light = -FLT_MAX;
	swarm->best = malloc(sizeof(float) * num_var);
	for(j = 0; j < num_var;j++){
		swarm->best[j] = 0.0f;
	}
	
	
	for(i = 0; i < num_fly;i++){
		swarm->light[i] = 0.0f;
		swarm->moved[i] = 0;
		swarm->fly[i] = malloc(sizeof(float) * num_var);
		for(j = 0; j < num_var;j++){
			swarm->fly[i][j] = 0.0f;
		}
	}
	
	return swarm;
}

void free_swarm(struct swarm* swarm )
{
	int i;
	for(i = 0; i < swarm->num_fly;i++){
		free(swarm->fly[i]);
	}
	free(swarm->best );
	free(swarm->moved);
	free(swarm->fly);
	free(swarm->light);
	free(swarm);
	
}







