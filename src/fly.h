//
//  fly.h
//  tagdust2
//
//  Created by lassmann on 5/21/13.
//  Copyright (c) 2013 lassmann. All rights reserved.
//

#ifndef tagdust2_fly_h
#define tagdust2_fly_h


#define EVD_lambda 1 
#define EVD_mu 0

#include <assert.h>
#include <float.h>

struct swarm{
	float** fly;
	float* light;
	int* moved;
	int num_fly;
	float best_light;
	float* best;
};


struct fly_thread_data{
	struct swarm* swarm;
	float* data; //model_bag* mb;
	int n;
	int num_var;

	float (*fp)(float*,int ,  float* ,int );
	int thread_num;
	unsigned int seed;
	//float* best;
};


void flytest();
float* firefly(float* data, int n,int num_var,float* max,float* min, float (*fp)(float*,int ,  float* ,int ));
float*  run_firefly_thread(float* data,int n,int num_var,float* max, float* min, float (*fp)(float*,int ,  float* ,int ) ,int num_threads);
void* do_fly(void *threadarg);
float extreme_value_distribution_eval(float*data,int n, float* var, int num_var);


struct swarm* init_swarm(int num_fly, int num_var);
void free_swarm(struct swarm* swarm );

#endif
