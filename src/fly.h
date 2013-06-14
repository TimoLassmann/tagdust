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
	double** fly;
	double* light;
	int* moved;
	int num_fly;
	double best_light;
	double* best;
};


struct fly_thread_data{
	struct swarm* swarm;
	double* data; //model_bag* mb;
	int n;
	int num_var;

	double (*fp)(double*,int ,  double* ,int );
	int thread_num;
	unsigned int seed;
	//float* best;
};


void flytest();
double* firefly(double* data, int n,int num_var,double* max,double* min, double (*fp)(double*,int ,  double* ,int ));
double*  run_firefly_thread(double* data,int n,int num_var,double* max, double* min, double (*fp)(double*,int ,  double* ,int ) ,int num_threads);
void* do_fly(void *threadarg);
double extreme_value_distribution_eval(double*data,int n, double* var, int num_var);


struct swarm* init_swarm(int num_fly, int num_var);
void free_swarm(struct swarm* swarm );

#endif
