//
//  sim.h
//  tagdust2
//
//  Created by lassmann on 5/14/13.
//  Copyright (c) 2013 lassmann. All rights reserved.
//


/*! \file sim.h
 \brief Functions to simulate sequences.
 Initializes nucleotide arrays.
 \author Timo Lassmann
 \bug No known bugs.
 */

#ifndef tagdust2_sim_h
#define tagdust2_sim_h


struct eval_results{
	double average_read_similarity;
	int num_extracted;
	int num_rand_extracted;
	int num_wrong_bc;
};

void simulate(struct parameters* param);
void simulation_for_benchmark(struct parameters* param);
struct eval_results* get_results(struct eval_results* eval,struct parameters* param, char* filename,char* program );


#endif
