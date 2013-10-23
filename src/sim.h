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
 along with Tagdust.  If not, see <http://www.gnu.org/licenses/>.
 
 */


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
