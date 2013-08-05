//
//  exact.h
//  tagdust2
//
//  Created by lassmann on 5/14/13.
//  Copyright (c) 2013 lassmann. All rights reserved.
//


/*! \file exact.h
 \brief Functions to extract sequences without using errors.

 \author Timo Lassmann
 \bug No known bugs.
 */

#ifndef tagdust2_exact_h
#define tagdust2_exact_h
void exact_controller(struct parameters* param,int (*fp)(struct read_info** ,struct parameters*,FILE* ),int file_num);

int byg_end_barcode(const char* pattern,const char*text, int m, int n);// m is pattern length , n is text length....;
void print_seq_from_position_x(struct read_info* ri,FILE* out,int x);
#endif
