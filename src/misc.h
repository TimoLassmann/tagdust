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
 along with Delve.  If not, see <http://www.gnu.org/licenses/>.
 
 */


#ifndef tagdust2_misc_h
#define tagdust2_misc_h

#include <stdio.h>
#include <string.h>
#include "math.h"

#define LOGSUM_SIZE 16000
#define SCALE 1000.0f

#define HMMER3_MIN(a,b)          (((a)<(b))?(a):(b))
#define HMMER3_MAX(a,b)          (((a)>(b))?(a):(b))
#ifdef HUGE_VAL
#define SCALEINFTY HUGE_VAL
#endif

#endif

int byg_end(const char* pattern,const char*text);



void init_logsum();
float logsum(float a,float b);

float logsum_print(float a,float b);


float prob2scaledprob(float p);
float scaledprob2prob(float p);

double binomial_distribution(double p , int n, int k);
double gammln(const double xx);

int qsort_string_cmp(const void *a, const void *b);




int binsearch_down(const char*p,const char** suffix,int h,int len);
int binsearch_up(const char*p,const char** suffix,int h,int len);

int count_string(const char*p,const char** suffix,int h,int len);
