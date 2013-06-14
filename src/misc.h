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


#ifndef tagdust2_misc_h
#define tagdust2_misc_h

#include <stdio.h>

#include <stdlib.h>
#include <string.h>
#include "math.h"


#define ALPHABET_LEN 255
#define max(a, b) ((a < b) ? b : a)

#define LOGSUM_SIZE 16000
#define SCALE 1000.0f

#define HMMER3_MIN(a,b)          (((a)<(b))?(a):(b))
#define HMMER3_MAX(a,b)          (((a)>(b))?(a):(b))
#ifdef HUGE_VAL
#define SCALEINFTY HUGE_VAL
#endif

#define SQRT2M_PI 2.506628274631
#define INV_SQRT_2PI 0.3989422804014327


#endif

int byg_end(const char* pattern,const char*text);

int byg_count(char* pattern,char*text);


void init_logsum();
float logsum(float a,float b);

float logsum_print(float a,float b);
int bindoublesearch_up(double x ,double* y,int h);
unsigned int pop(int x);
float prob2scaledprob(float p);
float scaledprob2prob(float p);

int qsort_string_cmp(const void *a, const void *b);
int qsort_flt_cmp(const void * a, const void * b);

int bpm(const char* t,const char* p,int n,int m);


int binsearch_down(const char*p,const char** suffix,int h,int len);
int binsearch_up(const char*p,const char** suffix,int h,int len);

int count_string(const char*p,const char** suffix,int h,int len);

double log_pdf(double x, double mean,double stdev);
double gaussian_pdf(double x, double m,double s);





