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


#define PI 3.14159265
#define MAXIT 1000
#define EPS 3.0e-7
#define FPMIN 1.0e-30
#define ITMAX 100


#define SQRT2M_PI 2.506628274631
#define MIN_STDEV 0.3989422804014

#define INV_SQRT_2PI 0.3989422804014327


#endif

int byg_end(const char* pattern,const char*text);
int byg_count(char* pattern,char*text);


void init_logsum();
float logsum(float a,float b);

float logsum_print(float a,float b);

unsigned int pop(int x);
float prob2scaledprob(float p);
float scaledprob2prob(float p);

double binomial_distribution(double p , int n, int k);
double gammln(const double xx);

int qsort_string_cmp(const void *a, const void *b);




int binsearch_down(const char*p,const char** suffix,int h,int len);
int binsearch_up(const char*p,const char** suffix,int h,int len);

int count_string(const char*p,const char** suffix,int h,int len);

double cdf(double x, double mean,double stdev);
double gammp(double a, double x);
double gammq(double a, double x);
double erffc(double x);

void gser(double *gamser, double a, double x, double *gln);
void gcf(double *gammcf, double a, double x, double *gln);

double log_pdf(double x, double mean,double stdev);
double log_truncated_pdf(double x, double mean,double stdev,double a, double b);
double gaussian_pdf(double x, double m,double s);





