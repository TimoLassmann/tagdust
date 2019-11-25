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

/*! \file misc.h
 \brief Header for misc.c 
 */



#ifndef tagdust2_misc_h
#define tagdust2_misc_h

#include <stdio.h>
#include <emmintrin.h>
#include <mmintrin.h>
#include <stdlib.h>
#include <string.h>
#include "math.h"



/** \def LOGSUM_SIZE
 \brief Precision of the lookup table.
 
 */
#define LOGSUM_SIZE 16000

/** \def SCALE
 \brief defines the precision of the calculation; the
 default of 1000.0 means rounding differences to the nearest 0.001
  nat
 
 */
#define SCALE 1000.0f

/** \def HMMER3_MIN(a,b) 
 \brief Computes the minimum of \a a and \a b.
 
 */

#define HMMER3_MIN(a,b)          (((a)<(b))?(a):(b))

/** \def HMMER3_MAX(a,b)
 \brief Computes the maximum of \a a and \a b.
 
 */
#define HMMER3_MAX(a,b)          (((a)>(b))?(a):(b))
#ifdef HUGE_VAL
#define SCALEINFTY HUGE_VAL
#endif

/** \def SQRT2M_PI
 \brief Square root of 2* pi .
 
 */
#define SQRT2M_PI 2.506628274631


#define INV_SQRT_2PI 0.3989422804014327


#ifndef _MM_ALIGN16
#ifdef __GNUC__
#define _MM_ALIGN16 __attribute__((aligned (16)))
#endif
#ifdef __MSVC__
#define _MM_ALIGN16 __declspec(align(16))
#endif
#endif


#endif

#ifdef RTEST
#define rand() myrand()
#define srand(x) mysrand(x)

#endif



int myrand(void);
void mysrand(unsigned seed);

int byg_end(const char* pattern,const char*text);

int byg_count(char* pattern,char*text);

char* append_message(char* old_message, char* new_message);
void init_logsum();
float logsum(float a,float b);

float logsum_print(float a,float b);
unsigned int pop(int x);
float prob2scaledprob(float p);
float scaledprob2prob(float p);

int qsort_string_cmp(const void *a, const void *b);
int qsort_flt_cmp(const void * a, const void * b);

int bpm(const char* t,const char* p,int n,int m);
int bpm_global(const  char* t,const  char* p,int n,int m);


int binsearch_down(const char*p,const char** suffix,int h,int len);
int binsearch_up(const char*p,const char** suffix,int h,int len);

int count_string(const char*p,const char** suffix,int h,int len);

double gaussian_pdf(double x, double m,double s);

int validate_bpm_sse(unsigned char**  query, int* query_lengths,unsigned char* t,int n,int num);
int bpm_check_error(const unsigned char* t,const unsigned char* p,int n,int m);
int bpm_check_error_global(const unsigned char* t,const unsigned char* p,int n,int m);
unsigned char* reverse_complement(unsigned char* p,int len);
void reverse_sequence(char* p,int len);
char* shorten_pathname(char* p);
int file_exists(char* name);
int bitcount64(long long int i);
int highest_bit(long long int n);
