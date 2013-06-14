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




#include "misc.h"
#include <ctype.h>

float logsum_lookup[LOGSUM_SIZE];

void init_logsum()
{
	int i;
	for(i = 0; i < LOGSUM_SIZE;i++){
		//logsum_lookup[i] = SCALE * 1.442695041f * (log(1.0f + exp(0.693147181f*(float)-i/SCALE)));
		logsum_lookup[i] =  log(1. + exp((double) -i / SCALE));
	}
}

float logsum(float a,float b)
{
	
	const float max = HMMER3_MAX(a, b);
	const float min = HMMER3_MIN(a, b);
	return (min == -SCALEINFTY || (max-min) >= 15.7f) ? max : max + logsum_lookup[(int)((max-min)*SCALE)];
}



float prob2scaledprob(float p)
{
	if(p == 0.0){
		return -SCALEINFTY;
	}else{
		return  log(p);
		//return SCALE * sreLOG2(p);
	}
}


float scaledprob2prob(float p)
{
	if(p == -SCALEINFTY){
		return 0.0;
	}else{
		return exp(p);
		//return sreEXP2(p / SCALE);
	}
}


int byg_count(char* pattern,char*text)
{
	int Tc;
	int i  = 0;
	int s = 0;
	int T[256];
	for (i = 0;i < 256;i++){
		T[i] = 0;
	}
	
	int m = (int) strlen(pattern);
	int n = (int) strlen(text);
	int count = 0;
	
	if(m > n){
		return -1;
	}
	int mb = (1 << (m-1));
	
	for (i= 0;i < m;i++){
		T[(int)toupper(pattern[i])] |= (1 << i);
	}
	
	for (i = 0;i < n;i++){
		s <<= 1;
		s |= 1;
		Tc = T[(int)toupper(text[i])];
		s &= Tc;
		if(s & mb){
			count++;
		}
	}
	return count;
}



int byg_end(const char* pattern,const char*text)
{
	const char* tmp = 0;
	int Tc;
	int i  = 0;
	int s = 0;
	int T[256];
	for (i = 0;i < 256;i++){
		T[i] = 0;
	}
	
	int m = (int)strlen(pattern);
	int n = (int)strlen(text);
	if (m > n){
		i = m;
		m = n;
		n = i;
		tmp = text;
		text = pattern;
		pattern = tmp;
	}
	
	int mb = (1 << (m-1));
	
	for (i= 0;i < m;i++){
		T[(int)pattern[i]] |= (1 << i);
	}
	
	for (i = 0;i < n;i++){
		s <<= 1;
		s |= 1;
		if(!text[i]){
			return -1;
		}
		Tc = T[(int)text[i]];
		s &= Tc;
		if(s & mb){
			return i+1;
		}
	}
	return 0;
}



int count_string(const char*p,const char** suffix,int h,int len)
{
	int a,b;
	//for(i = 0; i < 1000000;i++){
	a = binsearch_down(p,suffix,h,len);
	b = binsearch_up(p,suffix,h,len);
	return b-a;
}


int binsearch_down(const char*p,const char** suffix,int h,int len)
{
	int m = 0;
	int l = 0;
	/*if (t_long_strncmp(p,text+suffix[l],len)<= 0){
	 l = l;
	 }else */
	if(strncmp(p,suffix[h],len) >  0){
		return h;
	}else{
		while(h-l > 1){
			//m = (l+h)/2;
			m = (l + h) >> 1;
			if(strncmp(p,suffix[m],len) <= 0){
				h = m;
			}else{
				l = m;
			}
		}
	}
	return l+1;
}

int binsearch_up(const char*p,const char** suffix,int h,int len)
{
	int m = 0;
	int l = 0;
	/*if (t_long_strncmp(p,text+suffix[l],len)<= 0){
	 l = l;
	 }else*/
	if(strncmp(p,suffix[h],len) >  0){
		return h;
	}else{
		while(h-l > 1){
			//m = (l+h)/2;
			m = (l + h) >> 1;
			if(strncmp(p,suffix[m],len) < 0){
				h = m;
			}else{
				l = m;
			}
		}
	}
	return l+1;
}



int bindoublesearch_up(double x ,double* y,int h)
{
	int m = 0;
	int l = 0;
	/*if (t_long_strncmp(p,text+suffix[l],len)<= 0){
	 l = l;
	 }else*/
	if(x < y[h]){
		return h+1;
	}else{
		while(h-l > 1){
			//m = (l+h)/2;
			m = (l + h) >> 1;
			
			if(x <= y[m]){
				l =m;
			}else{
				h = m;
			}
		}
	}
	return l+1;
}



int qsort_string_cmp(const void *a, const void *b)
{
	const char **one = (const char **)a;
	const char **two = (const char **)b;
	return strcmp(*one, *two);
}

int qsort_flt_cmp(const void * a, const void * b)
{
	//const float a  = (float) *elem1;
	if(*(const float*)a > *(const float*)b)
		return -1;
	return *(const float*)a < *(const float*)b;
	
	//return (*(float*) a) - (*(float*) b);
}



double gaussian_pdf(double x, double m,double s)
{
	double a = (x-m) / s;
	return INV_SQRT_2PI / s *exp(-0.5 * a * a);
}

double log_pdf(double x, double mean,double stdev)
{
	double out;
	
	//out = 1.0 / sqrt (2 * M_PI * pow( stdev,2) ) * exp(-1.0 *( (pow(x - mean,2)) / ( 2* pow( stdev,2)   )));
	
	
	out = log(1.0 / (stdev * SQRT2M_PI)) +  (-1.0 *( (x - mean) * (x - mean) / ( 2.0 * stdev * stdev)   ));
	//if(out < FLT_MIN){
	//	return FLT_MIN;
	//}
	
	
	if(isnan(out)){
		fprintf(stderr,"logpdf problem.... \n");
	}
	
	return out;
}


unsigned int pop(int x)
{
	unsigned int n = 0;
	while(x != 0){
		n = n +1;
		x = x &(x-1);
	}
	return n;
}

int bpm(const  char* t,const  char* p,int n,int m)
{
	register unsigned long int i;//,c;
	unsigned long int diff;
	unsigned long int B[255];
	if(m > 31){
		m = 31;
	}
	
	unsigned long int k = m;
	//static int counter = 0;
	register unsigned long int VP,VN,D0,HN,HP,X;
	
	long int MASK = 0;
	//c = 0;
	
	diff = m;
	
	for(i = 0; i < 255;i++){
		B[i] = 0;
	}
	
	for(i = 0; i < m;i++){
		B[(int)(p[i] )] |= (1ul << i);
	}
	
	//c = 0;
	VP = 0xFFFFFFFFFFFFFFFFul;
	VN = 0ul;
	m--;
	MASK = 1ul << m;
	for(i = 0; i < n;i++){
		X = (B[(int)(t[i])  ] | VN);
		D0 = ((VP+(X&VP)) ^ VP) | X ;
		HN = VP & D0;
		HP = VN | ~(VP | D0);
		X = HP << 1ul;
		VN = X & D0;
		VP = (HN << 1ul) | ~(X | D0);
		diff += (HP & MASK) >> m;
		diff -= (HN & MASK) >> m;
		if(diff < k){
			k = diff;
			//fprintf(stderr,"%ld	%ld\n",i,k);
			//if(k <= limit){
			//	return (int)k;
			//}
			
		}
	}
	return (int)k;
}
