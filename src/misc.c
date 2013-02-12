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



#include "misc.h"


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




double binomial_distribution(double p , int n, int k)
{
	double fac = gammln((double)n + 1.0);
	if(k < 0){
		//fprintf(stderr,"Bad k in binomialdist\n");
	}
	if(k > n){
		return 0;
	}
	return exp(k*log(p) +(n-k)*log(1.0 - p) + fac -gammln(k + 1.0) - gammln(n-k+1.0));
	
}



double gammln(const double xx)
{
	int j;
	double x,tmp,y,ser;
	static const double cof[14] = {57.1562356658629235,-59.5979603554754912,14.1360979747417471, -0.491913816097620199,0.339946499848118887e-4,0.465236289270485756e-4, -0.983744753048795646e-4, 0.0158088703224912494e-3,-0.210264441724104883e-3, 0.217439618115212643e-3, -0.164318106536763890e-3,0.888182239838527433e-4, -0.261908384015814087e-4, 0.368991826595316234e-5};
	if(xx <= 0.0){
		//fprintf(stderr,"bar arg in gammln");
	}
	y = xx;
	x = xx;
	tmp = x+5.24218750000000000;
	tmp = (x + 0.5) * log(tmp) - tmp;
	ser = 0.999999999999997092;
	for(j = 0; j < 14;j++){
		ser += cof[j] / ++y;
	}
	return tmp + log(2.5066282746310005*ser/x);
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



int qsort_string_cmp(const void *a, const void *b)
{
	const char **one = (const char **)a;
	const char **two = (const char **)b;
	return strcmp(*one, *two);
}










