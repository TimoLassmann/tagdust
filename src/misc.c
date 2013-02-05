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



