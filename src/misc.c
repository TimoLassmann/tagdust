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


#include "tagdust2.h"
#include "nuc_code.h"
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



int bpm_check_error(const unsigned char* t,const unsigned char* p,int n,int m,int limit)
{
	register unsigned long int i;//,c;
	unsigned long int diff;
	unsigned long int B[5];
	if(m > 31){
		m = 31;
	}
	
	unsigned long int k = m;
	//static int counter = 0;
	register unsigned long int VP,VN,D0,HN,HP,X;
	
	long int MASK = 0;
	//c = 0;
	
	diff = m;
	
	for(i = 0; i < 5;i++){
		B[i] = 0;
	}
	
	for(i = 0; i < m;i++){
		B[(int)(p[i] & 0x3)] |= (1ul << i);
	}
	
	//c = 0;
	VP = 0xFFFFFFFFFFFFFFFFul;
	VN = 0ul;
	m--;
	MASK = 1ul << m;
	for(i = 0; i < n;i++){
		X = (B[(int)(t[i] &0x3)  ] | VN);
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





int validate_bpm_sse(unsigned char**  query, int* query_lengths,unsigned char* t,int n,int num)
{
	int i,j;
	int len = 0;
	//struct seq_info* si = 0;
	//unsigned long int new;
	unsigned int _MM_ALIGN16 nuc[16];
	//int _MM_ALIGN16 positions[4];
	//int _MM_ALIGN16 errors[4];
	
	unsigned int _MM_ALIGN16 lengths[4];
	
	__m128i VP,VN,D0,HN,HP,X,MASK,K,NOTONE,POS,diff,zero,one;
	__m128i* nuc_p;
	__m128i xmm1,xmm2;
	//long int MASK = 0;
	
	for(i = 0; i < 16;i++){
		nuc[i] = 0ul;
	}
	
	for(i = 0; i < num;i++){
		//si = qs->seq_info[abs(assignment[i])];
		len = query_lengths[i];
		if(len > 31){
			len = 31;
		}
		
		lengths[i] = len;
		
		
		
		for(j = 0; j < len;j++){
			nuc[((int)(query[i][j] & 0x3u) << 2) + i] |=  (1ul << j);// (unsigned long int)(len-1-j));
			//			fprintf(stderr,"%c",si->seq[j]+65);
		}
		//}else{
		//	for(j = 0; j < len;j++){
		//		nuc[((int)(si->reverse_seq[j] & 0x3u) << 2) + i] |=  (1ul << (unsigned long int)(len-1-j));
		//		//			fprintf(stderr,"%c",si->reverse_seq[j] + 65);
		//	}
		//}
		//	fprintf(stderr,"\n" );
	}
	nuc_p = (__m128i*) nuc;
	//S = _mm_set1_epi32(0);
	zero = _mm_set1_epi32(0);
	one = _mm_set1_epi32(1);
	diff =  _mm_load_si128 ( (__m128i*) lengths );  // _mm_set1_epi32(m);
	VP =  _mm_set1_epi32(0xFFFFFFFFu);
	VN =  _mm_set1_epi32(0);
	NOTONE =  _mm_set1_epi32(0xFFFFFFFF);
	K =  _mm_set1_epi32(0x7FFFFFFF);
	POS = _mm_set1_epi32(0);
	
	//VP = 0xFFFFFFFFFFFFFFFFul;
	//VN = 0ul;
	
	for(i = 0; i< 4;i++){
		lengths[i]--;
		lengths[i] = 1 << lengths[i];
		//	fprintf(stderr,"%d	%d	LEN:%d\n",i,num,lengths[i]);
	}
	
	// m--;
	MASK =  _mm_load_si128 ( (__m128i*) lengths ); //  _mm_set1_epi32(1ul << m);
	for(i = 0; i < n ;i++){
		//fprintf(stderr,"%c",*t + 65);
		X = _mm_or_si128 (*(nuc_p +( (int)(*t)  & 0x3u) ) , VN);
		//X = (B[(int) *t] | VN);
		xmm1 = _mm_and_si128(X, VP);
		xmm2 = _mm_add_epi32(VP ,xmm1);
		xmm1 = _mm_xor_si128 (xmm2, VP);
		D0 = _mm_or_si128(xmm1, X);
		//D0 = ((VP+(X&VP)) ^ VP) | X ;
		HN = _mm_and_si128(VP, D0);
		//HN = VP & D0;
		xmm1 = _mm_or_si128(VP, D0);
		xmm2 = _mm_andnot_si128 (xmm1,NOTONE);
		HP = _mm_or_si128(VN, xmm2);
		//HP = VN | ~(VP | D0);
		X = _mm_slli_epi32(HP,1);
		//X = HP << 1ul;
		VN = _mm_and_si128(X, D0);
		//VN = X & D0;
		xmm1 = _mm_slli_epi32(HN,1);
		xmm2 = _mm_or_si128(X, D0);
		xmm2 = _mm_andnot_si128 (xmm2,NOTONE);
		VP = _mm_or_si128(xmm1, xmm2);
		//VP = (HN << 1ul) | ~(X | D0);
		xmm1 = _mm_and_si128(HP, MASK);
		xmm2 = _mm_cmpgt_epi32(xmm1, zero);
		diff = _mm_add_epi32(diff , _mm_and_si128( xmm2, one));
		//diff += (HP & MASK) >> m;
		xmm1 = _mm_and_si128(HN, MASK);
		xmm2 = _mm_cmpgt_epi32(xmm1, zero);
		diff = _mm_sub_epi32(diff,  _mm_and_si128( xmm2, one));
		//diff -= (HN & MASK) >> m;
		xmm1 = _mm_cmplt_epi32(diff, K);
		xmm2 = _mm_and_si128(xmm1, diff);
		K = _mm_or_si128(xmm2, _mm_andnot_si128  (xmm1,K));
		//xmm2 = _mm_and_si128(xmm1, _mm_set1_epi32(i));
		//POS = _mm_or_si128(xmm2, _mm_andnot_si128  (xmm1,POS));
		t++;
	}
	//fprintf(stderr,"\n");
	//_mm_store_si128 ((__m128i*) positions, POS);
	//fprintf(stderr,"%d	%d	%d	%d	",positions[0],positions[1],positions[2],positions[3]);
	_mm_store_si128 ((__m128i*) query_lengths, K);
	//fprintf(stderr,"%d	%d	%d	%d\n",query_lengths[0],query_lengths[1],query_lengths[2],query_lengths[3]);
	
	/*for(i = 0; i < num;i++){
	 si = qs->seq_info[abs(assignment[i])];
	 //len = si->len;
	 //lengths[i] = len;
	 //	fprintf(stderr,"%d	%d	%d	%d	%d\n",  assignment[i] , si->len,num,out2[i],out1[i] + offset);
	 if(assignment[i] >= 0){
	 new = ((unsigned long int)(si->len - errors[i])) << 56ul;
	 new |= ((unsigned long int) ((unsigned long int)offset - (unsigned long int) positions[i]) << 1ul);
	 BinaryInsertionSortOne(si->match_units,LIST_STORE_SIZE ,new);
	 }else{
	 new = ((unsigned long int)(si->len - errors[i])) << 56ul;
	 new |= ((unsigned long int) ((unsigned long int) offset - (unsigned long int) positions[i]) << 1ul);
	 new |= 1ul ;
	 BinaryInsertionSortOne(si->match_units,LIST_STORE_SIZE ,new);
	 }
	 
	 }*/
	return 1;
}





unsigned char* reverse_complement2(unsigned char* p,int len)
{
	unsigned char* tmp = malloc(sizeof(unsigned char)*MAX_SEQ_LEN);
	int i,c;
	c = 0;
	for(i =len-1; i >= 0;i--){
		tmp[c] = rev_nuc_code[(int)p[i]];
		c++;
	}
	tmp[c] = 0;
	for(i= 0; i < len;i++){
		p[i] = tmp[i];
	}
	free(tmp);
	return p;
}

