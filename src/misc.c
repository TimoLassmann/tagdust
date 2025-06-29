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
/*! \file misc.c
 \brief Miscellaneous Functions
 
 A collection of functions used to:
 - parse strings
 - compare strings 
 - add probabilities in logspace
 */

#if HAVE_CONFIG_H
#include "config.h"
#endif
#include "kslib.h"

#include "tagdust2.h"
#include "nuc_code.h"
#include "misc.h"
#include <ctype.h>



static unsigned long next = 1;


/** \var float logsum_lookup
 \brief Lookup table.
 */
float logsum_lookup[LOGSUM_SIZE];


/** \fn void init_logsum()
 \brief Initializes lookup table.

 
 \warning Remember to call before using the logsum() function.
 */
void init_logsum()
{
	int i;
	for(i = 0; i < LOGSUM_SIZE;i++){
		logsum_lookup[i] =  log(1. + exp((double) -i / SCALE));
	}
}


/** \fn float logsum(float a,float b)
 \brief Sums two probabilities in log space.
 
   \return  log(exp(a) + exp(b)
 \warning Remember to call init_logsum() before using logsum().
 */
float logsum(float a,float b)
{
	
	const float max = HMMER3_MAX(a, b);
	const float min = HMMER3_MIN(a, b);
	return (min == -SCALEINFTY || (max-min) >= 15.7f) ? max : max + logsum_lookup[(int)((max-min)*SCALE)];
}


/** \fn float prob2scaledprob(float p)
\brief Returns log(p).
  \return log(p) 
*/
float prob2scaledprob(float p)
{
	if(p == 0.0){
		return -SCALEINFTY;
	}else{
		return  log(p);
	}
}

/** \fn float scaledprob2prob(float p)
 \brief Returns exp(p).
 \return exp(p)
 */
float scaledprob2prob(float p)
{
	if(p == -SCALEINFTY){
		return 0.0;
	}else{
		return exp(p);
	}
}


/** \fn int byg_count(char* pattern,char*text)
 \brief Counts occurance of pattern text. 
 \param pattern string containing pattern.
 \param text string containing text.
 \exception Input strings must be 0 terminated. 
 
 \return number of occurences.
 */
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


/** \fn int byg_end(const char* pattern,const char*text)
 \brief Finds pattern in text and returns the end coordinate of match.
 \param pattern string containing pattern.
 \param text string containing text.
 \exception Input strings must be 0 terminated.
 \return index.
 */

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

/** \fn int count_string(const char*p,const char** suffix,int h,int len)
 \brief Counts occurance of p in a suffix array.
 \param p string containing pattern.
 \param suffix suffix array. 
 \param h size of suffix array.
 \param len length of pattern.
 \return index.
 */

int count_string(const char*p,const char** suffix,int h,int len)
{
	int a,b;
	//for(i = 0; i < 1000000;i++){
	a = binsearch_down(p,suffix,h,len);
	b = binsearch_up(p,suffix,h,len);
	return b-a;
}

/** \fn int binsearch_down(const char*p,const char** suffix,int h,int len)
 \brief finds first occurance of p in suffix array.
 \param p string containing pattern.
 \param suffix suffix array.
 \param h size of suffix array.
 \param len length of pattern.
 \return index.
 */
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

/** \fn int binsearch_up(const char*p,const char** suffix,int h,int len)
 \brief finds last occurance of p in suffix array.
 \param p string containing pattern.
 \param suffix suffix array.
 \param h size of suffix array.
 \param len length of pattern.
 \return index.
 */
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



char* append_message(char* old_message, char* new_message)
{
	static size_t message_len = 0;
	struct tm *ptr;
	int status;
	
	char time_string[200];
	int hour;
	//char am_or_pm;
	time_t current = time(NULL);
	ptr = localtime(&current);
	hour = ptr->tm_hour;
	/*if (hour <= 11)
		am_or_pm = 'a';
	else {
		hour -= 12;
		am_or_pm = 'p';
	}*/
	if (hour == 0){
		hour = 12;
	}
	
	strftime(time_string, 200, "[%F %H:%M:%S]\t", ptr);
	fprintf(stderr,"%s%s",time_string,new_message );
	//%H:%M:%S.000
	//sprintf(time_string,"[%.2d-%.2d-%d %2d:%.2d%cm\t",ptr->tm_mon + 1,ptr->tm_mday, ptr->tm_year + 1900,hour,ptr->tm_min, am_or_pm );
	size_t time_len = strlen(time_string);
	
	size_t added_len = strlen(new_message);
	
	if(message_len == 0){
		MMALLOC(old_message,sizeof(char) *( time_len+added_len+1));
	}else{
		MREALLOC(old_message,  sizeof(char) *( message_len + time_len+added_len + 1));

	}
	
	
		//char *concat = (char*) malloc(len1 + len2 + 1);
	memcpy(old_message+message_len, time_string, time_len+1);
	
	memcpy(old_message+message_len+time_len, new_message, added_len+1);
	
	message_len =strlen(old_message);

	
	return old_message;
ERROR:
	KSLIB_MESSAGE(status,"Something wrong in append_message.\n");
	return NULL;
}


/** \fn int qsort_string_cmp(const void *a, const void *b)
 \brief Compares two strings. 
 Used to sort arrays of string using qsort.   
 \param a void pointer to first string. 
 \param b void pointer to second string. 
 */

int qsort_string_cmp(const void *a, const void *b)
{
	const char **one = (const char **)a;
	const char **two = (const char **)b;
	return strcmp(*one, *two);
}

/** \fn int qsort_flt_cmp(const void *a, const void *b)
 \brief Compares two floats.
 Used to sort arrays of floats.
 \param a void pointer to first float.
 \param b void pointer to second float.
 */
int qsort_flt_cmp(const void * a, const void * b)
{
	//const float a  = (float) *elem1;
	if(*(const float*)a > *(const float*)b)
		return -1;
	return *(const float*)a < *(const float*)b;
	
	//return (*(float*) a) - (*(float*) b);
}


/** \fn double gaussian_pdf(double x, double m,double s)
 Calculates the gaussian probability density function $P(X=x)$ for a normal distribution.
 \param x value.
 \param m mean.
 \param s standard deviation.
 */
double gaussian_pdf(double x, double m,double s)
{
	double a = (x-m) / s;
	return INV_SQRT_2PI / s *exp(-0.5 * a * a);
}


/** \fn unsigned int pop(int x)
 \brief Counts bits in x.
 \param x value.
 \return number of bits in \a x
 */
unsigned int pop(int x)
{
	unsigned int n = 0;
	while(x != 0){
		n = n +1;
		x = x &(x-1);
	}
	return n;
}

/** \fn int bpm(const  char* t,const  char* p,int n,int m)
 \brief Calculates edit distance between two strings. 
 \param t string1.
  \param p string 2.
  \param n length of t.
  \param m length of p.
 \return exit distance.
 */
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
			//if(k <= limit){
			//	return (int)k;
			//}
			
		}
		
	}
	return (int)k;
}



/** \fn int bpm(const  char* t,const  char* p,int n,int m)
 \brief Calculates edit distance between two strings.
 \param t string1.
 \param p string 2.
 \param n length of t.
 \param m length of p.
 \return exit distance.
 */
int bpm_global(const  char* t,const  char* p,int n,int m)
{
	register unsigned long int i;//,c;
	unsigned long int diff;
	unsigned long int B[255];
	int status;
	
	int c;
	char* p1= 0;
	char* p2 = 0;
	
	MMALLOC(p1, sizeof(char) * (n+11));
	MMALLOC(p2, sizeof(char) * (m+11));
	
	
	for(i = 0; i < 5;i++){
		p1[i] = 'F';
		p2[i] = 'F';
	}
	c = 5;
	for(i = 0; i < n;i++){
		p1[c] = t[i];
		c++;
	}
	for(i = 0; i < 5;i++){
		p1[c] = 'Q';
		c++;
	}
	p1[c] = 0;
	n = c;
	
	
	c = 5;
	for(i = 0; i < m;i++){
		p2[c] = p[i];
		c++;
	}
	for(i = 0; i < 5;i++){
		p2[c] = 'Q';
		c++;
	}
	p2[c] = 0;
	m = c;
	
	
	
	
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
		B[(int)(p2[i] )] |= (1ul << i);
	}
	
	//c = 0;
	VP = 0xFFFFFFFFFFFFFFFFul;
	VN = 0ul;
	m--;
	MASK = 1ul << m;
	
	for(i = 0; i < n;i++){
		X = (B[(int)(p1[i])  ] | VN);
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
			//if(k <= limit){
			//	return (int)k;
			//}
			
		}
		
	}
	MFREE(p1);
	MFREE(p2);
	
	return (int)k;
ERROR:
	KSLIB_MESSAGE(status,"Somethign wrong in bpm_global.\n");
	return kslFAIL;
}



/** \fn int bpm_check_error(const unsigned char* t,const unsigned char* p,int n,int m)
 \brief Calculates edit distance between two strings.
 \param t string1.
 \param p string 2.
 \param n length of t.
 \param m length of p.
 \return exit distance.
 */
int bpm_check_error(const unsigned char* t,const unsigned char* p,int n,int m)
{
	register unsigned long int i;//,c;
	unsigned long int diff;
	unsigned long int B[5];
	
	int new_len = 0;
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
		if(p[i] != 65){
		B[(int)(p[i] & 0x3)] |= (1ul << i);
		new_len++;
		}
		
	}
	if(new_len > 31){
		new_len = 31;
	}
	m = new_len;
	k =new_len;
	
	
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




/** \fn int bpm_check_error_global(const unsigned char* t,const unsigned char* p,int n,int m)
 \brief Calculates edit distance between two strings.
 \param t string1.
 \param p string 2.
 \param n length of t.
 \param m length of p.
 \return exit distance.
 */
int bpm_check_error_global(const unsigned char* t,const unsigned char* p,int n,int m)
{
	register unsigned long int i;//,c;
	unsigned long int diff;
	unsigned long int B[5];
	if(m > 63){
		m = 63;
	}
	
	//unsigned long int k = m;
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
		//if(diff < k){
		//	k = diff;
			//fprintf(stderr,"%ld	%ld\n",i,k);
			//if(k <= limit){
			//	return (int)k;
			//}
			
		//}
	}
	return (int)diff;
}


/** Helper macros for portable bit-parallel Myers algorithm */
#define DIV_CEIL(a,b) (a == 0 ? 1 : a/b+(a%b == 0 ? 0 : 1))

/** \fn int bmp_single(const unsigned char* t, const unsigned char* p, int n, int m)
 \brief Portable Myers bit-parallel algorithm for single sequence comparison.
 \param t target sequence.
 \param p pattern sequence.
 \param n length of target.
 \param m length of pattern.
 \return edit distance.
 */
static int bmp_single(const unsigned char* t, const unsigned char* p, int n, int m)
{
    register uint64_t VP,VN,D0,HN,HP,X;
    register uint64_t i;
    uint64_t MASK = 0;
    int64_t diff;
    uint64_t B[4];  // 4 nucleotides: A=0, C=1, G=2, T=3
    int k;

    if(m > 63){
        m = 63;
    }
    diff = m;
    k = m;
    
    for(i = 0; i < 4; i++){
        B[i] = 0;
    }

    // Build bit patterns for each nucleotide
    for(i = 0; i < (uint64_t)m; i++){
        if(p[i] != 65){  // Skip 'A' characters (65) like original
            B[p[i] & 0x3u] |= (1ul << i);
        }
    }

    VP = (1ul << (m))-1;
    VN = 0ul;
    m--;
    MASK = 1ul << (m);

    for(i = 0; i < (uint64_t) n; i++){
        X = (B[t[i] & 0x3u] | VN);
        D0 = ((VP+(X&VP)) ^ VP) | X;
        HN = VP & D0;
        HP = VN | (~(VP | D0));
        X = HP << 1ul;
        VN = X & D0;
        VP = (HN << 1ul) | (~(X | D0));

        diff += (HP & MASK)? 1 : 0;
        diff -= (HN & MASK)? 1 : 0;
        if(diff < k){
            k = diff;
        }
    }
    return k;
}

/** \fn int validate_bpm_sse(unsigned char**  query, int* query_lengths,unsigned char* t,int n,int num)
 \brief Calculates edit distance between four queries and one target using portable algorithm.
 \param query array of four strings.
 \param query_lengths lengths of query strings.
 \param t target sequence.
 \param n length of t.
 \param num number of queries (should be 4).
 */

int validate_bpm_sse(unsigned char**  query, int* query_lengths,unsigned char* t,int n,int num)
{
    int i;
    int results[4];
    
    // Process each query sequence individually using portable algorithm
    for(i = 0; i < num && i < 4; i++){
        if(query_lengths[i] > 0){
            results[i] = bmp_single(t, query[i], n, query_lengths[i]);
        } else {
            results[i] = n;  // Maximum possible distance for empty queries
        }
    }
    
    // Store results back in query_lengths array (same behavior as original)
    for(i = 0; i < num && i < 4; i++){
        query_lengths[i] = results[i];
    }
    
    return 1;
}


/** \fn char* shorten_pathname(char* p)
 \brief Moves string pointer to character after first backslash. 
 
 \param p file name including path. 
 \return pointer to character after last backslash.
 */


char* shorten_pathname(char* p)
{
	int i;
	char* tmp = p;
	for(i = 0; i< strlen(p);i++){
		if(p[i] == '/'){
			tmp = p+i +1;
		}
	}
	return tmp;
}


/** \fn unsigned char* reverse_complement(unsigned char* p,int len)
 \brief Reverses and complements nucleotide sequences. 
 
 \param p nucleotide sequence.
 \param len length.
 */



unsigned char* reverse_complement(unsigned char* p,int len)
{
	unsigned char* tmp = 0;
	int status;
	MMALLOC(tmp,sizeof(unsigned char)*(len +2));
	int i,c;
	c = 0;
	for(i =len-1; i >= 0;i--){
		if(p[i]== 65){
			tmp[c] = 65;
		}else{
			tmp[c] = rev_nuc_code[(int)p[i]];
		}
		c++;
	}
	tmp[c] = 0;
	for(i= 0; i < len;i++){
		p[i] = tmp[i];
	}
	MFREE(tmp);
	return p;
ERROR:
	KSLIB_MESSAGE(status,"Something wrong in reverse_complement.\n");
	return NULL;
}

/** \fn void reverse_sequence(char* p,int len)
 \brief Reverses sequences.
 
 \param p nucleotide sequence.
 \param len length.
 */


void reverse_sequence(char* p,int len)
{
	int c, i, j;
	
	for (i = 0, j = len - 1; i < j; i++, j--)
	{
		c = p[i];
		p[i] = p[j];
		p[j] = c;
	}
}



/* RAND_MAX assumed to be 32767 */
int myrand(void)
{
	next = next * 1103515245 + 12345;
	return((unsigned)(next/65536) % 32768);
}

void mysrand(unsigned seed)
{
	next = seed;
}

int file_exists(char* name)
{
	struct stat buf;
	int ret,local_ret;
	ret = 0;
	
	local_ret= stat ( name, &buf );
		/* File found */
	if ( local_ret == 0 )
	{
		ret++;
			//return 1;
	}
	return ret;
}

int bitcount64(long long int i)
{
	i = i - ((i >> 1) & 0x5555555555555555);
	i = (i & 0x3333333333333333) + ((i >> 2) & 0x3333333333333333);
	return (((i + (i >> 4)) & 0xF0F0F0F0F0F0F0F) * 0x101010101010101) >> 56;
}

/* CODE from Kaz Kylheku on Stackoverflow - works great - stays away from sign bit (i.e. works for 63 bits)*/
int highest_bit(long long int n)
{
	const long long mask[] = {
		0x000000007FFFFFFF,
		0x000000000000FFFF,
		0x00000000000000FF,
		0x000000000000000F,
		0x0000000000000003,
		0x0000000000000001
	};
	int hi = 64;
	int lo = 0;
	int i = 0;
	
	if (n == 0)
		return 0;
	
	for (i = 0; i < sizeof mask / sizeof mask[0]; i++) {
		int mi = lo + (hi - lo) / 2;
		
		if ((n >> mi) != 0)
			lo = mi;
		else if ((n & (mask[i] << lo)) != 0)
			hi = mi;
	}
	
	return lo + 1;
}

