//
//  suffic_drittel.c
//  tagdust2
// Simple Linear Work Suffix Array Construction⋆ Juha Ka ̈rkka ̈inen and Peter Sanders
// Max-Planck-Institut fu ̈r Informatik Stuhlsatzenhausweg 85, 66123 Saarbru ̈cken, Germany [juha, sanders]@mpi- sb.mpg.de.
//  Created by lassmann on 2/25/13.
//  Copyright (c) 2013 lassmann. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include "suffix_drittel.h"


int leq(int a1, int a2,   int b1, int b2)
{ // lexic. order for pairs
	return(a1 < b1 || (a1 == b1 && a2 <= b2));
	//return(a1 <= b1 && a2 <= b2);
	
}
// and triples
int leq3(int a1, int a2, int a3,   int b1, int b2, int b3)
{
	return(a1 < b1 || (a1 == b1 && leq(a2,a3, b2,b3)));
	//return(a1 <= b1 && leq(a2,a3, b2,b3));
}
// stably sort a[0..n-1] to b[0..n-1] with keys in 0..K from r
void radixPass(int* a, int* b, int* r, int n, int K)
{ // count occurrences
	//int* c = new int[K + 1];
	int* c = malloc (sizeof(int) * (K+1));
	
	// counter array
	for (int i = 0;  i <= K;  i++){
		c[i] = 0;         // reset counters

	}
	for (int i = 0;  i < n;  i++){
		c[r[a[i]]]++;    // count occurences
	}
	for (int i = 0, sum = 0;  i <= K;  i++) { // exclusive prefix sums
		int t = c[i];  c[i] = sum;  sum += t;
	}
	for (int i = 0;  i < n;  i++){
		b[c[r[a[i]]]++] = a[i];      // sort
	}
	free(c);
}

// find the suffix array SA of s[0..n-1] in {1..K}^n
// require s[n]=s[n+1]=s[n+2]=0, n>=2
void suffixArray(int* s, int* SA, int n, int K)
{
	int i,j,c0,c1,c2,name;
	int n0=(n+2)/3;
	int n1=(n+1)/3;
	int n2=n/3;
	int n02=n0+n2;
	
	int* s12  = malloc (sizeof(int) * (n02 +3));// )new int[n02 + 3];
	s12[n02] = 0;
	s12[n02+1] = 0;
	s12[n02+2] =0;
	
	
	int* SA12 = malloc (sizeof(int) * (n02 +3)); //new int[n02 + 3];
	
	SA12[n02] = 0;
	SA12[n02+1] = 0;
	SA12[n02+2] = 0;
	
	int* s0   = malloc (sizeof(int) * n0);// new int[n0];
	int* SA0  = malloc (sizeof(int) * n0);// new int[n0];
	
	// generate positions of mod 1 and mod  2 suffixes
	// the "+(n0-n1)" adds a dummy mod 1 suffix if n%3 == 1
	j = 0;
	for (i = 0; i < n+(n0-n1);  i++){
		if (i%3 != 0){
			s12[j++] = i;
		}
	}
	
	// lsb radix sort the mod 1 and mod 2 triples
	radixPass(s12 , SA12, s+2, n02, K);
	radixPass(SA12, s12 , s+1, n02, K);
	radixPass(s12 , SA12, s  , n02, K);
	
	// find lexicographic names of triples
	name = 0;
	c0 = -1;
	c1 = -1;
	c2 = -1;
	
	for (int i = 0;  i < n02;  i++) {
		if (s[SA12[i]] != c0 || s[SA12[i]+1] != c1 || s[SA12[i]+2] != c2) {
			name++;  c0 = s[SA12[i]];  c1 = s[SA12[i]+1];  c2 = s[SA12[i]+2];
		}
		if (SA12[i] % 3 == 1) {
			s12[SA12[i]/3]      = name;
		}else{//left half
			s12[SA12[i]/3 + n0] = name;
		} // right half
	}
	
	// recurse if names are not yet unique
	if (name < n02) {
		suffixArray(s12, SA12, n02, name);
		// store unique names in s12 using the suffix array
		for (int i = 0;  i < n02;  i++){
			s12[SA12[i]] = i + 1;
		}
	} else {// generate the suffix array of s12 directly
		for (int i = 0;  i < n02;  i++) SA12[s12[i] - 1] = i;
	}
	
	// stably sort the mod 0 suffixes from SA12 by their first character
	for (int i=0, j=0;  i < n02;  i++) if (SA12[i] < n0) s0[j++] = 3*SA12[i];
	radixPass(s0, SA0, s, n0, K);
	
	// merge sorted SA0 suffixes and sorted SA12 suffixes
	for (int p=0,  t=n0-n1,  k=0;  k < n;  k++) {

		int i = GetI(); // pos of current offset 12 suffix
		int j = SA0[p]; // pos of current offset 0  suffix
		if (SA12[t] < n0 ?
		    leq(s[i],       s12[SA12[t] + n0], s[j],       s12[j/3]) :
		    leq3(s[i],s[i+1],s12[SA12[t]-n0+1], s[j],s[j+1],s12[j/3+n0]))
		{ // suffix from SA12 is smaller
			SA[k] = i;  t++;
			if (t == n02) { // done --- only SA0 suffixes left
				for (k++;  p < n0;  p++, k++) SA[k] = SA0[p];
			}
		} else {
			SA[k] = j;  p++;
			if (p == n0)  { // done --- only SA12 suffixes left
				for (k++;  t < n02;  t++, k++) SA[k] = GetI();
			}
		}
	}
	free(s12);
	free(SA12);
	free(SA0);
	free(s0);
	//delete [] s12; delete [] SA12; delete [] SA0; delete [] s0;
}















