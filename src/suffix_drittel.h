//
//  suffix_drittel.h
//  tagdust2
//
//  Created by lassmann on 2/25/13.
//  Copyright (c) 2013 lassmann. All rights reserved.
//

#ifndef tagdust2_suffix_drittel_h
#define tagdust2_suffix_drittel_h

#define GetI() (SA12[t] < n0 ? SA12[t] * 3 + 1 : (SA12[t] - n0) * 3 + 2)


void suffixArray(int* s, int* SA, int n, int K);
int leq(int a1, int a2,   int b1, int b2);                                            // and triples
int leq3(int a1, int a2, int a3,   int b1, int b2, int b3);
void radixPass(int* a, int* b, int* r, int n, int K);

#endif
