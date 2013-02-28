//
//  dbscan.h
//  tagdust2
//
//  Created by lassmann on 2/28/13.
//  Copyright (c) 2013 lassmann. All rights reserved.
//

#ifndef tagdust2_dbscan_h
#define tagdust2_dbscan_h


#ifndef _MM_ALIGN16
#ifdef __GNUC__
#define _MM_ALIGN16 __attribute__((aligned (16)))
#endif
#ifdef __MSVC__
#define _MM_ALIGN16 __declspec(align(16))
#endif
#endif


float jaccard_coefficient(int* a, int* b , int n);
void cluster_reads_based_on_pst_patterns(struct pst_node** patterns, int num_patterns, int numseq);


int expand (int * cluster,int * visited , float eps, int min_size,int n,float** dm,int t,int ID);
#endif

