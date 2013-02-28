//
//  dbscan.c
//  tagdust2
//
//  Created by lassmann on 2/28/13.
//  Copyright (c) 2013 lassmann. All rights reserved.
//

#include <stdio.h>
#include <xmmintrin.h>

#include "interface.h"
#include "io.h"
#include "pst.h"
#include "misc.h"
#include "dbscan.h"


void cluster_reads_based_on_pst_patterns(struct pst_node** patterns, int num_patterns, int numseq)
{
	int i,j,c,max_len,representative;
	float** dm = 0;
	//num_patterns = 20;
	
	float eps = 0.1;
	float min_points = 2;
	
	int* closeby = malloc(sizeof(int) * num_patterns);
	int* closeby2 = malloc(sizeof(int) * num_patterns);
	int* clusters = malloc(sizeof(int) * num_patterns);
	int* visited = malloc(sizeof(int) * num_patterns);
		
	int clusterID = 0;
	for(i = 0 ; i  < num_patterns;i++){
		clusters[i] = -1;
		visited[i] = 0;
	}
	
	
	
	dm = malloc(sizeof(float*) * num_patterns );
	for(i =0 ; i < num_patterns;i++){
		patterns[i]->last_seen = 0;
		dm[i] =  malloc(sizeof(float) * num_patterns );
		for(j = 0 ; j < num_patterns;j++){
			dm[i][j] = 0.0f;
			
			dm[i][j] = 1.0f -  jaccard_coefficient(patterns[i]->bit_occ, patterns[j]->bit_occ ,1+ numseq  / BITSPERWORD);
			//if(strlen(patterns[i]->label) > 5 && strlen(patterns[j]->label)  ){
			//	fprintf(stderr,"%d	%d	%s	%s	%f\n",i,j,patterns[i]->label,  patterns[j]->label,  dm[i][j] );
			//}
		}
	}
	
	//dbscanning .. .
	
	clusterID = 1;
	for(i = 0 ; i  < num_patterns;i++){
		if(!visited[i]){
			
			
			
			///fprintf(stderr,"starting %d\n",i);
			if(expand(clusters, visited, eps, min_points,num_patterns ,dm,i,clusterID)){
				
				clusterID++;
			}
			//};/
		}
		/*for(c = 0;c < num_patterns;c++){
			fprintf(stderr,"%d ",c);// )
		}
		fprintf(stderr,"\n");// )
		for(c = 0;c < num_patterns;c++){
			fprintf(stderr,"%d ",visited[c]);// )
		}
		fprintf(stderr,"\n");// )
		for(c = 0;c < num_patterns;c++){
			fprintf(stderr,"%d ",clusters[c]);// )
		}
		fprintf(stderr,"\n");// )	*/	
		
	}
	
	
	
	for(c = 1;c < clusterID;c++){
		fprintf(stderr,"CLUSTER:%d\n", c);
		representative = -1;
		max_len = -1;
		for(i = 0 ; i  < num_patterns;i++){
			
			if(clusters[i] == c){
				if((int)strlen(patterns[i]->label) > max_len){
					max_len = (int)strlen(patterns[i]->label) ;
					representative = i;
				}
				//fprintf(stderr,"%s	%d\n",patterns[i]->label,(int)strlen(patterns[i]->label) );
			}
			
			
	
		}
		
		for(i = 0 ; i  < num_patterns;i++){
			
			if(clusters[i] == c){
				for(j = 0 ; j<1+ numseq  / BITSPERWORD ;j++){
					patterns[representative]->bit_occ[j] |= patterns[i]->bit_occ[j];
				}
				
			}
		}
		max_len = 0;
		//fprintf(stderr,"%d\n",j);
		for(j = 0 ; j<1+ numseq  / BITSPERWORD ;j++){
			max_len += pop( patterns[representative]->bit_occ[j]);// |= patterns[i]->bit_occ[j];
		}
		
		fprintf(stderr,"%s	%d	%d\n",patterns[representative]->label,  max_len,numseq );
		
	}
	/*fprintf(stderr,"Leftover:\n");
	for(i = 0 ; i  < num_patterns;i++){
		if(clusters[i] == -1){
			fprintf(stderr,"%s\n",patterns[i]->label);
		}
		
	}*/
	
	for(i =0 ; i < num_patterns;i++){
		free(dm[i]);/// =  malloc(sizeof(float) * num_patterns );
	}
	free(dm);
	free(closeby2);

	free(closeby);
	free(clusters);
	
}


int expand (int * cluster,int * visited , float eps, int min_size,int n,float** dm,int t,int ID)
{
	int j;
	int ex_size = 0;
	
	
	ex_size = 0;
	for(j = 0; j < n;j++){
		if(dm[t][j] < eps){
			ex_size++;
			//fprintf(stderr," %d near %d\n" ,  t,j );
		}
	}
	//fprintf(stderr,"%d	%d\n",t, ex_size);
	if(ex_size > min_size){
		cluster[t] = ID;
		for(j = 0; j < n;j++){
			
			if(dm[t][j] < eps){
				if(!visited[j]){
				//cluster[j] = 2;
				
					visited[j] = 1; //[j]->last_seen = 1;
					
					expand(cluster, visited, eps,min_size,n,dm,j,ID);
				}
				if(cluster[j] == -1){
					cluster[j]  = cluster[t];
					
				}
								
			}
		}
		return 1;
	}else{
		return 0;
	}
}




float jaccard_coefficient(int* a, int* b , int n)
{
	
	float j_intersection = 0.0;
	float j_union = 0.0;
	int i;//,x;
	
	/*for(i = 0; i < n;i++){
	 x = a[i] & b[i];
	 j_intersection += pop(x);
	 
	 x = a[i] | b[i];
	 j_union += pop(x);
	 }*/
	
	
	int j = n / 4;
	const __m128i* a_ptr = (__m128i*) a;
	const __m128i* b_ptr =  (__m128i*) b;
	__m128i xmm1;
	__m128i xmm2;
	__m128i xmm3;
	__m128i xmm4;
	
	const unsigned mu1 = 0x55555555;
	const unsigned mu2 = 0x33333333;
	const unsigned mu3 = 0x0F0F0F0F;
	const unsigned mu4 = 0x0000003F;
	
	__m128i m1 = _mm_set_epi32 (mu1, mu1, mu1, mu1);
	__m128i m2 = _mm_set_epi32 (mu2, mu2, mu2, mu2);
	__m128i m3 = _mm_set_epi32 (mu3, mu3, mu3, mu3);
	__m128i m4 = _mm_set_epi32 (mu4, mu4, mu4, mu4);
	__m128i mcnt_i;
	__m128i mcnt_u;
	mcnt_i =  _mm_set1_epi32(0);//   _mm_xor_si128(mcnt_i, mcnt_i); // cnt = 0
	mcnt_u = _mm_set1_epi32(0);//
	
	for(i = 0; i < j;i++){
		xmm1 = _mm_load_si128(a_ptr);
		xmm2 = _mm_load_si128(b_ptr);
		xmm4 = _mm_or_si128(xmm1, xmm2);
		xmm1 = _mm_and_si128(xmm1, xmm2);
		
		
		// b = (b & 0x55555555) + (b >> 1 & 0x55555555);
		xmm2 = _mm_srli_epi32(xmm1, 1);                    // tmp1 = (b >> 1 & 0x55555555)
		xmm2 = _mm_and_si128(xmm2, m1);
		xmm3 = _mm_and_si128(xmm1, m1);                    // tmp2 = (b & 0x55555555)
		xmm1    = _mm_add_epi32(xmm2, xmm3);               //  b = tmp1 + tmp2
		
		// b = (b & 0x33333333) + (b >> 2 & 0x33333333);
		xmm2 = _mm_srli_epi32(xmm1, 2);                    // (b >> 2 & 0x33333333)
		xmm2 = _mm_and_si128(xmm2, m2);
		xmm3 = _mm_and_si128(xmm1, m2);                    // (b & 0x33333333)
		xmm1    = _mm_add_epi32(xmm2, xmm3);               // b = tmp1 + tmp2
		
		// b = (b + (b >> 4)) & 0x0F0F0F0F;
		xmm2 = _mm_srli_epi32(xmm1, 4);                    // tmp1 = b >> 4
		xmm1 = _mm_add_epi32(xmm1, xmm2);                     // b = b + (b >> 4)
		xmm1 = _mm_and_si128(xmm1, m3);                       //           & 0x0F0F0F0F
		
		// b = b + (b >> 8);
		xmm2 = _mm_srli_epi32 (xmm1, 8);                   // tmp1 = b >> 8
		xmm1 = _mm_add_epi32(xmm1, xmm2);                     // b = b + (b >> 8)
		
		// b = (b + (b >> 16)) & 0x0000003F;
		xmm2 = _mm_srli_epi32 (xmm1, 16);                  // b >> 16
		xmm1 = _mm_add_epi32(xmm1, xmm2);                     // b + (b >> 16)
		xmm1 = _mm_and_si128(xmm1, m4);                       // (b >> 16) & 0x0000003F;
		
		mcnt_i = _mm_add_epi32(mcnt_i,xmm1);
		
		
		// b = (b & 0x55555555) + (b >> 1 & 0x55555555);
		xmm2 = _mm_srli_epi32(xmm4, 1);                    // tmp1 = (b >> 1 & 0x55555555)
		xmm2 = _mm_and_si128(xmm2, m1);
		xmm3 = _mm_and_si128(xmm4, m1);                    // tmp2 = (b & 0x55555555)
		xmm4    = _mm_add_epi32(xmm2, xmm3);               //  b = tmp1 + tmp2
		
		// b = (b & 0x33333333) + (b >> 2 & 0x33333333);
		xmm2 = _mm_srli_epi32(xmm4, 2);                    // (b >> 2 & 0x33333333)
		xmm2 = _mm_and_si128(xmm2, m2);
		xmm3 = _mm_and_si128(xmm4, m2);                    // (b & 0x33333333)
		xmm4    = _mm_add_epi32(xmm2, xmm3);               // b = tmp1 + tmp2
		
		// b = (b + (b >> 4)) & 0x0F0F0F0F;
		xmm2 = _mm_srli_epi32(xmm4, 4);                    // tmp1 = b >> 4
		xmm4 = _mm_add_epi32(xmm4, xmm2);                     // b = b + (b >> 4)
		xmm4 = _mm_and_si128(xmm4, m3);                       //           & 0x0F0F0F0F
		
		// b = b + (b >> 8);
		xmm2 = _mm_srli_epi32 (xmm4, 8);                   // tmp1 = b >> 8
		xmm4 = _mm_add_epi32(xmm4, xmm2);                     // b = b + (b >> 8)
		
		// b = (b + (b >> 16)) & 0x0000003F;
		xmm2 = _mm_srli_epi32 (xmm4, 16);                  // b >> 16
		xmm4 = _mm_add_epi32(xmm4, xmm2);                     // b + (b >> 16)
		xmm4 = _mm_and_si128(xmm4, m4);                       // (b >> 16) & 0x0000003F;
		
		mcnt_u = _mm_add_epi32(mcnt_u,xmm4);
		
		
		
		
		
		
		a_ptr++;
		b_ptr++;
		
	}
	_MM_ALIGN16 int result[4];
	_mm_storeu_si128((__m128i*)result, mcnt_i);
	//fprintf(stderr,"%f\t",j_intersection);
	j_intersection = (float)(result[0] + result[1] + result[2] + result[3]);
	//fprintf(stderr,"%f\t",j_intersection);
	_mm_storeu_si128((__m128i*)result, mcnt_u);
	//fprintf(stderr,"%f\t",j_union);
	j_union = (float)(result[0] + result[1] + result[2] + result[3]);
	//fprintf(stderr,"%f\n",j_union);
	
	//fprintf(stderr, "I:%f\tu:%f\n",j_intersection,j_union );
	return j_intersection / j_union;
}



