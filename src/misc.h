#ifndef MISC_H
#define MISC_H

#include "tldevel.h"

extern int byg_count(char* pattern,char*text);
int byg_end(const char* pattern,const char*text);
extern int bpm(const  char* t,const  char* p,int n,int m);
extern int bitcount64(long long int i);

extern int bpm_check_error(const unsigned char* t,const unsigned char* p,int n,int m);
extern int validate_bpm_sse(unsigned char**  query, int* query_lengths,unsigned char* t,int n,int num);

extern double gaussian_pdf(double x, double m,double s);

extern unsigned char* reverse_complement(unsigned char* p,int len);
#endif
