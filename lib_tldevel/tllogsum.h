#ifndef TLLOGSUM_H
#define TLLOGSUM_H

#ifdef TLLOGSUM_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN void init_logsum(void);
EXTERN float logsum(const float a,const float b);
EXTERN float prob2scaledprob(float p);
EXTERN float scaledprob2prob(float p);

#undef TLLOGSUM_IMPORT
#undef EXTERN
#endif
