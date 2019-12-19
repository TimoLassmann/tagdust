
#include <math.h>

#define TLLOGSUM_IMPORT
#include "tllogsum.h"

#define SCALE 100000.0
#define LOGSUM_SIZE 1600000

#define MACRO_MIN(a,b)          (((a)<(b))?(a):(b))
#define MACRO_MAX(a,b)          (((a)>(b))?(a):(b))

static float logsum_lookup[LOGSUM_SIZE];

void init_logsum()
{
        static int called = 0;
        int i;
        if(!called){
                called = 1;
                for(i = 0; i < LOGSUM_SIZE;i++){
                        logsum_lookup[i] = log(1.0 +exp((double) -i / SCALE));
                }
        }
}

float logsum(const float a,const float b)
{
        register const float max = MACRO_MAX(a, b);
        register const float min = MACRO_MIN(a, b);

        if(min == -INFINITY){
                return max;
        }
        if( (max-min) >= 15.7f){
                return max;
        }
        return  max+ logsum_lookup[(int)((max-min)*SCALE)];
}

float prob2scaledprob(float p)
{
        if(p == 0.0){
                return -INFINITY;
        }else{
                return  log(p);
        }
}

float scaledprob2prob(float p)
{
        if(p == -INFINITY){
                return 0.0;
        }else{
                return exp(p);
        }
}
