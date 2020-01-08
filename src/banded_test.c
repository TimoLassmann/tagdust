#include "tldevel.h"

#include "tlalphabet.h"

#include "tlrng.h"
#include <stdint.h>
#include "string.h"

uint8_t dyn_256_print(const uint8_t* t,const uint8_t* p,int n,int m);


int main(int argc, char *argv[])
{

        struct alphabet* alphabet = NULL;
        struct rng_state* rng;
        //long int r;
        int len = 0;
        int i,j,c;
        uint8_t* a = NULL;
        uint8_t* b = NULL;
        uint8_t score;


        RUNP(rng = init_rng(0));

        RUN(create_alphabet(&alphabet,rng,TLALPHABET_NOAMBIGUOUS_DNA));

        char seq[] = "GCTGACACGCTGTCCTCTGGCGACC";
        len = strlen(seq);
        if(len > 256){
                len = 256;
        }

        MMALLOC(a , sizeof(uint8_t) * len) ;
        MMALLOC(b , sizeof(uint8_t) * len) ;

        for(i = 0;i < len;i++){
                a[i] = tlalphabet_get_code(alphabet, seq[i]);
                b[i] = a[i];
                fprintf(stdout,"%d %d\n", a[i],b[i]);
        }

        score = dyn_256_print(a,b,len,10);

        score = dyn_256_print(a,b,12,10);

        score = dyn_256_print(a,b,10,12);

        return EXIT_SUCCESS;
ERROR:
        return EXIT_FAILURE;
}


uint8_t dyn_256_print(const uint8_t* query,const uint8_t* ref,int l_query,int l_ref)
{
        uint8_t* prev = NULL;
        uint8_t* cur = NULL;

        uint8_t* tmp = NULL;
        int i,j,c;
        int s,e,x;
        int bw = abs(l_query - l_ref)+2;
        fprintf(stdout,"band: %d\n",bw);


        bw = l_ref > l_query? l_ref : l_query;
        if (bw > 0) bw = 0;
        if (bw < abs(l_ref - l_query)) bw = abs(l_ref - l_query);

                //beg = (i - diff) > 0 )

        MMALLOC(prev, sizeof(uint8_t)* 257);
        MMALLOC(cur, sizeof(uint8_t)* 257);
        cur[0] = 0;
          fprintf(stdout,"%3d", cur[0]);

        s= 0, e = l_ref < bw + 1? l_ref : bw + 1;
        for(j = 1; j <= l_ref;j++){
                cur[j] = cur[j-1] +1;
                if(abs(0-j)> bw){
                        fprintf(stdout,"  N");
                }else{

                        fprintf(stdout,"%3d", cur[j]);
                }
        }

        fprintf(stdout,"   (%d %d)\n",s,e);
        tmp  = cur;
        cur = prev;
        prev = tmp;


        for(i = 1; i <= l_query;i++){
                s = 0;
                e = l_ref;

                x = i - bw; s = s> x? s : x;
                x = i + bw; e = e < x?e : x;

                cur[0] = prev[0]+1;
                if(abs(i-0)> bw){
                        fprintf(stdout,"  N");
                }else{
                        fprintf(stdout,"%3d", cur[0]);
                }
                for(j = 1; j < l_ref;j++){
                        c = 1;
                        if(query[i-1] == ref[j-1]){
                                c = 0;

                        }
                        cur[j] = prev[j-1] +c ;

                        cur[j] = MACRO_MIN(cur[j], prev[j]+1);
                        cur[j] = MACRO_MIN(cur[j], cur[j-1]+1);
                        if(abs(i-j)> bw){
                                fprintf(stdout,"  N");
                        }else{
                                fprintf(stdout,"%3d", cur[j]);
                        }

                }
                j = l_ref;
                c = 1;
                if(query[i-1] == ref[j-1]){
                        c = 0;

                }
                cur[j] = prev[j-1] +c ;

                cur[j] = MACRO_MIN(cur[j], prev[j]);
                cur[j] = MACRO_MIN(cur[j], cur[j-1]+1);
                if(abs(i-j)> bw){
                        fprintf(stdout,"  N");
                }else{
                        fprintf(stdout,"%3d", cur[j]);
                }

                fprintf(stdout,"   (%d %d)\n",s,e);
                tmp  = cur;
                cur = prev;
                prev = tmp;

        }
        c = prev[l_ref];
        MFREE(prev);
        MFREE(cur);
        return c;
ERROR:
        return 255;

}
