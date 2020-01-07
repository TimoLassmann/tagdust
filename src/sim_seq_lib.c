#define TAGDUST_SIM_SEQ_LIB_IMPORT
#include "sim_seq_lib.h"

#include "tldevel.h"
#include "tlrng.h"


static inline int nuc_to_internal(const char c);

int generate_random_seq(char** seq, int* l, struct rng_state* rng)
{
        char alphabet[4] = "ACGT";
        char* s;
        int len;
        int i;

        s = *seq;
        len = *l;
        if(!s){
                MMALLOC(s , sizeof(char) *(len+1));
        }

        for(i = 0; i < len;i++){
                s[i] = alphabet[ tl_random_int(rng,4)];
        }
        s[len] = 0;

        *seq = s;
        *l = len;
        return OK;
ERROR:
        return FAIL;
}

int insert_seq(char* r, int r_len,char* insert, int i_len, struct rng_state* rng)
{
        ASSERT(r_len >= i_len,"ref too short: %d -> %d", r_len,i_len);
        int t;
        int i;
        t = tl_random_int(rng, r_len-i_len);
        for(i = 0; i < i_len;i++){
                r[t+i] = insert[i];
                //LOG_MSG("Inserting %c at %d\n", insert[i],t+i);
        }
        return OK;
ERROR:
        return FAIL;
}


int mutate_seq(char* ref,char* target,int len, float error_rate, struct rng_state* rng)
{
        char alphabet[4] = "ACGT";
        int i;
        double r;
        for(i = 0; i < len;i++){
                r = tl_random_double(rng);
                if(r < error_rate){
                        //LOG_MSG("Has error");
                        target[i] = alphabet[ tl_random_int(rng,4)];
                }else{
                        target[i] = ref[i];
                }
                //target[i] = sb->sequences[0]->seq[i];
        }
        target[len]= 0;
        return OK;
}

int seq_to_internal(char* seq, int len, uint8_t** internal, int* i_len)
{
        uint8_t* o = NULL;
        int i;
        o = *internal;
        if(!o){
                MMALLOC(o, sizeof(uint8_t) * (len+1));
        }

        for(i = 0; i < len;i++){
                o[i] = nuc_to_internal(seq[i]);
        }

        *i_len= len;
        return OK;
}


static inline int nuc_to_internal(const char c)
{
        switch (c) {
        case 'A':
        case 'a':
                return 0;
                break;
        case 'C':
        case 'c':
                return 1;
                break;
        case 'G':
        case 'g':
                return 2;
                break;
        case 'T':
        case 't':
                return 3;
                break;
        case 'N':
        case 'n':
                return 0;
                break;
        default:
                return 0;
                break;
        }
        return -1;
}
