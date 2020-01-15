
#include <stdint.h>
#include <stdio.h>

#include "tldevel.h"


#define CORRECT_IMPORT
#include "correct.h"



static inline unsigned int nuc_to_internal(const char c);

static inline uint32_t compress(uint32_t x, uint32_t m);
/* do stuff  */

int search_1M(khash_t(exact) * h , uint32_t key , int len);
int fill_exact_hash(khash_t(exact) ** hash, struct tl_seq_buffer* sb)
{
        khash_t(exact) *h = NULL;
        khiter_t k;
        int i;
        uint32_t key;
        int ret;

        h = *hash;
        if(!h){
                h = kh_init(exact);
        }

        for(i = 0; i < sb->num_seq;i++){
                key = seq_to_code(sb->sequences[i]->seq, MACRO_MIN(16, sb->sequences[i]->len));
                k = kh_put(exact, h, key, &ret);
                //is_missing = (k == kh_end(h));
                //fprintf(stdout,"put: %d   ret:%d\n",k,ret);
                if (!ret){
                        //k = kh_get(32, h, r);
                        //fprintf(stdout,"Exists - increment\n");
                        kh_value(h, k) = 0;
                }else{
                        kh_value(h, k) = i;
                        //fprintf(stdout,"set to 1 \n");
                }

        }
        i = 0;
        for (k = kh_begin(h); k != kh_end(h); ++k){
                if (kh_exist(h, k)){
                        fprintf(stdout, "%d %d ",k,kh_value(h, k));
                        key = kh_key(h,k);
                        search_1M(h , key , 16);
                        code_to_seq(key, 16);
                        i++;
                        if(i == 3){
                                break;
                        }
                        //kh_value(h, k) = 1;
                }
        }

        *hash = h;
        return OK;
ERROR:
        return FAIL;
}

int search_1M(khash_t(exact) * h , uint32_t key , int len)
{
        khiter_t k;
        uint32_t i;
        uint32_t search;
        uint32_t mask[16] = {
                0x3FFFFFFFu,
                0xCFFFFFFFu,
                0xF3FFFFFFu,
                0xFCFFFFFFu,
                0xFF3FFFFFu,
                0xFFCFFFFFu,
                0xFFF3FFFFu,
                0xFFFCFFFFu,
                0xFFFF3FFFu,
                0xFFFFCFFFu,
                0xFFFFF3FFu,
                0xFFFFFCFFu,
                0xFFFFFF3Fu,
                0xFFFFFFCFu,
                0xFFFFFFF3u,
                0xFFFFFFFCu
        };
        LOG_MSG("Searching for: ");
        code_to_seq(key, 16);
        for(i = 0; i < len;i++){
                fprintf(stdout," ");
                search = key;

                search = compress(search, mask[i]);
                code_to_seq(search, 15);
                k = kh_get(exact, h, search);
                if(k != kh_end(h)){
                        LOG_MSG("Found!!!");
                }

        }
        LOG_MSG("DONE");
        return OK;
}



static inline uint32_t compress(uint32_t x, uint32_t m)
{
        uint32_t mk,mp,mv,t;
        int i;
        x = x & m;
        mk = ~m << 1;
        for(i = 0; i < 5;i++){
                mp = mk ^(mk << 1);
                mp = mp ^(mp << 2);
                mp = mp ^(mp << 4);
                mp = mp ^(mp << 8);
                mp = mp ^(mp << 16);
                mv = mp & m;
                m = m ^ mv | ( mv >> (1 << i));
                t = x & mv;
                x = x ^ t | (t >> (1 << i));
                mk = mk & ~mp;
        }
        return x;
}


/* core functions  */


uint32_t seq_to_code(char* seq, int len)
{
        uint32_t k;
        register uint32_t i;
        k = 0UL;
        i = len;
        while(i){
                i--;
                k |= nuc_to_internal(seq[i]) << (i << 1);
                //fprintf(stdout,"%d %lx %c %d\n",i, k, seq[i], i << 1);
        }
        return k;
}

void code_to_seq(uint32_t k,int len)
{

        register uint32_t i;
        i = len;
        while(i){
                i--;
                fprintf(stdout,"%c",  "ACGT"[k & 0x3u]);
                k = k >> 2;
                //fprintf(stdout,"%d %lx %c %d\n",i, k, seq[i], i << 1);
        }
        fprintf(stdout,"\n");
}



static inline unsigned int nuc_to_internal(const char c)
{
        switch (c) {
        case 'A':
        case 'a':
                return 0u;
                break;
        case 'C':
        case 'c':
                return 1u;
                break;
        case 'G':
        case 'g':
                return 2u;
                break;
        case 'T':
        case 't':
                return 3u;
                break;
        case 'N':
        case 'n':
                return 0u;
                break;
        default:
                break;
        }
}
