
#include <stdint.h>
#include <stdio.h>

#include "tldevel.h"

#include "tllogsum.h"

#include "kstring.h"



#define CORRECT_IMPORT
#include "correct.h"

static inline unsigned int nuc_to_internal(const char c);
static inline uint32_t compress(uint32_t x, uint32_t m);
/* do stuff  */
//static int code_to_seq_bit(struct seq_bit** sb, uint32_t k,int len);
//int search_1M(khash_t(exact) * h , uint32_t key , int len);

int fill_exact_hash(khash_t(exact) ** hash, struct tl_seq_buffer* sb)
{
        khash_t(exact) *h = NULL;
        khiter_t k;
        int i;
        uint32_t key;
        int ret;

        struct qsubscore* subm = NULL;
        struct seq_bit* sbit = NULL;

        calc_score_matrix(&subm, 0.02f, 0.1f);

        h = *hash;
        if(!h){
                h = kh_init(exact);
        }

        for(i = 0; i < sb->num_seq;i++){
                key = seq_to_code(sb->sequences[i]->seq, MACRO_MIN(16, sb->sequences[i]->len));
                k = kh_put(exact, h, key, &ret);
                //is_missing = (k == kh_end(h));
                //fprintf(stdout,"put: %d   ret:%d\n",k,ret);
                if(!ret){
                        //k = kh_get(32, h, r);
                        //fprintf(stdout,"Exists - increment\n");
                        kh_value(h, k)++;
                }else{
                        kh_value(h, k) =1;
                        //fprintf(stdout,"set to 1 \n");
                }
        }
        /*i = 0;
        for (k = kh_begin(h); k != kh_end(h); ++k){
                if (kh_exist(h, k)){
                        //fprintf(stdout, "%d %d ",k,kh_value(h, k));
                        key = kh_key(h,k);
                        //search_1M(h , key , 16);
                        //code_to_seq(key, 16);
                        code_t_seq_bit(&sbit, key, 16);
                        ref_correct(h, subm, sbit,33);
                        LOG_MSG("%s %s %s", sbit->p.s, sbit->p_corr.s, sbit->q.s);
                        sbit->p.l = 0;
                        sbit->p_corr.l = 0;
                        sbit->q.l = 0;
                        i++;
                        //if(i == 1000){
                        //break;
                        //}
                        //kh_value(h, k) = 1;
                        if(i % 100000 == 0){
                                LOG_MSG("done %d",i);
                        }
                }
        }
        */
        *hash = h;
        return OK;
ERROR:
        return FAIL;
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

/*int code_to_seq_bit(struct seq_bit** sb, uint32_t k,int len)
{
        struct seq_bit* s = NULL;
        uint8_t q;
        s = *sb;
        if(!s){
                MMALLOC(s, sizeof(struct seq_bit));
                s->p.l = 0;
                s->p.m = 0;
                s->p.s = NULL;

                s->p_corr.l = 0;
                s->p_corr.m = 0;
                s->p_corr.s = NULL;

                s->q.l = 0;
                s->q.m = 0;
                s->q.s = NULL;


                s->fail = 0;
                s->type = 0;

        }
        register uint32_t i;
        q = 40;
        i = len;
        while(i){
                i--;
                kputc("ACGT"[k & 0x3u], &s->p);
                kputc(q+33, &s->q);
                q-= 2;
                 //fprintf(stdout,"%c",  "ACGT"[k & 0x3u]);
                k = k >> 2;
                //fprintf(stdout,"%d %lx %c %d\n",i, k, seq[i], i << 1);
        }
//        fprintf(stdout,"%s\n%s\n", s->p.s, s->q.s);
        //fprintf(stdout,"\n");
        *sb = s;
        return OK;
ERROR:
        return FAIL;
        }*/



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
