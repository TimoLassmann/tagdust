#include "misc.h"
#include "nuc_code.h"


#include <ctype.h>
#include <string.h>
#include <emmintrin.h>
#include <mmintrin.h>
#include <math.h>


#define INV_SQRT_2PI 0.3989422804014327

#ifndef _MM_ALIGN16
#ifdef __GNUC__
#define _MM_ALIGN16 __attribute__((aligned (16)))
#endif
#ifdef __MSVC__
#define _MM_ALIGN16 __declspec(align(16))
#endif
#endif

int byg_count(char* pattern,char*text)
{
        int Tc;
        int i  = 0;
        int s = 0;
        int T[256];
        for (i = 0;i < 256;i++){
                T[i] = 0;
        }

        int m = (int) strlen(pattern);
        int n = (int) strlen(text);
        int count = 0;

        if(m > n){
                return -1;
        }
        int mb = (1 << (m-1));

        for (i= 0;i < m;i++){
                T[(int)toupper(pattern[i])] |= (1 << i);
        }

        for (i = 0;i < n;i++){
                s <<= 1;
                s |= 1;
                Tc = T[(int)toupper(text[i])];
                s &= Tc;
                if(s & mb){
                        count++;
                }
        }
        return count;
}

int byg_end(const char* pattern,const char*text)
{
        const char* tmp = 0;
        int Tc;
        int i  = 0;
        int s = 0;
        int T[256];
        for (i = 0;i < 256;i++){
                T[i] = 0;
        }

        int m = (int)strlen(pattern);
        int n = (int)strlen(text);
        if (m > n){
                i = m;
                m = n;
                n = i;
                tmp = text;
                text = pattern;
                pattern = tmp;
        }

        int mb = (1 << (m-1));

        for (i= 0;i < m;i++){
                T[(int)pattern[i]] |= (1 << i);
        }

        for (i = 0;i < n;i++){
                s <<= 1;
                s |= 1;
                if(!text[i]){
                        return -1;
                }
                Tc = T[(int)text[i]];
                s &= Tc;
                if(s & mb){
                        return i+1;
                }
        }
        return 0;
}


int bpm(const  char* t,const  char* p,int n,int m)
{
        register unsigned long int i;//,c;
        unsigned long int diff;
        unsigned long int B[255];
        if(m > 31){
                m = 31;
        }

        unsigned long int k = m;
        //static int counter = 0;
        register unsigned long int VP,VN,D0,HN,HP,X;

        long int MASK = 0;
        //c = 0;

        diff = m;

        for(i = 0; i < 255;i++){
                B[i] = 0;
        }

        for(i = 0; i < m;i++){
                B[(int)(p[i] )] |= (1ul << i);
        }

        //c = 0;
        VP = 0xFFFFFFFFFFFFFFFFul;
        VN = 0ul;
        m--;
        MASK = 1ul << m;

        for(i = 0; i < n;i++){
                X = (B[(int)(t[i])  ] | VN);
                D0 = ((VP+(X&VP)) ^ VP) | X ;
                HN = VP & D0;
                HP = VN | ~(VP | D0);
                X = HP << 1ul;
                VN = X & D0;
                VP = (HN << 1ul) | ~(X | D0);
                diff += (HP & MASK) >> m;
                diff -= (HN & MASK) >> m;
                if(diff < k){
                        k = diff;
                        //if(k <= limit){
                        //	return (int)k;
                        //}

                }

        }
        return (int)k;
}

int bitcount64(long long int i)
{
        i = i - ((i >> 1) & 0x5555555555555555);
        i = (i & 0x3333333333333333) + ((i >> 2) & 0x3333333333333333);
        return (((i + (i >> 4)) & 0xF0F0F0F0F0F0F0F) * 0x101010101010101) >> 56;
}


double gaussian_pdf(double x, double m,double s)
{
        double a = (x-m) / s;
        return INV_SQRT_2PI / s *exp(-0.5 * a * a);
}




int validate_bpm_sse(unsigned char**  query, int* query_lengths,unsigned char* t,int n,int num)
{
        int i,j;
        int len = 0;
        int new_len = 0;

        unsigned int _MM_ALIGN16 nuc[16];

        unsigned int _MM_ALIGN16 lengths[4];

        __m128i VP,VN,D0,HN,HP,X,MASK,K,NOTONE,diff,zero,one;
        __m128i* nuc_p;
        __m128i xmm1,xmm2;

        for(i = 0; i < 16;i++){
                nuc[i] = 0ul;
        }

        for(i = 0; i < num;i++){
                len = query_lengths[i];



                new_len = 0;
                for(j = 0; j < len;j++){
                        if(query[i][j] != 65){
                                nuc[((int)(query[i][j] & 0x3u) << 2) + i] |=  (1ul << j);// (unsigned long int)(len-1-j));
                                new_len++;
                        }
                }

                if(new_len > 31){
                        new_len = 31;
                }

                lengths[i] = new_len;


        }
        nuc_p = (__m128i*) nuc;
        zero = _mm_set1_epi32(0);
        one = _mm_set1_epi32(1);
        diff =  _mm_load_si128 ( (__m128i*) lengths );  // _mm_set1_epi32(m);
        VP =  _mm_set1_epi32(0xFFFFFFFFu);
        VN =  _mm_set1_epi32(0);
        NOTONE =  _mm_set1_epi32(0xFFFFFFFF);
        K =  _mm_set1_epi32(0x7FFFFFFF);

        for(i = 0; i< 4;i++){
                lengths[i]--;
                lengths[i] = 1 << lengths[i];
        }

        MASK =  _mm_load_si128 ( (__m128i*) lengths ); //  _mm_set1_epi32(1ul << m);

        for(i = 0; i < n ;i++){
                //fprintf(stderr,"%c",*t + 65);
                X = _mm_or_si128 (*(nuc_p +( (int)(*t)  & 0x3u) ) , VN);
                //X = (B[(int) *t] | VN);
                xmm1 = _mm_and_si128(X, VP);
                xmm2 = _mm_add_epi32(VP ,xmm1);
                xmm1 = _mm_xor_si128 (xmm2, VP);
                D0 = _mm_or_si128(xmm1, X);
                //D0 = ((VP+(X&VP)) ^ VP) | X ;
                HN = _mm_and_si128(VP, D0);
                //HN = VP & D0;
                xmm1 = _mm_or_si128(VP, D0);
                xmm2 = _mm_andnot_si128 (xmm1,NOTONE);
                HP = _mm_or_si128(VN, xmm2);
                //HP = VN | ~(VP | D0);
                X = _mm_slli_epi32(HP,1);
                //X = HP << 1ul;
                VN = _mm_and_si128(X, D0);
                //VN = X & D0;
                xmm1 = _mm_slli_epi32(HN,1);
                xmm2 = _mm_or_si128(X, D0);
                xmm2 = _mm_andnot_si128 (xmm2,NOTONE);
                VP = _mm_or_si128(xmm1, xmm2);
                //VP = (HN << 1ul) | ~(X | D0);
                xmm1 = _mm_and_si128(HP, MASK);
                xmm2 = _mm_cmpgt_epi32(xmm1, zero);
                diff = _mm_add_epi32(diff , _mm_and_si128( xmm2, one));
                //diff += (HP & MASK) >> m;
                xmm1 = _mm_and_si128(HN, MASK);
                xmm2 = _mm_cmpgt_epi32(xmm1, zero);
                diff = _mm_sub_epi32(diff,  _mm_and_si128( xmm2, one));
                //diff -= (HN & MASK) >> m;
                xmm1 = _mm_cmplt_epi32(diff, K);
                xmm2 = _mm_and_si128(xmm1, diff);
                K = _mm_or_si128(xmm2, _mm_andnot_si128  (xmm1,K));
                t++;
        }

        _mm_store_si128 ((__m128i*) query_lengths, K);

        return 1;
}


int bpm_check_error(const unsigned char* t,const unsigned char* p,int n,int m)
{
        register unsigned long int i;//,c;
        unsigned long int diff;
        unsigned long int B[5];

        int new_len = 0;
        unsigned long int k = m;
        //static int counter = 0;
        register unsigned long int VP,VN,D0,HN,HP,X;

        long int MASK = 0;
        //c = 0;

        diff = m;

        for(i = 0; i < 5;i++){
                B[i] = 0;
        }

        for(i = 0; i < m;i++){
                if(p[i] != 65){
                        B[(int)(p[i] & 0x3)] |= (1ul << i);
                        new_len++;
                }

        }
        if(new_len > 31){
                new_len = 31;
        }
        m = new_len;
        k =new_len;


        //c = 0;
        VP = 0xFFFFFFFFFFFFFFFFul;
        VN = 0ul;
        m--;
        MASK = 1ul << m;
        for(i = 0; i < n;i++){
                X = (B[(int)(t[i] &0x3)  ] | VN);
                D0 = ((VP+(X&VP)) ^ VP) | X ;
                HN = VP & D0;
                HP = VN | ~(VP | D0);
                X = HP << 1ul;
                VN = X & D0;
                VP = (HN << 1ul) | ~(X | D0);
                diff += (HP & MASK) >> m;
                diff -= (HN & MASK) >> m;
                if(diff < k){
                        k = diff;
                        //fprintf(stderr,"%ld	%ld\n",i,k);
                        //if(k <= limit){
                        //	return (int)k;
                        //}

                }
        }
        return (int)k;
}


unsigned char* reverse_complement(unsigned char* p,int len)
{
        unsigned char* tmp = 0;

        MMALLOC(tmp,sizeof(unsigned char)*(len +2));
        int i,c;
        c = 0;
        for(i =len-1; i >= 0;i--){
                if(p[i]== 65){
                        tmp[c] = 65;
                }else{
                        tmp[c] = rev_nuc_code[(int)p[i]];
                }
                c++;
        }
        tmp[c] = 0;
        for(i= 0; i < len;i++){
                p[i] = tmp[i];

        }
        MFREE(tmp);
        return p;
ERROR:
        //KSLIB_MESSAGE(status,"Something wrong in reverse_complement.\n");
        return NULL;
}
