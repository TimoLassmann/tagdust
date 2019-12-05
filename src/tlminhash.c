
#include <limits.h>
#include <stdint.h>

#include "tldevel.h"
#include "tlrng.h"
#include "tlbitvec.h"


#define TLMINHASH_IMPORT
#include "tlminhash.h"

struct Boolean_matrix{
        struct bitvec** m;
        int n_row;
        int n_column;
};


struct minhash{
        unsigned int** sig;
        int n_signatures;
        int n_columns;
        int n_samples;
};


static int shuffle_arr_minhash(int* arr,int n, struct rng_state* rng);

struct minhash* create_min_hash(struct Boolean_matrix* bm, int num_sig,long int seed)
{
        struct minhash* min_h = NULL;
        struct bitvec*  col = NULL;
        int* list = NULL;
        int i,j,c;
        int n,m;
        int ret;
        //struct drand48_data randBuffer;

        struct rng_state* rng = NULL;
        ASSERT(bm!= NULL, "No matrix");

        RUNP(rng = init_rng(seed));

        MMALLOC(min_h, sizeof(struct minhash));
        min_h->sig = NULL;
        //min_h->a = NULL;
        //min_h->b = NULL;


        min_h->n_signatures = num_sig;
        min_h->n_columns = bm->n_column;
        min_h->n_samples = bm->n_row;

        MMALLOC( min_h->sig, sizeof(uint32_t*) * bm->n_column);

        for(i = 0; i < min_h->n_columns;i++){
                min_h->sig[i] = NULL;
                MMALLOC(min_h->sig[i],sizeof(uint32_t) * min_h->n_signatures);
                for(j = 0; j < min_h->n_signatures;j++){
                        min_h->sig[i][j] = INT_MAX;
                }
        }
        //MMALLOC(min_h->a, sizeof(uint32_t)* min_h->n_signatures);
        //MMALLOC(min_h->b, sizeof(uint32_t)* min_h->n_signatures);



        m = bm->n_column;
        n = bm->n_row;

        //fprintf(stdout,"%d columns\n", m);
        MMALLOC(list,sizeof(int)* n);
        for(i = 0; i < n;i++){
                list[i] = i;
        }
        /* Apply has functions  */
        for(c = 0;c < min_h->n_signatures;c++){

                shuffle_arr_minhash(list, n ,rng);
                for(i = 0; i < m;i++){
                        col = bm->m[i];
                        for(j = 0; j < n;j++){
                                RUN(bit_test(col, j, &ret));
                                if(ret){
                                        if(list[j]+1 < min_h->sig[i][c]){
                                                min_h->sig[i][c] = list[j]+1;
                                        }
                                }
                        }
                }
        }
        MFREE(list);
        MFREE(rng);
        return min_h;
ERROR:
        return NULL;
}



int jaccard_sim_min_multihash(struct minhash* min_h , int* S, int n, double* jac_sim, double *p_S_in_X)
{
        int i,j;
        unsigned int** m = NULL;
        double set_intersection = 0.0;
        double c,min;
        int num_samples;

        double min_stuff = 0.0;
        ASSERT(min_h != NULL, "No minhash");
        m = min_h->sig;
        num_samples = min_h->n_samples;

        for(i = 0; i < min_h->n_signatures;i++){
                min =  m[S[0]][i];
                c = 1.0;
                for(j = 1; j < n;j++){
                        min = MACRO_MIN(min, m[S[j]][i]);
                        if(m[S[j]][i] != m[S[0]][i]){
                                c =0.0;
                                //break;
                        }
                }

                set_intersection += c;

                //fprintf(stdout,"%d %d : %f\n",col_a[i],col_b[i],(double)MACRO_MIN(col_a[i],col_b[i]));
                min_stuff += min;
        }
        //fprintf(stdout,"%f\n",min_stuff);
        min_stuff = min_stuff / (double)  min_h->n_signatures;
        *jac_sim = set_intersection / (double) (min_h->n_signatures);
        //fprintf(stdout,"%f\n",min_stuff);

        *p_S_in_X = *jac_sim *  (((double)(num_samples+1) /(double) num_samples) *(  1.0 /min_stuff) -(1.0/(double)(num_samples +1)));
        //fprintf(stdout,"%d %d %f %f\n",a,b,set_intersection,  *avg_min_sig_diff);
        //LOG_MSG("samples:%d %d", num_samples,min_h->n_samples);

        return OK;
ERROR:
        return FAIL;
}


int jaccard_sim_min_hash(struct minhash* min_h , int a, int b, double* jac_sim, double *avg_min_sig_diff)
{
        int i;
        double set_intersection = 0.0;


        unsigned int* col_a = NULL;
        unsigned int* col_b = NULL;
        double min_stuff = 0.0;
        ASSERT(min_h != NULL, "No minhash");
        col_a = min_h->sig[a];
        col_b = min_h->sig[b];
        for(i = 0; i < min_h->n_signatures;i++){
                if(col_a[i] == col_b[i]){
                        set_intersection += 1.0;
                }
                min_stuff += (double)MACRO_MIN(col_a[i],col_b[i]);
        }
        *avg_min_sig_diff = (min_stuff /= (double)  min_h->n_signatures);
        *jac_sim = set_intersection / (double) (min_h->n_signatures);

        return OK;
ERROR:
        return FAIL;
}

int jaccard_sim(struct Boolean_matrix* bm, int*S , int n, double* jac_sim)
{
        int i,j,c;
        int n_row;
        int ret;
        double set_intersection = 0.0;
        double set_union = 0.0;

        ASSERT(bm!= NULL, "No matrix");
        n_row = bm->n_row;

        for(i = 0; i < n_row;i++){
                c = 0;
                for(j = 0; j < n;j++){
                        RUN(bit_test(bm->m[S[j]],i,&ret));
                        c += ret;
                }
                if(c == n){
                        set_intersection += 1.0;
                }
                if(c){
                        set_union += 1.0;
                }
        }
        *jac_sim = 0.0;

        if(set_union){
                *jac_sim = set_intersection / set_union;
        }
        return OK;
ERROR:
        return FAIL;
}

struct Boolean_matrix* init_Bmatrix( int columns,int rows)
{
        struct Boolean_matrix* bm = NULL;
        int i;

        ASSERT(rows > 0, "No rows.");
        ASSERT(columns > 0, "No columns.");

        MMALLOC(bm, sizeof(struct Boolean_matrix));
        bm->m = NULL;
        bm->n_column = columns;
        bm->n_row = rows;

        MMALLOC(bm->m, sizeof(uint32_t*) * bm->n_column);
        for(i = 0; i < bm->n_column;i++){
                bm->m[i] = NULL;
                RUN(make_bitvector(&bm->m[i], bm->n_row));
        }
        return bm;
ERROR:
        return NULL;
}

void free_minhash(struct minhash* min_h)
{
        int i;
        if(min_h){
                if(min_h->sig){
                        for(i = 0; i < min_h->n_columns;i++){
                                MFREE(min_h->sig[i]);
                        }
                        MFREE(min_h->sig);
                }
                MFREE(min_h);
        }
}

void free_Boolean_matrix(struct Boolean_matrix* bm)
{
        int i;
        if(bm){
                if(bm->m){
                        for(i = 0; i < bm->n_column;i++){
                                free_bitvector(&bm->m[i]);
                        }
                        MFREE(bm->m);
                }
                MFREE(bm);
        }
}


int shuffle_arr_minhash(int* arr,int n, struct rng_state* rng)
{
        int i,j;
        int tmp;
        long int r;


        for (i = 0; i < n - 1; i++) {

                r = tl_random_int(rng, n-i);
                //RUN(lrand48_r(randBuffer,&r));
                //int lrand48_r(struct drand48_data *buffer, long int *result);
                j = i + r;// ((int) r % (int) (n-i));
                tmp = arr[j];
                arr[j] = arr[i];
                arr[i] = tmp;
        }
        return OK;
}


