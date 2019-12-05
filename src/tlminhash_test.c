
#include <stdio.h>
#include <math.h>
#include "tldevel.h"

#include "tlminhash.h"

#include "tlbitvec.h"

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



int jaccard_sim(struct Boolean_matrix* bm, int*S , int n, double* jac_sim);
int jaccard_sim_min_hash(struct minhash* min_h , int a, int b, double* jac_sim, double *avg_min_sig_diff);

static struct Boolean_matrix* init_random_Bmatrix(int columns,int rows,  double alpha,struct drand48_data *rd);
int print_Boolean_matrix(struct Boolean_matrix* bm);
int print_minhash_signatures(struct minhash* min_h);


struct Boolean_matrix* init_random_Bmatrix( int columns,int rows, double alpha,struct drand48_data* rd)
{
        struct Boolean_matrix* bm = NULL;
        double  r;
        int i,j;

        //long int seed = 0;

        ASSERT(rows > 0, "No rows.");
        ASSERT(columns > 0, "No columns.");
        ASSERT(alpha > 0.0, "No alpha - the matrix will be empty.");


        //seed  =  (long int) (time(NULL) * (42));

        //srand48(seed);
        RUNP( bm = init_Bmatrix(columns, rows));
        /*MMALLOC(bm, sizeof(struct Boolean_matrix));
        bm->m = NULL;
        bm->n_column = columns;
        bm->n_row = rows;

        MMALLOC(bm->m, sizeof(uint32_t*) * bm->n_column);

        for(i = 0; i < bm->n_column;i++){
                bm->m[i] = NULL;
                RUNP(bm->m[i] = make_bitvector(bm->n_row));
                }*/


        for(i = 0; i < bm->n_column;i++){
                for(j = 0; j < bm->n_row;j++){
                        //r = drand48();
                        RUN(drand48_r(rd, &r));
                        if(r <= alpha){
                                bit_set(bm->m[i], j);

                        }
                }
        }

        return bm;
ERROR:
        return NULL;
}

int print_Boolean_matrix(struct Boolean_matrix* bm)
{
        int i,j;
        int res;
        for(i = 0; i < bm->n_column;i++){
                for(j = 0; j < bm->n_row;j++){
                        RUN(bit_test(bm->m[i], j, &res));
                        fprintf(stdout,"%d",res);
                        //fprintf(stdout,"%d",bit_test(bm->m[i], j));
                }
                fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
        return OK;
ERROR:
        return FAIL;
}

int print_minhash_signatures(struct minhash* min_h)
{
        int i,j;
        fprintf(stdout,"\t");
        for(i = 0; i < min_h->n_columns  ;i++){
                fprintf(stdout,"Col:%d\t",i);
        }
        fprintf(stdout,"\n");
        for(j = 0; j < min_h->n_signatures;j++){
                fprintf(stdout,"h%d\t",j);
                for(i = 0; i < min_h->n_columns  ;i++){
                        fprintf(stdout,"%d ",  min_h->sig[i][j]);
                }
                fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
        return OK;
}

int main (int argc,char * const argv[])
{
        fprintf(stdout,"Hello world\n");

        struct Boolean_matrix* bm = NULL;
        struct minhash* min_h = NULL;
        double sim;
        double sim_min = 0.0;
        int num_samples = 1000;
        int i;
        int trials = 10;
        double s1 = 0.0;
        double s2 = 0.0;
        double s1_p = 0.0;
        double s2_p= 0.0;

        double diff;
        double diff_p;
        double p_S_in_X = 0.0;
        double alpha = 0.999;

        int num_hash_functions = 200;
        int iter;
        struct drand48_data randBuffer;
        int* index =NULL;
        int S_size = 6;

        MMALLOC(index, sizeof(int) * S_size);
        for(i = 0; i < S_size;i++){
                index[i] = i;
        }


        srand48_r(42, &randBuffer);

        for(iter = 0; iter < trials;iter++){
                /* S_size is just for simulation - this should be the number of variables... */
                RUNP(bm = init_random_Bmatrix(S_size,num_samples,alpha, &randBuffer));
                //RUN(print_Boolean_matrix(bm));
                min_h = create_min_hash(bm, num_hash_functions, 0);
                RUN(jaccard_sim(bm,index, S_size, &sim));
                jaccard_sim_min_multihash(min_h, index, S_size, &sim_min,&p_S_in_X);

                diff = fabs(sim-sim_min);
                s1 += diff;
                s2 += diff * diff;
                diff_p = fabs(p_S_in_X-  pow(alpha,(double)S_size));
                s1_p += diff_p;
                s2_p += diff_p * diff_p;
                fprintf(stdout,"%f %f delta: %f\t",sim,sim_min, diff);
                fprintf(stdout,"P seeing: %f (%f) delta: %f\n",p_S_in_X,   pow(alpha, (double)S_size),diff_p);
                free_minhash(min_h);
                free_Boolean_matrix(bm);

        }

        s2 = sqrt(((double) trials * s2 - s1 * s1)/ ((double) trials * ((double) trials -1.0)));
        s1 = s1 / (double) trials;
        fprintf(stdout,"mean: %f stdev:%f\n", s1,s2);

s2_p = sqrt(((double) trials * s2_p - s1_p * s1_p)/ ((double) trials * ((double) trials -1.0)));
        s1_p = s1_p / (double) trials;
        fprintf(stdout,"mean: %f stdev:%f\n", s1_p,s2_p);
//        free_bitvector(&bm);
        MFREE(index);
        return EXIT_SUCCESS;
ERROR:
        return EXIT_FAILURE;
}
