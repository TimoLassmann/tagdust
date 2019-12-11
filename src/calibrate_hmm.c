#include <math.h>

#include "tllogsum.h"
#include "tlrng.h"

#include "calibrate_hmm.h"


#include "arch_lib.h"
#include "seq_stats.h"
#include "hmm_model_bag.h"
#include "core_hmm_functions.h"


#define NUM_TEST_SEQ 1000


static int generate_test_sequences(struct read_info*** ri_b, struct model_bag* mb,int average_length);
static int run_scoring(struct model_bag* mb, struct read_info** ri);
static int qsort_ri_mapq_compare(const void *a, const void *b);
static int calibrate(struct arch_library* al, struct seq_stats* si,int i_file,int i_hmm);

int calibrate_architectures(struct arch_library* al, struct seq_stats* si)
{
        int i,j;
#pragma omp parallel default(shared)
#pragma omp for private(i)
        for(i = 0; i < al->num_file;i++){
                j = al->arch_to_read_assignment[i];
                LOG_MSG("File %d HMM: %d",i,j);
                calibrate(al, si, i, j);

        }

        for(i = 0; i < al->num_file;i++){
                LOG_MSG("File %d threshold:%f", i, al->confidence_thresholds[i]);
        }


        return OK;
}


int calibrate(struct arch_library* al, struct seq_stats* si,int i_file,int i_hmm)
{
        struct model_bag* mb = NULL;
        struct read_info** ri = NULL;

        int i,j;
        double TP,FP,TN,FN;
        TP = 0.0;
        FP = 0.0;
        TN = 0.0;
        FN = 0.0;

        RUNP(mb = init_model_bag(al->read_structure[i_hmm], si->ssi[i_file], i_hmm));

        for(i = 0; i < mb->num_models;i++){
                if(al->read_structure[i_hmm]->type[i] == 'B'){
                        for(j = 0 ; j < mb->model[i]->num_hmms-1;j++){
                                mb->model[i]->silent_to_M[j][0] = prob2scaledprob(1.0 / (float)( mb->model[i]->num_hmms-1));
                        }
                        mb->model[i]->silent_to_M[mb->model[i]->num_hmms-1][0] = prob2scaledprob(0.0);
                }
                if(al->read_structure[i_hmm]->type[i] == 'S'){
                        for(j = 0 ; j < mb->model[i]->num_hmms-1;j++){
                                mb->model[i]->silent_to_M[j][0] = prob2scaledprob(1.0 / (float)( mb->model[i]->num_hmms-1));
                        }
                        mb->model[i]->silent_to_M[mb->model[i]->num_hmms-1][0] = prob2scaledprob(0.0);
                }
        }

        RUN(generate_test_sequences(&ri, mb, si->ssi[i_file]->average_length));
        free_model_bag(mb);


        RUNP(mb = init_model_bag(al->read_structure[i_hmm], si->ssi[i_file], i_hmm));

        RUN(run_scoring(mb, ri));


        free_model_bag(mb);
        qsort(ri,NUM_TEST_SEQ, sizeof(struct read_info*), qsort_ri_mapq_compare);

        double kappa = 0.0;
        double tmp = 0.0;

        double P_e = 0.0;

        double P_o = 0.0;



        float thres[6];
        thres[0] = 1000.0;
        thres[1] = 1000.0;
        thres[2] = 1000.0;
        thres[3] = 0.0;
        thres[4] = 1000.0;
        thres[5] = 1000.0;

        float sensitivity, specificity;

        FN = NUM_TEST_SEQ / 2;
        TN = NUM_TEST_SEQ / 2;
        for(i = 0; i < NUM_TEST_SEQ;i++){
                if(ri[i]->read_type){
                        FP += 1.0;
                        TN -= 1.0;
                }else{
                        TP += 1.0;
                        FN -= 1.0;
                }
                //if(i < 100){
                //fprintf(stderr,"%d	%f	%d\n",i,ri[i]->mapq, ri[i]->read_type);
                //}
                sensitivity = TP/( TP + FN );
                specificity =  TN / ( TN + FP);

                if(FP /(FP+ TP) < 0.01){
                        thres[0] =  ri[i]->mapq;
                }else if(FP /(FP+ TP) < 0.05){
                        thres[1] = ri[i]->mapq;
                }else if(FP /(FP+ TP) < 0.1){
                        thres[2] = ri[i]->mapq;
                }

                if(sensitivity + specificity > thres[3]){
                        thres[3] = specificity + sensitivity;
                        thres[4] = ri[i]->mapq;
                }
                P_e = ((TP+FN) / (double) NUM_TEST_SEQ) * ((TP+FP) / (double)NUM_TEST_SEQ) +  ( ((FP+TN) / (double)NUM_TEST_SEQ ) * ((FN+TN) / (double)NUM_TEST_SEQ));
                P_o =(TP+TN)/(double)NUM_TEST_SEQ ;

                tmp = (P_o - P_e) / (1.0 - P_e);


                if(tmp > kappa ){

                        kappa = tmp;
                        thres[5] = ri[i]->mapq;
                        //fprintf(stderr,"%d	KAPPA:%f	%f\n",i,kappa,ri[i]->mapq);
                }

        }

        if(thres[4] < 20){
                al->confidence_thresholds[i_file] = thres[4];
        }else{
                al->confidence_thresholds[i_file] =  20;
        }

        fprintf(stderr,"FDR:0.01: %f\n", thres[0]);
        fprintf(stderr,"FDR:0.05: %f\n", thres[1]);
        fprintf(stderr,"FDR:0.1: %f\n", thres[2]);
        fprintf(stderr,"Sen+spe: %f\n", thres[4]);

        fprintf(stderr,"Kappa: %f at %f\n",  kappa,thres[5]   );

        fprintf(stderr,"Selected Threshold: %f\n", al->confidence_thresholds[i_file]);
        //sprintf(param->buffer,"Selected Threshold:: %f\n", param->confidence_threshold );

        //param->messages = append_message(param->messages, param->buffer);
        free_read_info(ri,NUM_TEST_SEQ);
        return OK;
ERROR:
        return FAIL;
}

int generate_test_sequences(struct read_info*** ri_b, struct model_bag* mb,int average_length)
{
        struct read_info** ri = NULL;
        struct rng_state* rng = NULL;
        int i,j;


        RUNP(rng = init_rng(42));
        RUNP(ri = malloc_read_info(ri,NUM_TEST_SEQ));
        RUNP(ri = clear_read_info(ri,NUM_TEST_SEQ));
        //LOG_MSG("target: %d", average_length);
        j = NUM_TEST_SEQ /2;
        for(i = 0; i < j;i++){
                RUN(emit_read_sequence(mb, ri[i], average_length, rng));
                ASSERT(ri[i]->len < average_length,"Emit failed! %d < %d",ri[i]->len,average_length);
                //print_seq(ri[i], stdout);
                //fprintf(stdout,"%d ", ri[i]->len);
                ri[i]->read_type = 0;
        }
        for(i = j; i < NUM_TEST_SEQ ;i++){
                RUN(emit_random_sequence(mb, ri[i], average_length, rng));
                ASSERT(ri[i]->len < average_length,"RAND emit failed!");
                //print_seq(ri[i], stdout);
                ri[i]->read_type = 1;
                //fprintf(stdout,"%d ", ri[i]->len);
        }
        //fprintf(stdout,"\n");
        //LOG_MSG("Emit Done");
        free_rng(rng);
        *ri_b = ri;
        return OK;
ERROR:
        return FAIL;
}

int run_scoring(struct model_bag* mb, struct read_info** ri)
{
        float pbest = 0.0f;
        float Q = 0.0f;
        int i;

        for(i = 0; i < NUM_TEST_SEQ;i++){
                RUN(backward(mb, ri[i]->seq,ri[i]->len));
                RUN(forward_max_posterior_decoding(mb, ri[i], ri[i]->seq ,ri[i]->len));

                pbest = prob2scaledprob(0.0f);
                //fprintf(stdout,"%d %f\n",i, pbest);
                pbest = logsum(pbest, mb->f_score);
                pbest = logsum(pbest, mb->r_score);

                pbest = 1.0 - scaledprob2prob(  (ri[i]->bar_prob + mb->f_score ) - pbest);

                if(!pbest){
                        Q = 40.0;
                }else if(pbest == 1.0){
                        Q = 0.0;
                }else{
                        Q = -10.0 * log10(pbest) ;
                }

                ri[i]->mapq = Q;
        }
        return OK;
ERROR:
        return FAIL;
}



int qsort_ri_mapq_compare(const void *a, const void *b)
{
        const struct read_info **elem1 = (const struct read_info**) a;
        const struct read_info **elem2 = (const struct read_info**) b;
        if ( (*elem1)->mapq > (*elem2)->mapq){
                return -1;
        }else if ((*elem1)->mapq < (*elem2)->mapq){
                return 1;
        }else{
                return 0;
        }
}
