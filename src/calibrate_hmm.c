#include <math.h>

#include "tllogsum.h"


#include "calibrate_hmm.h"


#include "arch_lib.h"
#include "seq_stats.h"
#include "hmm_model_bag.h"
#include "core_hmm_functions.h"


#define NUM_TEST_SEQ 100

struct calibrate_seq{
        uint8_t* seq;
        int len;
        uint8_t type;
        float Q;
};

struct calibrate_buffer{
        struct calibrate_seq** seq;
        int num_seq;
};



static int generate_test_sequences(struct  calibrate_buffer** ri_b, struct model_bag* mb,int average_length,struct rng_state* rng);


static int run_scoring(struct model_bag* mb, struct calibrate_buffer* cb);
static int qsort_cb_q_compare(const void *a, const void *b);
static int calibrate(struct arch_library* al, struct seq_stats* si,int* seeds,int i_file,int i_hmm);


int calibrate_architectures(struct arch_library* al, struct seq_stats* si,struct rng_state* main_rng)
{
        int i,j;

        int* seeds = NULL;

        MMALLOC(seeds, sizeof(int) * al->num_file);
        for(i = 0; i < al->num_file;i++){
                seeds[i] = tl_random_int(main_rng, INT32_MAX);
        }


#pragma omp parallel default(shared)
#pragma omp for private(i)
        for(i = 0; i < al->num_file;i++){
                j = al->arch_to_read_assignment[i];
                LOG_MSG("File %d HMM: %d",i,j);
                calibrate(al, si,seeds, i, j);

        }

        for(i = 0; i < al->num_file;i++){
                LOG_MSG("File %d threshold:%f", i, al->confidence_thresholds[i]);
        }

        MFREE(seeds);
        return OK;
ERROR:
        return FAIL;
}

int calibrate(struct arch_library* al, struct seq_stats* si,int* seeds,int i_file,int i_hmm)
{
        struct model_bag* mb = NULL;
        struct calibrate_buffer* cb = NULL;
        struct rng_state* local_rng = NULL;
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

        RUNP(local_rng = init_rng(seeds[i_file]));
        RUN(generate_test_sequences(&cb, mb, si->ssi[i_file]->average_length,local_rng));
        free_model_bag(mb);


        RUNP(mb = init_model_bag(al->read_structure[i_hmm], si->ssi[i_file], i_hmm));

        RUN(run_scoring(mb, cb));


        free_model_bag(mb);
        qsort(cb->seq,NUM_TEST_SEQ, sizeof(struct calibrate_seq*), qsort_cb_q_compare);

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
                if(cb->seq[i]->type){
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
                        thres[0] = cb->seq[i]->Q;
                }else if(FP /(FP+ TP) < 0.05){
                        thres[1] = cb->seq[i]->Q;
                }else if(FP /(FP+ TP) < 0.1){
                        thres[2] = cb->seq[i]->Q;
                }

                if(sensitivity + specificity > thres[3]){
                        thres[3] = specificity + sensitivity;
                        thres[4] = cb->seq[i]->Q;
                }
                P_e = ((TP+FN) / (double) NUM_TEST_SEQ) * ((TP+FP) / (double)NUM_TEST_SEQ) +  ( ((FP+TN) / (double)NUM_TEST_SEQ ) * ((FN+TN) / (double)NUM_TEST_SEQ));
                P_o =(TP+TN)/(double)NUM_TEST_SEQ ;

                tmp = (P_o - P_e) / (1.0 - P_e);


                if(tmp > kappa ){

                        kappa = tmp;
                        thres[5] = cb->seq[i]->Q;
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
        //free_read_info(ri,NUM_TEST_SEQ);
        for(i = 0; i < cb->num_seq;i++){
                MFREE(cb->seq[i]->seq);
                MFREE(cb->seq[i]);
        }
        MFREE(cb->seq);
        MFREE(cb);
        free_rng(local_rng);
        return OK;
ERROR:
        return FAIL;
}

int generate_test_sequences(struct  calibrate_buffer** ri_b, struct model_bag* mb,int average_length,struct rng_state* rng)
{

        struct calibrate_buffer* cb = NULL;

        int i,j;

        MMALLOC(cb, sizeof(struct calibrate_buffer));
        cb->num_seq = NUM_TEST_SEQ;
        cb->seq = NULL;

        MMALLOC(cb->seq, sizeof(struct calibrate_seq*) * cb->num_seq);
        for(i = 0; i < cb->num_seq;i++){
                cb->seq[i] = NULL;
                MMALLOC(cb->seq[i], sizeof(struct calibrate_seq));
                cb->seq[i]->len = 0;
                cb->seq[i]->seq = NULL;
                cb->seq[i]->type = 0;
        }

         //LOG_MSG("target: %d", average_length);
        j = NUM_TEST_SEQ /2;
        for(i = 0; i < j;i++){
                RUN(emit_read_sequence(mb, &cb->seq[i]->seq, &cb->seq[i]->len, average_length, rng));
                cb->seq[i]->type = 0;
        }
        for(i = j; i < NUM_TEST_SEQ ;i++){
                RUN(emit_random_sequence(mb, &cb->seq[i]->seq, &cb->seq[i]->len, average_length, rng));
                cb->seq[i]->type = 1;

                //fprintf(stdout,"%d ", ri[i]->len);
        }

        //free_rng(rng);
        *ri_b = cb;
        return OK;
ERROR:
        return FAIL;
}

int run_scoring(struct model_bag* mb, struct calibrate_buffer* cb)
{
        float pbest = 0.0f;
        float Q = 0.0f;
        int i;
        uint8_t* seq = NULL;
        int len;
        for(i = 0; i < NUM_TEST_SEQ;i++){

                seq = cb->seq[i]->seq;
                len = cb->seq[i]->len;
                RUN(backward(mb, seq,len));
                RUN(forward_max_posterior_decoding(mb,seq,NULL,len));
                RUN(random_score(mb, seq, len));
                pbest = prob2scaledprob(0.0f);
                //fprintf(stdout,"%d %f\n",i, pbest);
                pbest = logsum(pbest, mb->f_score);
                pbest = logsum(pbest, mb->r_score);

                pbest = 1.0 - scaledprob2prob(  (mb->bar_score + mb->f_score ) - pbest);

                if(!pbest){
                        Q = 40.0;
                }else if(pbest == 1.0){
                        Q = 0.0;
                }else{
                        Q = -10.0 * log10(pbest) ;
                }
                cb->seq[i]->Q = Q;
        }
        return OK;
ERROR:
        return FAIL;
}



int qsort_cb_q_compare(const void *a, const void *b)
{
        const struct calibrate_seq **elem1 = (const struct calibrate_seq**) a;
        const struct calibrate_seq **elem2 = (const struct calibrate_seq**) b;
        if ( (*elem1)->Q > (*elem2)->Q){
                return -1;
        }else if ((*elem1)->Q < (*elem2)->Q){
                return 1;
        }else{
                return 0;
        }
}
