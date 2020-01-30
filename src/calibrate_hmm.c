#include <math.h>

#include "tllogsum.h"

#include "seq_stats.h"
#include "calibrate_hmm.h"


#include "arch_lib.h"
#include "arch_lib_sim.h"

//#include "hmm_model_bag.h"
//#include "core_hmm_functions.h"


#include "poahmm.h"
#include "poahmm_structs.h"
#include "init_poahmm.h"


#define NUM_TEST_SEQ 1000

struct calibrate_seq{
        uint8_t* seq;
        uint8_t* qual;
        int len;
        uint8_t type;
        float Q;
};

struct calibrate_buffer{
        struct calibrate_seq** seq;
        int num_seq;
        int max_len;
};

static int generate_test_sequences(struct  calibrate_buffer** ri_b, struct read_structure* rs  ,double mean,double stdev,struct rng_state* rng , struct alphabet* a);


static int run_scoring(struct poahmm* poahmm, struct calibrate_buffer* cb);

static int qsort_cb_q_compare(const void *a, const void *b);
static int calibrate(struct arch_library* al, struct sequence_stats_info* ssi, struct alphabet*a, int* seeds,int i_file,int i_hmm);


int calibrate_architectures(struct arch_library* al, struct read_ensembl* e,struct rng_state* main_rng)
{
        int i,j;

        int* seeds = NULL;

        MMALLOC(seeds, sizeof(int) * al->num_file);
        for(i = 0; i < al->num_file;i++){
                seeds[i] = tl_random_int(main_rng, INT32_MAX);
        }

        //LOG_MSG("Testmg:%d files", al->num_file);
#ifdef HAVE_OPENMP
#pragma omp parallel default(shared)
#pragma omp for private(i)
#endif
        for(i = 0; i < e->num_files;i++){
                j =  e->arch_to_read_assignment[i];
                //LOG_MSG("File %d HMM: %d",i,j);
                calibrate(al, e->ssi[i] ,e->a, seeds, i, j);
        }

        //for(i = 0; i < al->num_file;i++){
        //LOG_MSG("File %d threshold:%f", i, al->confidence_thresholds[i]);
        //}

        MFREE(seeds);
        return OK;
ERROR:
        return FAIL;
}

int calibrate(struct arch_library* al, struct sequence_stats_info* ssi, struct alphabet*a, int* seeds,int i_file,int i_hmm)
{
        struct global_poahmm_param* p = NULL;
        struct poahmm* poahmm = NULL;
        //struct model_bag* mb = NULL;
        struct calibrate_buffer* cb = NULL;
        struct rng_state* local_rng = NULL;
        int i;

        double TP,FP,TN,FN;
        double kappa = 0.0;
        double tmp = 0.0;

        double P_e = 0.0;

        double P_o = 0.0;

        float thres[6];

        TP = 0.0;
        FP = 0.0;
        TN = 0.0;
        FN = 0.0;

        /*RUN(init_model_bag(&mb,al->read_structure[i_hmm], si->ssi[i_file], si->a, i_hmm));
        ASSERT(mb!= NULL, "Could not init model...");
        for(i = 0; i < mb->num_models;i++){

                if(al->read_structure[i_hmm]->seg_spec[i]->extract == ARCH_ETYPE_SPLIT){
                        for(j = 1 ; j < mb->model[i]->num_hmms;j++){
                                mb->model[i]->silent_to_M[j][0] = prob2scaledprob(1.0 / (float)( mb->model[i]->num_hmms-1));
                        }
                        mb->model[i]->silent_to_M[0][0] = prob2scaledprob(0.0);
                }
        }*/

        RUNP(local_rng = init_rng(seeds[i_file]));
        RUN(generate_test_sequences(&cb, al->read_structure[i_hmm], ssi->mean_seq_len,ssi->stdev_seq_len,local_rng,a));


        //LOG_MSG("SIMLEN: %d", cb->max_len);

        //RUN(init_model_bag(&mb,al->read_structure[i_hmm], si->ssi[i_file],si->a, i_hmm));
        MMALLOC(p, sizeof(struct global_poahmm_param));
        p->min_seq_len = ssi->average_length;
        p->max_seq_len = MACRO_MAX(cb->max_len, ssi->max_seq_len);
        p->base_error = 0.05f;
        p->indel_freq = 0.1f;
        for(i =0; i < 5;i++){
                p->back[i] = ssi->background[i];
        }


        RUN(poahmm_from_read_structure(&poahmm, p, al->read_structure[i_hmm],  a));
        RUN(run_scoring(poahmm, cb));
        free_poahmm(poahmm);
        //free_model_bag(mb);

        /* Do I set thresholds at all???  */
        tmp=0.0;
        for(i = 0; i < NUM_TEST_SEQ / 2;i++){

                tmp += cb->seq[i]->Q;
        }
        tmp = tmp / (double)(NUM_TEST_SEQ /2);
        //LOG_MSG("Average Q: %f", tmp);

        if(tmp < 10){
                al->confidence_thresholds[i_file] = 0.0f;
                LOG_MSG("Arch: %s thres: %f", al->spec_line[i_hmm], al->confidence_thresholds[i_file]);
                goto SKIP;
        }

        qsort(cb->seq,NUM_TEST_SEQ, sizeof(struct calibrate_seq*), qsort_cb_q_compare);


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
                //fprintf(stdout,"%d %f %d %f\n",i, cb->seq[i]->Q, cb->seq[i]->type, powf(10.0f,-1.0f * cb->seq[i]->Q / 10));
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


        /*fprintf(stderr,"FDR:0.01: %f\n", thres[0]);
        fprintf(stderr,"FDR:0.05: %f\n", thres[1]);
        fprintf(stderr,"FDR:0.1: %f\n", thres[2]);
        fprintf(stderr,"Sen+spe: %f\n", thres[4]);

        fprintf(stderr,"Kappa: %f at %f\n",  kappa,thres[5]   );
        fprintf(stderr,"Selected Threshold: %f\n", al->confidence_thresholds[i_file]);*/
        //sprintf(param->buffer,"Selected Threshold:: %f\n", param->confidence_threshold );

        //param->messages = append_message(param->messages, param->buffer);
        //free_read_info(ri,NUM_TEST_SEQ);

        LOG_MSG("Arch: %s thres: %f", al->spec_line[i_hmm], al->confidence_thresholds[i_file]);
SKIP:
        for(i = 0; i < cb->num_seq;i++){
                MFREE(cb->seq[i]->qual);
                MFREE(cb->seq[i]->seq);
                MFREE(cb->seq[i]);
        }
        MFREE(cb->seq);
        MFREE(cb);
        free_rng(local_rng);
        MFREE(p);

        return OK;
ERROR:
        return FAIL;
}












static int generate_test_sequences(struct  calibrate_buffer** ri_b, struct read_structure* rs  ,double mean,double stdev,struct rng_state* rng , struct alphabet* a)
{

        struct calibrate_buffer* cb = NULL;

        int i,j,c;
        int t_len;

        MMALLOC(cb, sizeof(struct calibrate_buffer));
        cb->num_seq = NUM_TEST_SEQ;
        cb->seq = NULL;
        cb->max_len = 0;
        MMALLOC(cb->seq, sizeof(struct calibrate_seq*) * cb->num_seq);
        for(i = 0; i < cb->num_seq;i++){
                cb->seq[i] = NULL;
                MMALLOC(cb->seq[i], sizeof(struct calibrate_seq));
                cb->seq[i]->len = 0;
                cb->seq[i]->seq = NULL;
                cb->seq[i]->qual = NULL;
                cb->seq[i]->type = 0;
        }

        //fprintf(stdout,"Generating: %f %f", mean,stdev);

        //LOG_MSG("target: %d", average_length);
        j = NUM_TEST_SEQ /2;
        for(i = 0; i < j;i++){
                t_len = (int) round(tl_random_gaussian(rng, mean, stdev));
                RUN(emit_from_rs(rs, rng,a, &cb->seq[i]->seq,&cb->seq[i]->qual,  &cb->seq[i]->len, t_len));
                if(cb->seq[i]->len > cb->max_len){
                        cb->max_len = cb->seq[i]->len;
                }
//RUN(emit_read_sequence(mb, &cb->seq[i]->seq, &cb->seq[i]->len, t_len, rng));
                cb->seq[i]->type = 0;
        }
        for(i = j; i < NUM_TEST_SEQ ;i++){
                t_len = (int) round(tl_random_gaussian(rng, mean, stdev));
                /* makesure t_len is >= the minimum length the model can handle. */
                t_len = MACRO_MAX(10, t_len);

                MMALLOC(cb->seq[i]->seq, sizeof(uint8_t) * t_len);
                MMALLOC(cb->seq[i]->qual, sizeof(uint8_t) * t_len);
                for(c = 0; c < t_len;c++){
                        cb->seq[i]->seq[c] = tl_random_int(rng, 4);
                        cb->seq[i]->qual[c] = tl_random_gaussian(rng, 20.0, 2.0);
                        //s->label[l] = etype[spec->extract];

                }
                cb->seq[i]->len = t_len;
                cb->seq[i]->type = 1;
                if(cb->seq[i]->len > cb->max_len){
                        cb->max_len = cb->seq[i]->len;
                }
                //fprintf(stdout,"%d ", ri[i]->len);
        }

        //free_rng(rng);
        *ri_b = cb;
        return OK;
ERROR:
        return FAIL;
}

int run_scoring(struct poahmm* poahmm, struct calibrate_buffer* cb)
{
        float pbest = 0.0f;
        float Q = 0.0f;
        int i;
        uint8_t* seq = NULL;
        uint8_t* qual = NULL;
        int len;
        for(i = 0; i < NUM_TEST_SEQ;i++){

                seq = cb->seq[i]->seq;
                qual= cb->seq[i]->qual;
                len = cb->seq[i]->len;

                //LOG_MSG("LEN: %d    model: %d %d %d", len, poahmm->min_model_len, poahmm->max_model_len,poahmm->alloc_seq_len);
                viterbi_poahmm_banded(poahmm, seq, qual, len, NULL, 0);

                //RUN(backward(mb, seq,len));
                //LOG_MSG("F:%f\tR:%f", poahmm->f_score, poahmm->random_scores[len]);
                pbest = logsum(poahmm->f_score, poahmm->random_scores[len]);
                pbest = 1.0 - scaledprob2prob(poahmm->f_score - pbest);
                //RUN(forward_max_posterior_decoding(mb,seq,NULL,len));
                //RUN(random_score(mb, seq, len));


                //pbest = prob2scaledprob(0.0f);
                //fprintf(stdout,"%d %f\n",i, pbest);
                //pbest = logsum(mb->f_score + mb->bar_score,mb->r_score);
                        //pbest = logsum(pbest, mb->r_score);
                //LOG_MSG("%d %f %f %f %f sum:%f %d", len,mb->f_score, mb->b_score, mb->r_score,mb->bar_score, pbest, cb->seq[i]->type);
                //                pbest = 1.0 - scaledprob2prob(  (mb->bar_score + mb->f_score ) - pbest);
                //LOG_MSG("M:%f R:%f",1.0 - scaledprob2prob(  (mb->bar_score + mb->f_score ) - pbest),1.0 - scaledprob2prob(  (mb->r_score ) - pbest));

                //LOG_MSG("M:%f R:%f   best:%f",scaledprob2prob(mb->f_score+mb->bar_score  - pbest),scaledprob2prob(mb->r_score  - pbest),pbest);

                //pbest = 1.0 - scaledprob2prob(mb->f_score+mb->bar_score  - pbest);
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
