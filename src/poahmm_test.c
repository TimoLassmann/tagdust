
#include <stdint.h>
#include <string.h>
#include "tldevel.h"

#include "arch_lib.h"

#include "tlrng.h"
#include "tllogsum.h"
#include "tlalphabet.h"

#include "poahmm.h"
#include "poahmm_structs.h"

#include "init_poahmm.h"

#include "sim_seq_lib.h"

struct sim_seq{
        char* buffer;
        char* seq;
        char* label;
        uint8_t* qual;
        int len;
        int alloc_len;
};

struct shared_sim_data{
        struct alphabet* a;
        struct rng_state* main_rng;
        struct global_poahmm_param* gp;
        struct arch_library* al;
        struct poahmm* poahmm;
        struct sim_seq* seq;
};

static int init_sim_data(struct shared_sim_data** sim_data, int seed);
static void free_sim_data(struct shared_sim_data* sd);

static int sim_seq_from_read_struct(struct sim_seq** simseq, struct read_structure* rs,int sim_len, int base_q, struct rng_state* rng);
static void free_sim_seq(struct sim_seq* s);

int score_labelling(struct poahmm* poahmm, uint32_t* path,char* label, float* score);
/* for use in tagdust */
//int poahmm_from_read_structure(struct poahmm** poahmm,struct global_poahmm_param* p, struct read_structure* rs,struct alphabet* a);

int single_seq_test(void);
int arch_test(void);


//int init_nodes_from_read_structure(struct poahmm* poahmm, struct read_structure* rs, struct alphabet* a, int random);

int print_poa_graph(struct poahmm* poahmm);
int print_path(struct poahmm* poahmm, uint32_t* path,char* seq,char* label);

int print_poahmm_param(struct poahmm* poahmm);
int poahmm_to_dot(struct poahmm* poahmm,char* filename);

int test_model_init(struct shared_sim_data* sd);
int test_simple_N_plus_arch(struct shared_sim_data* sd);
int test_banded(struct shared_sim_data* sd);
int test_indel(struct shared_sim_data* sd);
int single_tests(struct shared_sim_data* sd);

int run_single_test(struct shared_sim_data* sd, char** arch,int n, char* seq,char* label, float mean_q, float stdev_q);

int main(void)
{
        struct shared_sim_data* sim_data=NULL;
        init_logsum();
        RUN(init_sim_data(&sim_data, 0));

        LOG_MSG("-----------------------------------");
        LOG_MSG("model init test");
        LOG_MSG("-----------------------------------");
        RUN(test_model_init(sim_data));
        LOG_MSG("-----------------------------------");
        LOG_MSG("banded test");
        LOG_MSG("-----------------------------------");
        RUN(test_banded(sim_data));
        LOG_MSG("-----------------------------------");
        LOG_MSG("Simple plus test");
        LOG_MSG("-----------------------------------");
        RUN(test_simple_N_plus_arch(sim_data));
        //RUN(test_indel(sim_data));
        //RUN(single_seq_test());
        //RUN(arch_test());
        single_tests(sim_data);


        free_sim_data(sim_data);
        return EXIT_SUCCESS;
ERROR:
        return EXIT_FAILURE;
}


int single_tests(struct shared_sim_data* sd)
{
        char** arch = NULL;
        int i;

        MMALLOC(arch, sizeof(char*) * 10);
        for(i = 0; i < 10;i++){
                arch[i] = NULL;
                MMALLOC(arch[i], sizeof(char) * 50);
        }


        snprintf(arch[0], 50, "E:NNNNNN");
        run_single_test(sd, arch, 1, "AAAAAA","EEEEEE", 40, 1);


        snprintf(arch[0], 50, "S:NNNNNN");
        snprintf(arch[1], 50, "R:N{2,4}");
        run_single_test(sd, arch, 2, "ACGTACGGG","SSSSSSEEE", 40, 1);



        for(i = 0; i < 10;i++){
                MFREE(arch[i]);
        }
        MFREE(arch);
        return OK;
ERROR:
        return FAIL;

}

int test_model_init(struct shared_sim_data* sd)
{
        struct poahmm* poahmm = NULL;

        char* in[] = {
                "E:N"
                //"BC:A:ACAGTG,ACTTGA,TTAGGC"
                /* "L:ACGT", */
                /* "R2:L:AACCTT", */
                /* "R2:E:N+", */
                /* "A:AAA,TTT", */
                /* "R:AAACCCGGGTT" */
        };

        int active_read_structure;
        int size;
        size = sizeof(in) / sizeof(char*);
        RUN(read_arch_into_lib(sd->al, in, size));

        active_read_structure = sd->al->num_arch-1;

        sd->gp->min_seq_len = 30;
        sd->gp->max_seq_len = 30;

        poahmm = sd->poahmm;

        LOG_MSG("Current poa size: nodes: %d seq_len:%d", poahmm->alloced_num_nodes , poahmm->alloc_seq_len);
        RUN(poahmm_from_read_structure(&poahmm, sd->gp, sd->al->read_structure[active_read_structure ], sd->a));

        if(poahmm){
                poahmm_to_dot(poahmm, "test_singleN.dot");
        }
        //free_poahmm(poahmm);
        //poahmm = NULL;
        return OK;
ERROR:
        free_poahmm(poahmm);
        return FAIL;
}
/* run various tests  */
/* in each test I want to vary the input length and insert errors and test if the sequence is still extracted  */
/* 1) R+ */


int run_single_test(struct shared_sim_data* sd, char** arch,int n, char* seq,char* label, float mean_q, float stdev_q)
{

        struct poahmm* poahmm = NULL;
        int i,j;
        int len;

        int active_read_structure;
        uint8_t* i_seq = NULL;
        uint8_t* i_qual =  NULL;
        uint32_t* path = NULL;

        float e;


        MMALLOC(path, sizeof(int)*1024);

        poahmm = sd->poahmm;

        LOG_MSG("Current poa size: nodes: %d seq_len:%d", poahmm->alloced_num_nodes , poahmm->alloc_seq_len);

        RUN(read_arch_into_lib(sd->al, arch, n));
        for(i = 0; i < sd->al->num_arch;i++){

                LOG_MSG("Arch %d:  %s",i,sd->al->spec_line[i]);
        }


        active_read_structure = sd->al->num_arch -1;
        LOG_MSG("Working with: %s",sd->al->spec_line[active_read_structure]);

        free_poahmm(poahmm);
        poahmm = NULL;
        sd->gp->min_seq_len = 10;
        sd->gp->max_seq_len = 16;
        RUN(poahmm_from_read_structure(&poahmm, sd->gp, sd->al->read_structure[active_read_structure ], sd->a));

        len = strlen(seq);
        MMALLOC(i_seq, sizeof(uint8_t) * (len+1));
        MMALLOC(i_qual, sizeof(uint8_t) * (len+1));
        for(i = 0; i < len; i++){
                i_seq[i] = tlalphabet_get_code(sd->a, seq[i]);

                i_qual[i] = tl_random_gaussian(sd->main_rng, mean_q, stdev_q);
                LOG_MSG("%c %d %d", seq[i],i_seq[i], i_qual[i]);
        }
        i_seq[len] = 0;
        RUN(viterbi_poahmm_banded(poahmm,i_seq, i_qual, len,path,2));


        e = poahmm->f_score - poahmm->random_scores[i] ;//  poahmm->random_scores[i];
        e = exp2f(e) / (1.0 + exp2f(e));
        LOG_MSG("Len: %d:", i);
        fprintf(stdout,"e:%f F:%5.3f\tR:%5.3f\n",e, poahmm->f_score, poahmm->random_scores[i]);
        RUN(print_path(poahmm,path ,seq,label));

        sd->poahmm = poahmm;
        MFREE(path);
        MFREE(i_seq);
        MFREE(i_qual);

        return OK;
ERROR:
        return FAIL;
}

int test_simple_N_plus_arch(struct shared_sim_data* sd)
{
        struct poahmm* poahmm = NULL;
        char* in[] = {
                //"I:NNNNNN"
                "R2:E:N{4}",
                "A:AA,TT",
                "E:N{4}"
        };

        struct stats{
                double s_0_correct[5];
                double s_1_correct[5];
                double s_2_correct[5];

                double s_0_false[5];
                double s_1_false[5];
                double s_2_false[5];
                double pass;
        } stats;

        int size;

        uint8_t* i_seq = NULL;
        uint32_t* path = NULL;

        int i_len;
        int i;
        int j;

        int base_q;
        float e;
        float accuracy;
        float error_rate;

        int num_error;

        int active_read_structure;

        int num_tests = 1;
        size = sizeof(in) / sizeof(char*);

        MMALLOC(path, sizeof(int)*1024);

        poahmm = sd->poahmm;

        LOG_MSG("Current poa size: nodes: %d seq_len:%d", poahmm->alloced_num_nodes , poahmm->alloc_seq_len);

        RUN(read_arch_into_lib(sd->al, in, size));
        for(i = 0; i < sd->al->num_arch;i++){

                LOG_MSG("Arch %d:  %s",i,sd->al->spec_line[i]);
        }


        active_read_structure = sd->al->num_arch -1;
        LOG_MSG("Working with: %s",sd->al->spec_line[active_read_structure]);

        sd->gp->min_seq_len = 10;
        sd->gp->max_seq_len = 16;
        RUN(poahmm_from_read_structure(&poahmm, sd->gp, sd->al->read_structure[active_read_structure ], sd->a));
        //print_poahmm_param(poahmm);
        poahmm_to_dot(poahmm, "test.dot");


        //LOG_MSG("MM N : %f %f", scaledprob2prob(poahmm->MM), (1.0/16.0));
        //LOG_MSG("YY B : %f %f", scaledprob2prob(poahmm->YY_boundary ), scaledprob2prob(poahmm->background[1]));

        //RUN(set_random_scores(poahmm, 50));
        LOG_MSG("Model_len: %d -%d ", poahmm->min_model_len,poahmm->max_model_len);
        for(i = poahmm->min_model_len ;i <=  poahmm->max_model_len;i++){
                for(base_q = 40; base_q >= 5;base_q -= 5){
                //RUN(sim_seq_from_read_struct(&sd->seq, sd->al->read_structure[active_read_structure ] , i, sd->main_rng));
                //RUN(generate_random_seq(sd->seq->seq, i, sd->main_rng));
                //RUN(generate_random_seq(&seq, &i, sd->main_rng ));
                //nif(i > 10){
                //RUN(insert_seq(seq, i, "AACCGGTT", 8, sd->main_rng));
                //}
                //RUN(seq_to_internal(sd->seq->seq, i, &i_seq, &i_len));
                //poahmm->r_len = 1.0f;//   MACRO_MAX(1, i - poahmm->max_rank);
                        LOG_MSG("Simulating reads of length %d, base quality %d",i, base_q);
                        for(j = 0; j < 5;j ++){

                                stats.s_0_correct[j] = 0.0;
                                stats.s_1_correct[j] = 0.0;
                                stats.s_2_correct[j] = 0.0;

                                stats.s_0_false[j] = 0.0;
                                stats.s_1_false[j] = 0.0;
                                stats.s_2_false[j] = 0.0;
                                stats.pass = 0.0;
                        }

                        //set_terminal_gap_prob(poahmm, i);
                        //LOG_MSG("%f %f", scaledprob2prob(poahmm->YY_boundary), scaledprob2prob(poahmm->YY_boundary_exit));
                        error_rate = 0.2f;
                        for(j = 0; j < num_tests; j++){
                                //error_rate += 0.001f;
                                RUN(sim_seq_from_read_struct(&sd->seq, sd->al->read_structure[active_read_structure ] , i,base_q,sd->main_rng));

                                mutate_seq(sd->seq->seq,sd->seq->buffer,  sd->seq->len, error_rate , sd->main_rng, &num_error);
                                RUN(seq_to_internal(sd->seq->buffer, i, &i_seq, &i_len));


                                RUN(viterbi_poahmm_banded(poahmm,i_seq, sd->seq->qual, i_len,path,2));
                                //random_poahmm(poahmm, i_seq, i_len);
                                RUN(score_labelling(poahmm, path,sd->seq->label, &accuracy));


                                e = poahmm->f_score - poahmm->random_scores[i] ;//  poahmm->random_scores[i];
                                e = exp2f(e) / (1.0 + exp2f(e));
                                //LOG_MSG("Len: %d:", i);
                                fprintf(stdout,"%d a:%f e:%f F:%5.3f\tR:%5.3f nume:%d\n",i,accuracy,e, poahmm->f_score, poahmm->random_scores[i],  num_error);
                                RUN(print_path(poahmm,path ,sd->seq->buffer, sd->seq->label ));
                                //e = scaledprob2prob(e);
                                if(e - 0.0001 >= 0.5f){
                                        stats.pass++;
                                }

                                if(num_error > 4){
                                        num_error = 4 ;
                                }
                                if(accuracy == 1.0){
                                        stats.s_0_correct[num_error]++;
                                        stats.s_1_correct[num_error] += e;
                                        stats.s_2_correct[num_error] += e * e;
                                }else{
                                        stats.s_0_false[num_error]++;
                                        stats.s_1_false[num_error] += e;
                                        stats.s_2_false[num_error] += e * e;
                                }

                                /*if(accuracy != 1.0 && num_error <= 5){
                                  LOG_MSG("Wrong assignment at error rate: %f (%d)", error_rate, num_error);
//fprintf(stdout,"score: %f\n",e);
//forward_poahmm(poahmm, i_seq, i_len);
//backward_poahmm(poahmm, i_seq, i_len);
//random_poahmm(poahmm, i_seq, i_len);
//e = poahmm->f_score - poahmm->random_scores[i];
//e = exp2f(e) / (1.0 + exp2f(e));
//forward_poahmm(poahmm, i_seq, i_len);
//backward_poahmm(poahmm, i_seq, i_len);
//LOG_MSG("%f", poahmm->random_scores[i]);
fprintf(stdout,"%d a:%f e:%f F:%5.3f\tR:%5.3f nume:%d\t",i,accuracy,e, poahmm->f_score, poahmm->random_scores[i], num_error);
fprintf(stdout,"%s\n",sd->seq->buffer);

if((accuracy != 1.0 && e >= 0.9)|| num_error < 3 ){
RUN(print_path(poahmm,path ,sd->seq->buffer, sd->seq->label ));
}
//fprintf(stdout,"%f\n",1.0 - scaledprob2prob(poahmm->f_score - poahmm->r_score));

//break;
}*/
                                //RUN(viterbi_poahmm(poahmm,i_seq, i_len,path));


                                //RUN(print_path(poahmm,path ,sd->seq->buffer, sd->seq->label ));
                        //exit(0);
                }

                LOG_MSG("------------------------------");
                double stdev;
                double mean;
                for(j = 0;j < 5;j++){
                        fprintf(stdout,"ERRORS: %d\t",j);

                        if(stats.s_0_correct[j] == 0.0){
                                mean = 0.0;
                                stdev = 0.0;

                        }else{
                                mean = stats.s_1_correct[j] / stats.s_0_correct[j];
                                stdev = sqrt ( (stats.s_0_correct[j] * stats.s_2_correct[j] -  pow(stats.s_1_correct[j], 2.0)) /  (stats.s_0_correct[j] * ( stats.s_0_correct[j] - 1.0)));
                        }
                        fprintf(stdout,"%f (+-%f) N=%d\t", mean,stdev, (int)stats.s_0_correct[j]);

                        if(stats.s_0_false[j] == 0.0){
                                mean = 0.0;
                                stdev = 0.0;
                        }else{
                                mean = stats.s_1_false[j] / stats.s_0_false[j];
                                stdev = sqrt ( (stats.s_0_false[j] * stats.s_2_false[j] -  pow(stats.s_1_false[j], 2.0)) /  (stats.s_0_false[j] * ( stats.s_0_false[j] - 1.0)));
                        }
                        fprintf(stdout,"%f (+-%f) N=%d\t", mean,stdev,(int) stats.s_0_false[j]);
                        //fprintf(stdout,"%f (+-%f)\t", mean,stdev);
                        fprintf(stdout,"\n");

                }
                fprintf(stdout,"Pass: %d (%3.1f)\n", (int) stats.pass, stats.pass / (double) num_tests * 100.0);

                LOG_MSG("------------------------------");

                MFREE(i_seq);

                i_seq = NULL;
                }
        }
        MFREE(path);


/* resport stats  */
        return OK;
ERROR:
        return FAIL;
}

int test_banded(struct shared_sim_data* sd)
{
        struct poahmm* poahmm = NULL;
        char* in[] = {
                "I:ACGTACGTACGT",
                "R2:E:N{15}",
                "S:AAAAAAAA,CCCCCCCC,GGGGGGGG,TTTTTTTT",
                "A:AAACCCGGGTTT,TTTGGGCCCAAA"
        };
        int size;


        uint8_t* i_seq = NULL;
        uint32_t* path = NULL;
        int runs = 10;
        int i_len;
        int i;
        int j;
        int c;
        float error_rate;

        int num_error;

        int active_read_structure;
        size = sizeof(in) / sizeof(char*);


        poahmm = sd->poahmm;

        LOG_MSG("Current poa size: nodes: %d seq_len:%d", poahmm->alloced_num_nodes , poahmm->alloc_seq_len);

        RUN(read_arch_into_lib(sd->al, in, size));

        active_read_structure = sd->al->num_arch-1;

        RUN(poahmm_from_read_structure(&poahmm, sd->gp, sd->al->read_structure[active_read_structure  ], sd->a));

        LOG_MSG("Model Len: %d %d", poahmm->min_model_len,poahmm->max_model_len);

        LOG_MSG("Seq Len: %d %d", poahmm->min_seq_len, poahmm->max_seq_len);
        //print_poahmm_param(poahmm);
        poahmm_to_dot(poahmm, "test2.dot");
        //exit(0);

        MMALLOC(path, sizeof(int)*poahmm->max_model_len + poahmm->alloc_seq_len);


        LOG_MSG("Starting run");
        DECLARE_TIMER(timer);
        //RUN(set_random_scores(poahmm, 50));
        for(i = poahmm->min_model_len ;i <= poahmm->max_model_len;i++){
                //RUN(sim_seq_from_read_struct(&sd->seq, sd->al->read_structure[1] , i, sd->main_rng));
                //RUN(generate_random_seq(sd->seq->seq, i, sd->main_rng));
                //RUN(generate_random_seq(&seq, &i, sd->main_rng ));
                //nif(i > 10){
                //RUN(insert_seq(seq, i, "AACCGGTT", 8, sd->main_rng));
                //}
                //RUN(seq_to_internal(sd->seq->seq, i, &i_seq, &i_len));
                //poahmm->r_len = 1.0f;//   MACRO_MAX(1, i - poahmm->max_rank);

                //set_terminal_gap_prob(poahmm, i);
                if(i >= poahmm->alloc_seq_len){
                        LOG_MSG("Resize");
                        resize_poahmm(poahmm, poahmm->num_nodes, i);
                        LOG_MSG("%d", poahmm->alloc_seq_len);
                }

                error_rate = 0.0;
                RUN(sim_seq_from_read_struct(&sd->seq, sd->al->read_structure[ active_read_structure] , i,49, sd->main_rng));
                mutate_seq(sd->seq->seq,sd->seq->buffer,  sd->seq->len, error_rate , sd->main_rng, &num_error);
                RUN(seq_to_internal(sd->seq->buffer, i, &i_seq, &i_len));


                START_TIMER(timer);
                for(j = 0; j < runs;j++){
                        RUN(viterbi_poahmm(poahmm,i_seq,i_len,path));

                }
                STOP_TIMER(timer);
                LOG_MSG("LEN:%d\t%f\tNormal Viterbi\ttook %f",i, poahmm->f_score,GET_TIMING(timer));
                poahmm->f_score = -2;
                for(c = 10; c > 0;c--){
                        START_TIMER(timer);
                        for(j = 0; j < runs;j++){
                                RUN(viterbi_poahmm_banded(poahmm,i_seq,sd->seq->qual, i_len,path,c));
                        }
                        STOP_TIMER(timer);
                        LOG_MSG("LEN:%d\t%f\tbanded %d\ttook: %f",i, poahmm->f_score,c,GET_TIMING(timer));


                }

                MFREE(i_seq);

                i_seq = NULL;
        }
        MFREE(path);
        return OK;
ERROR:
        return FAIL;
}

int test_indel(struct shared_sim_data* sd)
{
        struct poahmm* poahmm = NULL;
        char* in[] = {
                "I:CCCCCCCCCCCCCC",
                "I:AAAAAAAA"
        };
        int size;


        uint8_t* i_seq = NULL;
        uint32_t* path = NULL;

        int i_len;
        //int i;
        //int j;

        float e;




        int active_read_structure;

        size = sizeof(in) / sizeof(char*);


        poahmm = sd->poahmm;
        RUN(read_arch_into_lib(sd->al, in, size));

        active_read_structure = sd->al->num_arch-1;
        RUN(poahmm_from_read_structure(&poahmm, sd->gp, sd->al->read_structure[active_read_structure ], sd->a));
        //print_poahmm_param(poahmm);
        poahmm_to_dot(poahmm, "test3.dot");


        MMALLOC(path, sizeof(int)*poahmm->max_model_len + poahmm->alloc_seq_len);
        LOG_MSG("Starting run");


        //RUN(set_random_scores(poahmm, 50));
        snprintf(sd->seq->seq, sd->seq->alloc_len,   "%s%s","CCCCCCCCCCCCCC", "AAAAAAAA");
        snprintf(sd->seq->label, sd->seq->alloc_len, "%s%s","IIIIIIIIIIIIII", "IIIIIIII");
        sd->seq->len = strnlen(sd->seq->seq, sd->seq->alloc_len);

        LOG_MSG("LEN: %d",sd->seq->len);

        if(sd->seq->len >= poahmm->alloc_seq_len){
                LOG_MSG("Resize");
                resize_poahmm(poahmm, poahmm->num_nodes, sd->seq->len);
                LOG_MSG("%d", poahmm->alloc_seq_len);
        }

        RUN(seq_to_internal(sd->seq->seq, sd->seq->len, &i_seq, &i_len));
        RUN(viterbi_poahmm(poahmm,i_seq,i_len,path));
        e = poahmm->f_score - poahmm->random_scores[sd->seq->len];
        e = exp2f(e) / (1.0 + exp2f(e));
        LOG_MSG("Score: %f (%f,%f)",e, poahmm->f_score, poahmm->random_scores[sd->seq->len]);
        print_path(poahmm, path,sd->seq->seq, sd->seq->label);

        MFREE(i_seq);
        i_seq = NULL;

        snprintf(sd->seq->seq, sd->seq->alloc_len,   "%s%c%s","CCCCCCCCCCCCCC",'T', "AAAAAAAA");
        snprintf(sd->seq->label, sd->seq->alloc_len, "%s%c%s","IIIIIIIIIIIIII",'X', "IIIIIIII");
        sd->seq->len = strnlen(sd->seq->seq, sd->seq->alloc_len);
        LOG_MSG("LEN: %d",sd->seq->len);


        RUN(seq_to_internal(sd->seq->seq, sd->seq->len, &i_seq, &i_len));
        RUN(viterbi_poahmm(poahmm,i_seq,i_len,path));
        e = poahmm->f_score - poahmm->random_scores[sd->seq->len];
        e = exp2f(e) / (1.0 + exp2f(e));
        LOG_MSG("Score: %f (%f,%f)",e, poahmm->f_score, poahmm->random_scores[sd->seq->len]);
        print_path(poahmm, path,sd->seq->seq, sd->seq->label);
        MFREE(i_seq);
        i_seq = NULL;


        snprintf(sd->seq->seq, sd->seq->alloc_len,   "%s%s%s","CCCCCCCCCCCCCC","TTTT", "AAAAAAAA");
        snprintf(sd->seq->label, sd->seq->alloc_len, "%s%s%s","IIIIIIIIIIIIII","XXXX", "IIIIIIII");
        sd->seq->len = strnlen(sd->seq->seq, sd->seq->alloc_len);
        LOG_MSG("LEN: %d",sd->seq->len);

        RUN(seq_to_internal(sd->seq->seq, sd->seq->len, &i_seq, &i_len));
        RUN(viterbi_poahmm(poahmm,i_seq,i_len,path));
        e = poahmm->f_score - poahmm->random_scores[sd->seq->len];
        e = exp2f(e) / (1.0 + exp2f(e));
        LOG_MSG("Score: %f (%f,%f)",e, poahmm->f_score, poahmm->random_scores[sd->seq->len]);
        print_path(poahmm, path,sd->seq->seq, sd->seq->label);
        MFREE(i_seq);
        i_seq = NULL;

        MFREE(path);
        return OK;
ERROR:
        return FAIL;
}


int single_seq_test(void)
{
        struct rng_state* rng = NULL;
        struct poahmm* poahmm = NULL;
        char* seq = NULL;
        uint8_t* i_seq = NULL;
        uint32_t* path = NULL;
        int len = 32;
        int i_len;
        //int nuc_count[4];

        //int i;

        int malloced_path = 100;

        MMALLOC(path , sizeof( uint32_t) *malloced_path);



        RUNP(poahmm = init_poahmm(NULL));

        RUNP(rng = init_rng(0));

        RUN(generate_random_seq(&seq, &len, rng));
        RUN(seq_to_internal(seq, len, &i_seq, &i_len));


        RUN(init_nodes_from_single_sequence(poahmm, i_seq+8,i_len-16));

        forward_poahmm(poahmm, i_seq, i_len);
        backward_poahmm(poahmm, i_seq, i_len);
        fprintf(stdout,"F:%f\nB:%f\n", poahmm->f_score,poahmm->b_score);



        RUN(viterbi_poahmm(poahmm,i_seq,i_len,path));
        //update_poahmm_arch(poahmm, path, 10);
                //RUN(print_dyn_matrix(poahmm,poahmm_data->len[i] ),"print matrix failed.");


        //RUN(print_path(poahmm,path));
        MFREE(path);
        if(i_seq){
                MFREE(i_seq);
        }
        if(seq){
                MFREE(seq);
        }
        free_rng(rng);
        free_poahmm(poahmm);
        return OK;
ERROR:
        return FAIL;

}


int arch_test(void)
{
        struct poahmm* poahmm = NULL;
        struct arch_library* al = NULL;

        struct alphabet* a = NULL;
        struct rng_state* rng = NULL;

        char* seq = NULL;
        uint8_t* i_seq = NULL;
        uint32_t* path = NULL;
        int len = 32;
        int i_len;

        char* in[] = {
                 "I:A{0,3}",
                 "S:GTA,AAC",
                 "MYREAD:E:NNN",
                 "I:CCTTAA",
                 "S:ACAGTG,ACTTGA,TTAGGC",
                "R2:E:N{5,20}"
        };
        //int nuc_count[4];

        int size;
        //int i;
        int malloced_path = 100;

        MMALLOC(path , sizeof(uint32_t ) *malloced_path);


        RUN(alloc_arch_lib(&al));

        size= sizeof(in) / sizeof(char*);

        RUN(read_arch_into_lib(al, in, size));

        LOG_MSG("Read in %d architectures.",al->num_arch);
        RUNP(rng= init_rng(0));
        RUNP(poahmm = init_poahmm(NULL));
        RUN(create_alphabet(&a, rng, TLALPHABET_DEFAULT_DNA ));
        RUN(init_nodes_from_read_structure(poahmm, al->read_structure[1],a,0,10,10));
        poahmm_to_dot(poahmm, "test.dot");
        //print_poahmm_param(poahmm);
        //exit(0);
        RUN(generate_random_seq(&seq, &len, rng));
        //exit(0);
        char match_seq[] = "TTTTTTCTTTCCTTAATTATGGCACGTACGT";
        len = strlen(match_seq);

        RUN(seq_to_internal(match_seq, len, &i_seq, &i_len));

        forward_poahmm(poahmm, i_seq, i_len);
        backward_poahmm(poahmm, i_seq, i_len);
        random_poahmm(poahmm, i_seq, i_len);
        fprintf(stdout,"F:%f\nB:%f\nR:%f\n", poahmm->f_score,poahmm->b_score, poahmm->r_score);
        RUN(viterbi_poahmm(poahmm,i_seq,i_len,path));
        fprintf(stdout,"%s\n",match_seq);
        //RUN(print_path(poahmm,path));
        fprintf(stdout,"F:%f\n", poahmm->f_score);

        if(seq){
                MFREE(seq);
        }
        if(i_seq){
                MFREE(i_seq);

        }

        MFREE(path);
        free_arch_lib(al);
        free_alphabet(a);
        free_rng(rng);
        al = NULL;
        free_poahmm(poahmm);
        return OK;
ERROR:
        return FAIL;

}


int print_path(struct poahmm* poahmm, uint32_t* path,char* seq,char* label)
{
        uint32_t i;
        char alphabet[5] = "ACGTN";
        char etype[8] = "_EASIPLR";

        uint32_t seq_pos;
        uint32_t node_pos;
        fprintf(stdout, "PATH:\n");

        for(i = 1; i < path[0];i++){
                seq_pos = path[i] >> 16u;
                node_pos = path[i] & 0xFFFFu;
                //fprintf(stdout,"Position %d: %d %d\n",i,seq_pos,node_pos);
                if(seq_pos != 0xFFFFu){
                        fprintf(stdout, " %3c", seq[seq_pos]);
                }else{
                        fprintf(stdout, " %3c", '-');
                }
        }
        fprintf(stdout,"\n");
        for(i = 1; i < path[0];i++){
                seq_pos = path[i] >> 16u;
                node_pos = path[i] & 0xFFFFu;

                if(node_pos!= 0xFFFFu){
                        fprintf(stdout, " %3c", alphabet[poahmm->nodes[node_pos]->nuc]);
                }else{
                        fprintf(stdout, " %3c", '-');
                }
        }
        fprintf(stdout,"\n");


                for(i = 1; i < path[0];i++){
                seq_pos = path[i] >> 16u;
                node_pos = path[i] & 0xFFFFu;
                //fprintf(stdout,"Position %d: %d %d\n",i,seq_pos,node_pos);
                if(seq_pos != 0xFFFFu){
                        fprintf(stdout, " %3c", label[seq_pos]);
                }else{
                        fprintf(stdout, " %3c", '-');
                }
        }
        fprintf(stdout,"\n");
        for(i = 1; i < path[0];i++){
                seq_pos = path[i] >> 16u;
                node_pos = path[i] & 0xFFFFu;

                if(node_pos!= 0xFFFFu){
                        fprintf(stdout, " %3c", etype[  poahmm->nodes[node_pos]->type]);
                }else{
                        fprintf(stdout, " %3c", '-');
                }
        }
        fprintf(stdout,"\n");
        for(i = 1; i < path[0];i++){
                seq_pos = path[i] >> 16u;
                node_pos = path[i] & 0xFFFFu;

                if(node_pos!= 0xFFFFu){
                        fprintf(stdout, " %3d",poahmm->nodes[node_pos]->alt);
                }else{
                        fprintf(stdout, " %3c", '-');
                }
        }
        fprintf(stdout,"\n");

        return OK;
}

int score_labelling(struct poahmm* poahmm, uint32_t* path,char* label, float* score)
{

        char etype[6] = "_EASIP";
        float len;

        uint32_t i;

        uint32_t seq_pos;
        uint32_t node_pos;

        len = 0.0;
        *score= 0.0;

        for(i = 1; i < path[0];i++){
                seq_pos = path[i] >> 16u;
                node_pos = path[i] & 0xFFFFu;
                //fprintf(stdout,"Position %d: %d %d\n",i,seq_pos,node_pos);
                if(seq_pos != 0xFFFFu && node_pos!= 0xFFFFu){
                        if (label[seq_pos]  == etype[  poahmm->nodes[node_pos]->type]){
                                        *score = *score + 1.0;


                        }
                        len++;
                }
        }
        *score = *score / len;

        return OK;
//ERROR:
        //return FAIL;
}





int print_poa_graph(struct poahmm* poahmm)
{
        int i,j;
        fprintf(stdout,"POA graph:\n");
        fprintf(stdout,"%3d ",0);
        for(i = 0; i < poahmm->num_nodes;i++){
                fprintf(stdout,"%3d ",i);
        }
        fprintf(stdout,"\n");
        fprintf(stdout,"%3d ",0);
        for(i = 0; i < poahmm->num_nodes;i++){
                fprintf(stdout,"%f ",poahmm->entry_probabilities[i]);
        }
        fprintf(stdout,"\n");



        fprintf(stdout,"%3d ",0);
        for(i = 0; i < poahmm->num_nodes;i++){
                fprintf(stdout,"%3d ",i);
        }
        fprintf(stdout,"\n");
        for(i = 0; i <  poahmm->alloced_num_nodes;i++){
                fprintf(stdout,"%3d ",i);
                for(j = 0; j <  poahmm->alloced_num_nodes;j++){
                        //if(poahmm->poa_graph[i][j] != prob2scaledprob(0.0)){
                        fprintf(stdout,"%f ",poahmm->e_poa_graph[i][j] );
                        //}else{
                        //	fprintf(stdout,"%f ",0);
                        //}
                }
                fprintf(stdout,"	%c\n","ACGTN"[poahmm->nodes[i]->nuc]);
        }
        fprintf(stdout,"\n");


        fprintf(stdout,"To: graph:\n");
        for(i = 0; i < poahmm->num_nodes;i++){
                fprintf(stdout,"to:%d :",i);
                for(j = 1; j <  poahmm->to_tindex[i][0];j++){
                        fprintf(stdout,"%d ",poahmm->to_tindex[i][j]);
                }
                fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");

        fprintf(stdout,"From: graph:\n");
        for(i = 0; i < poahmm->num_nodes;i++){
                fprintf(stdout,"to:%d :",i);
                for(j = 1; j <  poahmm->from_tindex[i][0];j++){
                        fprintf(stdout,"%d ",poahmm->from_tindex[i][j]);
                }
                fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");


        return OK;
}






int print_poahmm_param(struct poahmm* poahmm)
{


        int i,j;
        float sum;


        fprintf(stdout,"emission_M:\n");
        for(i = 0; i < 5;i++){
                sum = prob2scaledprob(0.0);
                for(j = 0; j < 5;j++){
                        fprintf(stdout,"%0.4f ", scaledprob2prob( poahmm->emission_M[i][j]));
                        sum = logsum(sum, poahmm->emission_M[i][j]);
                }
                fprintf(stdout,"	%f\n",scaledprob2prob(sum));
        }
        fprintf(stdout,"\n");

        fprintf(stdout,"emission_X:\n");
        sum = prob2scaledprob(0.0);
        for(i = 0; i < 4;i++){

                fprintf(stdout,"%0.4f ", scaledprob2prob( poahmm->emission_X[i]));
                sum = logsum(sum, poahmm->emission_X[i]);

        }
        fprintf(stdout,"	%f\n",scaledprob2prob(sum));
        sum = prob2scaledprob(0.0);
        fprintf(stdout,"emission_Y:\n");
        for(i = 0; i < 4;i++){

                fprintf(stdout,"%0.4f ", scaledprob2prob( poahmm->emission_Y[i]));
                sum = logsum(sum, poahmm->emission_Y[i]);

        }
        fprintf(stdout,"	%f\n",scaledprob2prob(sum));


        fprintf(stdout,"\n\n");


        fprintf(stdout,"e_emission_M:\n");
        for(i = 0; i < 4;i++){
                sum = prob2scaledprob(0.0);
                for(j = 0; j < 4;j++){
                        fprintf(stdout,"%0.4f ", scaledprob2prob( poahmm->e_emission_M[i][j]));
                        sum = logsum(sum, poahmm->e_emission_M[i][j]);
                }
                fprintf(stdout,"	%f\n",scaledprob2prob(sum));
        }
        fprintf(stdout,"\n");

        fprintf(stdout,"e_emission_X:\n");
        sum = prob2scaledprob(0.0);
        for(i = 0; i < 4;i++){

                fprintf(stdout,"%0.4f ", scaledprob2prob( poahmm->e_emission_X[i]));
                sum = logsum(sum, poahmm->e_emission_X[i]);

        }
        fprintf(stdout,"	%f\n",scaledprob2prob(sum));
        sum = prob2scaledprob(0.0);
        fprintf(stdout,"e_emission_Y:\n");
        for(i = 0; i < 4;i++){

                fprintf(stdout,"%0.4f ", scaledprob2prob( poahmm->e_emission_Y[i]));
                sum = logsum(sum, poahmm->e_emission_Y[i]);

        }
        fprintf(stdout,"	%f\n",scaledprob2prob(sum));

        sum = prob2scaledprob(0.0);

        fprintf(stdout,"e_random:\n");
        for(i = 0; i < 4;i++){

                fprintf(stdout,"%0.4f ", scaledprob2prob( poahmm->background[i]));
                sum = logsum(sum, poahmm->background[i]);

        }
        fprintf(stdout,"	%f\n",scaledprob2prob(sum));




        fprintf(stdout,"M: %0.4f %0.4f %0.4f\teM:%0.4f %0.4f %0.4f\n",scaledprob2prob(poahmm->MM), scaledprob2prob(poahmm->MX),scaledprob2prob(poahmm->MY),scaledprob2prob(poahmm->e_MM), scaledprob2prob(poahmm->e_MX),scaledprob2prob(poahmm->e_MY));

        fprintf(stdout,"X: %0.4f %0.4f %0.4f\teX:%0.4f %0.4f %0.4f\n",scaledprob2prob(poahmm->XM), scaledprob2prob(poahmm->XX),scaledprob2prob(poahmm->XY),scaledprob2prob(poahmm->e_XM), scaledprob2prob(poahmm->e_XX),scaledprob2prob(poahmm->e_XY));

        fprintf(stdout,"Y: %0.4f %0.4f %0.4f\teY:%0.4f %0.4f %0.4f\n",scaledprob2prob(poahmm->YM), scaledprob2prob(poahmm->YX),scaledprob2prob(poahmm->YY),scaledprob2prob(poahmm->e_YM), scaledprob2prob(poahmm->e_YX),scaledprob2prob(poahmm->e_YY));
        fprintf(stdout,"entry\nexit\n");
        for(i = 0; i < poahmm->num_nodes;i++){
                fprintf(stdout,"%d %f %f\n",i, scaledprob2prob(poahmm->entry_probabilities[i]),scaledprob2prob(poahmm->exit_probabilities[i]));
        }
        return OK;
}









int init_sim_data(struct shared_sim_data** sim_data, int seed)
{
        struct shared_sim_data* sd = NULL;

        int i;
        MMALLOC(sd, sizeof(struct shared_sim_data));
        sd->a = NULL;
        sd->main_rng = NULL;
        sd->gp = NULL;
        sd->al = NULL;
        sd->seq = NULL;

        RUNP(sd->main_rng = init_rng(seed));

        RUN(create_alphabet(&sd->a, sd->main_rng, TLALPHABET_DEFAULT_DNA));

        MMALLOC(sd->gp, sizeof(struct global_poahmm_param));
        sd->gp->base_error = 0.05f;
        sd->gp->indel_freq = 0.01f;
        for(i = 0; i < 4;i++){
                sd->gp->back[i] = prob2scaledprob(0.25f);
        }
        sd->gp->back[4] = prob2scaledprob(1.0);
        sd->gp->average_seq_length = 128;
        sd->gp->min_seq_len = 1;
        sd->gp->max_seq_len = 256;

        RUNP(sd->poahmm = init_poahmm( sd->gp));

        RUN(sim_seq_from_read_struct(&sd->seq, NULL,10,49,sd->main_rng));
        RUN(alloc_arch_lib(&sd->al));

        *sim_data = sd;
        return OK;
ERROR:
        return FAIL;
}

void free_sim_data(struct shared_sim_data* sd)
{
        if(sd){
                free_sim_seq(sd->seq);
                free_poahmm(sd->poahmm);
                free_rng(sd->main_rng);
                free_alphabet(sd->a);
                free_arch_lib(sd->al);
                MFREE(sd->gp);
                MFREE(sd);
        }

}

int sim_seq_from_read_struct(struct sim_seq** simseq, struct read_structure* rs,int sim_len, int base_q, struct rng_state* rng)
{
        struct sim_seq* s = NULL;
        struct segment_specs* spec = NULL;
        char alphabet[4] = "ACGT";
        char etype[6] = "_EASIP";
        //int total_len;
        int i;
        int j;
        int c;
        int g;
        int b;
        s= *simseq;
        if(s == NULL){
                MMALLOC(s, sizeof(struct sim_seq));

                s->alloc_len = 1024;
                s->len = 0;
                s->label = NULL;
                s->seq = NULL;
                s->buffer = NULL;
                s->qual = NULL;
                MMALLOC(s->seq, sizeof(char)* s->alloc_len);
                MMALLOC(s->label, sizeof(char)* s->alloc_len);
                MMALLOC(s->buffer, sizeof(char)* s->alloc_len);
                MMALLOC(s->qual, sizeof(uint8_t)* s->alloc_len);
        }
        if(sim_len >= s->alloc_len){
                s->alloc_len = s->alloc_len + s->alloc_len /2;
                MREALLOC(s->seq, sizeof(char)* s->alloc_len);
                MREALLOC(s->label, sizeof(char)* s->alloc_len);
                MREALLOC(s->buffer, sizeof(char)* s->alloc_len);
                MREALLOC(s->qual, sizeof(uint8_t)* s->alloc_len);
        }

        /* step one: add random residues  */
        RUN(generate_random_seq(&s->seq,&sim_len,rng));
        for(i = 0; i < sim_len;i++){
                s->label[i] = etype[0];
                s->qual[i] = base_q;
        }

        s->label[sim_len] = 0;
        s->len = sim_len;

        if(rs){
                int plus_min_len, plus_max_len;
                RUN(set_len_of_unknown_poa(rs,&plus_min_len,&plus_max_len, sim_len,sim_len));
                //total_len = 0;
                c = 0;
                for(i = 0; i < rs->num_segments;i++){
                        spec = rs->seg_spec[i];
                        //print_segment_spec(spec);
                        if(spec->max_len == INT32_MAX){
                                g = tl_random_int(rng, plus_max_len - plus_min_len) + plus_min_len;
                        }else{
                                g = tl_random_int(rng, spec->max_len- spec->min_len) + spec->min_len;
                        }
                        //LOG_MSG("G: %d", g);
                        b = 0;
                        if(spec->num_seq >1){
                                b = tl_random_int(rng, spec->num_seq-1) + 1;
                        }
                        for(j = 0; j < g;j++){
                                //LOG_MSG("j:%d",j);
                                if(spec->max_len == INT32_MAX){
                                        s->buffer[c] = alphabet[tl_random_int(rng, 4)];
                                }else{

                                        if( spec->seq[b][j] == 'N'){
                                                s->buffer[c] = alphabet[tl_random_int(rng, 4)];
                                        }else{
                                                s->buffer[c] = spec->seq[b][j];
                                        }
                                }
                                s->label[c] = etype[spec->extract];
                                c++;
                        }
                }
                s->buffer[c] = 0;
                s->label[c] = 0;
                if(c > s->len){
                        c = s->len;
                        s->label[c] = 0;
                }
                RUN(insert_seq(s->seq, s->len, s->buffer,c,rng,&g));
                for(i = 0; i < c;i++){
                        s->buffer[i] = s->label[i];
                }
                for(i = 0; i < sim_len;i++){
                        s->label[i] = etype[0];
                }

                //memcpy(s->label+g, s->buffer, c);
                for(i = 0; i < c;i++){
                        s->label[i+g] = s->buffer[i];
                }
        }
        *simseq = s;



        return OK;
ERROR:
        return FAIL;
}


void free_sim_seq(struct sim_seq* s)
{

        if(s){
                MFREE(s->buffer);
                MFREE(s->seq);
                MFREE(s->label);
                MFREE(s->qual);
                MFREE(s);
        }
}
