
#include <stdint.h>
#include <string.h>
#include "tldevel.h"

#include "arch_lib.h"

#include "tlrng.h"
#include "tllogsum.h"
#include "tlalphabet.h"

#include "poahmm.h"



#include "sim_seq_lib.h"

struct sim_seq{
        char* buffer;
        char* seq;
        char* label;
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

static int sim_seq_from_read_struct(struct sim_seq** simseq, struct read_structure* rs,int sim_len,struct rng_state* rng);
static void free_sim_seq(struct sim_seq* s);

int score_labelling(struct poahmm* poahmm, uint32_t* path,char* label, float* score);
int test_simple_N_plus_arch(struct shared_sim_data* sd);
int test_banded(struct shared_sim_data* sd);
/* for use in tagdust */
int poahmm_from_read_structure(struct poahmm** poahmm,struct global_poahmm_param* p, struct read_structure* rs,struct alphabet* a);

int single_seq_test(void);
int arch_test(void);
int init_nodes_from_read_structure(struct poahmm* poahmm, struct read_structure* rs, struct alphabet* a, int random);

int print_poa_graph(struct poahmm* poahmm);
int print_path(struct poahmm* poahmm, uint32_t* path,char* seq,char* label);

int print_poahmm_param(struct poahmm* poahmm);
int poahmm_to_dot(struct poahmm* poahmm,char* filename);

int poahmm_from_read_structure(struct poahmm** poahmm,struct global_poahmm_param* p, struct read_structure* rs,struct alphabet* a)
{
        struct poahmm* ph = NULL;
        uint8_t* nnn = NULL;
        uint32_t* path = NULL;
        int random;
        int i;
        int j;
        int c;
        struct rng_state* rng = NULL;
        ph = *poahmm;
        if(!ph){
                RUNP(ph = init_poahmm(p));
        }

        random = 1;
        RUN(init_nodes_from_read_structure(ph, rs,a,random));
        MMALLOC(path, sizeof(int)* p->max_seq_len *2);
        MMALLOC(nnn, sizeof(uint8_t) * (p->max_seq_len+1));
        rng = init_rng(0);

        //ph->max_rank
        for(i = ph->max_rank; i < p->max_seq_len;i++){
                //LOG_MSG("%d", i);
                //for(c = 0; c < 10;c++){
                for(j = 0; j < p->max_seq_len+1;j++){
                        nnn[j] = tl_random_int(rng, 4);
                }
                RUN(viterbi_poahmm(ph, nnn, i, path));
                ph->random_scores[i] = ph->f_score;
                //fprintf(stdout,"%f ", ph->f_score);
                //}
                //fprintf(stdout,"\n");
                //exit(0);
        }

        //exit(0);
        free_rng(rng);
        MFREE(path);
        MFREE(nnn);


        RUN(init_nodes_from_read_structure(ph, rs,a,0));
        *poahmm = ph;
        return OK;
ERROR:
        return FAIL;
}


int main(int argc, char *argv[])
{
        struct shared_sim_data* sim_data=NULL;
        init_logsum();
        RUN(init_sim_data(&sim_data, 0));
        RUN(test_banded(sim_data));
        RUN(test_simple_N_plus_arch(sim_data));
        //RUN(single_seq_test());
        //RUN(arch_test());

        free_sim_data(sim_data);
        return EXIT_SUCCESS;
ERROR:
        return EXIT_FAILURE;
}


/* run various tests  */
/* in each test I want to vary the input length and insert errors and test if the sequence is still extracted  */
/* 1) R+ */


int test_simple_N_plus_arch(struct shared_sim_data* sd)
{
        struct poahmm* poahmm = NULL;
        char* in[] = {
                "I:ACGT",
                "R2:E:N{10}",
                "A:AAA,TTT"
        };
        int size;

        uint8_t* i_seq = NULL;
        uint32_t* path = NULL;

        int i_len;
        int i;
        int j;
        float e;
        float accuracy;
        float error_rate;

        int num_error;

        int active_read_structure;
        size = sizeof(in) / sizeof(char*);

        MMALLOC(path, sizeof(int)*1024);
        poahmm = sd->poahmm;
        RUN(read_arch_into_lib(sd->al, in, size));

        active_read_structure = sd->al->num_arch-1;

        RUN(poahmm_from_read_structure(&poahmm, sd->gp, sd->al->read_structure[active_read_structure ], sd->a));
        //print_poahmm_param(poahmm);
        poahmm_to_dot(poahmm, "test.dot");

        LOG_MSG("%d", poahmm->max_rank);
        //RUN(set_random_scores(poahmm, 50));
        for(i = poahmm->max_rank+1 ;i < 30;i++){
                //RUN(sim_seq_from_read_struct(&sd->seq, sd->al->read_structure[active_read_structure ] , i, sd->main_rng));
                //RUN(generate_random_seq(sd->seq->seq, i, sd->main_rng));
                //RUN(generate_random_seq(&seq, &i, sd->main_rng ));
                //nif(i > 10){
                //RUN(insert_seq(seq, i, "AACCGGTT", 8, sd->main_rng));
                //}
                //RUN(seq_to_internal(sd->seq->seq, i, &i_seq, &i_len));
                //poahmm->r_len = 1.0f;//   MACRO_MAX(1, i - poahmm->max_rank);
                if(i < poahmm->max_rank){
                        LOG_MSG("Sequence too short to match");
                }else if(poahmm->max_rank < i){
                        j = i - poahmm->max_rank;
                        poahmm->YY_boundary = (float)j / (float)(j + 2);
                        poahmm->YY_boundary_exit = 1.0 - poahmm->YY_boundary;
                        poahmm->YY_boundary = prob2scaledprob(poahmm->YY_boundary);
                        poahmm->YY_boundary_exit = prob2scaledprob(poahmm->YY_boundary_exit);

                }


                error_rate = 0.0f;
                for(j = 0; j < 1000; j++){
                        error_rate += 0.001f;
                        RUN(sim_seq_from_read_struct(&sd->seq, sd->al->read_structure[active_read_structure ] , i, sd->main_rng));

                        mutate_seq(sd->seq->seq,sd->seq->buffer,  sd->seq->len, error_rate , sd->main_rng, &num_error);
                        RUN(seq_to_internal(sd->seq->buffer, i, &i_seq, &i_len));

                        RUN(viterbi_poahmm_banded(poahmm,i_seq,i_len,path,2));

                        RUN(score_labelling(poahmm, path,sd->seq->label, &accuracy));



                        if(accuracy != 1.0 && num_error <= 5){
                                LOG_MSG("Wrong assignment at error rate: %f (%d)", error_rate, num_error);
//fprintf(stdout,"score: %f\n",e);
                                //forward_poahmm(poahmm, i_seq, i_len);
                                //backward_poahmm(poahmm, i_seq, i_len);
                                //random_poahmm(poahmm, i_seq, i_len);
                                e = poahmm->f_score - poahmm->random_scores[i];
                                e = exp2f(e) / (1.0 + exp2f(e));
                                forward_poahmm(poahmm, i_seq, i_len);
                                //backward_poahmm(poahmm, i_seq, i_len);

                                fprintf(stdout,"%d a:%f e:%f F:%5.3f\tB:%5.3f\tR:%5.3f nume:%d\t",i,accuracy,e, poahmm->f_score,poahmm->b_score, poahmm->random_scores[i], num_error);
                                fprintf(stdout,"%s\n",sd->seq->buffer);

                                if(accuracy != 1.0 && e >= 0.9){
                                        RUN(print_path(poahmm,path ,sd->seq->buffer, sd->seq->label ));
                                }
                                //fprintf(stdout,"%f\n",1.0 - scaledprob2prob(poahmm->f_score - poahmm->r_score));

                                //break;
                        }
                }
                MFREE(i_seq);

                i_seq = NULL;
        }
        MFREE(path);
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
        float e;
        float accuracy;
        float error_rate;

        int num_error;
        size = sizeof(in) / sizeof(char*);


        poahmm = sd->poahmm;
        RUN(read_arch_into_lib(sd->al, in, size));
        RUN(poahmm_from_read_structure(&poahmm, sd->gp, sd->al->read_structure[sd->al->num_arch-1  ], sd->a));
        //print_poahmm_param(poahmm);
        poahmm_to_dot(poahmm, "test2.dot");


        MMALLOC(path, sizeof(int)*poahmm->max_rank + poahmm->alloc_seq_len);
        LOG_MSG("Starting run");
        DECLARE_TIMER(timer);
        //RUN(set_random_scores(poahmm, 50));
        for(i = poahmm->max_rank ;i < poahmm->max_rank +1;i++){
                //RUN(sim_seq_from_read_struct(&sd->seq, sd->al->read_structure[1] , i, sd->main_rng));
                //RUN(generate_random_seq(sd->seq->seq, i, sd->main_rng));
                //RUN(generate_random_seq(&seq, &i, sd->main_rng ));
                //nif(i > 10){
                //RUN(insert_seq(seq, i, "AACCGGTT", 8, sd->main_rng));
                //}
                //RUN(seq_to_internal(sd->seq->seq, i, &i_seq, &i_len));
                //poahmm->r_len = 1.0f;//   MACRO_MAX(1, i - poahmm->max_rank);
                LOG_MSG("LEN: %d",i+3);
                if(i < poahmm->max_rank){
                        LOG_MSG("Sequence too short to match");
                }else if(poahmm->max_rank < i){
                        j = i - poahmm->max_rank;
                        poahmm->YY_boundary = (float)j / (float)(j + 2);
                        poahmm->YY_boundary_exit = 1.0 - poahmm->YY_boundary;
                        poahmm->YY_boundary = prob2scaledprob(poahmm->YY_boundary);
                        poahmm->YY_boundary_exit = prob2scaledprob(poahmm->YY_boundary_exit);

                }
                if(i >= poahmm->alloc_seq_len){
                        LOG_MSG("Resize");
                        resize_poahmm(poahmm, poahmm->num_nodes, i);
                        LOG_MSG("%d", poahmm->alloc_seq_len);
                }

                error_rate = 0.0;
                RUN(sim_seq_from_read_struct(&sd->seq, sd->al->read_structure[1] , i, sd->main_rng));
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

                                RUN(viterbi_poahmm_banded(poahmm,i_seq,i_len,path,c));
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
        RUN(init_nodes_from_read_structure(poahmm, al->read_structure[1],a,0));
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
        char etype[6] = "_EASIP";

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



int poahmm_to_dot(struct poahmm* poahmm,char* filename)
{
        FILE* fp_write;
        static int node_number = 0;
        int i,j;
        int nuc = 0;

        char color[10];

        int min_color[5][3] = {
                {5 ,255,5 }, // A green
                {255 ,5,5}, // C red
                {5 ,5,255}, // G blue
                {255 ,255,5}, // T yellow
                {125 ,125,125} // N black
        };
        int max_color[5][3] = {
                {255 ,255 ,255}, // A green
                {255 ,255 ,255}, // C red
                {255 ,255 ,255}, // G blue
                {255 ,255 ,255}, // T yellow
                {255 ,255 ,255} // N
        };
        double max_val = 0.0;
        for(i = 0 ;i < poahmm->num_nodes;i++){
                if(poahmm->nodes[i]->total_signal > max_val){
                        max_val = poahmm->nodes[i]->total_signal;
                }
        }

        //fillcolor="#0000ff"

        //if(!node_number){
                RUNP(fp_write = fopen(filename, "w"));
                //}else{
                //RUNP(fp_write = fopen(filename, "a"));
                //}
        fprintf(fp_write, " // The graph name and the semicolons are optional\n");
        fprintf(fp_write, "digraph graphname {\n");
        //fprintf(fp_write, "graph [fontname=\"Monospace\"]\n");
        fprintf(fp_write, "\"p\" [shape=record, label=\"Emission");

        for(i = 0; i < 4;i++){
                fprintf(fp_write, "| {%c|%7.2f|%7.2f|%7.2f|%7.2f} ", "ACGTN"[i], poahmm->emission_M[i][0], poahmm->emission_M[i][1], poahmm->emission_M[i][2], poahmm->emission_M[i][3]);

        }

        fprintf(fp_write, "| {%c| %3.2f | %3.2f | %3.2f | %3.2f } ", 'X', poahmm->emission_X[0], poahmm->emission_X[1], poahmm->emission_X[2], poahmm->emission_X[3]);
        fprintf(fp_write, "| {%c| %3.2f | %3.2f | %3.2f | %3.2f } ", 'Y', poahmm->emission_Y[0], poahmm->emission_Y[1], poahmm->emission_Y[2], poahmm->emission_Y[3]);
        fprintf(fp_write, "| {%c| %3.2f | %3.2f | %3.2f | %3.2f } ", 'B', poahmm->background[0], poahmm->background[1], poahmm->background[2], poahmm->background[3]);

        fprintf(fp_write, "\"]\n");

        fprintf(fp_write, "rankdir=\"LR\";\n");
        for(i = 0 ;i < poahmm->num_nodes;i++){
                nuc =poahmm->nodes[i]->nuc;
                snprintf(color, 10, "#%02X%02X%02X",	 min_color[nuc][0] + (int)((max_color[nuc][0] - min_color[nuc][0]) *( 1.0 - ((double)poahmm->nodes[i]->total_signal  / max_val))),
                         min_color[nuc][1] + (int)((max_color[nuc][1] - min_color[nuc][1]) * ( 1.0 - ((double)poahmm->nodes[i]->total_signal  / max_val))),
                         min_color[nuc][2] + (int)((max_color[nuc][2] - min_color[nuc][2]) *( 1.0 - ((double)poahmm->nodes[i]->total_signal  / max_val))));



                fprintf(fp_write, "n%d [label=\"%c,%d\n%d\nS:%f;E:%f;\" ,style=filled  fillcolor=\"%s\"];\n", i+node_number,"ACGTN"[poahmm->nodes[i]->nuc], poahmm->nodes[i]->rank, poahmm->nodes[i]->identifier,poahmm->entry_probabilities[i],poahmm->exit_probabilities[i],color);



        }

        for(i = 0 ;i < poahmm->num_nodes;i++){
                for(j = 0 ;j < poahmm->num_nodes;j++){
                        if(poahmm->poa_graph[i][j] != prob2scaledprob(0.0)){
                                fprintf(fp_write, "n%d -> n%d [ label=\"p = %0.2f\" ];\n",i+node_number,j+node_number, scaledprob2prob(poahmm->poa_graph[i][j]));
                        }
                }
        }
        node_number +=poahmm->num_nodes;
        fprintf(fp_write, "}\n");
        fclose(fp_write);
        return OK;
ERROR:
        return FAIL;
}

int init_nodes_from_read_structure(struct poahmm* poahmm, struct read_structure* rs, struct alphabet* a, int random)
{
        int i,j,c;
        int* b_nodes = NULL;
        int* e_nodes_new = NULL;
        int* e_nodes = NULL;
        int* tmp;
        int n_index;
        int num_e;
        int skipN;
        //int num_b;
        int num_e_new;
        int num_nodes;
        struct poahmm_node* node_ptr = NULL;
        struct segment_specs* s;

        num_nodes = 0;

        for(i= 0; i < rs->num_segments;i++){
                s = rs->seg_spec[i];
                for(j = 0; j < s->num_seq;j++){
                        if(s->max_len == INT32_MAX  ){
                                ERROR_MSG("N+ segments not len assigned");
                        }
                }
        }

        for(i= 0; i < rs->num_segments;i++){
                s = rs->seg_spec[i];
                for(j = 0; j < s->num_seq;j++){
                        if(s->max_len == INT32_MAX  ){
                                num_nodes++;
                        }else{
                                num_nodes+= s->max_len;
                        }
                        fprintf(stdout,"%s\n",s->seq[j]);
                }
                //if(j = )
        }
        //num_nodes+10;
        if(poahmm->alloced_num_nodes <= num_nodes){
                RUN(resize_poahmm(poahmm, num_nodes, 0));
                //ERROR_MSG("No enough space in poahmm: want %d have %d", num_nodes,poahmm->alloced_num_nodes);
        }
        LOG_MSG("Need  %d nodes", num_nodes);
        //exit(0);
        for(i = 0; i< poahmm->alloced_num_nodes;i++){
                for(j = 0 ; j < poahmm->alloced_num_nodes;j++){
                        poahmm->poa_graph[i][j] = prob2scaledprob(0.0f);
                        poahmm->e_poa_graph[i][j] = 0.0;
                }
                poahmm->entry_probabilities[i] = prob2scaledprob(0.0f);
                poahmm->e_entry_probabilities[i] = 0.0;

                poahmm->exit_probabilities[i] = prob2scaledprob(0.0f);
                poahmm->e_exit_probabilities[i] = 0.0;
        }

        poahmm->YY_boundary = prob2scaledprob(0.0f);
        poahmm->YY_boundary_exit = prob2scaledprob(1.0f);
        MMALLOC(e_nodes, sizeof(int) * num_nodes);
        MMALLOC(b_nodes, sizeof(int) * num_nodes);
        MMALLOC(e_nodes_new, sizeof(int) * num_nodes);
        num_e = 0;
        //num_b = 0;
        num_e_new = 0;
        n_index = 0;

        e_nodes[num_e] = 0;
        num_e++;
        for(i= 0; i < rs->num_segments;i++){
                s = rs->seg_spec[i];
                skipN = 0;
                if(s->num_seq> 1){
                        skipN =1;
                }
                for(j = skipN;j < s->num_seq;j++){
                        for(c = 0; c < num_e;c++){
                                if(!i){
                                        poahmm->e_entry_probabilities[n_index] = 1.0;
                                }else{
                                        poahmm->e_poa_graph[e_nodes[c]][n_index] = 1.0;
                                }
                        }
                        for(c = 0; c < s->max_len;c++){
                                if(c >= s->min_len-1){
                                        e_nodes_new[num_e_new] = n_index;
                                        num_e_new++;
                                }
                                if(c == s->max_len -1){
                                        e_nodes_new[num_e_new] = n_index;
                                        num_e_new++;
                                }else{
                                        poahmm->e_poa_graph[n_index][n_index+1] = 1.0;
                                }
                                node_ptr = poahmm->nodes[n_index];
                                if(random){
                                        node_ptr->nuc = 4;
                                }else{
                                        node_ptr->nuc = tlalphabet_get_code(a,s->seq[j][c]);
                                }
                                node_ptr->rank = UINT32_MAX;
                                node_ptr->type = s->extract;
                                node_ptr->identifier = n_index;
                                n_index++;

                        }
                }
                num_e = num_e_new;
                num_e_new = 0;
                tmp = e_nodes;
                e_nodes = e_nodes_new;
                e_nodes_new = tmp;
        }

        for(i = 0; i < num_e;i++){
                poahmm->e_exit_probabilities[e_nodes[i]] = 1.0;
        }
        //viterbi_poahmm(struct poahmm *poahmm, uint8_t *seq, int len, uint32_t *path)
//stop state
        //node_ptr = poahmm->nodes[POAHMM_ENDSTATE];
        //node_ptr->nuc = 4;
        //node_ptr->identifier = POAHMM_ENDSTATE;
        //exit(0);
        //poahmm->num_nodes = len;
        LOG_MSG("assigned %d",n_index);
        poahmm->num_nodes = n_index;
        for(i = 0; i < poahmm->alloced_num_nodes;i++){
                poahmm->nodes[i]->rank = UINT32_MAX;
        }
        //poahmm->alloced_num_nodes
        MFREE(e_nodes);
        MFREE(e_nodes_new);
        MFREE(b_nodes);

        RUN(set_rank_transition_poahmm(poahmm));
        poahmm->max_rank = 0;
        for(i = 0; i< poahmm->num_nodes;i++){
                if(poahmm->max_rank < poahmm->nodes[i]->rank){
                        poahmm->max_rank = poahmm->nodes[i]->rank;

                }
                /*fprintf(stdout,"RANK: %d %d", i, poahmm->nodes[i]->rank);
                fprintf(stdout,"ENTRY: %f EXIT:%f\t", poahmm->entry_probabilities[i],poahmm->exit_probabilities[i]);
                for(j = 0 ; j < poahmm->num_nodes;j++){
                        fprintf(stdout,"%f ",poahmm->poa_graph[i][j]);
                        //poahmm->e_poa_graph[i][j] = 0.0;
                }
                fprintf(stdout,"\n");*/
        }
        //exit(0);
        return OK;
ERROR:
        return FAIL;
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
        sd->gp->max_seq_len = 2048;

        RUNP(sd->poahmm = init_poahmm( sd->gp));

        RUN(sim_seq_from_read_struct(&sd->seq, NULL,10,sd->main_rng));
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

int sim_seq_from_read_struct(struct sim_seq** simseq, struct read_structure* rs,int sim_len,struct rng_state* rng)
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
                MMALLOC(s->seq, sizeof(char)* s->alloc_len);
                MMALLOC(s->label, sizeof(char)* s->alloc_len);
                MMALLOC(s->buffer, sizeof(char)* s->alloc_len);
        }
        if(sim_len >= s->alloc_len){
                s->alloc_len = s->alloc_len + s->alloc_len /2;
                MREALLOC(s->seq, sizeof(char)* s->alloc_len);
                MREALLOC(s->label, sizeof(char)* s->alloc_len);
                MREALLOC(s->buffer, sizeof(char)* s->alloc_len);
        }

        /* step one: add random residues  */
        RUN(generate_random_seq(&s->seq,&sim_len,rng));
        for(i = 0; i < sim_len;i++){
                s->label[i] = etype[0];
        }

        s->label[sim_len] = 0;
        s->len = sim_len;

        if(rs){
                //total_len = 0;
                c = 0;
                for(i = 0; i < rs->num_segments;i++){
                        spec = rs->seg_spec[i];
                        g = tl_random_int(rng, spec->max_len- spec->min_len) + spec->min_len;
                        b = 0;
                        if(spec->num_seq >1){
                                b = tl_random_int(rng, spec->num_seq-1) + 1;
                        }
                        for(j = 0; j < g;j++){
                                if( spec->seq[b][j] == 'N'){
                                        s->buffer[c] = alphabet[tl_random_int(rng, 4)];
                                }else{
                                        s->buffer[c] = spec->seq[b][j];
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
                MFREE(s);
        }
}
