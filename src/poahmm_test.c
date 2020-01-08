
#include <stdint.h>
#include <string.h>
#include "tldevel.h"

#include "arch_lib.h"

#include "tlrng.h"
#include "tllogsum.h"
#include "tlalphabet.h"

#include "poahmm.h"



#include "sim_seq_lib.h"
int single_seq_test(void);
int arch_test(void);

int init_nodes_from_read_structure(struct poahmm* poahmm, struct read_structure* rs, struct alphabet* a);

int print_poa_graph(struct poahmm* poahmm);

int print_path(struct poahmm* poahmm, int* path);
int print_poahmm_param(struct poahmm* poahmm);
int poahmm_to_dot(struct poahmm* poahmm,char* filename);



int main(int argc, char *argv[])
{
        RUN(single_seq_test());
        RUN(arch_test());
        return EXIT_SUCCESS;
ERROR:
        return EXIT_FAILURE;
}

int single_seq_test(void)
{
        struct rng_state* rng = NULL;
        struct poahmm* poahmm = NULL;
        char* seq = NULL;
        uint8_t* i_seq = NULL;
        int* path = NULL;
        int len = 32;
        int i_len;
        int nuc_count[4];

        int i;

        int malloced_path = 100;

        MMALLOC(path , sizeof(int) *malloced_path);


        for(i = 0; i < 4;i++){
                nuc_count[i] = 1000;
        }
        RUNP(poahmm = init_poahmm(128,nuc_count,500.0));

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


        RUN(print_path(poahmm,path));
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
        int* path = NULL;
        int len = 32;
        int i_len;

        char* in[] = {
                "I:A{0,3}",
                "S:GTA,AAC",
                "MYREAD:E:NNN",
                "I:CCTTAA",
                "S:ACAGTG,ACTTGA,TTAGGC" ,
                "R2:E:NNNNNNNNNNNNNNNNNNNNNN"
        };
        int nuc_count[4];

        int size;
        int i;
        int malloced_path = 100;

        MMALLOC(path , sizeof(int) *malloced_path);


        RUN(alloc_arch_lib(&al));

        size= sizeof(in) / sizeof(char*);

        RUN(read_arch_into_lib(al, in, size));

        LOG_MSG("Read in %d architectures.",al->num_arch);
        for(i = 0; i < 4;i++){
                nuc_count[i] = 1000;
        }
        RUNP(rng= init_rng(0));
        RUNP(poahmm = init_poahmm(128,nuc_count,500.0));
        RUN(create_alphabet(&a, rng, TLALPHABET_DEFAULT_DNA ));
        RUN(init_nodes_from_read_structure(poahmm, al->read_structure[1],a));

        //print_poahmm_param(poahmm);
        //exit(0);
        RUN(generate_random_seq(&seq, &len, rng));

        char match_seq[] = "TTTTTTCTTTCCTTAATTATGGCACGTACGT";
        len = strlen(match_seq);
                
        RUN(seq_to_internal(match_seq, len, &i_seq, &i_len));

        forward_poahmm(poahmm, i_seq, i_len);
        backward_poahmm(poahmm, i_seq, i_len);
        random_poahmm(poahmm, i_seq, i_len);
        fprintf(stdout,"F:%f\nB:%f\nR:%f\n", poahmm->f_score,poahmm->b_score, poahmm->r_score);
        RUN(viterbi_poahmm(poahmm,i_seq,i_len,path));
        fprintf(stdout,"%s\n",match_seq);
        RUN(print_path(poahmm,path));
        fprintf(stdout,"F:%f\n", poahmm->f_score);


        poahmm_to_dot(poahmm, "test.dot");



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


int print_path(struct poahmm* poahmm, int* path)
{
        int i;
        char alphabet[5] = "ACGTN";
        char etype[5] = "_EASIP";
        fprintf(stdout, "PATH:\n");

        for(i = 1; i < path[0];i++){
                fprintf(stdout, " %3c", alphabet[path[i] >> 28]);
        }
        fprintf(stdout,"\n");
        for(i = 1; i < path[0];i++){

                if((path[i] & 0xFFFFFFF) == 0xFFFFFFF){
                        fprintf(stdout, " %3c",'-');
                }else{
                        fprintf(stdout, " %3c", alphabet[poahmm->nodes[path[i] & 0xFFFFFFF]->nuc]);
                }

        }
        fprintf(stdout,"\n");
        for(i = 1; i < path[0];i++){
                if((path[i] & 0xFFFFFFF) == 0xFFFFFFF){
                        fprintf(stdout, " %3c",'-');
                }else{
                        fprintf(stdout, " %3d", path[i] & 0xFFFFFFF);
                }
        }
        fprintf(stdout,"\n");
        for(i = 1; i < path[0];i++){
                if((path[i] & 0xFFFFFFF) == 0xFFFFFFF){
                        fprintf(stdout, " %3c",'-');
                }else{
                        fprintf(stdout, " %3c",  etype[ poahmm->nodes[path[i] & 0xFFFFFFF]->total_signal]);
                }
        }
        fprintf(stdout,"\n");

        return OK;
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

int init_nodes_from_read_structure(struct poahmm* poahmm, struct read_structure* rs, struct alphabet* a)
{
        int i,j,c;
        int* b_nodes = NULL;
        int* e_nodes_new = NULL;
        int* e_nodes = NULL;
        int* tmp;
        int n_index;
        int num_e;
        int num_b;
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
                }
        }
        if(poahmm->alloced_num_nodes <= num_nodes){
                ERROR_MSG("No enough space in poahmm");
        }
        LOG_MSG("Need  %d nodes", num_nodes);

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

        MMALLOC(e_nodes, sizeof(int) * 1024);
        MMALLOC(b_nodes, sizeof(int) * 1024);
        MMALLOC(e_nodes_new, sizeof(int) * 1024);
        num_e = 0;
        num_b = 0;
        num_e_new = 0;
        n_index = 0;

        e_nodes[num_e] = 0;
        num_e++;
        for(i= 0; i < rs->num_segments;i++){
                s = rs->seg_spec[i];
                for(j = 0;j < s->num_seq;j++){
                        for(c = 0; c < num_e;c++){
                                if(!i){
                                        poahmm->e_entry_probabilities[n_index] = 1.0;
                                }else{
                                        poahmm->e_poa_graph[e_nodes[c]][n_index] = 1.0;
                                }

                        }
                        for(c = 0; c < s->max_len;c++){
                                if(i == rs->num_segments-1){
                                        if(c > s->max_len -15 ){
                                                e_nodes_new[num_e_new] = n_index;
                                                num_e_new++;
                                        }else{
                                                poahmm->e_poa_graph[n_index][n_index+1] = 1.0;
                                        }
                                }else if(c == s->max_len -1){
                                        e_nodes_new[num_e_new] = n_index;
                                        num_e_new++;
                                }else{
                                        poahmm->e_poa_graph[n_index][n_index+1] = 1.0;
                                }
                                node_ptr = poahmm->nodes[n_index];
                                node_ptr->nuc = tlalphabet_get_code(a,s->seq[j][c]);
                                node_ptr->rank = UINT32_MAX;
                                //node_ptr->part = (s->extract == ARCH_ETYPE_PARTIAL) ? 1: 0;
                                //node_ptr->signal[poahmm_data->seq_id[index]] = 1;
                                node_ptr->total_signal = s->extract;
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

//stop state
        //node_ptr = poahmm->nodes[POAHMM_ENDSTATE];
        //node_ptr->nuc = 4;
        //node_ptr->identifier = POAHMM_ENDSTATE;

        //poahmm->num_nodes = len;
        poahmm->num_nodes = n_index;
        for(i = 0; i < poahmm->num_nodes;i++){
                poahmm->nodes[i]->rank = UINT32_MAX;
        }
        MFREE(e_nodes);
        MFREE(e_nodes_new);
        MFREE(b_nodes);

        RUN(set_rank_transition_poahmm(poahmm));

        /*for(i = 0; i< poahmm->num_nodes;i++){

                fprintf(stdout,"ENTRY: %f EXIT:%f\t", poahmm->entry_probabilities[i],poahmm->exit_probabilities[i]);
                for(j = 0 ; j < poahmm->num_nodes;j++){
                        fprintf(stdout,"%f ",poahmm->poa_graph[i][j]);
                        //poahmm->e_poa_graph[i][j] = 0.0;
                }
                fprintf(stdout,"\n");
                }*/


//poahmm =  set_rank(poahmm, 0, 0);

        //RUN(reset_to_from_index(poahmm));
        //RUN(reset_poa_graph_transitions_based_on_counts(poahmm));

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
        for(i = 0; i < 4;i++){
                sum = prob2scaledprob(0.0);
                for(j = 0; j < 4;j++){
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

        return OK;
}
