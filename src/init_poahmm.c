


#include <stdint.h>
#include <string.h>
#include "tldevel.h"

#include "arch_lib.h"

//#include "tlrng.h"
#include "tllogsum.h"
#include "tlalphabet.h"

#include "poahmm.h"
#include "poahmm_structs.h"

#include "init_poahmm.h"



static int set_entry_prob(struct poahmm* poahmm, int* e_nodes, int num, int at_start,int index);

int poahmm_from_read_structure(struct poahmm** poahmm,struct global_poahmm_param* p, struct read_structure* rs,struct alphabet* a)
{
        struct poahmm* ph = NULL;
        uint8_t* nnn = NULL;
        uint8_t* qqq = NULL;
        uint32_t* path = NULL;
        int random;
        int i;
        int j;
        int c;
        int plus_min_len;
        int plus_max_len;
        ph = *poahmm;
        if(!ph){
                RUNP(ph = init_poahmm(p));
        }
        ph->max_seq_len = p->max_seq_len;
        ph->min_seq_len = p->min_seq_len;

        random = 1;

        RUN(set_len_of_unknown_poa(rs,&plus_min_len,&plus_max_len, p->min_seq_len, p->max_seq_len));

        if(plus_max_len == -1 || plus_min_len == -1){
                WARNING_MSG("Model too long from sequences.");
                *poahmm = ph;
                return OK;
        }
        RUN(init_nodes_from_read_structure(ph, rs,a,random, plus_min_len, plus_max_len));



        //LOG_MSG("Model Len: %d %d", ph->min_model_len,ph->max_model_len);

        //LOG_MSG("Seq Len: %d %d alloc: %d", ph->min_seq_len, ph->max_seq_len, ph->alloc_seq_len);

        //LOG_MSG("Current poa size: nodes: %d seq_len:%d", ph->alloced_num_nodes , ph->alloc_seq_len);

        //LOG_MSG("Writnign ");
        //poahmm_to_dot(ph, "test_singleN.dot");
        //exit(0);


        c = MACRO_MAX(p->max_seq_len, ph->max_model_len) +2;
        //LOG_MSG("path will be %d long", c*2);
        MMALLOC(path, sizeof(int)* c * 2);
        MMALLOC(nnn, sizeof(uint8_t) * c);
        MMALLOC(qqq, sizeof(uint8_t) * c);

        for(j = 0; j < c;j++){
                nnn[j] = j %4 ;//tl_random_int(rng, 4);
                qqq[j] = 0;
        }

        //ph->max_rank
        for(i = p->min_seq_len; i < p->max_seq_len;i++){
                ph->random_scores[i] = prob2scaledprob(0.0f);
        }
        for(i = ph->min_model_len; i <=  ph->max_seq_len ;i++){
                //LOG_MSG("tsting length : %d (seq_len: %d)", i, c);
                //for(c = 0; c < 10;c++){
                //set_terminal_gap_prob(ph, i);
                RUN(viterbi_poahmm_banded(ph, nnn,qqq, i, NULL,2));
                ph->random_scores[i] = ph->f_score;

        }

        MFREE(path);
        MFREE(nnn);
        MFREE(qqq);

        RUN(init_nodes_from_read_structure(ph, rs,a,0, plus_min_len,plus_max_len));

        *poahmm = ph;
        return OK;
ERROR:
        return FAIL;
}



int init_nodes_from_read_structure(struct poahmm* poahmm, struct read_structure* rs, struct alphabet* a, int random, int plus_min_len, int plus_max_len)
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
        int min,max;
        int is_plus;
        int is_L;
        int is_R;
        struct poahmm_node* node_ptr = NULL;
        struct segment_specs* s;

        num_nodes = 0;

        poahmm->max_model_len = 0;
        poahmm->min_model_len = 0;
        for(i= 0; i < rs->num_segments;i++){
                s = rs->seg_spec[i];
                if(s->max_len == INT32_MAX  ){
                        poahmm->min_model_len += plus_min_len;
                        poahmm->max_model_len += plus_max_len;
                }else{
                        poahmm->min_model_len += s->min_len;
                        poahmm->max_model_len += s->max_len;
                }
                //print_segment_spec(s);
        }

        for(i= 0; i < rs->num_segments;i++){
                s = rs->seg_spec[i];
                for(j = 0; j < s->num_seq;j++){
                        if(s->max_len == INT32_MAX  ){
                                num_nodes+= plus_max_len;
                        }else{
                                num_nodes+= s->max_len;
                        }
                        //fprintf(stdout,"%s\n",s->seq[j]);
                }
                //fprintf(stdout,"\n");
        }
        //num_nodes+10;
        if(poahmm->alloced_num_nodes <= num_nodes+10){
                RUN(resize_poahmm(poahmm, num_nodes+10, 0));
                //ERROR_MSG("No enough space in poahmm: want %d have %d", num_nodes,poahmm->alloced_num_nodes);
        }
        //LOG_MSG("Need  %d nodes", num_nodes);
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

        MMALLOC(e_nodes, sizeof(int) * num_nodes*2);
        MMALLOC(b_nodes, sizeof(int) * num_nodes*2);
        MMALLOC(e_nodes_new, sizeof(int) * num_nodes*2);
        num_e = 0;
        //num_b = 0;
        num_e_new = 0;
        n_index = 0;

        e_nodes[num_e] = 0;
        num_e++;
        for(i= 0; i < rs->num_segments;i++){
                s = rs->seg_spec[i];
                min = s->min_len;
                max = s->max_len;
                is_plus = 0;
                is_L = 0;
                is_R = 0;

                if(s->extract == ARCH_ETYPE_WOBBLE_LEFT){
                        is_L = 1;
                }
                if(s->extract == ARCH_ETYPE_WOBBLE_RIGHT){
                        is_R = 1;
                }

                if(max == INT32_MAX){
                        min = plus_min_len;
                        max = plus_max_len;
                        is_plus = 1;
                }


                skipN = 0;
                if(s->num_seq> 1){
                        skipN =1;
                }
                for(j = skipN;j < s->num_seq;j++){
                        //if(!is_L){
                        //RUN(set_entry_prob(poahmm, e_nodes, num_e, !i,n_index));
                        //}
                        /*for(c = 0; c < num_e;c++){

                                if(!i){
                                        poahmm->e_entry_probabilities[n_index] = 1.0;
                                }else{
                                        poahmm->e_poa_graph[e_nodes[c]][n_index] = 1.0;
                                }
                                }*/
                        for(c = 0; c < max;c++){
                                if(is_L){
                                        //LOG_MSG("Setting L)");
                                        RUN(set_entry_prob(poahmm, e_nodes, num_e, !i,n_index));
                                }else{
                                        if(c == 0){
                                                RUN(set_entry_prob(poahmm, e_nodes, num_e, !i,n_index));
                                        }
                                }
                                if(is_R){
                                        e_nodes_new[num_e_new] = n_index;
                                        num_e_new++;
                                }else{
                                        if(c >= min-1){
                                                //LOG_MSG("%d", num_e_new);
                                                e_nodes_new[num_e_new] = n_index;
                                                num_e_new++;
                                        }
                                }
                                if(c == max -1){
                                        e_nodes_new[num_e_new] = n_index;
                                        num_e_new++;
                                }else{
                                        poahmm->e_poa_graph[n_index][n_index+1] = 1.0;
                                }
                                node_ptr = poahmm->nodes[n_index];
                                if(random){
                                        node_ptr->nuc = 4;
                                }else{
                                        if(is_plus){
                                                node_ptr->nuc = tlalphabet_get_code(a,s->seq[j][0]);
                                        }else{
                                                node_ptr->nuc = tlalphabet_get_code(a,s->seq[j][c]);
                                        }
                                }
                                node_ptr->rank = -1 * INT32_MAX;
                                node_ptr->type = s->extract;
                                node_ptr->identifier = n_index;
                                node_ptr->alt = j;
                                node_ptr->segment = i;
                                n_index++;

                        }
                }
                /*for(c = 0; c < num_e_new;c++){
                        fprintf(stdout,"%d ",e_nodes_new[c]);
                }
                fprintf(stdout,"\n");*/
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
        //LOG_MSG("assigned %d",n_index);
        poahmm->num_nodes = n_index;
        /*for(i = 0; i < poahmm->alloced_num_nodes;i++){
                poahmm->nodes[i]->rank = INT32_MAX;
                }*/

        MFREE(e_nodes);
        MFREE(e_nodes_new);
        MFREE(b_nodes);


        RUN(set_rank_transition_poahmm(poahmm));
        //LOG_MSG("Model len: %d -> %d", poahmm->min_model_len, poahmm->max_model_len);
        //LOG_MSG("Seq len: %d -> %d", poahmm->min_seq_len, poahmm->max_seq_len);
        /*poahmm->max_rank = 0;
        for(i = 0; i< poahmm->num_nodes;i++){
                if(poahmm->max_rank < poahmm->nodes[i]->rank){
                        poahmm->max_rank = poahmm->nodes[i]->rank;
                        }*/
                //fprintf(stdout,"RANK: %d %d\n", i, poahmm->nodes[i]->rank);
                /*fprintf(stdout,"ENTRY: %f EXIT:%f\t", poahmm->entry_probabilities[i],poahmm->exit_probabilities[i]);
                for(j = 0 ; j < poahmm->num_nodes;j++){
                        fprintf(stdout,"%f ",poahmm->poa_graph[i][j]);
                        //poahmm->e_poa_graph[i][j] = 0.0;
                }
                fprintf(stdout,"\n");*/
        //}
        //LOG_MSG("MAXRANK: %d %d", poahmm->max_rank,poahmm->num_nodes);
        //exit(0);
        return OK;
ERROR:
        return FAIL;
}

int set_terminal_gap_prob(struct poahmm* poahmm, int seq_len)
{

        poahmm->YY_boundary = prob2scaledprob(0.0f);
        poahmm->YY_boundary_exit = prob2scaledprob(1.0f);
        if(poahmm->min_model_len < seq_len){
                float target_len = (float)seq_len - ( (float) poahmm->max_model_len - (float) poahmm->min_model_len) /4.0;
                //j = i - poahmm->min_model_len
                poahmm->YY_boundary = target_len / (target_len+1.0);
                poahmm->YY_boundary_exit = 1.0 - poahmm->YY_boundary;
                poahmm->YY_boundary = prob2scaledprob(poahmm->YY_boundary);
                poahmm->YY_boundary_exit = prob2scaledprob(poahmm->YY_boundary_exit);
        }
        return OK;

}


int set_entry_prob(struct poahmm* poahmm, int* e_nodes, int num, int at_start,int index)
{
        int i;
        for(i = 0; i < num;i++){
                if(at_start){
                        poahmm->e_entry_probabilities[index] = 1.0;
                }else{
                        poahmm->e_poa_graph[e_nodes[i]][index] = 1.0;
                }
        }
        return OK;
}

int set_len_of_unknown_poa(const struct read_structure*rs, int* min_plus_len,int* max_plus_len, int min_seq_len, int max_seq_len)
{
        int i;
        int known_len;
        int num_unknown;
        int min;

        known_len = 0;
        num_unknown = 0;

        for(i = 0; i < rs->num_segments  ;i++){
                //LOG_MSG("%d: %d->%d %s", i, rs->seg_spec[i]->min_len, rs->seg_spec[i]->max_len, rs->seg_spec[i]->seq[0]);
                if(rs->seg_spec[i]->max_len == INT32_MAX){
                        num_unknown++;
                }else{
                        known_len += rs->seg_spec[i]->max_len;
                }
        }
        *min_plus_len = 0;
        *max_plus_len = 0;
        //LOG_MSG("NumUnknown: %d len known:%d    %d %d ", num_unknown,known_len, min_seq_len, max_seq_len);


        if(num_unknown){
                if(known_len >= max_seq_len){
                        /* oh dear sequence too short to fit this models  */
                        *max_plus_len = -1;
                        *min_plus_len = -1;
                }else{
                        *max_plus_len = (int) ceilf (((float) max_seq_len - (float) known_len) / (float) num_unknown);
                        *min_plus_len = (int) floorf(((float) min_seq_len - (float) known_len) / (float) num_unknown);
                }
        }
        return OK;
ERROR:
        return FAIL;
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
        double max_val = 1.0;
        /*for(i = 0 ;i < poahmm->num_nodes;i++){
                if(poahmm->nodes[i]->total_signal > max_val){
                        max_val = poahmm->nodes[i]->total_signal;
                }
                }*/

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
                //dLOG_MSG("i:%d", i);
                fprintf(fp_write, "| {%c|%7.2f|%7.2f|%7.2f|%7.2f} ", "ACGTN"[i], poahmm->emission_M[i][0], poahmm->emission_M[i][1], poahmm->emission_M[i][2], poahmm->emission_M[i][3]);

        }

        fprintf(fp_write, "| {%c| %3.2f | %3.2f | %3.2f | %3.2f } ", 'X', poahmm->emission_X[0], poahmm->emission_X[1], poahmm->emission_X[2], poahmm->emission_X[3]);
        fprintf(fp_write, "| {%c| %3.2f | %3.2f | %3.2f | %3.2f } ", 'Y', poahmm->emission_Y[0], poahmm->emission_Y[1], poahmm->emission_Y[2], poahmm->emission_Y[3]);
        fprintf(fp_write, "| {%c| %3.2f | %3.2f | %3.2f | %3.2f } ", 'B', poahmm->background[0], poahmm->background[1], poahmm->background[2], poahmm->background[3]);

        fprintf(fp_write, "\"]\n");

        fprintf(fp_write, "rankdir=\"LR\";\n");
        node_number = 0;
        /* start  */
        fprintf(fp_write, "n%d [label=\"START\"];\n",node_number);
        node_number++;

        /* end */
        fprintf(fp_write, "n%d [label=\"END\"];\n",node_number);

        node_number++;
        for(i = 0 ;i < poahmm->num_nodes;i++){
                nuc =poahmm->nodes[i]->nuc;
                snprintf(color, 10, "#%02X%02X%02X",	 min_color[nuc][0] + (int)((max_color[nuc][0] - min_color[nuc][0]) *( 1.0 - ((double)1.0  / max_val))),
                         min_color[nuc][1] + (int)((max_color[nuc][1] - min_color[nuc][1]) * ( 1.0 - ((double)1.0  / max_val))),
                         min_color[nuc][2] + (int)((max_color[nuc][2] - min_color[nuc][2]) *( 1.0 - ((double)1.0  / max_val))));

                snprintf(color, 10, "#%02X%02X%02X",	 min_color[nuc][0] + (int)((max_color[nuc][0] - min_color[nuc][0]) *0.1),
                         min_color[nuc][1] + (int)((max_color[nuc][1] - min_color[nuc][1])  * 0.1),
                         min_color[nuc][2] + (int)((max_color[nuc][2] - min_color[nuc][2]) * 0.1));



                fprintf(fp_write, "n%d [label=\"%c,%d\n%d\nS:%f;E:%f;\" ,style=filled  fillcolor=\"%s\"];\n", i+node_number,"ACGTN"[poahmm->nodes[i]->nuc], poahmm->nodes[i]->rank, poahmm->nodes[i]->identifier,poahmm->entry_probabilities[i],poahmm->exit_probabilities[i],color);



        }
        for(i = 0 ;i < poahmm->num_nodes;i++){
                if(poahmm->entry_probabilities[i] != prob2scaledprob(0.0)){
                        fprintf(fp_write, "n%d -> n%d [ label=\"p = %0.2f\" ];\n",0,i+node_number, scaledprob2prob(poahmm->entry_probabilities[i]));
                }

                if(poahmm->exit_probabilities[i] != prob2scaledprob(0.0)){
                        fprintf(fp_write, "n%d -> n%d [ label=\"p = %0.2f\" ];\n",i+node_number,1, scaledprob2prob(poahmm->exit_probabilities[i]));
                }
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
