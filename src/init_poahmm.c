


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

static int set_len_of_unknown_poa(const struct read_structure*rs, int* min_plus_len,int* max_plus_len, int min_seq_len, int max_seq_len);

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
        LOG_MSG("Plus : %d -  %d", plus_min_len,plus_max_len);
        RUN(init_nodes_from_read_structure(ph, rs,a,random, plus_min_len, plus_max_len));

        MMALLOC(path, sizeof(int)* p->max_seq_len *2);
        MMALLOC(nnn, sizeof(uint8_t) * (p->max_seq_len+1));
        MMALLOC(qqq, sizeof(uint8_t) * (p->max_seq_len+1));


        //ph->max_rank
        for(i = p->min_seq_len; i < p->max_seq_len;i++){
                ph->random_scores[i] = prob2scaledprob(0.0f);
        }
        for(i = ph->max_rank+1; i < p->max_seq_len;i++){
                //LOG_MSG("%d", i);
                //for(c = 0; c < 10;c++){
                for(j = 0; j < p->max_seq_len+1;j++){
                        nnn[j] = j %4 ;//tl_random_int(rng, 4);
                        qqq[j] = 0;
                }
                if(i < ph->max_rank){
                        ERROR_MSG("Sequence too short to match");
                }else if(ph->max_rank < i){
                        j = i - ph->max_rank;
                        ph->YY_boundary = (float)j / (float)(j + 2);
                        ph->YY_boundary_exit = 1.0 - ph->YY_boundary;
                        ph->YY_boundary = prob2scaledprob(ph->YY_boundary);
                        ph->YY_boundary_exit = prob2scaledprob(ph->YY_boundary_exit);

                }


                //RUN(viterbi_poahmm(ph, nnn, i, path));
                //fprintf(stdout,"%f ", ph->f_score);
                RUN(viterbi_poahmm_banded(ph, nnn,qqq, i, path,1));
                //fprintf(stdout,"%f ", ph->f_score);
                //print_path(poahmm, path,nnn,nnn);
                ph->random_scores[i] = ph->f_score;
                //fprintf(stdout,"%f ", ph->f_score);
                //}

                  //fprintf(stdout,"\n");
                //exit(0);
        }
        //exit(0);
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




        LOG_MSG("Model len: %d -> %d", poahmm->min_model_len, poahmm->max_model_len);
        LOG_MSG("Seq len: %d -> %d", poahmm->min_seq_len, poahmm->max_seq_len);

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

        //LOG_MSG("MAXRANK: %d %d", poahmm->max_rank,poahmm->num_nodes);


        //exit(0);
        return OK;
ERROR:
        return FAIL;
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

static int set_len_of_unknown_poa(const struct read_structure*rs, int* min_plus_len,int* max_plus_len, int min_seq_len, int max_seq_len)
{
        int i;

        int known_len;
        int num_unknown;
        int min;

        known_len = 0;
        num_unknown = 0;

        for(i = 0; i < rs->num_segments  ;i++){
                if(rs->seg_spec[i]->max_len == INT32_MAX){
                        num_unknown++;
                }else{
                        known_len += rs->seg_spec[i]->max_len;
                }
        }
        *min_plus_len = 0;
        *max_plus_len = 0;
        if(num_unknown){
                *max_plus_len = (int) ceilf (((float) max_seq_len - (float) known_len) / (float) num_unknown);
                *min_plus_len = (int) floorf(((float) min_seq_len - (float) known_len) / (float) num_unknown);
        }
        return OK;
ERROR:
        return FAIL;
}
