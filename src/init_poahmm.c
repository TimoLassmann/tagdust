
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
        ph = *poahmm;
        if(!ph){
                RUNP(ph = init_poahmm(p));
        }

        random = 1;
        RUN(init_nodes_from_read_structure(ph, rs,a,random));
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

        RUN(init_nodes_from_read_structure(ph, rs,a,0));
        *poahmm = ph;
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
                fprintf(stdout,"\n");
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
