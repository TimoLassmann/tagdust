
#include <stdint.h>
#include "tldevel.h"
#include "tllogsum.h"
#include "poahmm.h"





static int reestimate_param_poahmm(struct poahmm* poahmm);



static int reset_to_from_index(struct poahmm* poahmm);
static int reset_poa_graph_transitions_based_on_counts(struct poahmm* poahmm);

static int set_pseudocount(struct poahmm* poahmm);


static int check_and_extend_poahmm(struct poahmm* poahmm,int num_states, int new_maxlen);
//static struct poahmm_node* malloc_a_node(int maxseq_len,int num_samples );
static struct poahmm_node* malloc_a_node(int maxseq_len);
static void free_a_node(struct poahmm_node* node);
static struct  poahmm_boundary_node* malloc_a_boundary_node(int maxseq_len);
static void free_a_boundary_node(struct poahmm_boundary_node* node);

static int cmp_node_rank_low_to_high(const void * a, const void * b);
static int cmp_node_rank_high_to_low(const void * a, const void * b);
static struct poahmm* set_rank(struct poahmm* poahmm,int index,int rank);
static double gamma_fraction(unsigned int* seed , double a);



int forward_poahmm(struct poahmm* poahmm, uint8_t* seq, int len)
{
        int i,j,c,n;
        int node_id;

        struct cell* cells;
        struct cell* prev_cells;
        //struct cell* tmp;

        uint8_t node_nuc;

        seq = seq-1;

        //uint8_t* seq = poahmm_data->seq[index] -1; // so that seq[1]  is the first letter...
        //int len = poahmm_data->len[index];

        float MM = poahmm->MM;
        float MX = poahmm->MX;
        float MY = poahmm->MY;

        float XM = poahmm->XM;
        float XX = poahmm->XX;
        float XY = poahmm->XY;

        float YM = poahmm->YM;
        float YX = poahmm->YX;
        float YY = poahmm->YY;

        float* eY =poahmm->emission_Y;
        float* eX = poahmm->emission_X;

        float* eR = poahmm->background;
        float* entry = poahmm->entry_probabilities;

        float* exit = poahmm->exit_probabilities;

        float** eM = poahmm->emission_M;
        float tmp;

        float YY_boundary = prob2scaledprob((float)len / ((float)len+2.0));

        float YY_boundary_exit = prob2scaledprob(1.0 -(float)len / ((float)len+2.0) );

        RUN(check_and_extend_poahmm(poahmm,poahmm->num_nodes,len));

        qsort(poahmm->rank_sorted_nodes,poahmm->num_nodes,sizeof(struct poahmm_node*), cmp_node_rank_low_to_high);

        cells = poahmm->begin->cells;

        //poahmm->begin->fY[0] = prob2scaledprob(1.0f);
        //FIRST COLUMN
        cells[0].fY = prob2scaledprob(1.0f);
        for(i = 1; i < len;i++){
                cells[i].fY = cells[i-1].fY + YY_boundary + eR[seq[i]];
        }
        cells[len].fY = prob2scaledprob(0.0f);
        cells[len+1].fY = prob2scaledprob(0.0f);

        // start proper DP.
        for(j = 0;  j < poahmm->num_nodes;j++){
                node_id = poahmm->rank_sorted_nodes[j]->identifier;
                //LOG_MSG("Start with node: %d (rank: %d)",node_id,poahmm->rank_sorted_nodes[j]->rank);
                node_nuc = poahmm->nodes[node_id]->nuc;
                cells = poahmm->nodes[node_id]->cells;
                //LOG_MSG("NUC:%d",node_nuc);
                i = 0;
                cells[i].fM= prob2scaledprob(0.0);
                cells[i].fY = prob2scaledprob(0.0);
                cells[i].fX = poahmm->begin->cells[i].fY + MX +YY_boundary_exit+ entry[node_id];



                for(c = 1; c < poahmm->to_tindex[node_id][0];c++){
                        n =poahmm->to_tindex[node_id][c];
                        prev_cells = poahmm->nodes[n]->cells;

                        cells[i].fX = logsum(cells[i].fX ,prev_cells[i].fX + YY_boundary +  poahmm->poa_graph[n][node_id] );
                }
                //LOG_MSG("%f ", eR[node_nuc]);
                cells[i].fX += eR[node_nuc];

                for(i = 1; i < len;i++){
                        //index_dp = &cells[i].f;
                        prev_cells = poahmm->begin->cells;

                        tmp = YY_boundary_exit+entry[node_id];

                        cells[i].fM = prev_cells[i-1].fY + MM + tmp;//YY_boundary_exit+entry[node_id];
                        cells[i].fX =  prev_cells[i].fY + MX + tmp;//YY_boundary_exit+entry[node_id];
                        cells[i].fY =  prev_cells[i-1].fY + MY + tmp;//YY_boundary_exit+entry[node_id];

                        cells[i].fY = logsum(cells[i].fY,cells[i-1].fY+ YY);
                        cells[i].fY = logsum(cells[i].fY,cells[i-1].fM + MY);
                        cells[i].fY = logsum(cells[i].fY,cells[i-1].fX + XY);


                        for(c = 1; c < poahmm->to_tindex[node_id][0];c++){
                                n = poahmm->to_tindex[node_id][c];
                                prev_cells = poahmm->nodes[n]->cells;
                                tmp = poahmm->poa_graph[n][node_id];

                                cells[i].fM = logsum(cells[i].fM, prev_cells[i-1].fM + MM + tmp);// poahmm->poa_graph[n][node_id]);
                                cells[i].fM = logsum(cells[i].fM, prev_cells[i-1].fX  + XM + tmp);// poahmm->poa_graph[n][node_id]);
                                cells[i].fM = logsum(cells[i].fM, prev_cells[i-1].fY  + YM + tmp);// poahmm->poa_graph[n][node_id]);

                                cells[i].fX = logsum(cells[i].fX, prev_cells[i].fM + MX+ tmp);// poahmm->poa_graph[n][node_id]);
                                cells[i].fX = logsum(cells[i].fX, prev_cells[i].fX  + XX + tmp);//poahmm->poa_graph[n][node_id]  );
                                cells[i].fX = logsum(cells[i].fX, prev_cells[i].fY  + YX + tmp);//poahmm->poa_graph[n][node_id]  );
                        }

                        cells[i].fM += eM[node_nuc][seq[i]];
                        cells[i].fX +=  eX[node_nuc];
                        cells[i].fY += eY[seq[i]];
                }

                i = len;
                prev_cells = poahmm->begin->cells;
                tmp =YY_boundary_exit+entry[node_id];
                cells[i].fM = prev_cells[i-1].fY + MM + tmp;//YY_boundary_exit+entry[node_id];
                cells[i].fX =  prev_cells[i].fY + MX + tmp;//YY_boundary_exit+entry[node_id];
                cells[i].fY =  prev_cells[i-1].fY + MY + tmp;//YY_boundary_exit+entry[node_id];

                cells[i].fY = logsum(cells[i].fY,cells[i-1].fY+ YY);
                cells[i].fY = logsum(cells[i].fY,cells[i-1].fM + MY);
                cells[i].fY = logsum(cells[i].fY,cells[i-1].fX + XY);


                for(c = 1; c < poahmm->to_tindex[node_id][0];c++){
                        n = poahmm->to_tindex[node_id][c];
                        prev_cells = poahmm->nodes[n]->cells;
                        tmp = poahmm->poa_graph[n][node_id];
                        cells[i].fM = logsum(cells[i].fM, prev_cells[i-1].fM + MM + tmp);// poahmm->poa_graph[n][node_id]);
                        cells[i].fM = logsum(cells[i].fM, prev_cells[i-1].fX  + XM + tmp);// poahmm->poa_graph[n][node_id]);
                        cells[i].fM = logsum(cells[i].fM, prev_cells[i-1].fY  + YM + tmp);// poahmm->poa_graph[n][node_id]);

                        cells[i].fX = logsum(cells[i].fX, prev_cells[i].fM  + MX + tmp);// poahmm->poa_graph[n][node_id]);
                        cells[i].fX = logsum(cells[i].fX, prev_cells[i].fX + YY_boundary + tmp);// poahmm->poa_graph[n][node_id]  );
                        cells[i].fX = logsum(cells[i].fX, prev_cells[i].fY + YX + tmp);// poahmm->poa_graph[n][node_id]  );
                }
                cells[i].fM += eM[node_nuc][seq[i]];
                cells[i].fX +=  eR[node_nuc];
                cells[i].fY += eY[seq[i]];
        }

        cells = poahmm->end->cells;


        cells[0].fY = prob2scaledprob(0.0f);
        cells[1].fY = prob2scaledprob(0.0f);

        for(i = 2; i <= len;i++){
                cells[i].fY = prob2scaledprob(0.0f);
                for(j = 0;  j < poahmm->num_nodes;j++){
                        node_id = poahmm->rank_sorted_nodes[j]->identifier;
                        prev_cells =poahmm->nodes[node_id]->cells;
                        tmp =YY_boundary+ exit[node_id];
                        cells[i].fY = logsum(cells[i].fY, prev_cells[i-1].fM + MM + tmp);//YY_boundary+ exit[node_id];
                        cells[i].fY = logsum(cells[i].fY, prev_cells[i-1].fX + XM + tmp);//YY_boundary+ exit[node_id]);
                        cells[i].fY = logsum(cells[i].fY, prev_cells[i-1].fY + YM  + tmp);//YY_boundary+ exit[node_id]);
                }
                cells[i].fY = logsum(cells[i].fY, cells[i-1].fY + YY_boundary);
                cells[i].fY += eR[ seq[i]];
        }
        i = len+1;
        cells[i].fY = prob2scaledprob(0.0f);

        for(j = 0;  j < poahmm->num_nodes;j++){
                node_id = poahmm->rank_sorted_nodes[j]->identifier;
                prev_cells =poahmm->nodes[node_id]->cells;
                tmp = YY_boundary_exit+ exit[node_id];
                cells[i].fY = logsum(cells[i].fY, prev_cells[i-1].fM + MM + tmp);//YY_boundary_exit+ exit[node_id];
                cells[i].fY = logsum(cells[i].fY, prev_cells[i-1].fX + XM + tmp);//+ YY_boundary_exit+ exit[node_id]);
                cells[i].fY = logsum(cells[i].fY, prev_cells[i-1].fY + YM  + tmp);//+ YY_boundary_exit+ exit[node_id]);
        }

        cells[i].fY = logsum(cells[i].fY, cells[i-1].fY + YY_boundary_exit);
        poahmm->f_score = cells[i].fY;



        return OK;
ERROR:
        return FAIL;
}



int backward_poahmm(struct poahmm* poahmm, uint8_t* seq, int len)
{
        int i,j,c,n;
        int node_id;

        struct cell* cells;
        struct cell* next_cells;
        uint8_t node_nuc;

        seq = seq -1;
        //uint8_t* seq = poahmm_data->seq[index] -1; // so that seq[1]  is the first letter...
        //int len = poahmm_data->len[index];

        float tmp;

        float MM = poahmm->MM;
        float MX = poahmm->MX;
        float MY = poahmm->MY;

        float XM = poahmm->XM;
        float XX = poahmm->XX;
        float XY = poahmm->XY;

        float YM = poahmm->YM;
        float YX = poahmm->YX;
        float YY = poahmm->YY;

        float* eY = poahmm->emission_Y;
        float* eX = poahmm->emission_X;
        float* eR = poahmm->background;
        float* entry = poahmm->entry_probabilities;
        float* exit = poahmm->exit_probabilities;
        float** eM = poahmm->emission_M;

        float incoming_M;
        float incoming_X;

        float YY_boundary = prob2scaledprob((float)len / ((float)len+2.0));
        float YY_boundary_exit = prob2scaledprob(1.0 -(float)len / ((float)len+2.0) );

        qsort(poahmm->rank_sorted_nodes,poahmm->num_nodes,sizeof(struct poahmm_node*), cmp_node_rank_high_to_low);

        cells = poahmm->end->cells;
        i = len+1;
        cells[i].bY = prob2scaledprob(1.0);

        i = len;
        cells[i].bY = cells[i+1].bY +   YY_boundary_exit ;
        for(i = len-1;i >= 2; i--){ // the last two cells can never be used...
                cells[i].bY = cells[i+1].bY +YY_boundary + eR[seq[i+1]];

        }

        cells[1].bY = prob2scaledprob(0.0);

        cells[0].bY = prob2scaledprob(0.0);


        for(j = 0;  j < poahmm->num_nodes;j++){
                node_id = poahmm->rank_sorted_nodes[j]->identifier;
                cells = poahmm->nodes[node_id]->cells;
                next_cells = poahmm->end->cells;

                i = len+1;
                cells[i].bM = prob2scaledprob(0.0);
                cells[i].bX = prob2scaledprob(0.0);
                cells[i].bY = prob2scaledprob(0.0);

                i = len;

                tmp = next_cells[i+1].bY  +YY_boundary_exit + exit[node_id];
                cells[i].bM = tmp + MM;// +  YY_boundary_exit + exit[node_id];
                cells[i].bX = tmp + XM;// +  YY_boundary_exit + exit[node_id];
                cells[i].bY = tmp +YM;// +  YY_boundary_exit + exit[node_id];

                for(c = 1; c < poahmm->from_tindex[node_id][0];c++){
                        n = poahmm->from_tindex[node_id][c];
                        next_cells = poahmm->nodes[n]->cells;
                        node_nuc = poahmm->nodes[n]->nuc;
                        incoming_X =  next_cells[i].bX +  poahmm->poa_graph[node_id][n]  + eR[node_nuc];


                        cells[i].bM = logsum(cells[i].bM , MX  + incoming_X);
                        cells[i].bX = logsum(cells[i].bX ,YY_boundary + incoming_X);
                }

                for(i = len-1;i >= 1; i--){

                        next_cells = poahmm->end->cells;
                        tmp =next_cells[i+1].bY +YY_boundary + eR[seq[i+1]]+ exit[node_id];

                        cells[i].bM = tmp + MM;//next_cells[i+1].bY + MM + YY_boundary + eR[seq[i+1]]+ exit[node_id] ;
                        cells[i].bX = tmp + XM;//next_cells[i+1].bY + XM + YY_boundary + eR[seq[i+1]]+ exit[node_id] ;
                        cells[i].bY = tmp + YM;//next_cells[i+1].bY + YM + YY_boundary + eR[seq[i+1]]+ exit[node_id];
                        for(c = 1; c < poahmm->from_tindex[node_id][0];c++){
                                n = poahmm->from_tindex[node_id][c];
                                next_cells = poahmm->nodes[n]->cells;
                                node_nuc = poahmm->nodes[n]->nuc;
                                incoming_M = next_cells[i+1].bM +  poahmm->poa_graph[node_id][n] + eM[node_nuc][seq[i+1]] ;
                                incoming_X = next_cells[i].bX +  poahmm->poa_graph[node_id][n] + eX[node_nuc];

                                cells[i].bM = logsum(cells[i].bM, MM + incoming_M);
                                cells[i].bM = logsum(cells[i].bM, MX  + incoming_X);

                                cells[i].bX = logsum(cells[i].bX, XM + incoming_M );
                                cells[i].bX = logsum(cells[i].bX, XX + incoming_X );

                                cells[i].bY = logsum(cells[i].bY, YM + incoming_M);
                                cells[i].bY = logsum(cells[i].bY, YX + incoming_X);
                        }
                        tmp =cells[i+1].bY +eY[seq[i+1]];
                        cells[i].bM = logsum(cells[i].bM, tmp+ MY);//cells[i+1].bY + MY + eY[seq[i+1]]);
                        cells[i].bX = logsum(cells[i].bX , tmp+ XY);//cells[i+1].bY + XY + eY[seq[i+1]]);
                        cells[i].bY = logsum(cells[i].bY , tmp+ YY);//cells[i+1].bY + YY + eY[seq[i+1]]);
                }
                i = 0;
                next_cells = poahmm->end->cells;
                tmp =next_cells[i+1].bY +  YY_boundary + eR[seq[i+1]]+ exit[node_id];
                cells[i].bM = tmp + MM;//next_cells[i+1].bY +MM +  YY_boundary + eR[seq[i+1]]+ exit[node_id];
                cells[i].bX = tmp + XM;//next_cells[i+1].bY + XM +  YY_boundary+ eR[seq[i+1]]+ exit[node_id];
                cells[i].bY = tmp + YM;//next_cells[i+1].bY + YM +  YY_boundary + eR[seq[i+1]]+ exit[node_id];

                for(c = 1; c < poahmm->from_tindex[node_id][0];c++){
                        //	DPRINTF3("Looking at %d (from %d)",poahmm->to_tindex[node_id][c], node_id  );
                        n =poahmm->from_tindex[node_id][c];
                        //ancestor_ptr = poahmm->nodes[n];
                        next_cells = poahmm->nodes[n]->cells;
                        node_nuc = poahmm->nodes[n]->nuc;

                        incoming_M = next_cells[i+1].bM +  poahmm->poa_graph[node_id][n] + eM[node_nuc][seq[i+1]] ;
                        incoming_X = next_cells[i].bX + poahmm->poa_graph[node_id][n]  + eR[node_nuc];// eX[ancestor_ptr->nuc];

                        cells[i].bX = logsum(cells[i].bX , XM + incoming_M);
                        cells[i].bX = logsum(cells[i].bX , YY_boundary + incoming_X);
                }
                cells[i].bX = logsum(cells[i].bX ,  cells[i+1].bY + XY + eY[seq[i+1]]);
        }


        cells = poahmm->begin->cells;

        cells[len+1].bY =prob2scaledprob(0.0);
        cells[len].bY =prob2scaledprob(0.0);

        //poahmm->begin->bY[len+1] = prob2scaledprob(0.0);
        //poahmm->begin->bY[len] = prob2scaledprob(0.0);

        for(i = len-1;i > 0; i--){
                cells[i].bY = prob2scaledprob(0.0);
                for(j = 0;  j < poahmm->num_nodes;j++){
                        node_id = poahmm->rank_sorted_nodes[j]->identifier;
                        //node_ptr = poahmm->nodes[node_id];
                        next_cells = poahmm->nodes[node_id]->cells;
                        node_nuc =  poahmm->nodes[node_id]->nuc;
                        //if(poahmm->to_tindex[node_id][0] == 1){
                        tmp =YY_boundary_exit+entry[node_id];
                        cells[i].bY = logsum(cells[i].bY, next_cells[i+1].bM + eM[node_nuc][seq[i+1]] +  MM + tmp);
                        cells[i].bY = logsum(cells[i].bY, next_cells[i].bX + eX[node_nuc] +  MX + tmp);
                        cells[i].bY = logsum(cells[i].bY, next_cells[i+1].bY + eY[seq[i+1]] +  MY + tmp);
                }
                cells[i].bY = logsum(cells[i].bY,cells[i+1].bY+ YY_boundary + eR[seq[i+1]]);
        }
        i = 0;
        cells[i].bY = prob2scaledprob(0.0);
        for(j = 0;  j < poahmm->num_nodes;j++){
                node_id = poahmm->rank_sorted_nodes[j]->identifier;
                next_cells = poahmm->nodes[node_id]->cells;
                node_nuc =  poahmm->nodes[node_id]->nuc;
                tmp =YY_boundary_exit+entry[node_id] ;
                cells[i].bY = logsum(cells[i].bY, next_cells[i+1].bM + eM[node_nuc][seq[i+1]] +  MM +tmp );
                cells[i].bY = logsum(cells[i].bY, next_cells[i].bX + eR[node_nuc] +  MX + tmp);
                cells[i].bY = logsum(cells[i].bY, next_cells[i+1].bY + eY[seq[i+1]] +  MY + tmp);
        }
        cells[i].bY = logsum(cells[i].bY,cells[i+1].bY+ YY_boundary + eR[seq[i+1]]);
        poahmm->b_score = cells[i].bY;//poahmm->begin->bY[0];

        return OK;
}

#define TO_M 0u
#define TO_X 1u
#define TO_Y 2u
#define TO_B 3u
#define TO_E 4u
#define START_STATE 0x1FFFFFFFu
#define END_STATE 0x1FFFFFFEu


int viterbi_poahmm(struct poahmm* poahmm, uint8_t* seq, int len, int* path)
{
        int i,j,c,n;
        int node_id;

        struct cell* cells;
        struct cell* prev_cells;
        //struct cell* tmp;

        uint8_t node_nuc;

        seq = seq -1;


        //uint8_t* seq = poahmm_data->seq[index] -1; // so that seq[1]  is the first letter...
        //int len = poahmm_data->len[index];

        float MM = poahmm->MM;
        float MX = poahmm->MX;
        float MY = poahmm->MY;

        float XM = poahmm->XM;
        float XX = poahmm->XX;
        float XY = poahmm->XY;

        float YM = poahmm->YM;
        float YX = poahmm->YX;
        float YY = poahmm->YY;

        float* eY =poahmm->emission_Y;
        float* eX = poahmm->emission_X;

        float* eR = poahmm->background;
        float* entry = poahmm->entry_probabilities;

        float* exit = poahmm->exit_probabilities;

        float** eM = poahmm->emission_M;
        float tmp;
        //float max;
        float new_max;

        float YY_boundary = prob2scaledprob((float)len / ((float)len+2.0));

        float YY_boundary_exit = prob2scaledprob(1.0 -(float)len / ((float)len+2.0) );



        RUN(check_and_extend_poahmm(poahmm,poahmm->num_nodes,len));

        qsort(poahmm->rank_sorted_nodes,poahmm->num_nodes,sizeof(struct poahmm_node*), cmp_node_rank_low_to_high);

        cells = poahmm->begin->cells;

        //poahmm->begin->fY[0] = prob2scaledprob(1.0f);
        //FIRST COLUMN
        cells[0].fY = prob2scaledprob(1.0f);
        for(i = 1; i < len;i++){
                cells[i].fY = cells[i-1].fY + YY_boundary + eR[seq[i]];
                cells[i].Y_to_state = START_STATE;
                cells[i].Y_trans = TO_B;
        }
        cells[len].fY = prob2scaledprob(0.0f);
        cells[len+1].fY = prob2scaledprob(0.0f);

        // start proper DP.
        for(j = 0;  j < poahmm->num_nodes;j++){
                node_id = poahmm->rank_sorted_nodes[j]->identifier;
                //DPRINTF3("Start with node: %d (rank: %d)",node_id,poahmm->rank_sorted_nodes[j]->rank);
                node_nuc = poahmm->nodes[node_id]->nuc;
                cells = poahmm->nodes[node_id]->cells;

                i = 0;

                cells[i].fM= prob2scaledprob(0.0);
                cells[i].fY = prob2scaledprob(0.0);
                cells[i].fX =poahmm->begin->cells[i].fY + MX +YY_boundary_exit+ entry[node_id];
                cells[i].X_to_state = START_STATE;
                cells[i].X_trans = TO_B;
                //cells[i]


                for(c = 1; c < poahmm->to_tindex[node_id][0];c++){
                        n =poahmm->to_tindex[node_id][c];
                        prev_cells = poahmm->nodes[n]->cells;
                        if(cells[i].fX < (new_max = prev_cells[i].fX + YY_boundary +  poahmm->poa_graph[n][node_id]) ){
                                //max = new_max;
                                cells[i].fX = new_max;
                                cells[i].X_to_state = n;
                                cells[i].X_trans = TO_X;
                        }
                }

                cells[i].fX += eR[node_nuc];

                for(i = 1; i < len;i++){
                        //index_dp = &cells[i].f;
                        prev_cells = poahmm->begin->cells;

                        tmp = YY_boundary_exit+entry[node_id];

                        cells[i].fM = prev_cells[i-1].fY + MM + tmp;//YY_boundary_exit+entry[node_id];
                        cells[i].M_to_state = START_STATE;
                        cells[i].M_trans = TO_B;
                        cells[i].fX =  prev_cells[i].fY + MX + tmp;//YY_boundary_exit+entry[node_id];
                        cells[i].X_to_state = START_STATE;
                        cells[i].X_trans = TO_B;
                        cells[i].fY =  prev_cells[i-1].fY + MY + tmp;//YY_boundary_exit+entry[node_id];
                        cells[i].Y_to_state = START_STATE;
                        cells[i].Y_trans = TO_B;

                        if((new_max = cells[i-1].fY+ YY) > cells[i].fY){
                                cells[i].fY = new_max;
                                cells[i].Y_to_state = node_id;
                                cells[i].Y_trans = TO_Y;

                        }
                        //cells[i].fY = logsum(cells[i].fY,cells[i-1].fY+ YY);
                        if((new_max = cells[i-1].fM+ MY) > cells[i].fY){
                                cells[i].fY = new_max;
                                cells[i].Y_to_state = node_id;
                                cells[i].Y_trans = TO_M;

                        }
                        //cells[i].fY = logsum(cells[i].fY,cells[i-1].fM + MY);

                        if((new_max = cells[i-1].fX+ XY) > cells[i].fY){
                                cells[i].fY = new_max;
                                cells[i].Y_to_state = node_id;
                                cells[i].Y_trans = TO_X;

                        }
                        //cells[i].fY = logsum(cells[i].fY,cells[i-1].fX + XY);


                        for(c = 1; c < poahmm->to_tindex[node_id][0];c++){

                                n = poahmm->to_tindex[node_id][c];
                                prev_cells = poahmm->nodes[n]->cells;
                                tmp = poahmm->poa_graph[n][node_id];


                                if((new_max = prev_cells[i-1].fM + MM + tmp) > cells[i].fM){
                                        cells[i].fM = new_max;
                                        cells[i].M_to_state = n;
                                        cells[i].M_trans = TO_M;
                                }
                                //cells[i].fM = logsum(cells[i].fM, prev_cells[i-1].fM + MM + tmp);// poahmm->poa_graph[n][node_id]);

                                if((new_max = prev_cells[i-1].fX  + XM + tmp) > cells[i].fM){
                                        cells[i].fM = new_max;
                                        cells[i].M_to_state = n;
                                        cells[i].M_trans = TO_X;
                                }

                                //cells[i].fM = logsum(cells[i].fM, prev_cells[i-1].fX  + XM + tmp);// poahmm->poa_graph[n][node_id]);

                                if((new_max = prev_cells[i-1].fY  + YM + tmp) > cells[i].fM){
                                        cells[i].fM = new_max;
                                        cells[i].M_to_state = n;
                                        cells[i].M_trans = TO_Y;
                                }
                                //cells[i].fM = logsum(cells[i].fM, prev_cells[i-1].fY  + YM + tmp);// poahmm->poa_graph[n][node_id]);


                                if((new_max  = prev_cells[i].fM + MX+ tmp)  > cells[i].fX){
                                        cells[i].fX = new_max;
                                        cells[i].X_to_state = n;
                                        cells[i].X_trans = TO_M;

                                }

                                //cells[i].fX = logsum(cells[i].fX, prev_cells[i].fM + MX+ tmp);// poahmm->poa_graph[n][node_id]);
                                if((new_max  = prev_cells[i].fX  + XX + tmp)  > cells[i].fX){
                                        cells[i].fX = new_max;
                                        cells[i].X_to_state = n;
                                        cells[i].X_trans = TO_X;

                                }
                                //cells[i].fX = logsum(cells[i].fX, prev_cells[i].fX  + XX + tmp);//poahmm->poa_graph[n][node_id]  );

                                if((new_max  = prev_cells[i].fY  + YX + tmp)  > cells[i].fX){
                                        cells[i].fX = new_max;
                                        cells[i].X_to_state = n;
                                        cells[i].X_trans = TO_Y;

                                }
                                //cells[i].fX = logsum(cells[i].fX, prev_cells[i].fY  + YX + tmp);//poahmm->poa_graph[n][node_id]  );
                        }

                        cells[i].fM += eM[node_nuc][seq[i]];
                        cells[i].fX +=  eX[node_nuc];
                        cells[i].fY += eY[seq[i]];
                }

                i = len;
                prev_cells = poahmm->begin->cells;
                tmp =YY_boundary_exit+entry[node_id];
                cells[i].fM = prev_cells[i-1].fY + MM + tmp;//YY_boundary_exit+entry[node_id];
                cells[i].M_to_state = START_STATE;
                cells[i].M_trans = TO_B;

                cells[i].fX =  prev_cells[i].fY + MX + tmp;//YY_boundary_exit+entry[node_id];
                cells[i].X_to_state = START_STATE;
                cells[i].X_trans = TO_B;

                cells[i].fY =  prev_cells[i-1].fY + MY + tmp;//YY_boundary_exit+entry[node_id];
                cells[i].Y_to_state = START_STATE;
                cells[i].Y_trans = TO_B;


                if((new_max =cells[i-1].fY+ YY) > cells[i].fY ){
                        cells[i].fY = new_max;
                        cells[i].Y_to_state = node_id;
                        cells[i].Y_trans = TO_Y;


                }
                //cells[i].fY = logsum(cells[i].fY,cells[i-1].fY+ YY);

                if((new_max = cells[i-1].fM + MY) > cells[i].fY ){
                        cells[i].fY = new_max;
                        cells[i].Y_to_state = node_id;
                        cells[i].Y_trans = TO_M;


                }
                //cells[i].fY = logsum(cells[i].fY,cells[i-1].fM + MY);

                if((new_max = cells[i-1].fX + XY) > cells[i].fY ){
                        cells[i].fY = new_max;
                        cells[i].Y_to_state = node_id;
                        cells[i].Y_trans = TO_X;


                }
                //cells[i].fY = logsum(cells[i].fY,cells[i-1].fX + XY);


                for(c = 1; c < poahmm->to_tindex[node_id][0];c++){
                        n = poahmm->to_tindex[node_id][c];
                        prev_cells = poahmm->nodes[n]->cells;
                        tmp = poahmm->poa_graph[n][node_id];


                        if((new_max =  prev_cells[i-1].fM + MM + tmp) > cells[i].fM){
                                cells[i].fM  = new_max;
                                cells[i].M_to_state = n;
                                cells[i].M_trans = TO_M;
                        }
                        //cells[i].fM = logsum(cells[i].fM, prev_cells[i-1].fM + MM + tmp);// poahmm->poa_graph[n][node_id]);
                        if((new_max = prev_cells[i-1].fX  + XM + tmp) > cells[i].fM){
                                cells[i].fM = new_max;
                                cells[i].M_to_state = n;
                                cells[i].M_trans = TO_X;
                        }
                        //cells[i].fM = logsum(cells[i].fM, prev_cells[i-1].fX  + XM + tmp);// poahmm->poa_graph[n][node_id]);
                        if((new_max = prev_cells[i-1].fY  + YM + tmp) > cells[i].fM){
                                cells[i].fM = new_max;
                                cells[i].M_to_state = n;
                                cells[i].M_trans = TO_Y;
                        }
                        //cells[i].fM = logsum(cells[i].fM, prev_cells[i-1].fY  + YM + tmp);// poahmm->poa_graph[n][node_id]);
                        if((new_max  = prev_cells[i].fM + MX+ tmp)  > cells[i].fX){
                                cells[i].fX = new_max;
                                cells[i].X_to_state = n;
                                cells[i].X_trans = TO_M;

                        }

                        //cells[i].fX = logsum(cells[i].fX, prev_cells[i].fM  + MX + tmp);// poahmm->poa_graph[n][node_id]);
                        if((new_max  = prev_cells[i].fX  + YY_boundary + tmp)  > cells[i].fX){
                                cells[i].fX = new_max;
                                cells[i].X_to_state = n;
                                cells[i].X_trans = TO_X;

                        }
                        //cells[i].fX = logsum(cells[i].fX, prev_cells[i].fX + YY_boundary + tmp);// poahmm->poa_graph[n][node_id]  );
                        if((new_max  = prev_cells[i].fY  + YX + tmp)  > cells[i].fX){
                                cells[i].fX = new_max;
                                cells[i].X_to_state = n;
                                cells[i].X_trans = TO_Y;

                        }
                        //cells[i].fX = logsum(cells[i].fX, prev_cells[i].fY + YX + tmp);// poahmm->poa_graph[n][node_id]  );
                }
                cells[i].fM += eM[node_nuc][seq[i]];
                cells[i].fX +=  eR[node_nuc];
                cells[i].fY += eY[seq[i]];
        }

        cells = poahmm->end->cells;


        cells[0].fY = prob2scaledprob(0.0f);
        cells[1].fY = prob2scaledprob(0.0f);

        for(i = 2; i <= len;i++){
                //cells[i].fY = prob2scaledprob(0.0f);

                cells[i].fY = cells[i-1].fY + YY_boundary;
                cells[i].Y_to_state = END_STATE;
                cells[i].Y_trans = TO_E;
                //}

                for(j = 0;  j < poahmm->num_nodes;j++){
                        node_id = poahmm->rank_sorted_nodes[j]->identifier;
                        prev_cells =poahmm->nodes[node_id]->cells;
                        tmp =YY_boundary+ exit[node_id];
                        if((new_max = prev_cells[i-1].fM + MM + tmp) > cells[i].fY){
                                cells[i].fY = new_max;
                                cells[i].Y_to_state = node_id;
                                cells[i].Y_trans = TO_M;
                                //			DPRINTF3("Picking M:%d	%d\n",node_id,i);
                        }
                        //cells[i].fY = logsum(cells[i].fY, prev_cells[i-1].fM + MM + tmp);//YY_boundary+ exit[node_id];

                        if((new_max = prev_cells[i-1].fX + XM + tmp) > cells[i].fY){
                                cells[i].fY = new_max;
                                cells[i].Y_to_state = node_id;
                                cells[i].Y_trans = TO_X;
                                //			DPRINTF3("Picking X:%d	%d\n",node_id,i);
                        }
                        //cells[i].fY = logsum(cells[i].fY, prev_cells[i-1].fX + XM + tmp);//YY_boundary+ exit[node_id]);
                        if((new_max = prev_cells[i-1].fY + YM  + tmp) > cells[i].fY){
                                cells[i].fY = new_max;
                                cells[i].Y_to_state = node_id;
                                cells[i].Y_trans = TO_Y;
                                //			DPRINTF3("Picking Y:%d	%d\n",node_id,i);
                        }
                        //cells[i].fY = logsum(cells[i].fY, prev_cells[i-1].fY + YM  + tmp);//YY_boundary+ exit[node_id]);
                }
                //cells[i].fY = logsum(cells[i].fY, cells[i-1].fY + YY_boundary);
                cells[i].fY += eR[ seq[i]];
        }
        i = len+1;
        //cells[i].fY = prob2scaledprob(0.0f);
        cells[i].fY = cells[i-1].fY + YY_boundary_exit;
        cells[i].Y_to_state = END_STATE;
        cells[i].Y_trans = TO_E;

        for(j = 0;  j < poahmm->num_nodes;j++){
                node_id = poahmm->rank_sorted_nodes[j]->identifier;
                prev_cells =poahmm->nodes[node_id]->cells;
                tmp = YY_boundary_exit+ exit[node_id];
                if((new_max = prev_cells[i-1].fM + MM + tmp) > cells[i].fY){
                        cells[i].fY = new_max;
                        cells[i].Y_to_state = node_id;
                        cells[i].Y_trans = TO_M;
                        //		DPRINTF3("Picking M:%d	%d set: %d\n",node_id,i,(int) cells[i].Y_to_state);
                }
                //cells[i].fY = logsum(cells[i].fY, prev_cells[i-1].fM + MM + tmp);//YY_boundary_exit+ exit[node_id];
                if((new_max = prev_cells[i-1].fX + XM + tmp) > cells[i].fY){
                        cells[i].fY = new_max;
                        cells[i].Y_to_state = node_id;
                        cells[i].Y_trans = TO_X;
                        //		DPRINTF3("Picking X:%d	%d\n",node_id,i);
                }
                //cells[i].fY = logsum(cells[i].fY, prev_cells[i-1].fX + XM + tmp);//+ YY_boundary_exit+ exit[node_id]);
                if((new_max = prev_cells[i-1].fY + YM  + tmp) > cells[i].fY){
                        cells[i].fY = new_max;
                        cells[i].Y_to_state = node_id;
                        cells[i].Y_trans = TO_Y;
                        //		DPRINTF3("Picking Y:%d	%d\n",node_id,i);
                }
                //cells[i].fY = logsum(cells[i].fY, prev_cells[i-1].fY + YM  + tmp);//+ YY_boundary_exit+ exit[node_id]);
        }

        //cells[i].fY = logsum(cells[i].fY, cells[i-1].fY + YY_boundary_exit);
        poahmm->f_score = cells[i].fY;

        //DPRINTF3("NOW to: %d move: %d   \n",poahmm->end->cells[len+1].Y_to_state, poahmm->end->cells[len+1].Y_trans);

        //print_viterbi_matrix(poahmm,  len);


        //Traceback!


        unsigned int next_state;
        unsigned int in_node;
        unsigned int state;
        //unsigned int mode;
        i = len+1;
        c = 1;
        //mode = poahmm->end->cells[i].Y_trans;
        in_node =poahmm->end->cells[i].Y_to_state;

        state = poahmm->end->cells[i].Y_trans;


        i = len;
        while(1){
                next_state = -1;
                if(state == TO_M){
                        path[c] = ((int) seq[i] << 28) | in_node;
                        //		DPRINTF3("%d - %d	%d",seq[i], poahmm->nodes[in_node]->nuc, in_node  );
                        c++;

                        next_state = poahmm->nodes[in_node]->cells[i].M_trans;
                        in_node = poahmm->nodes[in_node]->cells[i].M_to_state;
                        i--;
                }

                if(state == TO_B){
                        if(!i){
                                break;
                        }
                        path[c] = ((int) seq[i] << 28) | 0xFFFFFFF;
                        //		DPRINTF3("%d - -",seq[i]  );
                        c++;

                        next_state = poahmm->begin->cells[i].Y_trans;
                        in_node = poahmm->begin->cells[i].Y_to_state;
                        i--;
                }
                if(state == TO_Y){
                        path[c] = ((int) seq[i] << 28) | 0xFFFFFFF;
                        //		DPRINTF3("%d - -",seq[i]  );
                        c++;

                        //i--;
                        next_state = poahmm->nodes[in_node]->cells[i].Y_trans;
                        in_node = poahmm->nodes[in_node]->cells[i].Y_to_state;
                        i--;
                }

                if(state == TO_X){
                        path[c] = (4 << 28) | in_node;
                        //path[c] = ((int) seq[i] << 28) | 0xFFFFFFF;
                        //		DPRINTF3("- - %d",poahmm->nodes[in_node]->nuc );
                        c++;
                        //i--;
                        next_state = poahmm->nodes[in_node]->cells[i].X_trans;
                        in_node = poahmm->nodes[in_node]->cells[i].X_to_state;
                }
                if(state == TO_E){
                        path[c] = ((int) seq[i] << 28) | 0xFFFFFFF;
                        //		DPRINTF3("%d - -",seq[i]  );
                        c++;
                        //i--;
                        next_state = poahmm->end->cells[i].Y_trans;
                        in_node = poahmm->end->cells[i].Y_to_state;
                        i--;
                }


                state = next_state;
        }

        len = c;

        path[0] = c;

        for (i = 1, j = len - 1; i < j; i++, j--){
                c = path[i];
                path[i] = path[j];
                path[j] = c;
        }

        return OK;
ERROR:
        return FAIL;
}


#undef END_STATE
#undef START_STATE
#undef TO_E
#undef TO_B
#undef TO_Y
#undef TO_X
#undef TO_M





struct poahmm*  init_poahmm(int max_len, int* nuc_counts,float weight)
{
        struct poahmm* poahmm = NULL;
        int i;
        float sum;
        MMALLOC(poahmm, sizeof(struct poahmm));

        int maxmodel_len = max_len;

        poahmm->maxseq_len = 128;
        poahmm->seed =  (unsigned int) (time(NULL) * (42));

        poahmm->pseudo_weight = weight;

        poahmm->emission_M = NULL;
        poahmm->emission_X = NULL;
        poahmm->emission_Y = NULL;

        poahmm->e_emission_M = NULL;
        poahmm->e_emission_X = NULL;
        poahmm->e_emission_Y = NULL;

        poahmm->poa_graph = NULL;
        poahmm->e_poa_graph = NULL;

        RUN(galloc(&poahmm->emission_M,5,5));
        //poahmm->emission_M = malloc_2d_float(poahmm->emission_M, 5, 5,  prob2scaledprob(0.0));

        MMALLOC(poahmm->emission_X, sizeof(float) * 5);
        MMALLOC(poahmm->emission_Y, sizeof(float) * 5);

        RUN(galloc(&poahmm->e_emission_M,5,5));
        //poahmm->e_emission_M = malloc_2d_float(poahmm->e_emission_M, 5, 5,  prob2scaledprob(0.0));
        MMALLOC(poahmm->e_emission_X, sizeof(float) * 5);
        MMALLOC(poahmm->e_emission_Y, sizeof(float) * 5);

        poahmm->background = NULL;
        MMALLOC(poahmm->background, sizeof(float) * 5);

        if(nuc_counts){
                sum = 0.0;
                for(i = 0; i < 4;i++){
                        sum += nuc_counts[i];
                }

                for(i = 0; i < 4;i++){
                        poahmm->background[i] = prob2scaledprob( (float) nuc_counts[i] / (float) sum);
                }
                //poahmm->background[4] = prob2scaledprob(1.0);
        }

        RUN(set_pseudocount(poahmm));

        RUN(reestimate_param_poahmm(poahmm));

        poahmm->begin = NULL;
        poahmm->end = NULL;

        poahmm->nodes = NULL;
        poahmm->rank_sorted_nodes = NULL;

        poahmm->to_tindex = NULL;
        poahmm->from_tindex = NULL;
        poahmm->entry_probabilities = NULL;
        poahmm->e_entry_probabilities = NULL;

        poahmm->exit_probabilities = NULL;
        poahmm->e_exit_probabilities = NULL;


        poahmm->num_nodes = maxmodel_len;
        poahmm->alloced_num_nodes = maxmodel_len;

        MMALLOC(poahmm->entry_probabilities, sizeof(float)* poahmm->alloced_num_nodes);
        MMALLOC(poahmm->e_entry_probabilities, sizeof(float)* poahmm->alloced_num_nodes);

        MMALLOC(poahmm->exit_probabilities, sizeof(float)* poahmm->alloced_num_nodes);
        MMALLOC(poahmm->e_exit_probabilities, sizeof(float)* poahmm->alloced_num_nodes);




        RUN(galloc(&poahmm->poa_graph, poahmm->alloced_num_nodes, poahmm->alloced_num_nodes));
        RUN(galloc(&poahmm->e_poa_graph, poahmm->alloced_num_nodes, poahmm->alloced_num_nodes));
        RUN(galloc(&poahmm->to_tindex, poahmm->alloced_num_nodes, poahmm->alloced_num_nodes+1));
        RUN(galloc(&poahmm->from_tindex , poahmm->alloced_num_nodes, poahmm->alloced_num_nodes+1));

        MMALLOC(poahmm->nodes, sizeof(struct poahmm_node*) * poahmm->alloced_num_nodes);
        MMALLOC(poahmm->rank_sorted_nodes, sizeof(struct poahmm_node*) * poahmm->alloced_num_nodes);

        RUNP(poahmm->begin  = malloc_a_boundary_node(poahmm->maxseq_len));
        RUNP(poahmm->end  = malloc_a_boundary_node(poahmm->maxseq_len));

        for(i = 0; i < poahmm->alloced_num_nodes;i++){
                RUNP(poahmm->nodes[i] = malloc_a_node(poahmm->maxseq_len));
                poahmm->rank_sorted_nodes[i] =poahmm->nodes[i];
        }
        return poahmm;

ERROR:
        free_poahmm(poahmm);
        return NULL;
}

int reestimate_param_poahmm(struct poahmm* poahmm)
{
        int i,j;
        float sum;



        sum = prob2scaledprob(0.0);

        for(i = 0;i < 4;i++){

                for(j = 0;j < 4;j++){
                        sum = logsum(sum,poahmm->e_emission_M[i][j] );
                }
        }
        for(i = 0;i < 4;i++){
                sum = prob2scaledprob(0.0);
                for(j = 0;j < 4;j++){
                        sum = logsum(sum,poahmm->e_emission_M[i][j] );
                }
                for(j = 0;j < 4;j++){

                        poahmm->emission_M[i][j]  = poahmm->e_emission_M[i][j] - sum;
                        //sum = logsum(sum,poahmm->e_emission_M[i][j] );
                }
        }


        sum = prob2scaledprob(0.0);

        for(j = 0;j < 4;j++){
                sum = logsum(sum,poahmm->e_emission_X[j]);
        }
        for(j = 0;j < 4;j++){
                poahmm->emission_X[j] =    poahmm->e_emission_X[j]-sum;
        }

        sum = prob2scaledprob(0.0);

        for(j = 0;j < 4;j++){
                sum = logsum(sum,poahmm->e_emission_Y[j]);
        }
        for(j = 0;j < 4;j++){
                poahmm->emission_Y[j] =   poahmm->e_emission_Y[j] - sum;
        }


        sum = prob2scaledprob(0.0);
        sum = logsum(sum, poahmm->e_MM);
        sum = logsum(sum, poahmm->e_MX);
        sum = logsum(sum, poahmm->e_MY);

        //fprintf(stdout,"Resetting M-> %f %f %f (sum: %f)\n",poahmm->e_MM,poahmm->e_MX,poahmm->e_MY,sum);

        poahmm->MM = poahmm->e_MM - sum;
        poahmm->MX = poahmm->e_MX - sum;
        poahmm->MY = poahmm->e_MY - sum;


        sum = prob2scaledprob(0.0);
        sum = logsum(sum, poahmm->e_XM);
        sum = logsum(sum, poahmm->e_XX);
        sum = logsum(sum, poahmm->e_XY);


        poahmm->XM = poahmm->e_XM - sum;
        poahmm->XX = poahmm->e_XX - sum;
        poahmm->XY = poahmm->e_XY - sum;

        sum = prob2scaledprob(0.0);
        sum = logsum(sum, poahmm->e_YM);
        sum = logsum(sum, poahmm->e_YX);
        sum = logsum(sum, poahmm->e_YY);


        poahmm->YM = poahmm->e_YM - sum;
        poahmm->YY = poahmm->e_YY - sum;
        poahmm->YX = poahmm->e_YX - sum;

        RUN(set_pseudocount(poahmm));
        return OK;
ERROR:
        return FAIL;
}

int set_pseudocount(struct poahmm* poahmm)
{
        int i,j;

        float sum;
        sum = 0.0;

        // MM emission probabilities;
        for(i = 0;i < 4;i++){
                for(j = 0;j < 4;j++){
                        poahmm->e_emission_M[i][j] = gamma_fraction(&poahmm->seed, 0.99)+ 1.0;

                        if(i ==j){
                                poahmm->e_emission_M[i][j] += 100.0;//prob2scaledprob(0.95f);
                        }
                        sum +=poahmm->e_emission_M[i][j] ;
                }
                poahmm->e_emission_M[i][4] =prob2scaledprob(0.0f); // buffer for start and stop
                poahmm->e_emission_M[4][i] =prob2scaledprob(0.0f);// buffer for start and stop
        }

        for(i = 0;i < 4;i++){
                for(j = 0;j < 4;j++){
                        poahmm->e_emission_M[i][j] = prob2scaledprob(poahmm->e_emission_M[i][j]  / sum)  + prob2scaledprob(poahmm->pseudo_weight);
                }
        }
        poahmm->e_emission_M[4][4] = prob2scaledprob(0.0f);// buffer for start and stop

        for(i = 0;i < 4;i++){
                poahmm->e_emission_X[i] = gamma_fraction(&poahmm->seed, 0.99) + 100.0;
                poahmm->e_emission_Y[i] = gamma_fraction(&poahmm->seed, 0.99) + 100.0;
        }

        sum = 0.0;

        for(j = 0;j < 4;j++){
                sum += poahmm->e_emission_X[j];
        }
        for(j = 0;j < 4;j++){
                poahmm->e_emission_X[j] =prob2scaledprob( poahmm->e_emission_X[j]/ sum)+ prob2scaledprob(poahmm->pseudo_weight);
        }

        sum = 0.0;

        for(j = 0;j < 4;j++){
                sum += poahmm->e_emission_Y[j];
        }
        for(j = 0;j < 4;j++){
                poahmm->e_emission_Y[j] =  prob2scaledprob( poahmm->e_emission_Y[j]/ sum)+ prob2scaledprob(poahmm->pseudo_weight);
        }


        poahmm->emission_X[4] = prob2scaledprob(1.0f);// buffer for start and stop
        poahmm->emission_Y[4] = prob2scaledprob(1.0f);// buffer for start and stop


        poahmm->e_MM = prob2scaledprob(0.98)+ prob2scaledprob(poahmm->pseudo_weight);
        poahmm->e_MX = prob2scaledprob(0.01)+ prob2scaledprob(poahmm->pseudo_weight);
        poahmm->e_MY = prob2scaledprob(0.01)+ prob2scaledprob(poahmm->pseudo_weight);

        //poahmm->e_MM = prob2scaledprob( (float) 50.0 / (float)(50.0 +1))+ prob2scaledprob(poahmm->pseudo_weight);
        //poahmm->e_MX = prob2scaledprob((1.0- (float) 50.0 / (float)(50.0 +1))  / 2.0)+ prob2scaledprob(poahmm->pseudo_weight);
        //poahmm->e_MY = prob2scaledprob((1.0-  (float) 50.0 / (float)(50.0 +1)) / 2.0)+ prob2scaledprob(poahmm->pseudo_weight);


        poahmm->e_XM = prob2scaledprob(0.99f)+ prob2scaledprob(poahmm->pseudo_weight);
        poahmm->e_XX = prob2scaledprob(0.01f)+ prob2scaledprob(poahmm->pseudo_weight);
        poahmm->e_XY = prob2scaledprob(0.0f)+ prob2scaledprob(poahmm->pseudo_weight);

        poahmm->e_YM = prob2scaledprob(0.99f)+ prob2scaledprob(poahmm->pseudo_weight);
        poahmm->e_YX = prob2scaledprob(0.0f)+ prob2scaledprob(poahmm->pseudo_weight);
        poahmm->e_YY = prob2scaledprob(0.01f)+ prob2scaledprob(poahmm->pseudo_weight);



        return OK;
}

int check_and_extend_poahmm(struct poahmm* poahmm,int num_states, int new_maxlen)
{
        int i;
        int old = 0;

        if(new_maxlen+2 > poahmm->maxseq_len){
                poahmm->maxseq_len = new_maxlen+2;
                for(i = 0; i < poahmm->num_nodes;i++){
                        /*MREALLOC(poahmm->nodes[i]->fM, sizeof(float) * poahmm->maxseq_len );
                        MREALLOC(poahmm->nodes[i]->fX, sizeof(float) * poahmm->maxseq_len );
                        MREALLOC(poahmm->nodes[i]->fY, sizeof(float) * poahmm->maxseq_len );

                        MREALLOC(poahmm->nodes[i]->bM, sizeof(float) * poahmm->maxseq_len );
                        MREALLOC(poahmm->nodes[i]->bX, sizeof(float) * poahmm->maxseq_len );
                        MREALLOC(poahmm->nodes[i]->bY, sizeof(float) * poahmm->maxseq_len );
                        */
                        MREALLOC(poahmm->nodes[i]->cells, sizeof(struct cell) * poahmm->maxseq_len);

                }
                MREALLOC(poahmm->begin->fY , sizeof(float) * poahmm->maxseq_len );
                MREALLOC(poahmm->begin->bY , sizeof(float) * poahmm->maxseq_len );

                MREALLOC(poahmm->begin->cells ,sizeof(struct cell)* poahmm->maxseq_len );

        }

        if(num_states > poahmm->alloced_num_nodes){
                old =poahmm->alloced_num_nodes;
                poahmm->alloced_num_nodes = num_states;
                MREALLOC(poahmm->nodes, sizeof(struct poahmm_node*) *poahmm->alloced_num_nodes );
                MREALLOC(poahmm->rank_sorted_nodes, sizeof(struct poahmm_node*) *poahmm->alloced_num_nodes );

                //DPRINTF3("%d -> %d", poahmm->alloced_num_nodes, old);
                MREALLOC(poahmm->entry_probabilities, sizeof(float)* poahmm->alloced_num_nodes);
                MREALLOC(poahmm->e_entry_probabilities, sizeof(float)* poahmm->alloced_num_nodes);

                MREALLOC(poahmm->exit_probabilities, sizeof(float)* poahmm->alloced_num_nodes);
                MREALLOC(poahmm->e_exit_probabilities, sizeof(float)* poahmm->alloced_num_nodes);

                RUN(galloc(&poahmm->poa_graph, poahmm->alloced_num_nodes, poahmm->alloced_num_nodes));
                //poahmm->poa_graph = malloc_2d_float(poahmm->poa_graph, poahmm->alloced_num_nodes, poahmm->alloced_num_nodes, prob2scaledprob(0.0));
                RUN(galloc(&poahmm->e_poa_graph, poahmm->alloced_num_nodes, poahmm->alloced_num_nodes));
                //poahmm->e_poa_graph = malloc_2d_float(poahmm->e_poa_graph, poahmm->alloced_num_nodes, poahmm->alloced_num_nodes, 0.0);
                RUN(galloc(&poahmm->to_tindex, poahmm->alloced_num_nodes, poahmm->alloced_num_nodes+1));
                //poahmm->to_tindex = malloc_2d_int(poahmm->to_tindex, poahmm->alloced_num_nodes, poahmm->alloced_num_nodes+1, 0);
                RUN(galloc(&poahmm->from_tindex , poahmm->alloced_num_nodes, poahmm->alloced_num_nodes+1));
                //poahmm->from_tindex = malloc_2d_int(poahmm->from_tindex, poahmm->alloced_num_nodes, poahmm->alloced_num_nodes+1, 0);
                for(i = old; i < poahmm->alloced_num_nodes;i++){
                        RUNP(poahmm->nodes[i] = malloc_a_node(poahmm->maxseq_len));
                        poahmm->e_entry_probabilities[i] =0;
                        poahmm->entry_probabilities[i] = prob2scaledprob(0.0);
                        poahmm->e_exit_probabilities[i] =0;
                        poahmm->exit_probabilities[i] = prob2scaledprob(0.0);

                }
        }

        /*if(new_num_samples > poahmm->num_samples){
                poahmm->num_samples = new_num_samples;
                for(i = 0; i < poahmm->alloced_num_nodes;i++){
                        MREALLOC(poahmm->nodes[i]->signal, sizeof(uint32_t) *poahmm->num_samples );
                }
                }*/

        return OK;
ERROR:
        return FAIL;
}

void free_poahmm (struct poahmm* poahmm)
{
        int i;
        if(poahmm){

                if(poahmm->begin){
                        free_a_boundary_node(poahmm->begin);
                }

                if(poahmm->end){
                        free_a_boundary_node(poahmm->end);
                }

                if( poahmm->num_nodes){
                        for(i = 0; i < poahmm->alloced_num_nodes;i++){
                                free_a_node(poahmm->nodes[i]);
                        }
                        MFREE(poahmm->nodes);
                        MFREE(poahmm->rank_sorted_nodes);
                }

                if(poahmm->poa_graph){
                        gfree( poahmm->poa_graph);
                }

                if(poahmm->e_poa_graph){
                        gfree( poahmm->e_poa_graph);
                }

                if(poahmm->to_tindex){
                        gfree( poahmm->to_tindex);
                }

                if(poahmm->from_tindex){
                        gfree( poahmm->from_tindex);
                }

                if(poahmm->emission_M){
                        gfree( poahmm->emission_M);
                }

                if(poahmm->e_emission_M){
                        gfree( poahmm->e_emission_M);
                }
                if(poahmm->entry_probabilities){
                        MFREE(poahmm->entry_probabilities);
                }
                if(poahmm->e_entry_probabilities){
                        MFREE(poahmm->e_entry_probabilities);
                }

                if(poahmm->exit_probabilities){
                        MFREE(poahmm->exit_probabilities);
                }
                if(poahmm->e_exit_probabilities){
                        MFREE(poahmm->e_exit_probabilities);
                }

                if(poahmm->emission_X){
                        MFREE(poahmm->emission_X);//, sizeof(float) * 4);
                }
                if(poahmm->emission_Y){
                        MFREE(poahmm->emission_Y);//, sizeof(float) *4);
                }

                if(poahmm->e_emission_X){
                        MFREE(poahmm->e_emission_X);//, sizeof(float) *4);
                }

                if(poahmm->e_emission_Y){
                        MFREE(poahmm->e_emission_Y);//, sizeof(float) *4);
                }

                if(poahmm->background){
                        MFREE(poahmm->background);
                }

                MFREE(poahmm);
        }
}



struct poahmm_node* malloc_a_node(int maxseq_len)
{
        struct poahmm_node* node = NULL;
        MMALLOC(node , sizeof(struct poahmm_node));
        node->cells = NULL;
        node->signal = NULL; /// leave this for now...

        node->identifier = -1;
        node->nuc = 0;
        node->total_signal = 0;
        node->rank = 0;

        MMALLOC(node->cells, sizeof(struct cell) * maxseq_len);
        return node;
ERROR:
        free_a_node(node);
        return NULL;
}


void free_a_node(struct poahmm_node* node)
{
        if(node){
                if(node->cells){
                        MFREE(node->cells);
                }
                if(node->signal){
                        MFREE(node->signal);
                }
                MFREE(node);
        }
}

struct  poahmm_boundary_node* malloc_a_boundary_node(int maxseq_len)
{
        struct  poahmm_boundary_node* node  = NULL;
        MMALLOC(node, sizeof(struct poahmm_boundary_node));
        node->fY = NULL;
        node->bY = NULL;
        node->cells = NULL;
        MMALLOC(node->fY, sizeof(float) * maxseq_len);
        MMALLOC(node->bY, sizeof(float) * maxseq_len);

        MMALLOC(node->cells, sizeof(struct cell) * maxseq_len);

        node->identifier = 0;

        node->rank = 0;
        return node;
ERROR:
        free_a_boundary_node(node);
        return NULL;
}

void free_a_boundary_node(struct poahmm_boundary_node* node)
{
        if(node){
                if(node->fY){
                        MFREE(node->fY);
                }

                if(node->bY){
                        MFREE(node->bY);
                }

                if(node->cells){
                        MFREE(node->cells);
                }

                MFREE(node);
        }
}


//Stolen from the easel library!!!
static double gamma_fraction(unsigned int* seed , double a)	/* for fractional a, 0 < a < 1 */
{				/* Knuth 3.4.1, exercise 16, pp. 586-587 */
        double p, U, V, X, q;
        //static unsigned seed =  (unsigned int) (time(NULL) * (42));
        p = 2.71828182845904523536028747135 / (a + 2.71828182845904523536028747135);
        do {
                U = (double) rand_r(seed) / (double) RAND_MAX ;
                V = (double) rand_r(seed) / (double) RAND_MAX ;
                while( V == 0.0){
                        V = (double) rand_r(seed) / (double) RAND_MAX ;
                }

                if (U < p) {
                        X = pow(V, 1.0/a);
                        q = exp(-X);
                } else {
                        X = 1. - log(V);
                        q = pow(X, a-1.0);
                }
                U = (double) rand_r(seed) / (double) RAND_MAX ;
        } while (U >= q);
        return X;
}




int cmp_node_rank_low_to_high(const void * a, const void * b)
{
        struct poahmm_node* const *one = a;
        struct poahmm_node* const *two = b;
        if((*one)->rank <  (*two)->rank){
                return -1;
        }else{
                return 1;
        }
}

int cmp_node_rank_high_to_low(const void * a, const void * b)
{
        struct poahmm_node* const *one = a;
        struct poahmm_node* const *two = b;
        if((*one)->rank >  (*two)->rank){
                return -1;
        }else{
                return 1;
        }
}

struct poahmm* set_rank(struct poahmm* poahmm,int index,int rank)
{
        int i;

        LOG_MSG("At %d current rank: %d rank: %d nodes:%d", index, poahmm->nodes[index]->rank ,rank,poahmm->num_nodes);
        if(poahmm->nodes[index]->rank == UINT32_MAX){
                poahmm->nodes[index]->rank = rank;
        }
        rank++;

        for(i = 0; i < poahmm->num_nodes;i++){
                if(poahmm->e_poa_graph[index][i] != 0.0){
                        if(poahmm->nodes[i]->rank == UINT32_MAX){
                                poahmm = set_rank(poahmm, i,rank);
                        }
                }
        }
        return poahmm;
}

int reset_poa_graph_transitions_based_on_counts(struct poahmm* poahmm)
{
        float sum;
        int i,j;

        //for(i = 0; i < poahmm->num_nodes;i++){
        //	fprintf(stdout,"%d %d in: %d out:%d\n", i, poahmm->nodes[i]->rank, poahmm->to_tindex[i][0], poahmm->from_tindex[i][0]);
        //}

        for(i = 0; i < poahmm->num_nodes;i++){
                sum = 0.0f;
                for(j = 1;j < poahmm->from_tindex[i][0];j++){
                        //DPRINTF3("%d -> %d  (%d) ", i,poahmm->from_tindex[i][j],  poahmm->nodes[poahmm->from_tindex[i][j]]->total_signal );
                        //sum += poahmm->nodes[poahmm->from_tindex[i][j]]->total_signal;
                        sum += poahmm->e_poa_graph[i][poahmm->from_tindex[i][j]];
                }

                for(j = 1;j < poahmm->from_tindex[i][0];j++){
                        poahmm->poa_graph[i][poahmm->from_tindex[i][j]] = prob2scaledprob(1.0);// prob2scaledprob( poahmm->e_poa_graph[i][poahmm->from_tindex[i][j]] / sum);
                }

        }
        sum = 0.0f;
        for(i = 0; i < poahmm->num_nodes;i++){
                //if(poahmm->e_entry_probabilities[i]){
                if(poahmm->to_tindex[i][0] == 1){
                        sum += poahmm->e_entry_probabilities[i];
                }


        }
        for(i = 0; i < poahmm->num_nodes;i++){
                poahmm->entry_probabilities[i] = prob2scaledprob(1.0);
                //if(poahmm->e_entry_probabilities[i]){
                if(poahmm->to_tindex[i][0] == 1){
                        poahmm->entry_probabilities[i] = prob2scaledprob(1.0);//prob2scaledprob(poahmm->e_entry_probabilities[i] / sum);
                        //sum += poahmm->e_entry_probabilities[i];
                }
        }

        sum = 0.0f;
        for(i = 0; i < poahmm->num_nodes;i++){
                //if(poahmm->e_exit_probabilities[i]){
                if(poahmm->from_tindex[i][0] == 1){
                        sum += poahmm->e_exit_probabilities[i];
                }


        }
        for(i = 0; i < poahmm->num_nodes;i++){
                poahmm->exit_probabilities[i] = prob2scaledprob(1.0);
                //if(poahmm->e_exit_probabilities[i]){
                if(poahmm->from_tindex[i][0] == 1){
                        poahmm->exit_probabilities[i] = prob2scaledprob(1.0);// prob2scaledprob(poahmm->e_exit_probabilities[i] / sum);
                        //sum += poahmm->e_entry_probabilities[i];
                }
        }
        return OK;
}


int reset_to_from_index(struct poahmm* poahmm)
{
        int i,j,c;

        for(i = 0; i < poahmm->num_nodes;i++){
                c = 0;
                //poahmm->nodes[i]->n_nuc = 0;
                for(j = 0; j < poahmm->num_nodes;j++){
                        if(poahmm->e_poa_graph[i][j]  != 0.0){
                                //		poahmm->nodes[i]->n_nuc |= (1 << poahmm->nodes[j]->n_nuc);
                                poahmm->from_tindex[i][c+1] = j;
                                c++;
                        }
                }
                poahmm->from_tindex[i][0] = c+1;

                for(j = 1; j < poahmm->from_tindex[i][0];j++){
                        fprintf(stdout,"from %d -> %d\n",i, poahmm->from_tindex[i][j]);
                }


        }

        for(j = 0; j < poahmm->num_nodes;j++){
                c = 0;
                //poahmm->nodes[j]->p_nuc = 0;
                for(i = 0; i < poahmm->num_nodes;i++){
                        if(poahmm->e_poa_graph[i][j]  != 0.0){
                                //		poahmm->nodes[j]->p_nuc |= (1 << poahmm->nodes[i]->n_nuc);
                                poahmm->to_tindex[j][c+1] = i;
                                c++;
                        }
                }
                poahmm->to_tindex[j][0] = c+1;
                for(i = 1; i < poahmm->to_tindex[j][0];i++){
                        fprintf(stdout,"to:L %d -> %d\n",j, poahmm->to_tindex[j][i]);
                }
        }
        return OK;
}



/* for testing !!  */


int init_nodes_from_single_sequence(struct poahmm* poahmm, uint8_t* seq, int len)
{
        int i,j;
        struct poahmm_node* node_ptr = NULL;

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
        //start state
        //node_ptr = poahmm->nodes[POAHMM_STARTSTATE];
        //node_ptr->nuc = 4;
        //node_ptr->identifier = POAHMM_STARTSTATE;
        //poahmm->poa_graph[POAHMM_STARTSTATE][2] = prob2scaledprob(1.0);
        //poahmm->poa_graph[POAHMM_STARTSTATE][3] = prob2scaledprob(0.2);
        //poahmm->poa_graph[POAHMM_STARTSTATE][4] = prob2scaledprob(0.1);
        poahmm->e_entry_probabilities[0] = 1.0;

        poahmm->e_exit_probabilities[len-1] = 1.0;

        for(i = 0; i < len-1;i++){
                //poahmm->poa_graph[i][i+1] = prob2scaledprob(1.0);
                poahmm->e_poa_graph[i][i+1] = 1.0;
                //poahmm->poa_graph[i][POAHMM_ENDSTATE] = prob2scaledprob(0.1);
        }


        for(i = 0; i < poahmm->num_nodes;i++){
                poahmm->e_poa_graph[len-1][i] = 0.0;
                //poahmm->poa_graph[len-1][i] = prob2scaledprob(0.0);
        }

        //poahmm->poa_graph[len+1][POAHMM_ENDSTATE] = prob2scaledprob(1.0);

        for(i = 0; i < len;i++){
                //poahmm->poa_graph[n][n +1] = prob2scaledprob(1.0);//
                node_ptr = poahmm->nodes[i];
                node_ptr->nuc = seq[i];
                //node_ptr->signal[poahmm_data->seq_id[index]] = 1;
                node_ptr->total_signal = 1;
                node_ptr->identifier = i;
        }
        //stop state
        //node_ptr = poahmm->nodes[POAHMM_ENDSTATE];
        //node_ptr->nuc = 4;
        //node_ptr->identifier = POAHMM_ENDSTATE;

        poahmm->num_nodes = len;

        for(i = 0; i < poahmm->num_nodes;i++){
                poahmm->nodes[i]->rank = UINT32_MAX;
        }


        poahmm =  set_rank(poahmm, 0, 0);

        RUN(reset_to_from_index(poahmm));
        RUN(reset_poa_graph_transitions_based_on_counts(poahmm));

        for(i = 0; i< poahmm->num_nodes;i++){

                fprintf(stdout,"ENTRY: %f EXIT:%f\t", poahmm->e_entry_probabilities[i],poahmm->e_exit_probabilities[i]);
                for(j = 0 ; j < poahmm->num_nodes;j++){
                        fprintf(stdout,"%f ",poahmm->poa_graph[i][j]);
                        //poahmm->e_poa_graph[i][j] = 0.0;
                }
                fprintf(stdout,"\n");
        }

        return OK;
ERROR:
        return FAIL;
}


int set_rank_transition_poahmm(struct poahmm* poahmm)
{
        poahmm =  set_rank(poahmm, 0, 0);

        RUN(reset_to_from_index(poahmm));
        RUN(reset_poa_graph_transitions_based_on_counts(poahmm));
        return OK;
ERROR:
        return FAIL;
}



int update_poahmm_arch(struct poahmm* poahmm, int* path,int seq_id)
{
        int i,j;
        int c;
        int prev_state = -1;

        c = poahmm->num_nodes;

        //RUN(print_poa_graph(poahmm),"print graph failed.");

        //fprintf(stdout, "new nodes?:%d %d %d\n",poahmm->alloced_num_nodes,poahmm->num_nodes,c);
        for(i = 1; i < path[0];i++){

                if((path[i] & 0xFFFFFFF) == 0xFFFFFFF){
                        if(poahmm->alloced_num_nodes <= c){
                                RUN(check_and_extend_poahmm(poahmm,poahmm->alloced_num_nodes  << 1, 128));
                        }
                        if(prev_state != -1){
                                //			fprintf(stdout,"add edge between: %d and %d\n",prev_state,c);
                                //poahmm->poa_graph[prev_state][c] = -5;
                                poahmm->e_poa_graph[prev_state][c] += 1.0;
                        }
                        poahmm->nodes[c]->identifier = c;
                        poahmm->nodes[c]->signal[seq_id] = 1;
                        poahmm->nodes[c]->total_signal = 1;

                        poahmm->nodes[c]->nuc =path[i] >> 28;

                        if(i == 1){
                                poahmm->e_entry_probabilities[c] = poahmm->nodes[c]->total_signal;
                        }

                        if(i == (path[0]-1)){
                                poahmm->e_exit_probabilities[c] =poahmm->nodes[c]->total_signal;
                        }

                        prev_state = c;
                        c++;


                }else{
                        if((path[i] >> 28) == poahmm->nodes[path[i] & 0xFFFFFFF]->nuc){
                                if(prev_state != -1){
                                        //				fprintf(stdout,"add edge between: %d and %d\n",prev_state,path[i] & 0xFFFFFFF);
                                        //poahmm->poa_graph[prev_state][path[i] & 0xFFFFFFF] = -5;
                                        poahmm->e_poa_graph[prev_state][path[i] & 0xFFFFFFF] += 1.0;


                                }
                                //poahmm->nodes[path[i] & 0xFFFFFFF]->signal[seq_id] += 1;
                                poahmm->nodes[path[i] & 0xFFFFFFF]->total_signal  += 1;
                                if(i == 1){
                                        poahmm->e_entry_probabilities[path[i] & 0xFFFFFFF] = poahmm->nodes[path[i] & 0xFFFFFFF]->total_signal;
                                }
                                if(i == (path[0]-1)){
                                        poahmm->e_exit_probabilities[path[i] & 0xFFFFFFF] = poahmm->nodes[path[i] & 0xFFFFFFF]->total_signal;
                                }

                                prev_state = path[i] & 0xFFFFFFF;
                        }else{
                                if(path[i] >> 28 == 4){

                                }else{
                                        //DPRINTF3("old: %d	%d\n",poahmm->alloced_num_nodes,c);
                                        if(poahmm->alloced_num_nodes <= c){
                                                RUN(check_and_extend_poahmm(poahmm,poahmm->alloced_num_nodes  << 1, 128));
                                        }
                                        if(prev_state != -1){
                                                //fprintf(stdout,"add edge between: %d and %d\n",prev_state,c);
                                                //poahmm->poa_graph[prev_state][c] = -5;
                                                poahmm->e_poa_graph[prev_state][c] += 1.0;
                                        }
                                        poahmm->nodes[c]->nuc = path[i] >> 28;
//                                        poahmm->nodes[c]->signal[seq_id] = 1;
                                        poahmm->nodes[c]->total_signal = 1;
                                        poahmm->nodes[c]->identifier = c;
                                        if(i == 1){
                                                poahmm->e_entry_probabilities[c] = poahmm->nodes[c]->total_signal;
                                        }
                                        if(i == (path[0]-1)){
                                                poahmm->e_exit_probabilities[c] =poahmm->nodes[c]->total_signal;
                                        }

                                        prev_state = c;//path[i] & 0xFFFFFFF;
                                        c++;

                                }
                        }
                }
        }
        //poahmm->poa_graph[prev_state][POAHMM_ENDSTATE] = -5;

        poahmm->num_nodes = c;

        for(i = 0; i <poahmm->num_nodes;i++){
                poahmm->nodes[i]->rank = UINT32_MAX; // reset rank...
                poahmm->rank_sorted_nodes[i] =poahmm->nodes[i];
        }


        //exit(0);

        for(j = 0; j < poahmm->num_nodes;j++){
                c = 0;
                for(i = 0; i <poahmm->num_nodes; i++){
                        //		fprintf(stdout,"%f ",poahmm->poa_graph[i][j]);
                        if(poahmm->e_poa_graph[i][j] != 0.0){
                                c++;
                        }
                }
                if(!c){
                        //		fprintf(stdout,"NEWRUN:%d\n",j);
                        poahmm =  set_rank(poahmm, j, 0);
                }
                //	fprintf(stdout,"\n");
        }

        /*for(i = 0; i <poahmm->num_nodes; i++){
          sum = 0.0;
          for(j = 0; j < poahmm->num_nodes;j++){
          if(poahmm->poa_graph[i][j] != prob2scaledprob(0.0)){
          sum += 1.0;
          }
          }

          for(j = 0; j < poahmm->num_nodes;j++){
          if(poahmm->poa_graph[i][j] != prob2scaledprob(0.0)){
          poahmm->poa_graph[i][j] = prob2scaledprob(1.0 / sum);
          }
          }
          }*/

        RUN(reset_to_from_index(poahmm));

        //set counts on starting probs.

        for(i = 0; i < poahmm->num_nodes;i++){
                //if(poahmm->to_tindex[i][0] == 1){
                //	poahmm->e_entry_probabilities[i] = poahmm->nodes[i]->total_signal;
                //}else{
                //	poahmm->e_entry_probabilities[i] = 0;
                //}
                poahmm->entry_probabilities[i] = prob2scaledprob(0.0);
                poahmm->exit_probabilities[i] = prob2scaledprob(0.0);

        }


        RUN(reset_poa_graph_transitions_based_on_counts(poahmm));
        //exit(0);



        //print_poa_graph(poahmm);
        //exit(0);
        return OK;
ERROR:
        return FAIL;
}
