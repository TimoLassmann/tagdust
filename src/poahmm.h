#ifndef POAHMM_H
#define POAHMM_H

struct cell{
        float fM;
        float fX;
        float fY;
        union {
                float bM;
                struct{
                        unsigned int M_trans : 3;
                        unsigned int M_to_state : 29;
                };
        };
        union{
                float bX;
                struct{
                        unsigned int X_trans : 3;
                        unsigned int X_to_state : 29;
                };
        };
        union{
                float bY;
                struct{
                        unsigned int Y_trans : 3;
                        unsigned int Y_to_state : 29;
                };
        };
};

struct poahmm_boundary_node{
        float* fY;
        float* bY;
        struct cell* cells;
        uint32_t rank;
        int identifier;
};


struct poahmm_node{
        struct cell* cells;
        uint32_t* signal;
        uint32_t rank;
        uint8_t nuc;
        uint8_t type;
        int total_signal;
        int identifier;
};

struct poahmm{

        struct poahmm_node** nodes;
        struct poahmm_node** rank_sorted_nodes;

        struct poahmm_boundary_node* begin;
        struct poahmm_boundary_node* end;

        float* entry_probabilities;
        float* e_entry_probabilities;

        float* exit_probabilities;
        float* e_exit_probabilities;


        float** poa_graph;
        float** e_poa_graph;

        int** to_tindex;
        int** from_tindex;

        float** emission_M;
        float* emission_X;
        float* emission_Y;
        float** e_emission_M;
        float* e_emission_X;
        float* e_emission_Y;

        float* background;

        unsigned int seed;
        float pseudo_weight;
        int num_nodes;
        int num_samples;
        int alloced_num_nodes;
        int maxseq_len;

        float f_score;
        float b_score;
        float r_score;

        float MM;
        float MX;
        float MY;

        float XX;
        float XM;
        float XY;

        float YY;
        float YM;
        float YX;

        float e_MM;
        float e_MX;
        float e_MY;

        float e_XX;
        float e_XM;
        float e_XY;

        float e_YY;
        float e_YM;
        float e_YX;

};

struct global_poahmm_param{
        float back[5];
        float base_error;
        float indel_freq;
};

extern int random_poahmm(struct poahmm* poahmm, uint8_t* seq, int len);
extern int forward_poahmm(struct poahmm* poahmm, uint8_t* seq, int len);
extern int backward_poahmm(struct poahmm* poahmm, uint8_t* seq, int len);

extern int viterbi_poahmm(struct poahmm* poahmm, uint8_t* seq, int len,  uint32_t* path);

extern struct poahmm*  init_poahmm(struct global_poahmm_param* param, int max_len);
extern void free_poahmm (struct poahmm* poahmm);

extern int init_nodes_from_single_sequence(struct poahmm* poahmm, uint8_t* seq, int len);
extern int set_rank_transition_poahmm(struct poahmm* poahmm);
#endif
