#ifndef POAHMM_H
#define POAHMM_H

struct poahmm;

struct global_poahmm_param{
        float back[5];
        float base_error;
        float indel_freq;
        int average_seq_length;
        int max_seq_len;
};

extern int random_poahmm(struct poahmm* poahmm, uint8_t* seq, int len);
extern int forward_poahmm(struct poahmm* poahmm, uint8_t* seq, int len);
extern int backward_poahmm(struct poahmm* poahmm, uint8_t* seq, int len);

extern int viterbi_poahmm(struct poahmm* poahmm, uint8_t* seq, int len,  uint32_t* path);
extern int viterbi_poahmm_banded(struct poahmm* poahmm,const uint8_t* seq,const  int len,  uint32_t* path,const int band);

//int viterbi_poahmm_banded(struct poahmm* poahmm, uint8_t* seq, int len,  uint32_t* path,int band);
extern struct poahmm*  init_poahmm(struct global_poahmm_param* param);
extern int resize_poahmm(struct poahmm* poahmm,int num_states, int new_maxlen);
extern void free_poahmm (struct poahmm* poahmm);

extern int init_nodes_from_single_sequence(struct poahmm* poahmm, uint8_t* seq, int len);
extern int set_rank_transition_poahmm(struct poahmm* poahmm);
#endif
