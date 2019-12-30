#ifndef PST_H
#define PST_H


#include <stdint.h>

struct pst_node{
        struct pst_node* next[5];
        float nuc_probability[5];
        char* label;
};

struct pst {
        struct pst_node* pst_root;
        struct pst_node* ppt_root;
        char** suffix_array;

        uint32_t** counts;
        //char** suffix_array_local;
        int* seq_id_in_suffix;

        struct suffix_node** sn;

        int total_len;

        float p_min;
        float gamma_min;
        //float alpha;
        //float lamba;
        float r;
        int L;

        float numseq;
        float mean_length;
        int suffix_len;
        //int suffix_len_local;
        //int current_suffix_size;

};

extern int nuc_to_internal(char c);

#endif
