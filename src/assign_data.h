#ifndef ASSIGN_DATA_H
#define ASSIGN_DATA_H

#define READ_TYPE 1
#define BAR_TYPE 2
#define UMI_TYPE 3


#include "arch_lib.h"
#include <stdint.h>

#include "tlrbtree.h"

struct seq_bit{
        char* p;
        char* q;
        uint16_t len;
        uint8_t type;
        uint8_t file;
};

struct seq_bit_vec{
        struct seq_bit** bits;
        char* name;
        //char* bc;
        char* umi;
        float* Q;
        int sample_group;
        uint8_t fail;
        uint8_t num_bit;
};

struct demux_struct{
        char* name;
        int id;
        int count;
};

struct assign_struct{
        struct seq_bit_vec** bits;
        struct rbtree_root* demux_names;
        struct rbtree_root* file_names;
        int block_size;
        int max_seq_len;
        int max_bar_len;
        int num_files;
        int out_reads;
        int num_bits;
        int num_reads;
        int alloc_total;
};



extern int init_assign_structure(struct assign_struct** assign,struct arch_library* al, int total);

extern int sort_as_by_file_type(struct assign_struct* as);
//extern int alloc_assign_structure(struct assign_struct** assign,int num_files);
extern int reset_assign_structute(struct assign_struct* as);
//aextern int set_up_assign_structure(struct arch_library* al,struct assign_struct* as);
extern void free_assign_structure(struct assign_struct* as);


extern int post_process_assign(struct assign_struct* as);



#endif
