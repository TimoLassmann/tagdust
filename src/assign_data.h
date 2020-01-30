#ifndef ASSIGN_DATA_H
#define ASSIGN_DATA_H

#define READ_TYPE 1
#define BAR_TYPE 2
#define UMI_TYPE 3

#include "arch_lib.h"
#include <stdint.h>

#include "tlrbtree.h"
#include "kstring.h"
#include "khash.h"

#include "correct.h"

#include "read_groups.h"

#define READ_FAILQ 1
#define READ_FAILR 2
#define READ_FAILP 4
#define READ_NBAR 8
#define READ_FAILC 16


struct seq_bit{
//        struct seq_bit* next;
        kstring_t p;
        kstring_t p_corr;
        kstring_t q;
        char code;
        uint8_t type;
        uint8_t correct_index;
        uint8_t fail;
};

struct seq_bit_vec{
        struct seq_bit** bits;
        int* out_file_id;

        char* name;
        kstring_t append;
        int sample_group;
        uint8_t fail;
        uint8_t num_bit;
};

struct demux_struct{
        struct file_handler* f_hand;
        char* key;
        char* out_filename;
        int id;
        int count;
};

struct bit_annotation{
        khash_t(exact)* bar_hash;
        kstring_t name;
        kstring_t c_name;
        kstring_t q_name;
};

struct assign_struct{
        struct seq_bit_vec** bit_vec;
        struct bit_annotation** bit_ann;
        struct rbtree_root* demux_names;
        struct qsubscore* subm;
        //khash_t(exact)** exact;
        //struct rbtree_root* file_names;
        int* loc_out_reads;
        int* file_index;
        int block_size;
        int max_seq_len;
        int max_bar_len;
        //int append_len;
        int num_files;
        int out_reads;
        int n_exact_hash;
        int num_bits;
        int num_reads;
        int alloc_total;
};



extern int init_assign_structure(struct assign_struct** assign,struct arch_library* al, struct read_groups* rg, char* prefix, int total,int bam);


extern int sort_as_by_file_type(struct assign_struct* as);
//extern int alloc_assign_structure(struct assign_struct** assign,int num_files);
extern int reset_assign_structute(struct assign_struct* as);
//aextern int set_up_assign_structure(struct arch_library* al,struct assign_struct* as);
extern void free_assign_structure(struct assign_struct* as);


extern int post_process_assign(struct assign_struct* as);
extern int ref_correct(khash_t(exact) * h ,struct qsubscore* subm, struct seq_bit* sb, int q_offset);

#endif
