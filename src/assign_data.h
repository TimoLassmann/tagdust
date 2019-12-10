#ifndef ASSIGN_DATA_H
#define ASSIGN_DATA_H

#define READ_TYPE 1
#define BAR_TYPE 2
#define UMI_TYPE 3


#include "arch_lib.h"
#include <stdint.h>

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
        char* bc;
        char* umi;
        float* Q;
        int pass;
        uint8_t num_bit;
};

struct assign_struct{
        struct seq_bit_vec** bits;
        int num_files;
        int num_bits;
        int total;
};

extern int init_assign_structure(struct assign_struct** assign,struct arch_library* al, int total);
//extern int alloc_assign_structure(struct assign_struct** assign,int num_files);
extern int set_up_assign_structure(struct arch_library* al,struct assign_struct* as);
extern void free_assign_structure(struct assign_struct* as);



#endif
