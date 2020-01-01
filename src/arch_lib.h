#ifndef ARCH_LIB_H
#define ARCH_LIB_H

#include <stdint.h>

#define ARCH_ETYPE_EXTRACT 1
#define ARCH_ETYPE_APPEND 2
#define ARCH_ETYPE_SPLIT 3
#define ARCH_ETYPE_IGNORE 4

struct read_structure{
        char*** sequence_matrix;
        int* segment_length;
        int* numseq_in_segment;
        char* type;
        char** segment_name;
        uint8_t* extract;
        int* min_len;
        int* max_len;
        int num_segments;
        int alloc_num_seqments;
        //uint8_t assignment_to_read;
};

struct arch_library{
        struct read_structure** read_structure;
        char** spec_line;
        float** arch_posteriors;
        float* confidence_thresholds;
        int* arch_to_read_assignment;
        int num_arch;

        int alloc_num_arch;
        int num_file;
};


extern int read_architecture_files(struct arch_library* al, char* filename);
extern int read_arch_into_lib(struct arch_library* al, char** list, int len);

extern int alloc_arch_lib(struct arch_library** arch);
extern void free_arch_lib(struct arch_library* arch);

#endif
