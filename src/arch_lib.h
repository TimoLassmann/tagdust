#ifndef ARCH_LIB_H
#define ARCH_LIB_H

#include <stdint.h>

#define ARCH_ETYPE_EXTRACT 1
#define ARCH_ETYPE_APPEND 2
#define ARCH_ETYPE_SPLIT 3
#define ARCH_ETYPE_IGNORE 4
#define ARCH_ETYPE_PARTIAL 5
#define ARCH_ETYPE_WOBBLE_LEFT 6
#define ARCH_ETYPE_WOBBLE_RIGHT 7


struct segment_specs{
        char* name;
        char** seq;
        int num_seq;
        int max_len;
        int min_len;
        int alloc_len;
        uint8_t extract;
};

struct read_structure{
        struct segment_specs** seg_spec;
        int num_segments;
        int alloc_num_seqments;
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

extern int print_segment_spec(const struct segment_specs* spec);
/* emit sequences directly from read structure  */




#endif
