#ifndef ARCH_LIB_H
#define ARCH_LIB_H

struct read_structure{
        char*** sequence_matrix;
        int* numseq_in_segment;
        char* type;
        int num_segments;
        //uint8_t assignment_to_read;
};

struct arch_library{
        struct read_structure** read_structure;
        char** spec_line;
        int* arch_to_read_assignment;
        int num_arch;
        int num_alloc_arch;
};


extern int read_architecture_files(struct arch_library* al, char* filename);
extern int read_arch_into_lib(struct arch_library* al, char** list, int len);

extern int alloc_arch_lib(struct arch_library** arch);
extern void free_arch_lib(struct arch_library* arch);

#endif
