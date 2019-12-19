#ifndef TLSEQIO_H
#define TLSEQIO_H

#include <stdint.h>
#include <stdio.h>

#ifdef TLSEQIO_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif


#include <zlib.h>



#define TLSEQIO_READ 0

#define TLSEQIO_WRITE 1
#define TLSEQIO_WRITE_GZIPPED 2

typedef struct file_handler file_handler;

struct tl_seq{
        uint8_t* seq;
        char* name;
        char* qual;
        int malloc_len;
        int len;
};

struct tl_seq_buffer{
        struct tl_seq** sequences;
        int malloc_num;
        int num_seq;
        int max_len;
        int is_fastq;
        int L;
};

EXTERN int open_fasta_fastq_file(struct file_handler** fh,char* filename, int mode);
EXTERN int read_fasta_fastq_file(struct file_handler* fh, struct tl_seq_buffer** seq_buf, int num);
EXTERN int close_fasta_fastq_file(struct file_handler** fh);

EXTERN int write_fasta_fastq(struct tl_seq_buffer* sb, struct file_handler* fh);

EXTERN void free_tl_seq_buffer(struct tl_seq_buffer* sb);


EXTERN int echo_file(struct file_handler* f);
#undef TLSEQIO_IMPORT
#undef EXTERN

#endif
