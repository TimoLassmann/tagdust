#include "extract_reads.h"

#include <omp.h>

#include "arch_lib.h"
#include "seq_stats.h"
#include "hmm_model_bag.h"
#include "core_hmm_functions.h"

#define READ_CHUNK_SIZE 1000000

static int sanity_check_inputs(struct read_info_buffer** rb, int num_files);

int extract_reads(struct arch_library* al, struct seq_stats* si,struct parameters* param)
{
        struct read_info_buffer** rb = NULL;
        struct file_handler** f_hand = NULL;
        int (*fp)(struct read_info_buffer* rb, struct file_handler* f_handle) = NULL;
        int i,j;
        int total_read;

        omp_set_num_threads(param->num_threads);

        /* here I combine: architectures with sequence parameters of individual input files  */
        /* to: check which architecture belongs to which read */
        ASSERT(param != NULL, "no parameters");
        ASSERT(al != NULL, "No arch library");
        ASSERT(si != NULL, "no seqstats");


        MMALLOC(rb, sizeof(struct read_info_buffer*) * param->num_infiles);
        MMALLOC(f_hand, sizeof(struct file_handler*) * param->num_infiles);
        for(i = 0; i < param->num_infiles;i++){
                rb[i] = NULL;
                RUN(alloc_read_info_buffer(&rb[i], READ_CHUNK_SIZE));
                RUN(io_handler(&f_hand[i], param->infile[i]));
        }

        while(1){
                /* read everything in  */
                total_read = 0;

                for(i = 0; i < param->num_infiles;i++){
                        if(f_hand[i]->sam == 0){
                                fp = &read_fasta_fastq;
                        }else {
                                fp = &read_sam_chunk;
                        }
                        RUN(fp(rb[i],f_hand[i]));//  param,file,&numseq));
                        total_read += rb[i]->num_seq;
                }
                if(!total_read){
                        break;
                }


                /* sanity check input - are the files correctly sorted */
                RUN(sanity_check_inputs(rb,param->num_infiles));

                /* extract reads  */

        }
        for(i = 0; i < param->num_infiles;i++){
                pclose(f_hand[i]->f_ptr);
                MFREE(f_hand[i]);
        }

        MFREE(f_hand);

        return OK;

ERROR:
        return FAIL;
}

int sanity_check_inputs(struct read_info_buffer** rb, int num_infiles)
{
        int i,j,c;

        ASSERT(rb != NULL, "No reads");


        for(i = 0; i < num_infiles;i++){
                for(j = i +1; j < num_infiles;j++){
                        if(rb[i]->num_seq != rb[j]->num_seq){
                                ERROR_MSG("Input File:%d and %d differ in number of entries.\n",i,j);

                        }
                }
        }
        for(i = 0; i < num_infiles;i++){
                for(j = i +1; j < num_infiles;j++){
                        for(c = 0;c < MACRO_MIN(1000, rb[i]->num_seq);c++){
                                if(compare_read_names(rb[i]->ri[c]->name,rb[j]->ri[c]->name)){
                                        ERROR_MSG("Files seem to contain reads in different order:\n%s\n%s\n",rb[i]->ri[c]->name,rb[j]->ri[c]->name);
                                }
                        }
                }
        }


        return OK;
ERROR:
        return FAIL;
}
