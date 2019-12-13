#include "extract_reads.h"
#include <omp.h>
#include <math.h>
#include <string.h>

#include "hmm_model_bag.h"
#include "core_hmm_functions.h"
#include "assign_data.h"

#include "tllogsum.h"

#define READ_CHUNK_SIZE 100000
#define CHUNKS 10


static int process_read(struct read_info* ri, int* label, struct read_structure* rs , struct seq_bit_vec* b , int i_file);

static int sanity_check_inputs(struct read_info_buffer** rb, int num_files);

//static int run_extract( struct assign_struct* as,  struct read_info_buffer** rb, struct arch_library* al, strucnt seq_stats* si,int i_file,int i_hmm);

static int run_extract( struct assign_struct* as,  struct read_info_buffer** rb, struct arch_library* al, struct seq_stats* si,int i_file,int i_chunk);

int extract_reads(struct arch_library* al, struct seq_stats* si,struct parameters* param)
{
        struct read_info_buffer** rb = NULL;
        struct file_handler** f_hand = NULL;
        struct assign_struct* as = NULL;


        int (*fp)(struct read_info_buffer* rb, struct file_handler* f_handle) = NULL;
        int i,j,c;
        int total_read;

        DECLARE_TIMER(t1);


        /* figure out how many barcodes etc */
        RUN(init_assign_structure(&as, al,   CHUNKS* READ_CHUNK_SIZE));
        //RUN(galloc(&as->assignment, as->total, as->num_barcodes));
        /* not sure if this is required  */
        //exit(0);
        omp_set_num_threads(param->num_threads);

        /* here I combine: architectures with sequence parameters of individual input files  */
        /* to: check which architecture belongs to which read */
        ASSERT(param != NULL, "no parameters");
        ASSERT(al != NULL, "No arch library");
        ASSERT(si != NULL, "no seq stats");

        //LOG_MSG("Got here");
        MMALLOC(rb, sizeof(struct read_info_buffer*) * param->num_infiles * CHUNKS);
        MMALLOC(f_hand, sizeof(struct file_handler*) * param->num_infiles);

        for(i = 0; i < param->num_infiles;i++){
                f_hand[i] = NULL;
                RUN(io_handler(&f_hand[i], param->infile[i]));
        }
        for(i = 0; i < param->num_infiles * CHUNKS;i++){
                rb[i] = NULL;
                RUN(alloc_read_info_buffer(&rb[i], READ_CHUNK_SIZE));

        }

        while(1){
                /* read everything in  */
                total_read = 0;
                for(i = 0; i < param->num_infiles * CHUNKS;i++){
                        rb[i]->offset = 0;
                        rb[i]->num_seq = 0;
                }
                for(i = 0; i < param->num_infiles;i++){

                        if(f_hand[i]->sam == 0){
                                fp = &read_fasta_fastq;
                        }else {
                                fp = &read_sam_chunk;
                        }
                        /* set offset of first chunk to be the off ser  */
                        rb[i * CHUNKS]->offset =
                                rb[i* CHUNKS + CHUNKS-1]->offset
                                + rb[i* CHUNKS + CHUNKS-1]->num_seq;
                        for(j = 0; j < CHUNKS;j++){
                                //LOG_MSG("Reading %d chunk %d ->%d", i,j, i + j * param->num_infiles);
                                RUN(fp(rb[i * CHUNKS + j],f_hand[i]));//  param,file,&numseq));
                                total_read += rb[i* CHUNKS + j]->num_seq;
                        }
                        for(j =1; j < CHUNKS;j++){
                                rb[i* CHUNKS + j]->offset = rb[i* CHUNKS + j-1]->offset +  rb[i* CHUNKS + j-1]->num_seq;
                        }
                }
                /* Assign name to assign struct  */
                //for(i = 0; i < param->num_infiles;i++){
                i = 0;
                for(j = 0; j < CHUNKS;j++){
                        for(c = 0;c < rb[j]->num_seq;c++){
                                as->bits[i]->name = rb[j]->ri[c]->name;
                                i++;
                        }
                }
                as->num_reads = i;
                        //}
                //if(!  )
                if(!total_read){
                        LOG_MSG("Done");
                        break;
                }
                //LOG_MSG("Pausing here");
                //sleep(100);

                /* sanity check input - are the files correctly sorted */
                RUN(sanity_check_inputs(rb,param->num_infiles));
                fflush(stdout);
                /* extract reads  */

                START_TIMER(t1);

#pragma omp parallel default(shared)
#pragma omp for collapse(2) private(i, j)
                for(i = 0; i < param->num_infiles;i++){
                        for(j = 0; j < CHUNKS;j++){
                                run_extract(as, rb,al,si,i,j);
                        }
                }
                /* not necessary I think ... */
#pragma omp barrier
                RUN(post_process_assign(as));

                STOP_TIMER(t1);
                LOG_MSG("Took %f ",GET_TIMING(t1));
                RUN(write_all(as, param->outfile));
                RUN(reset_assign_structute(as));

        }
        free_assign_structure(as);
        for(i = 0; i < param->num_infiles* CHUNKS;i++){
                free_read_info_buffer(rb[i]);
        }
        MFREE(rb);

        for(i = 0; i < param->num_infiles;i++){
                pclose(f_hand[i]->f_ptr);
                MFREE(f_hand[i]);
        }
        MFREE(f_hand);

        return OK;
ERROR:
        return FAIL;
}



int run_extract( struct assign_struct* as,  struct read_info_buffer** rb, struct arch_library* al, struct seq_stats* si,int i_file,int i_chunk)
{
        struct model_bag* mb = NULL;
        struct read_info** ri = NULL;
        //int* tmp_bar;
        float pbest = 0.0f;
        float Q = 0.0f;

        int num_seq;
        int c;
        int i;
        int i_hmm;
        int seq_offset;

        c = i_file * CHUNKS + i_chunk;

        i_hmm = al->arch_to_read_assignment[i_file];
        mb = init_model_bag(al->read_structure[i_hmm], si->ssi[i_file], i_hmm);
        num_seq = rb[c]->num_seq;
        ri = rb[c]->ri;

        seq_offset = rb[c]->offset;
        //int tid = omp_get_thread_num();
        //LOG_MSG("Offset = %d",
        for(i = 0; i < num_seq;i++){
                //LOG_MSG("%d f:%d  on %d", tid,i_file,i);
                RUN(backward(mb, ri[i]->seq,ri[i]->len));
                RUN(forward_max_posterior_decoding(mb, ri[i], ri[i]->seq ,ri[i]->len));

                //pbest = ri[i]->mapq;
                pbest = prob2scaledprob(0.0f);
                pbest = logsum(pbest, mb->f_score);
                pbest = logsum(pbest, mb->r_score);

                pbest = 1.0 - scaledprob2prob(  (ri[i]->bar_prob + mb->f_score ) - pbest);

                if(!pbest){
                        Q = 40.0;
                }else if(pbest == 1.0){
                        Q = 0.0;
                }else{
                        Q = -10.0 * log10(pbest) ;
                }

                ri[i]->mapq = Q;
                ri[i]->bar_prob = 100;


                if(Q < al->confidence_thresholds[i_file]){
                        as->bits[seq_offset+i]->Q[i_file] = -Q;
                }else{
                        as->bits[seq_offset+i]->Q[i_file] = Q;
                }
                //fprintf(stdout,"File:%d %f %f\n",i_file, Q, al->confidence_thresholds[i_file]);
                //as->bits[seq_offset+i]->pass = 0;
                //}
                //tmp_bar = as->assignment[seq_offset + i] + assign_offset;
                //print_labelled_reads(mb,data->param ,ri[i]);
                //RUN(extract_reads(mb,data->param,ri[i]));
                //RUN(analyze_and_extract_reads(as,mb, al->read_structure[i_hmm], ri[i], i_file, al->confidence_thresholds[i_file]));
                //LOG_MSG("i:%d, seqoff = %d",i, seq_offset);
                //LOG_MSG("t: %d Working on file: %d  %d (%d seq) offset: %d",tid, i_file,i,num_seq,rb[c]->offset);

                RUN(process_read(ri[i], mb->label,al->read_structure[i_hmm], as->bits[seq_offset+i],i_file));
        }
        free_model_bag(mb);
        return OK;
ERROR:
        return FAIL;
}

int process_read(struct read_info* ri, int* label, struct read_structure* rs , struct seq_bit_vec* b , int i_file)
{
        struct seq_bit* sb = NULL;
        char* type;
        char* read_label;
        char c;
        int j;
        int segment;
        int hmm_in_segment;
        int c1;
        //int umi;
        //int umi_len;
        //int s_pos = 0;
        int len = ri->len;
        char old = '?';
        int local_bit_index;

        read_label = ri->labels+1;
        type = rs->type;

        local_bit_index = 0;
        for(j = 0; j < b->num_bit;j++){
                sb = b->bits[j];
                if(sb->file == i_file){
                        local_bit_index = j;
                        break;
                }
        }
        //assign_offset = as->num_per_file[i_file];
        for(j = 0; j < len;j++){
                c1 = label[(int)read_label[j]];
                segment = c1 & 0xFFFF; //which segment
                hmm_in_segment = (c1 >> 16) & 0x7FFF; // which HMM in a particular segment....
                c = type[segment];
                switch (c) {
                case 'F':
                        if(c != old){
                                sb = b->bits[local_bit_index];
                                sb->len = 0;
                                ASSERT(i_file == sb->file, "Oh dear: want %d got %d",i_file,sb->file);
                                //sb->file = i_file;
                                sb->type = UMI_TYPE;
                                sb->p = ri->seq+j;
                                local_bit_index++;
                        }

                        //umi_len++;
                        //umi = (umi << 2 )|  (ri->seq[j] & 0x3);
                        //ri->seq[s_pos] = 65; // 65 is the spacer! nucleotides are 0 -5....
                        //ri->qual[s_pos] = 65;
                        break;
                case 'B':
                        if(c != old){
                                sb = b->bits[local_bit_index];
                                sb->len = 0;
                                ASSERT(i_file == sb->file, "Oh dear: want %d got %d",i_file,sb->file);
                                //sb->file = i_file;
                                sb->type = BAR_TYPE;
                                sb->p = rs->sequence_matrix[segment][hmm_in_segment];
                                sb->len = rs->segment_length[segment];
                                local_bit_index++;
                        }
                        break;
                case 'R':
                        if(c != old){
                                sb = b->bits[local_bit_index];
                                sb->len = 0;
                                ASSERT(i_file == sb->file, "Oh dear: want %d got %d",i_file,sb->file);
                                //sb->file = i_file;
                                sb->type = READ_TYPE;
                                sb->p = ri->seq + j;
                                sb->q = ri->qual + j;
                                local_bit_index++;
                        }
                        sb->len++;
                        break;
                default:
                        break;
                }
                old = c;
        }
        //ri->len = s_pos;
        return OK;
ERROR:
        return FAIL;
}


int sanity_check_inputs(struct read_info_buffer** rb, int num_infiles)
{
        int i,j,c,g;

        ASSERT(rb != NULL, "No reads");

        for(g = 0; g < CHUNKS;g++){
                //offset = g * num_infiles;
                //rb += offset;
                for(i = 0; i < num_infiles;i++){
                        for(j = i +1; j < num_infiles;j++){
                                if(rb[i*CHUNKS + g ]->num_seq != rb[j*CHUNKS + g]->num_seq){
                                        ERROR_MSG("Input File:%d and %d differ in number of entries.\n",i,j);

                                }
                        }
                }
                for(i = 0; i < num_infiles;i++){
                        for(j = i +1; j < num_infiles;j++){
                                for(c = 0;c < MACRO_MIN(1000, rb[i*CHUNKS + g]->num_seq);c++){
                                        if(compare_read_names(rb[i*CHUNKS + g]->ri[c]->name,rb[j*CHUNKS + g]->ri[c]->name)){
                                                ERROR_MSG("Files seem to contain reads in different order:\n%s\n%s\n",rb[i*CHUNKS + g]->ri[c]->name,rb[j*CHUNKS + g]->ri[c]->name);
                                        }
                                }
                        }
                }
        }
        return OK;
ERROR:
        return FAIL;
}


