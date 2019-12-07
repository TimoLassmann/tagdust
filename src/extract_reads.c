#include "extract_reads.h"

#include <omp.h>
#include <math.h>
#include <string.h>

#include "hmm_model_bag.h"
#include "core_hmm_functions.h"


#include "tllogsum.h"

#define READ_CHUNK_SIZE 100000
#define CHUNKS 5

static int sanity_check_inputs(struct read_info_buffer** rb, int num_files);

static int run_extract(struct read_info_buffer** rb, struct arch_library* al, struct seq_stats* si,int i_file,int i_hmm);

static int analyze_and_extract_reads(struct model_bag* mb, struct read_structure* read_structure,  struct read_info* ri, float threshold);


static int make_extracted_read(struct model_bag* mb, struct read_structure* read_structure,  struct read_info* ri);

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
        ASSERT(si != NULL, "no seq stats");


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

                for(i = 0; i < param->num_infiles;i++){
                        LOG_MSG("Reading %d chunk %d ->%d", i,j, i + j * param->num_infiles);
                        if(f_hand[i]->sam == 0){
                                fp = &read_fasta_fastq;
                        }else {
                                fp = &read_sam_chunk;
                        }
                        for(j = 0; j < CHUNKS;j++){
                                RUN(fp(rb[i * CHUNKS + j],f_hand[i]));//  param,file,&numseq));
                                total_read += rb[i* CHUNKS + j]->num_seq;
                        }
                }

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

#pragma omp parallel default(shared)
#pragma omp for collapse(2) private(i, j)
                for(i = 0; i < param->num_infiles;i++){
                        for(j = 0; j < CHUNKS;j++){
                                run_extract(rb,al,si,i,j);
                        }
                }
        }

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

int run_extract(struct read_info_buffer** rb, struct arch_library* al, struct seq_stats* si,int i_file,int i_chunk)
{
        struct model_bag* mb = NULL;
        struct read_info** ri = NULL;
        float pbest = 0.0f;
        float Q = 0.0f;

        int num_seq;
        int c;
        int i;
        int i_hmm;

        c = i_file * CHUNKS + i_chunk;


        i_hmm = al->arch_to_read_assignment[i_file];
        mb = init_model_bag(al->read_structure[i_hmm], si->ssi[i_file], i_hmm);
        num_seq = rb[c]->num_seq;
        ri = rb[c]->ri;
        LOG_MSG("Working on %d (%d seq)",c,num_seq);
        for(i = 0; i < num_seq;i++){
                RUN(backward(mb, ri[i]->seq,ri[i]->len));
                RUN(forward_max_posterior_decoding(mb, ri[i], ri[i]->seq ,ri[i]->len));

                pbest = ri[i]->mapq;

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
        }

        for(i = 0; i < num_seq;i++){

                ri[i]->bar_prob = 100;
                //print_labelled_reads(mb,data->param ,ri[i]);
                //RUN(extract_reads(mb,data->param,ri[i]));
                RUN(analyze_and_extract_reads(mb, al->read_structure[i_hmm], ri[i], al->confidence_thresholds[i_file]));
        }

        free_model_bag(mb);
        return OK;
ERROR:
        return FAIL;
}




int analyze_and_extract_reads(struct model_bag* mb, struct read_structure* read_structure,  struct read_info* ri, float threshold)
{
        int j,c1,c2,c3,key,bar,mem,fingerlen,required_finger_len;

        int s_pos = 0;
        key = 0;
        bar = -1;
        mem = -1;
        //ret = 0;
        int offset = 0;
        int len;
        int hmm_has_barcode = 0;
        int too_short = 0;
        int in_read = 0;

        len = ri->len;
        /*if(param->matchstart != -1 || param->matchend != -1){
                offset = param->matchstart;
                len = param->matchend - param->matchstart;
                }*/
        required_finger_len = 0;
        for(j = 0; j < read_structure->num_segments;j++){
                if(read_structure->type[j] == 'F'){
                        required_finger_len += (int) strlen(read_structure->sequence_matrix[j][0]);
                }
        }
        //KSL_DPRINTF3(("Requiured_len: %d\n",required_finger_len ));

        if(threshold <=  ri->mapq){
                fingerlen = 0;
                for(j = 0; j < len;j++){
                        c1 = mb->label[(int)ri->labels[j+1]];
                        c2 = c1 & 0xFFFF; //which segment
                        c3 = (c1 >> 16) & 0x7FFF; // which HMM in a particular segment....
                        //KSL_DPRINTF3(("%c",   read_structure->type[c2] ));
                        if(read_structure->type[c2] == 'F'){
                                //	required_finger_len += (int) strlen(read_structure->sequence_matrix[c2][0]);
                                fingerlen++;

                                key = (key << 2 )|  (ri->seq[j+offset] & 0x3);
                        }
                        if(read_structure->type[c2] == 'B'){
                                hmm_has_barcode = 1;
                                bar = c3;

                                if(bar == read_structure->numseq_in_segment[c2]-1){
                                        //fprintf(stderr,"EXTRACTING N!!!!\n");
                                        hmm_has_barcode = -1;
                                }
                                mem = c2;
                        }

                        if(read_structure->type[c2] == 'R'){
                                s_pos++;
                                if(!in_read){
                                        in_read= 1;
                                }
                        }else{
                                if(in_read){
                                        if(s_pos < 10){ /* TOO SHORT - CONSTANT FIXME  */
                                                too_short = 1;
                                                break;
                                        }
                                }
                                in_read = 0;
                                s_pos = 0;
                        }
                }
                if(in_read){
                        if(s_pos < 10){ /* TOO SHORT - CONSTANT FIXME  */
                                too_short = 1;
                        }
                }
                //KSL_DPRINTF3(("\n"));
                //KSL_DPRINTF3(("len: %d\n",fingerlen ));
                if(!too_short){
                        if(hmm_has_barcode == -1){
                                ri->read_type  = EXTRACT_FAIL_BAR_FINGER_NOT_FOUND;
                        }else if(hmm_has_barcode && required_finger_len){
                                if(fingerlen == required_finger_len && bar != -1){
                                        RUN(make_extracted_read(mb, read_structure,ri));
                                        ri->barcode =  (mem << 16) |   bar;

                                        //ri->barcode_string = read_structure->sequence_matrix[mem][bar];
                                        if(required_finger_len <= 255){
                                                ri->fingerprint = (key <<  8) | required_finger_len;
                                        }else{
                                                ri->fingerprint = (key <<  8) | 255;
                                        }
                                        ri->read_type = EXTRACT_SUCCESS;
                                }else{
                                        ri->read_type  = EXTRACT_FAIL_BAR_FINGER_NOT_FOUND; // something wrong with the architecture
                                }
                        }else if(hmm_has_barcode){
                                if(bar != -1){
                                        RUN(make_extracted_read(mb, read_structure,ri));
                                        ri->barcode =  (mem << 16) |   bar;
                                        //ri->barcode_string = read_structure->sequence_matrix[mem][bar];
                                        ri->read_type = EXTRACT_SUCCESS;
                                }else{
                                        ri->read_type  = EXTRACT_FAIL_BAR_FINGER_NOT_FOUND;
                                }

                        }else if(required_finger_len){
                                if(fingerlen == required_finger_len){
                                        RUN(make_extracted_read(mb, read_structure,ri));
                                        if(required_finger_len <= 255){
                                                ri->fingerprint = (key <<  8) | required_finger_len;
                                        }else{
                                                ri->fingerprint = (key <<  8) | 255;
                                        }


                                        ri->read_type = EXTRACT_SUCCESS;
                                }else{
                                        ri->read_type  = EXTRACT_FAIL_BAR_FINGER_NOT_FOUND;
                                }
                        }else{
                                RUN(make_extracted_read(mb, read_structure,ri));


                                ri->read_type = EXTRACT_SUCCESS;
                        }
                }else{
                        ri->read_type = EXTRACT_FAIL_READ_TOO_SHORT;
                }
        }else{
                //fprintf(stderr,"UN\n");
                //print_labelled_reads(mb, param,ri);
                ri->read_type = EXTRACT_FAIL_ARCHITECTURE_MISMATCH;

        }

        ri->qual[ri->len] = 0;

        return OK;
ERROR:
        return FAIL;
}

int make_extracted_read(struct model_bag* mb, struct read_structure* read_structure,  struct read_info* ri)
{

        //print_labelled_reads(mb, param,ri);

        int s_pos,j,c1,c2;
        //int multireadread = 0;
        s_pos = 0;
        for(j = 0; j < ri->len;j++){
                c1 = mb->label[(int)ri->labels[j+1]];
                c2 = c1 & 0xFFFF; //which segment

                //fprintf(stderr,"%d	%c	\n", ri->seq[j],param->read_structure->type[c2]);

                if(read_structure->type[c2] == 'R'){
                        //if(multireadread == 0){
                        //	multireadread = 1;
                        //}
                        ri->seq[s_pos] = ri->seq[j];
                        ri->qual[s_pos] = ri->qual[j];
                        s_pos++;
                        //}else if (param->multiread){
                }else{
                        ri->seq[s_pos] = 65; // 65 is the spacer! nucleotides are 0 -5....
                        ri->qual[s_pos] = 65;
                        s_pos++;
                }
        }
        ri->len = s_pos;
        //exit(0);
        return OK;

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
