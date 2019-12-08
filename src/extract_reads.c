#include "extract_reads.h"

#include <omp.h>
#include <math.h>
#include <string.h>

#include "hmm_model_bag.h"
#include "core_hmm_functions.h"


#include "tllogsum.h"

#define READ_CHUNK_SIZE 30000
#define CHUNKS 5

struct assign_struct{
        int** assignment;
        int* num_per_file;
        int num_barcodes;
        int num_files;
        int total;
};

static int alloc_assign_structure(struct assign_struct** assign,int num_files);
static void free_assign_structure(struct assign_struct* as);

static int set_up_assign_structure(struct arch_library* al,struct assign_struct* as);

static int sanity_check_inputs(struct read_info_buffer** rb, int num_files);

//static int run_extract( struct assign_struct* as,  struct read_info_buffer** rb, struct arch_library* al, strucnt seq_stats* si,int i_file,int i_hmm);

static int run_extract( struct assign_struct* as,  struct read_info_buffer** rb, struct arch_library* al, struct seq_stats* si,int i_file,int i_chunk);

int analyze_and_extract_reads(struct assign_struct*as,  struct model_bag* mb, struct read_structure* read_structure,  struct read_info* ri, int i_file, float threshold);


int process_read(struct read_info* ri, int* label,char* type,int* bar_assign, int i_file);
//static int analyze_and_extract_reads(struct model_bag* mb, struct read_structure* read_structure,  struct read_info* ri, int i_file, float threshold);


static int make_extracted_read(struct model_bag* mb, struct read_structure* read_structure,  struct read_info* ri);



int extract_reads(struct arch_library* al, struct seq_stats* si,struct parameters* param)
{
        struct read_info_buffer** rb = NULL;
        struct file_handler** f_hand = NULL;
        struct assign_struct* as = NULL;
        int (*fp)(struct read_info_buffer* rb, struct file_handler* f_handle) = NULL;
        int i,j,c;
        int total_read;

        /* figure out how many barcodes etc */
        RUN(alloc_assign_structure(&as, al->num_file));
        RUN(set_up_assign_structure(al,as));
        as->total = si->ssi[0]->total_num_seq;
        RUN(galloc(&as->assignment, as->total, as->num_barcodes));
        /* not sure if this is required  */
        for(i = 0;i < as->total;i++){
                for(j = 0;j < as->num_barcodes;j++){
                        as->assignment[i][j] = 0;
                }
        }
        //exit(0);
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
                                LOG_MSG("Reading %d chunk %d ->%d", i,j, i + j * param->num_infiles);
                                RUN(fp(rb[i * CHUNKS + j],f_hand[i]));//  param,file,&numseq));
                                total_read += rb[i* CHUNKS + j]->num_seq;
                        }
                        for(j =1; j < CHUNKS;j++){
                                rb[i* CHUNKS + j]->offset = rb[i* CHUNKS + j-1]->offset +  rb[i* CHUNKS + j-1]->num_seq;
                        }
                }
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

#pragma omp parallel default(shared)
#pragma omp for collapse(2) private(i, j)
                for(i = 0; i < param->num_infiles;i++){
                        for(j = 0; j < CHUNKS;j++){
                                run_extract(as, rb,al,si,i,j);
                        }
                }

                //for(j = 0; j < CHUNKS;j++){
                for(c = 0; c < MACRO_MIN(10, READ_CHUNK_SIZE);c++){
                        for(i = 0; i < param->num_infiles;i++){
                                LOG_MSG("FILE: %s",param->infile[i]);
                                fprintf(stdout,"BARCODE:\t");
                                for(j = 0; j < as->num_barcodes;j++){
                                        fprintf(stdout,"%d ", as->assignment[rb[i*CHUNKS]->offset+c][j]);
                                }
                                fprintf(stdout,"\n");
                                print_seq(rb[i*CHUNKS]->ri[c], stdout);
                        }
                                //}
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

int run_extract( struct assign_struct* as,  struct read_info_buffer** rb, struct arch_library* al, struct seq_stats* si,int i_file,int i_chunk)
{
        struct model_bag* mb = NULL;
        struct read_info** ri = NULL;
        int* tmp_bar;
        float pbest = 0.0f;
        float Q = 0.0f;

        int num_seq;
        int c;
        int i;
        int i_hmm;

        int assign_offset;
        int seq_offset;
        assign_offset = as->num_per_file[i_file];


        c = i_file * CHUNKS + i_chunk;


        i_hmm = al->arch_to_read_assignment[i_file];
        mb = init_model_bag(al->read_structure[i_hmm], si->ssi[i_file], i_hmm);
        num_seq = rb[c]->num_seq;
        ri = rb[c]->ri;

        seq_offset = rb[c]->offset;
        LOG_MSG("Working on file: %d  (%d seq) offset: %d",i_file,num_seq,rb[c]->offset);

        //LOG_MSG("Offset = %d",
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
                //    }

                //for(i = 0; i < num_seq;i++){

                ri[i]->bar_prob = 100;

                tmp_bar = as->assignment[seq_offset + i] + assign_offset;
                //print_labelled_reads(mb,data->param ,ri[i]);
                //RUN(extract_reads(mb,data->param,ri[i]));
                //RUN(analyze_and_extract_reads(as,mb, al->read_structure[i_hmm], ri[i], i_file, al->confidence_thresholds[i_file]));
                RUN(process_read(ri[i], mb->label,al->read_structure[i_hmm]->type,tmp_bar,i_file));
        }
        free_model_bag(mb);
        return OK;
ERROR:
        return FAIL;
}

int process_read(struct read_info* ri, int* label,char* type,int* bar_assign, int i_file)
{
        char c;
        int j;
        int bar_segment = 0;
        int bar;
        int segment;
        int hmm_in_segment;
        int c1;
        int umi;
        int umi_len;
        int mem = -1;
        int s_pos = 0;
        int len = ri->len;


        //assign_offset = as->num_per_file[i_file];
        for(j = 0; j < len;j++){

                c1 = label[(int)ri->labels[j+1]];
                segment = c1 & 0xFFFF; //which segment
                hmm_in_segment = (c1 >> 16) & 0x7FFF; // which HMM in a particular segment....
                c = type[segment];
                switch (c) {
                case 'F':
                        umi_len++;

                        umi = (umi << 2 )|  (ri->seq[j] & 0x3);
                        ri->seq[s_pos] = 65; // 65 is the spacer! nucleotides are 0 -5....
                        ri->qual[s_pos] = 65;

                        break;
                case 'B':
                        bar = hmm_in_segment;

                        if(segment != mem){
                                bar_segment++;
                        }
                        bar_assign[ bar_segment] = bar;
                        mem = segment;
                        ri->seq[s_pos] = 65; // 65 is the spacer! nucleotides are 0 -5....
                        ri->qual[s_pos] = 65;

                        break;
                case 'R':
                        ri->seq[s_pos] = ri->seq[j];
                        ri->qual[s_pos] = ri->qual[j];
                        s_pos++;
                        break;
                default:
                        break;
                }
        }
        ri->len = s_pos;


        return OK;
ERROR:
        return FAIL;

}



int analyze_and_extract_reads(struct assign_struct*as,  struct model_bag* mb, struct read_structure* read_structure,  struct read_info* ri, int i_file, float threshold)
{
        int j,c1,c2,c3,key,bar,mem,fingerlen,required_finger_len;

        int s_pos = 0;
        //ret = 0;
        int offset = 0;
        int len;
        int hmm_has_barcode = 0;
        int too_short = 0;
        int in_read = 0;

        int assign_offset;

        assign_offset = as->num_per_file[i_file];

        key = 0;
        bar = -1;
        mem = -1;

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

                                //as->assignment[]
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


int set_up_assign_structure(struct arch_library* al,struct assign_struct* as)
{
        struct read_structure* read_structure = NULL;
        int i,j;
        ASSERT(al != NULL,"No archlib");
        ASSERT(as != NULL,"No assign struct ");

        for(i = 0; i < al->num_file;i++){
                read_structure = al->read_structure[al->arch_to_read_assignment[i]];
                fprintf(stdout,"Read %d: ",i);
                for(j = 0; j < read_structure->num_segments;j++){

                        if(read_structure->type[j]  == 'B'){
                                as->num_per_file[i]++;
                                as->num_barcodes++;
                                fprintf(stdout,"B(%d), ", j);
                        }
                        if(read_structure->type[j]  == 'F'){
                                fprintf(stdout,"F(%d), ", j);
                        }
                }
                fprintf(stdout,"\n");
        }
        /* create offsets  */
        for(i = 1; i < al->num_file;i++){
                as->num_per_file[i] = as->num_per_file[i-1];
        }


        return OK;
ERROR:
        return FAIL;
}


int alloc_assign_structure(struct assign_struct** assign,int num_files)
{

        struct assign_struct* as = NULL;
        int i;

        ASSERT(num_files >= 1,"no infiles");
        MMALLOC(as, sizeof(struct assign_struct));
        as->num_files = num_files;
        as->assignment = NULL;
        as->num_per_file = NULL;
        as->num_barcodes =0;
        MMALLOC(as->num_per_file, sizeof(int) * as->num_files);
        for(i = 0; i < num_files;i++){
                as->num_per_file[i] = 0;
        }
        *assign = as;
        return OK;
ERROR:
        free_assign_structure(as);
        return FAIL;
}

void free_assign_structure(struct assign_struct* as)
{
        if(as){
                if(as->assignment){
                        gfree(as->assignment);
                }
                if(as->num_per_file){
                        MFREE(as->num_per_file);
                }
                MFREE(as);
        }
}
