#include "extract_reads.h"
#include <omp.h>
#include <math.h>
#include <string.h>


#include "tllogsum.h"

#include "tlseqio.h"

#include "poahmm.h"
#include "poahmm_structs.h"
#include "init_poahmm.h"

//#include "hmm_model_bag.h"
//#include "core_hmm_functions.h"
#include "assign_data.h"

#include "filter.h"
//#include "filter_pst.h"
#include "pst.h"

#define READ_CHUNK_SIZE 100000
#define CHUNKS 10

#define MAX_OUTREADNAME 128

struct collect_read{
        char* seq;
        char* qual;
        uint32_t* path;
        char* label;
        int* hmm_label;
        int len;
        uint8_t f;
};

//static int process_read(struct collect_read* ri, struct read_structure* rs , struct seq_bit_vec* b , int i_file);
static int process_read(struct collect_read* ri,struct poahmm* poahmm, struct read_structure* rs , struct seq_bit_vec* b , int local_bit_index);

static int sanity_check_inputs(struct tl_seq_buffer** rb, int num_files);

static int run_extract( struct assign_struct* as,  struct tl_seq_buffer** rb, struct arch_library* al, struct seq_stats* si,int i_file,int i_chunk);

static int write_all(const struct assign_struct* as, struct tl_seq_buffer** wb, int bam);

/* just 4 debugging;.. */

static int print_path(struct poahmm* poahmm, uint32_t* path,char* seq);

int extract_reads(struct arch_library* al, struct seq_stats* si,struct parameters* param,struct rng_state* rng)
{
        struct kmer_counts* k = NULL;
        struct tl_seq_buffer**  rb = NULL;
        struct file_handler** f_hand = NULL;
        struct assign_struct* as = NULL;
        struct tl_seq_buffer*  wb = NULL;

        struct file_handler* r_fh = NULL;

        struct ref* ref = NULL;

        struct pst* pst = NULL;

        int i,j,c;
        int total_read;

        DECLARE_TIMER(t1);

        if(param->reference_fasta){
                //RUN(read_reference_sequences(&ref, param->reference_fasta,param->seed));
                RUN(open_fasta_fastq_file(&r_fh, param->reference_fasta, TLSEQIO_READ));

                RUN(read_fasta_fastq_file(r_fh, &wb, 65536));
                RUN(close_seq_file(&r_fh));
                RUN(init_ref(&ref, wb, rng));
                /* not necessary: */
                RUN(alloc_kmer_counts(&k, 12));
                RUN(add_counts(k, wb));



                RUN(run_build_pst(&pst, k));


//RUN(run_build_pst(&p, sb));
                free_kmer_counts(k);
                free_tl_seq_buffer(wb);
                wb = NULL;
        }

        //param->bam = 1;
        /* figure out how many barcodes etc */
        RUN(init_assign_structure(&as, al, param->outfile, CHUNKS* READ_CHUNK_SIZE, param->bam));

        as->block_size = READ_CHUNK_SIZE;
        //RUN(galloc(&as->assignment, as->total, as->num_barcodes));
        /* not sure if this is required  */
        //omp_set_num_threads(param->num_threads);

        /* here I combine: architectures with sequence parameters of individual input files  */
        /* to: check which architecture belongs to which read */
        ASSERT(param != NULL, "no parameters");
        ASSERT(al != NULL, "No arch library");
        ASSERT(si != NULL, "no seq stats");

        //LOG_MSG("Got here");
        MMALLOC(rb, sizeof(struct tl_seq_buffer*) * param->num_infiles * CHUNKS);
        MMALLOC(f_hand, sizeof(struct file_handler*) * param->num_infiles);

        for(i = 0; i < param->num_infiles;i++){
                f_hand[i] = NULL;
                RUN(open_fasta_fastq_file(&f_hand[i], param->infile[i], TLSEQIO_READ));
                //RUN(io_handler(&f_hand[i], param->infile[i]));
        }
        for(i = 0; i < param->num_infiles * CHUNKS;i++){
                rb[i] = NULL;
                //RUN(alloc_read_info_buffer(&rb[i], READ_CHUNK_SIZE));
        }


        while(1){
                /* read everything in  */
                total_read = 0;
                for(i = 0; i < param->num_infiles * CHUNKS;i++){
                        if(rb[i]){
                                rb[i]->offset = 0;
                                rb[i]->num_seq = 0;
                        }
                }
                for(i = 0; i < param->num_infiles;i++){

                        /* set offset of first chunk to be the off ser  */
                        if(rb[i*CHUNKS]){

                                rb[i * CHUNKS]->offset =
                                        rb[i* CHUNKS + CHUNKS-1]->offset
                                        + rb[i* CHUNKS + CHUNKS-1]->num_seq;
                        }
                        for(j = 0; j < CHUNKS;j++){
                                RUN(read_fasta_fastq_file(f_hand[i], &rb[i * CHUNKS + j],READ_CHUNK_SIZE));
                                //RUN(read_fasta_fastq(rb[i * CHUNKS + j],f_hand[i]));
                                //RUN(fp(rb[i * CHUNKS + j],f_hand[i]));//  param,file,&numseq));
                                total_read += rb[i* CHUNKS + j]->num_seq;
                        }
                        for(j =1; j < CHUNKS;j++){
                                if(rb[i* CHUNKS + j]){
                                        rb[i* CHUNKS + j]->offset = rb[i* CHUNKS + j-1]->offset +  rb[i* CHUNKS + j-1]->num_seq;
                                }
                        }
                }
                /* Assign name to assign struct  */
                //for(i = 0; i < param->num_infiles;i++){
                i = 0;
                for(j = 0; j < CHUNKS;j++){
                        for(c = 0;c < rb[j]->num_seq;c++){
                                as->bit_vec[i]->name = rb[j]->sequences[c]->name;
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
                //RUN(sanity_check_inputs(rb,param->num_infiles));
                //fflush(stdout);
                /* extract reads  */

                START_TIMER(t1);
#ifdef HAVE_OPENMP
                LOG_MSG("Run parallel");
#pragma omp parallel default(shared)
#pragma omp for collapse(2) private(i, j)
#endif
                for(i = 0; i < param->num_infiles;i++){
                        for(j = 0; j < CHUNKS;j++){
                                run_extract(as, rb,al,si,i,j);
                        }
                }

                STOP_TIMER(t1);
                LOG_MSG("extract Took %f ",GET_TIMING(t1));
                //RUN(sort_as_by_file_type(as));

                if(ref){
                        START_TIMER(t1);
#ifdef HAVE_OPENMP
#pragma omp parallel default(shared)

#pragma omp for private(i)
#endif
                        for(i = 0; i < as->num_reads;i++){
                                run_filter_exact(as,ref, i, param->filter_error);
                        }
                        STOP_TIMER(t1);
                        LOG_MSG("filter Took %f ",GET_TIMING(t1));
                }

                if(ref){
                        START_TIMER(t1);
#ifdef HAVE_OPENMP
#pragma omp parallel default(shared)
#pragma omp for private(i)
#endif
                        for(i = 0; i < as->num_reads;i++){
                                run_filter_pst(as,pst,i,0.5f);
                        }
                        STOP_TIMER(t1);
                        LOG_MSG("filter Took %f ",GET_TIMING(t1));
                }

                START_TIMER(t1);
                RUN(post_process_assign(as));

                STOP_TIMER(t1);
                LOG_MSG("Processing took: %f ",GET_TIMING(t1));
                //LOG_MSG("Write buff: %p",wb);
                START_TIMER(t1);
                RUN(write_all(as,&wb,param->bam));
                STOP_TIMER(t1);
                LOG_MSG("Write took: %f ",GET_TIMING(t1));

                RUN(reset_assign_structute(as));
                //exit(0);
        }
        /* FIXME */
        if(wb){
                free_tl_seq_buffer(wb);
        }
        if(param->reference_fasta){
                free_ref(&ref);
                free_pst(pst);
        }
        free_assign_structure(as);
        for(i = 0; i < param->num_infiles* CHUNKS;i++){
                if(rb[i]){
                        free_tl_seq_buffer(rb[i]);
                }
        }
        MFREE(rb);

        for(i = 0; i < param->num_infiles;i++){
                RUN(close_seq_file(&f_hand[i]));
        }
        MFREE(f_hand);

        return OK;
ERROR:
        return FAIL;
}

int run_extract( struct assign_struct* as,  struct tl_seq_buffer** rb, struct arch_library* al, struct seq_stats* si,int i_file,int i_chunk)
{
        struct global_poahmm_param* p = NULL;
        struct poahmm* poahmm = NULL;
        //struct model_bag* mb = NULL;
        struct alphabet* a = NULL;
        struct tl_seq** ri = NULL;
        struct collect_read cs;
        uint32_t* path = NULL;
        uint8_t* tmp_seq = NULL;
        uint8_t* tmp_qual = NULL;
        //char* label = NULL;

        float pbest = 0.0f;
        float Q = 0.0f;

        int num_seq;
        int c;
        int i;
        int j;
        int i_hmm;
        int base_q_offset = 0;
        int seq_offset;


        a = si->a;


        c = i_file * CHUNKS + i_chunk;

        i_hmm = al->arch_to_read_assignment[i_file];
        //RUN(init_model_bag(&mb,al->read_structure[i_hmm], si->ssi[i_file], si->a, i_hmm));

        MMALLOC(p, sizeof(struct global_poahmm_param));
        p->min_seq_len = si->ssi[i_file]->average_length;
        p->max_seq_len = MACRO_MAX(rb[c]->max_len, si->ssi[i_file]->max_seq_len);
        p->base_error = 0.05f;
        p->indel_freq = 0.1f;
        for(i =0; i < 5;i++){
                p->back[i] = si->ssi[i_file]->background[i];
        }


        RUN(poahmm_from_read_structure(&poahmm, p, al->read_structure[i_hmm],  si->a));


        num_seq = rb[c]->num_seq;
        ri = rb[c]->sequences;
        base_q_offset = rb[c]->base_quality_offset;
        seq_offset = rb[c]->offset;

        MMALLOC(tmp_seq, sizeof(uint8_t) * (si->ssi[i_file]->max_seq_len+1));
        MMALLOC(tmp_qual, sizeof(uint8_t) * (si->ssi[i_file]->max_seq_len+1));

        MMALLOC(path, sizeof(uint32_t) *(si->ssi[i_file]->max_seq_len+ poahmm->num_nodes +2));
        for(i = 0; i < num_seq;i++){
                for(j = 0; j < ri[i]->len;j++){
                        tmp_seq[j] = tlalphabet_get_code(a, ri[i]->seq[j]);
                }
                tmp_seq[ri[i]->len] = 0;
                for(j = 0; j < ri[i]->len;j++){
                        tmp_qual[j] = ri[i]->qual[j] - base_q_offset;
                }

                RUN(viterbi_poahmm_banded(poahmm, tmp_seq, tmp_qual,  ri[i]->len, path, 0));

                pbest = logsum(poahmm->f_score, poahmm->random_scores[ri[i]->len]);
                pbest = 1.0 - scaledprob2prob(poahmm->f_score - pbest);

                if(!pbest){
                        Q = 40.0;
                }else if(pbest == 1.0){
                        Q = 0.0;
                }else{
                        Q = -10.0 * log10(pbest) ;
                }



                cs.path = path;
                //cs.label = label;
                cs.qual = ri[i]->qual;
                cs.seq = ri[i]->seq;
                cs.len = ri[i]->len;
                //cs.hmm_label = mb->label;
                cs.f = 0;
                if(Q < al->confidence_thresholds[i_file]){
                        //LOG_MSG("Q: %f thres %f %d", Q, al->confidence_thresholds[i_file], i_file);
                        cs.f = READ_FAILQ;
                }

                RUN(process_read(&cs,poahmm, al->read_structure[i_hmm], as->bit_vec[seq_offset+i], as->file_index[i_file]));
        }
        //free_model_bag(mb);
        free_poahmm(poahmm);
        MFREE(p);
        MFREE(path);

        MFREE(tmp_seq);
        MFREE(tmp_qual);

        return OK;
ERROR:
        return FAIL;
}

int process_read(struct collect_read* ri,struct poahmm* poahmm, struct read_structure* rs , struct seq_bit_vec* b , int local_bit_index)
{
        struct seq_bit* sb = NULL;
        uint32_t* path = NULL;
        //char* type;
        char* read_label;
        int* label;
        uint8_t c;
        int j;
        int s_index;
        //int s_len;
        char* s_name;
        int segment;
        int hmm_in_segment;
        int c1;
        int len = ri->len;
        uint8_t old_c = 255;
        int read = 0;
        //      int local_bit_index;
        int seq_pos;
        int node_pos;
        read_label = ri->label+1;
        label = ri->hmm_label;

        path = ri->path;

        //print_path(poahmm, path, ri->seq);
        //type = rs->type;
        s_index = 0;
        //assign_offset = as->num_per_file[i_file];
        for(j = 1; j < ri->path[0];j++){
                seq_pos = (path[j] >> 16u );
                node_pos = path[j] & 0xFFFFu;
                //if(node_pos == 0xFFFFu){

                //}
                if(node_pos!= 0xFFFFu){
                        segment = poahmm->nodes[node_pos]->segment;
                        hmm_in_segment =  poahmm->nodes[node_pos]->alt;
                        c = poahmm->nodes[node_pos]->type;

                        //LOG_MSG("decode-ing: %d s:%c s:%d h:%d type:%d", j, ri->seq[seq_pos], segment,hmm_in_segment,c);
                        //c1 = label[(int)read_label[j]];
                        //segment = c1 & 0xFFFF; //which segment
                        //hmm_in_segment = (c1 >> 16) & 0x7FFF; // which HMM in a particular segment....
                        //c = type[segment];
                        //c = rs->seg_spec[segment]->extract;
                        //LOG_MSG("Decoding: %d read:%d hmmcode:%d  segment:%d seq:%d type: %d",j,read_label[j],label[read_label[j]], segment,hmm_in_segment, c);
                        switch (c) {
                        case ARCH_ETYPE_APPEND_CORRECT:
                                if(c != old_c){
                                        //s_index = 0;
                                        sb = b->bits[local_bit_index];
                                        //ASSERT(i_file == sb->file, "Oh dear: want %d got %d",i_file,sb->file);
                                        //sb->file = i_file;
                                        //sb->code = (char)(read + 33);
                                        sb->type = ARCH_ETYPE_APPEND_CORRECT;
                                        sb->fail = ri->f;
                                        local_bit_index++;
                                }
                                kputc(ri->seq[seq_pos],&sb->p);
                                kputc(ri->qual[seq_pos],&sb->q);
                                break;
                        case ARCH_ETYPE_APPEND:
                                if(c != old_c){
                                        s_index =0;
                                        sb = b->bits[local_bit_index];

                                        /*s_name = rs->seg_spec[segment]->name;
                                        if(b->append.l){
                                                kputc(' ', &b->append);
                                        }

                                        kputs(s_name, &b->append);

                                        kputs(":Z:", &b->append);*/

                                        //sb->len = 0;
                                        //ASSERT(i_file == sb->file, "Oh dear: want %d got %d",i_file,sb->file);
                                        //sb->file = i_file;
                                        sb->code = (char) ( hmm_in_segment + 33);
                                        sb->type = ARCH_ETYPE_APPEND;
                                        //sb->p = ri->seq+j;
                                        sb->fail = ri->f;
                                        if(hmm_in_segment == 0 && rs->seg_spec[segment]->num_seq > 1){
                                                //LOG_MSG("Read %s failrd", b->name);
                                                sb->fail |= READ_NBAR;
                                        }

                                        local_bit_index++;
                                }
                                if(rs->seg_spec[segment]->num_seq > 1){
                                        kputc(rs->seg_spec[segment]->seq[hmm_in_segment][s_index], &sb->p);
                                }else{
                                        kputc(ri->seq[seq_pos], &sb->p);
                                }
                                s_index++;
                                break;
                        case ARCH_ETYPE_SPLIT:
                                if(c != old_c){
                                        sb = b->bits[local_bit_index];

                                        //ASSERT(i_file == sb->file, "Oh dear: want %d got %d",i_file,sb->file);
                                        //sb->file = i_file;
                                        sb->type = ARCH_ETYPE_SPLIT;
                                        sb->code = (char) ( hmm_in_segment + 33);


                                        //kputs(rs->seg_spec[segment]->seq[hmm_in_segment], &sb->p);
                                        //sb->p = rs->seg_spec[segment]->seq[hmm_in_segment];
                                        //rs->sequence_matrix[segment][hmm_in_segment];
                                        //sb->len = rs->seg_spec[segment]->max_len;// should be identical to min_lena rs->segment_length[segment];
                                        sb->fail = ri->f;
                                        if(hmm_in_segment == 0){
                                                //LOG_MSG("Read %s failrd", b->name);
                                                sb->fail |= READ_NBAR;
                                        }
                                        local_bit_index++;
                                }
                                break;
                        case ARCH_ETYPE_EXTRACT:// case 'R':
                                if(c != old_c){
                                        s_index = 0;
                                        sb = b->bits[local_bit_index];
                                        //ASSERT(i_file == sb->file, "Oh dear: want %d got %d",i_file,sb->file);
                                        //sb->file = i_file;
                                        //sb->code = (char)(read + 33);
                                        sb->type = ARCH_ETYPE_EXTRACT;
                                        sb->fail = ri->f;
                                        local_bit_index++;
                                        read++;
                                }
                                kputc(ri->seq[seq_pos],&sb->p);
                                kputc(ri->qual[seq_pos],&sb->q);
                                break;
                        default:
                                break;
                        }
                        old_c = c;
                }
        }
        //b->append[b->a_len] = 0;
        //fprintf(stdout,"%s\n%s\n",sb->p,sb->q);
        //exit(0);
        //ri->len = s_pos;
        return OK;
ERROR:
        return FAIL;
}

int write_all(const struct assign_struct* as, struct tl_seq_buffer** wb, int bam)
{
        struct demux_struct** dm;
        struct demux_struct* tmp_ptr = NULL;
        struct seq_bit_vec* bv;
        struct seq_bit*sb;
        struct tl_seq_buffer* write_buf;
        struct file_handler* f_hand;

        int out_read,i;
        int file = -1;

        f_hand = NULL;

        if(*wb){
                write_buf = *wb;
        }else{
                write_buf = NULL;
                MMALLOC(write_buf,sizeof(struct tl_seq_buffer));
                write_buf->malloc_num = 1000000;
                write_buf->max_len = 0;
                write_buf->num_seq = 0;
                write_buf->L = 0;
                write_buf->is_fastq = 1;
                write_buf->offset = 0;
                write_buf->sequences = NULL;

                MMALLOC(write_buf->sequences, sizeof(struct tl_seq*) * write_buf->malloc_num);
                for(i = 0; i < write_buf->malloc_num;i++){
                        write_buf->sequences[i] = NULL;
                        MMALLOC(write_buf->sequences[i], sizeof(struct tl_seq));
                        write_buf->sequences[i]->name = NULL;
                        write_buf->sequences[i]->seq = NULL;
                        write_buf->sequences[i]->qual = NULL;
                        write_buf->sequences[i]->aux = NULL;
                        //MMALLOC(write_buf->sequences[i]->name, sizeof(char) * MAX_OUTREADNAME);
                }
                //RUN(alloc_tl_seq_buffer(&write_buf, 100000));

        }

        DECLARE_TIMER(t1);

        dm = (struct demux_struct**)as->demux_names->data_nodes;
//EXTERN int write_fasta_fastq(struct tl_seq_buffer* sb, struct file_handler* fh);
        START_TIMER(t1);

        for(out_read = 0; out_read < as->out_reads;out_read++){
                //LOG_MSG("OUT:%d",out_read);
                file = -1;
                for(i = 0; i < as->num_reads ;i++){


                        bv = as->bit_vec[i];
                        //fprintf(stdout,"outread: %d id: %d\n",out_read, dm[bv->out_file_id[out_read]]->id);
                        //LOG_MSG("Writing %d fail: %d",i, bv->fail);
                        //bv = as->bit_vec[as->loc_out_reads[i]];
                        if(bv->fail){
                                //LOG_MSG("Writing %d", write_buf->num_seq);
                                RUN(write_seq_buf(write_buf, f_hand));
                                write_buf->num_seq = 0;

                                //RUN(close_seq_file(&f_hand));
                                break;
                        }

                        /* check if we should write to new file */
                        if(dm[bv->out_file_id[out_read]]->id != file){
                                //LOG_MSG("New group: %d at %d", dm[bv->out_file_id[out_read]]->id,out_read);
                                //LOG_MSG("New group: %s", dm[bv->out_file_id[out_read]]->out_filename,out_read);
                                //LOG_MSG("New group: %p", dm[bv->out_file_id[out_read]]->f_hand,out_read);

                                if(file != -1){
                                        //      LOG_MSG("Writing %d", write_buf->num_seq);
                                        RUN(write_seq_buf(write_buf, f_hand));
                                        write_buf->num_seq = 0;
                                        //RUN(close_seq_file(&f_hand));
                                }

                                if(!dm[bv->out_file_id[out_read]]->f_hand){
                                        if(bam){
                                                LOG_MSG("Opening BAM : %s" , dm[bv->out_file_id[out_read]]->out_filename);
                                                RUN(open_sam_bam(&dm[bv->out_file_id[out_read]]->f_hand, dm[bv->out_file_id[out_read]]->out_filename, TLSEQIO_WRITE));
                                        }else{
                                                RUN(open_fasta_fastq_file(&dm[bv->out_file_id[out_read]]->f_hand, dm[bv->out_file_id[out_read]]->out_filename, TLSEQIO_WRITE));
                                        }
                                }
                                f_hand= dm[bv->out_file_id[out_read]]->f_hand;
                                file = dm[bv->out_file_id[out_read]]->id;
                                //LOG_MSG("New group: %p", dm[bv->out_file_id[out_read]]->f_hand,out_read);
                        }

                        sb = bv->bits[as->loc_out_reads[out_read]];
                        write_buf->sequences[write_buf->num_seq]->seq = sb->p.s;//  sb->p;
                        write_buf->sequences[write_buf->num_seq]->qual = sb->q.s;
                        write_buf->sequences[write_buf->num_seq]->len = sb->p.l;
                        write_buf->sequences[write_buf->num_seq]->name = bv->name;
                        //fprintf(stdout,"%s\n", bv->append.s);
                        write_buf->sequences[write_buf->num_seq]->aux = bv->append.s;
                        /*if(out_read == 1){
                        LOG_MSG("");
                        LOG_MSG("%s ", sb->p.s);
                        LOG_MSG("%s ", sb->q.s);
                        LOG_MSG("%s ", bv->name);
                        LOG_MSG("%s ", bv->append.s);
                        }*/
                        write_buf->num_seq++;

                        if(write_buf->num_seq == write_buf->malloc_num){
                                //LOG_MSG("Writing %d", write_buf->num_seq);
                                RUN(write_seq_buf(write_buf, f_hand));
                                write_buf->num_seq = 0;
                        }
                }
                if(write_buf->num_seq){
                        //LOG_MSG("Writing %d  %p", write_buf->num_seq,f_hand);
                        RUN(write_seq_buf(write_buf, f_hand));
                        write_buf->num_seq = 0;
                }
        }

        STOP_TIMER(t1);
        //LOG_MSG("took %f",GET_TIMING(t1));

        write_buf->num_seq = 0;
        for(i = 0; i < write_buf->malloc_num;i++){
                write_buf->sequences[i]->seq = NULL;
                write_buf->sequences[i]->qual = NULL;
                write_buf->sequences[i]->aux = NULL;
                write_buf->sequences[i]->name = NULL;

        }

        *wb = write_buf;
        return OK;
ERROR:
        return FAIL;

}


int sanity_check_inputs(struct tl_seq_buffer** rb, int num_infiles)
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
                                        //if(compare_read_names(rb[i*CHUNKS + g]->ri[c]->name,rb[j*CHUNKS + g]->ri[c]->name)){
                                        //ERROR_MSG("Files seem to contain reads in different order:\n%s\n%s\n",rb[i*CHUNKS + g]->ri[c]->name,rb[j*CHUNKS + g]->ri[c]->name);
                                                //}
                                }
                        }
                }
        }
        return OK;
ERROR:
        return FAIL;
}





int print_path(struct poahmm* poahmm, uint32_t* path,char* seq)
{
        uint32_t i;
        char alphabet[5] = "ACGTN";
        char etype[9] = "_EASIPLRC";

        uint32_t seq_pos;
        uint32_t node_pos;
        fprintf(stdout, "PATH:\n");

        for(i = 1; i < path[0];i++){
                seq_pos = path[i] >> 16u;
                node_pos = path[i] & 0xFFFFu;
                //fprintf(stdout,"Position %d: %d %d\n",i,seq_pos,node_pos);
                if(seq_pos != 0xFFFFu){
                        fprintf(stdout, " %3c", seq[seq_pos]);
                }else{
                        fprintf(stdout, " %3c", '-');
                }
        }
        fprintf(stdout,"\n");
        for(i = 1; i < path[0];i++){
                seq_pos = path[i] >> 16u;
                node_pos = path[i] & 0xFFFFu;

                if(node_pos!= 0xFFFFu){
                        fprintf(stdout, " %3c", alphabet[poahmm->nodes[node_pos]->nuc]);
                }else{
                        fprintf(stdout, " %3c", '-');
                }
        }
        fprintf(stdout,"\n");
        for(i = 1; i < path[0];i++){
                seq_pos = path[i] >> 16u;
                node_pos = path[i] & 0xFFFFu;

                if(node_pos!= 0xFFFFu){
                        fprintf(stdout, " %3c", etype[poahmm->nodes[node_pos]->type]);
                }else{
                        fprintf(stdout, " %3c", '-');
                }
        }
        fprintf(stdout,"\n");
}
