#include "extract_reads.h"
#include <omp.h>
#include <math.h>
#include <string.h>


#include "tllogsum.h"

#include "tlseqio.h"

#include "hmm_model_bag.h"
#include "core_hmm_functions.h"
#include "assign_data.h"

#include "filter.h"
//#include "filter_pst.h"
#include "pst.h"


#define READ_CHUNK_SIZE 10000
#define CHUNKS 10

#define MAX_OUTREADNAME 128

struct collect_read{
        char* seq;
        char* qual;
        char* label;
        int* hmm_label;
        int len;
        uint8_t f;
};

static int process_read(struct collect_read* ri, struct read_structure* rs , struct seq_bit_vec* b , int i_file);

static int sanity_check_inputs(struct tl_seq_buffer** rb, int num_files);

static int run_extract( struct assign_struct* as,  struct tl_seq_buffer** rb, struct arch_library* al, struct seq_stats* si,int i_file,int i_chunk);

static int write_all(const struct assign_struct* as, struct tl_seq_buffer** wb, int bam);

int extract_reads(struct arch_library* al, struct seq_stats* si,struct parameters* param,struct rng_state* rng)
{
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

                RUN(run_build_pst(&pst, wb));
//RUN(run_build_pst(&p, sb));
                free_tl_seq_buffer(wb);
                wb = NULL;
        }

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


                RUN(post_process_assign(as));
                START_TIMER(t1);
                STOP_TIMER(t1);
                LOG_MSG("Took %f ",GET_TIMING(t1));
                //LOG_MSG("Write buff: %p",wb);
                RUN(write_all(as,&wb,param->bam));
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
        struct model_bag* mb = NULL;
        struct alphabet* a = NULL;
        struct tl_seq** ri = NULL;
        struct collect_read cs;
        uint8_t* tmp_seq = NULL;
        char* label = NULL;
        float pbest = 0.0f;
        float Q = 0.0f;

        int num_seq;
        int c;
        int i;
        int j;
        int i_hmm;
        int seq_offset;

        a = si->a;

        MMALLOC(tmp_seq, sizeof(uint8_t) * (si->ssi[i_file]->max_seq_len+1));
        MMALLOC(label, sizeof(char) * (si->ssi[i_file]->max_seq_len+1));

        c = i_file * CHUNKS + i_chunk;

        i_hmm = al->arch_to_read_assignment[i_file];
        RUN(init_model_bag(&mb,al->read_structure[i_hmm], si->ssi[i_file], si->a, i_hmm));
        num_seq = rb[c]->num_seq;
        ri = rb[c]->sequences;

        seq_offset = rb[c]->offset;
        for(i = 0; i < num_seq;i++){
                for(j = 0; j < ri[i]->len;j++){
                        tmp_seq[j] = tlalphabet_get_code(a, ri[i]->seq[j]);
                }
                tmp_seq[ri[i]->len] = 0;


                //LOG_MSG("%d f:%d  on %d", tid,i_file,i);
                RUN(backward(mb, tmp_seq,ri[i]->len));

                /* FIXMEEEEE  */
                RUN(forward_max_posterior_decoding(mb, tmp_seq, &label, ri[i]->len));
                RUN(random_score(mb, tmp_seq, ri[i]->len));


                //pbest = ri[i]->mapq;
                //pbest = prob2scaledprob(0.0f);
                //pbest = logsum(pbest, mb->f_score);
                pbest = logsum(mb->bar_score + mb->f_score, mb->r_score);

                pbest = 1.0 - scaledprob2prob(  (mb->bar_score + mb->f_score ) - pbest);

                if(!pbest){
                        Q = 40.0;
                }else if(pbest == 1.0){
                        Q = 0.0;
                }else{
                        Q = -10.0 * log10(pbest) ;
                }

                //ri[i]->mapq = Q;
                //ri[i]->bar_prob = 100;


                        //as->bits[seq_offset+i]->Q[i_file] = Q;
                //}
                //fprintf(stdout,"File:%d %f %f\n",i_file, Q, al->confidence_thresholds[i_file]);
                //as->bits[seq_offset+i]->pass = 0;
                //}
                //tmp_bar = as->assignment[seq_offset + i] + assign_offset;
                //print_labelled_reads(mb,data->param ,ri[i]);
                //RUN(extract_reads(mb,data->param,ri[i]));
                //RUN(analyze_and_extract_reads(as,mb, al->read_structure[i_hmm], ri[i], i_file, al->confidence_thresholds[i_file]));
                //LOG_MSG("i:%d, seqoff = %d",i, seq_offset);
                //LOG_MSG("t: %d Working on file: %d  %d (%d seq) offset: %d",tid, i_file,i,num_seq,rb[c]->offset);
                cs.label = label;
                cs.qual = ri[i]->qual;
                cs.seq = ri[i]->seq;
                cs.len = ri[i]->len;
                cs.hmm_label = mb->label;
                cs.f = 0;
                if(Q < al->confidence_thresholds[i_file]){
                        LOG_MSG("Q: %f thres %f %d", Q, al->confidence_thresholds[i_file], i_file);
                        cs.f = READ_FAILQ;
                }

                RUN(process_read(&cs,al->read_structure[i_hmm], as->bit_vec[seq_offset+i], as->file_index[i_file]));
        }
        free_model_bag(mb);
        MFREE(label);
        MFREE(tmp_seq);
        return OK;
ERROR:
        return FAIL;
}

int process_read(struct collect_read* ri, struct read_structure* rs , struct seq_bit_vec* b , int local_bit_index)
{
        struct seq_bit* sb = NULL;
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

        read_label = ri->label+1;
        label = ri->hmm_label;
        //type = rs->type;
        s_index = 0;
        //assign_offset = as->num_per_file[i_file];
        for(j = 0; j < len;j++){
                c1 = label[(int)read_label[j]];
                segment = c1 & 0xFFFF; //which segment
                hmm_in_segment = (c1 >> 16) & 0x7FFF; // which HMM in a particular segment....
                //c = type[segment];
                c = rs->seg_spec[segment]->extract;
                //LOG_MSG("Decoding: %d read:%d hmmcode:%d  segment:%d seq:%d type: %d",j,read_label[j],label[read_label[j]], segment,hmm_in_segment, c);
                switch (c) {
                case ARCH_ETYPE_APPEND: //case 'F':
                        if(c != old_c){
                                s_index =0;
                                sb = b->bits[local_bit_index];

                                s_name = rs->seg_spec[segment]->name;
                                if(b->append.l){
                                        kputc(' ', &b->append);
                                }

                                kputs(s_name, &b->append);

                                kputs(":Z:", &b->append);


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
                        //kputc(ri->seq[j], &sb->p);

                        if(rs->seg_spec[segment]->num_seq > 1){
                                kputc(rs->seg_spec[segment]->seq[hmm_in_segment][s_index], &b->append);
                        }else{
                                kputc(ri->seq[j], &b->append);
                        }

                        s_index++;
                        break;
                case ARCH_ETYPE_SPLIT:// case 'B':
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
                        kputc(ri->seq[j],&sb->p);
                        kputc(ri->qual[j],&sb->q);
                        break;
                default:
                        break;
                }
                old_c = c;
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
                write_buf->malloc_num = 100000;
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
                                LOG_MSG("Writing %d", write_buf->num_seq);
                                RUN(write_seq_buf(write_buf, f_hand));
                                write_buf->num_seq = 0;

                                //RUN(close_seq_file(&f_hand));
                                break;
                        }

                        /* check if we should write to new file */
                        if(dm[bv->out_file_id[out_read]]->id != file){
                                LOG_MSG("New group: %d at %d", dm[bv->out_file_id[out_read]]->id,out_read);
                                LOG_MSG("New group: %s", dm[bv->out_file_id[out_read]]->out_filename,out_read);
                                LOG_MSG("New group: %p", dm[bv->out_file_id[out_read]]->f_hand,out_read);

                                if(file != -1){
                                        LOG_MSG("Writing %d", write_buf->num_seq);
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
                                LOG_MSG("New group: %p", dm[bv->out_file_id[out_read]]->f_hand,out_read);
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
                                LOG_MSG("Writing %d", write_buf->num_seq);
                                RUN(write_seq_buf(write_buf, f_hand));
                                write_buf->num_seq = 0;
                        }
                }
                if(write_buf->num_seq){
                        LOG_MSG("Writing %d  %p", write_buf->num_seq,f_hand);
                        RUN(write_seq_buf(write_buf, f_hand));
                        write_buf->num_seq = 0;
                }
        }

        STOP_TIMER(t1);
        LOG_MSG("took %f",GET_TIMING(t1));

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


