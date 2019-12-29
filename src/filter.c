#include <omp.h>
#include <math.h>

#include "tldevel.h"

#include "filter.h"
#include "bpm.h"


#include "tlseqio.h"
#include "tlrng.h"
#include "tlalphabet.h"

#include "tlmisc.h"

#define CHUNKS 10


static int init_ref(struct ref** r, struct tl_seq_buffer* sb,struct rng_state* rng);


int run_filter(struct assign_struct* as, struct ref* ref, int index, int thres)
{
        struct seq_bit_vec* bv;
        struct seq_bit* sb;

        uint8_t* seq_a = NULL;
        uint8_t* seq_b = NULL;
        uint8_t dist;
        int len_a;
        int len_b;
        char* tmp_seq;
        int i,j;

        MMALLOC(seq_a, sizeof(uint8_t)* as->max_seq_len);

        bv = as->bit_vec[index];
        for(i = 0; i < as->out_reads;i++){
                sb = bv->bits[i];
                //ASSERT(sb->type == READ_TYPE, "NO READ TYPE!!! ");
                tmp_seq = sb->p;
                len_a = sb->len;
                //fprintf(stdout,"%s\n", tmp_seq);
                for(j = 0; j < len_a;j++){
                        seq_a[j] = tlalphabet_get_code(ref->a, tmp_seq[j]);
                }
                for(j = 0; j < ref->num_seq;j++){
                        seq_b = ref->seq[j];
                        len_b = ref->len[j];
                        if(len_a > len_b){
                                dist = bpm_256(seq_a, seq_b, len_a, len_b);
                        }else{
                                dist = bpm_256(seq_b, seq_a, len_b, len_a);
                        }
                        //fprintf(stdout,"%d ", dist);

                        if(dist < thres){
                                //bv->Q[sb->file] = bv->Q[sb->file] - 1.0f;
#pragma omp critical
                                {
                                        bv->fail |= 2;
                                        ref->hits[j]++;
                                }
                                goto EARLY;

                        }
                }

        }
EARLY:

        MFREE(seq_a);
        return OK;
ERROR:
        return FAIL;
}

int read_reference_sequences(struct ref** r, char* filename,int seed)
{
        struct rng_state* rng = NULL;
        struct ref* ref = NULL;
        struct tl_seq_buffer* sb_fasta = NULL;
        struct file_handler* f_fasta = NULL;

        ASSERT(filename!= NULL, "No filename");
        ASSERT(my_file_exists(filename), "File %d not found", filename);

        RUNP(rng = init_rng(seed));

        if(filename){
                RUN(open_fasta_fastq_file(&f_fasta, filename, TLSEQIO_READ));

                RUN(read_fasta_fastq_file(f_fasta, &sb_fasta, 65536));
                RUN(close_seq_file(&f_fasta));

                /*for(i = 0; i < sb_fasta->num_seq;i++){
                        fprintf(stdout,">%s\n", sb_fasta->sequences[i]->name);
                        for(j = 0; j < sb_fasta->sequences[i]->len;j++){
                                fprintf(stdout,"%c", sb_fasta->sequences[i]->seq[j]);
                        }
                        fprintf(stdout,"\n");
                }
                LOG_MSG("Max_len: %d", sb_fasta->max_len);*/

                RUN(init_ref(&ref, sb_fasta,rng));
                free_tl_seq_buffer(sb_fasta);

        }else{
                ERROR_MSG("No reference sequences to map against!");
        }

#if HAVE_AVX2
        set_broadcast_mask();
#endif

        *r = ref;
        return OK;
ERROR:
        return FAIL;

}

//static int run_filter(struct read_info_buffer** rb, struct read_info_buffer* con, int i_file, int  i_chunk);



/*

int filter.




if(len_a > len_b){
                dist = bpm_256(seq_a, seq_b, len_a, len_b);
        }else{
                dist = bpm_256(seq_b, seq_a, len_b, len_a);
        }

*/



int init_ref(struct ref** r, struct tl_seq_buffer* sb, struct rng_state* rng)
{
        struct ref* ref = NULL;
        uint8_t* tmp_seq = NULL;
        char* org_seq = NULL;
        int i,j;
        ASSERT(sb != NULL, "No seqbuffer");

        MMALLOC(ref, sizeof(struct ref));
        ref->a = NULL;
        ref->seq = NULL;
        ref->num_seq = sb->num_seq;
        ref->len = NULL;
        ref->hits = NULL;

        RUN(create_alphabet(&ref->a, rng, TLALPHABET_NOAMBIGUOUS_DNA));

        MMALLOC(ref->len, sizeof(int) * ref->num_seq);
        MMALLOC(ref->seq, sizeof(uint8_t*) * ref->num_seq);
        MMALLOC(ref->hits, sizeof(int) * ref->num_seq);

        for(i = 0; i < ref->num_seq;i++){
                ref->len[i] = sb->sequences[i]->len;
                ref->seq[i] = NULL;
                ref->hits[i] = 0;
                MMALLOC(ref->seq[i], sizeof(uint8_t) * ref->len[i]);
                tmp_seq = ref->seq[i];
                org_seq = sb->sequences[i]->seq;
                for(j = 0; j < ref->len[i];j++){
                        tmp_seq[j] = tlalphabet_get_code(ref->a, org_seq[j]);
                        //fprintf(stdout,"%d", tmp_seq[j]);
                }
                //fprintf(stdout,"\n");
        }

        *r = ref;
        return OK;
ERROR:
        return FAIL;
}

int free_ref(struct ref** r)
{
        struct ref* ref = NULL;
        int i;

        ref = *r;
        if(ref){

                for(i = 0; i < ref->num_seq;i++){
                        MFREE(ref->seq[i]);
                }
                MFREE(ref->seq);
                MFREE(ref->len);
                MFREE(ref->hits);
                free_alphabet(ref->a);
                MFREE(ref);
                ref = NULL;

        }
        *r = ref;
        return OK;
}
