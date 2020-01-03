#include <omp.h>
#include <math.h>

#include "tldevel.h"

#include "filter.h"
#include "bpm.h"

#include "pst.h"

#include "tlseqio.h"
#include "tlrng.h"
#include "tlalphabet.h"

#include "tlmisc.h"

#define CHUNKS 10

static inline int nuc_to_internal(const char c);

int run_filter_exact(struct assign_struct* as, struct ref* ref, int index, int thres)
{
        struct seq_bit_vec* bv;
        struct seq_bit* sb;

        uint8_t seq_a[256];
        uint8_t* seq_b = NULL;
        uint8_t dist;
        int len_a;
        int len_b;
        char* tmp_seq;
        int i,j;


        bv = as->bit_vec[index];

        for(i = 0; i < as->out_reads;i++){

                sb = bv->bits[as->loc_out_reads[i]];

//ASSERT(sb->type == READ_TYPE, "NO READ TYPE!!! ");
                tmp_seq = sb->p;
                len_a = MACRO_MIN(256, sb->len);
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
                                sb->fail |= READ_FAILR;
                                //goto EARLY;

                        }
                }
                //RUN(scan_read_with_pst(pst,  sb->p , sb->len,&out));
                //fprintf(stdout,"ERR_%03d,%f\n",d,out);
                //fflush(stdout);

        }

        return OK;
}




int run_filter_pst(struct assign_struct* as, struct pst* pst, int index, float thres)

{
        struct seq_bit_vec* bv;
        struct seq_bit* sb;


        int i;
        float out;

        bv = as->bit_vec[index];

        for(i = 0; i < as->out_reads;i++){
                sb = bv->bits[as->loc_out_reads[i]];
                //sb = bv->bits[i];
                RUN(scan_read_with_pst(pst,  sb->p , sb->len,&out));
                if(out >= thres){
                        sb->fail |= READ_FAILP;
                }
        }
        return OK;
ERROR:
        return FAIL;
}

/*int read_reference_sequences(struct ref** r, struct tl_seq_buffer** ref,char* filename)
{
        //struct rng_state* rng = NULL;
        //struct ref* ref = NULL;
        struct tl_seq_buffer* sb_fasta = NULL;
        struct file_handler* f_fasta = NULL;

        ASSERT(filename!= NULL, "No filename");
        ASSERT(my_file_exists(filename), "File %d not found", filename);

        //RUNP(rng = init_rng(seed));

        if(filename){
                RUN(open_fasta_fastq_file(&f_fasta, filename, TLSEQIO_READ));

                RUN(read_fasta_fastq_file(f_fasta, &sb_fasta, 65536));
                RUN(close_seq_file(&f_fasta));
                //RUN(init_ref(&ref, sb_fasta,rng));
                //free_tl_seq_buffer(sb_fasta);
                *ref = sb_fasta;
        }else{
                ERROR_MSG("No reference sequences to map against!");
        }

        *r = ref;
        return OK;
ERROR:
        return FAIL;

}*/

int init_ref(struct ref** r, struct tl_seq_buffer* sb, struct rng_state* rng)
{
        struct ref* ref = NULL;
        uint8_t* tmp_seq = NULL;
        char* org_seq = NULL;
        int i,j;
        ASSERT(sb != NULL, "No seqbuffer");

#if HAVE_AVX2
        set_broadcast_mask();
#endif

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


static inline int nuc_to_internal(const char c)
{
        switch (c) {
        case 'A':
        case 'a':
                return 0;
                break;
        case 'C':
        case 'c':
                return 1;
                break;
        case 'G':
        case 'g':
                return 2;
                break;
        case 'T':
        case 't':
                return 3;
                break;
        case 'N':
        case 'n':
                return 0;
                break;
        default:
                break;
        }
        return -1;
}
