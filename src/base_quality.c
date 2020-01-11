

#include "base_quality.h"

#include "tllogsum.h"

#include "tldevel.h"
/* test code to incorporate base qualities into alignment scores. */
/* based on:

Frith MC, Wan R, Horton P. Incorporating sequence quality data into alignment improves DNA read mapping. Nucleic acids research. 2010 Jan 27;38(7):e100-.
 */

/* P (y|d) = 1 − 10 −Q y /10 */

struct qsubscore{
        float m[50][16];
        float d;
};

static float fill_qsubscore(float subm[4][4], int t, int q , float p,float n);

int calc_score_matrix(struct qsubscore** mat,  float base_error, float indel_freq)
{
        struct qsubscore* b = NULL;


        float subm[4][4];
        float diag;
        float off;
        float q,n;

        int i,j,c;
        int key;

        init_logsum();

        diag = prob2scaledprob((1.0  - base_error* (1.0- indel_freq)) / 4.0);
        off=  prob2scaledprob( base_error* (1.0- indel_freq)/ 12.0);
        for(i = 0;i < 4;i++){
                for(j = 0;j < 4;j++){
                        if(i ==j){
                                subm[i][j] = diag;
                        }else{
                                subm[i][j] = off;
                        }
                }
        }
        MMALLOC(b, sizeof(struct  qsubscore));
        b->d = diag;
        for(i = 0; i < 50;i++){
                q = prob2scaledprob(1.0 - powf(10.0f, (-1.0 * (float)i) / 10.0f));
                n = prob2scaledprob((powf(10.0f, (-1.0 * (float)i) / 10.0f) )/ 3.0f);
                for(j = 0; j < 4; j++){
                        for(c = 0; c < 4;c++){
                                key = (j << 2) | c;
                                b->m[i][key] = fill_qsubscore(subm,j,c,q,n);
                        }
                }
                /*fprintf(stdout,"%f %f\t", scaledprob2prob(q),scaledprob2prob(n));
                for(j = 0; j < 16; j++){
                        fprintf(stdout,"%f ", scaledprob2prob(b->m[i][j]));
                }
                fprintf(stdout,"\n");*/
        }
        //exit(0);
        *mat = b;
        return OK;
ERROR:
        return FAIL;

}

float get_qsubscore(struct qsubscore* subm, uint8_t a, uint8_t b, uint8_t q)
{

        register uint8_t k;
        if(a == 4 || b == 4){
                return subm->d;
        }
        k = a << 2 | b;

        return subm->m[q][k];
}


float fill_qsubscore(float subm[4][4], int t, int q , float p,float n)
{

        float sum = prob2scaledprob(0.0);

        int i;

        for(i = 0; i < 4; i++){
                if(i == q){
                        sum = logsum(sum, subm[t][i]+ p);
                }else{
                        sum = logsum(sum, subm[t][i]+ n);
                }
        }

        return sum;
}


#ifdef BASEQTEST

#include "tlseqio.h"


#include "tlalphabet.h"
#include "tlrng.h"


struct bq_score{
        float qual[128];

        float nqual[128];

        float m[5][5];
        float diag;
        float off_d;
};

int init_bq_score(struct bq_score** b, float base_error, float indel_freq);
void free_bq_score(struct bq_score*b);

int calc_score_simple(struct bq_score*b, int t, int q , int l, float* ret);

int test_bq(char* filename);
int test_timing(void);
int main(int argc, char *argv[])
{

        char* filename = NULL;
        if(argc == 2){
                filename = argv[1];
                test_bq(filename);
        }else{
                LOG_MSG("Run: %s <>fasta/q file", argv[0]);
        }
        //GGTTTACT
        RUN(test_timing());
        return EXIT_SUCCESS;
ERROR:
        return EXIT_FAILURE;
}

int test_timing(void)
{
        struct bq_score* score = NULL;
        struct qsubscore* m = NULL;
        struct rng_state* rng = NULL;
        const int seq_len = 1000000;
        uint8_t* seq = NULL;
        uint8_t* qual = NULL;
        float a,b;
        int i,j;
        RUN(init_bq_score(&score, 0.05f, 0.01f));

        RUN(calc_score_matrix(&m,0.05f,0.01f));

        rng = init_rng(0);
        MMALLOC(seq, sizeof(uint8_t) * seq_len);
        MMALLOC(qual, sizeof(uint8_t) * seq_len);

        for(i = 0; i < seq_len;i++){
                seq[i] = tl_random_int(rng, 5);
                qual[i] = tl_random_int(rng, 40);
        }
        /* test correctness */

        for(i = 0; i < MACRO_MIN(seq_len,1000);i++){
                for(j = 0; j < 4;j++){
                        RUN(calc_score_simple(score, j, seq[i],qual[i],&a));

                        b = get_qsubscore(m, j,seq[i],qual[i]);
                        ASSERT(a == b,"Scores differ ! %f %f", scaledprob2prob(a), scaledprob2prob(b));
                }


        }
        DECLARE_TIMER(timer);
        START_TIMER(timer);
        for(i = 0; i < seq_len;i++){
                for(j = 0; j < 4;j++){
                        RUN(calc_score_simple(score, j, seq[i],qual[i],&a));
                }
        }
        STOP_TIMER(timer);
        LOG_MSG("Default took %f seconds", GET_TIMING(timer));
        START_TIMER(timer);
        for(i = 0; i < seq_len;i++){
                for(j = 0; j < 4;j++){
                        b = get_qsubscore(m, j,seq[i],qual[i]);
                }
        }
        STOP_TIMER(timer);
        LOG_MSG("Default took %f seconds", GET_TIMING(timer));

        MFREE(seq);
        MFREE(qual);
        free_bq_score(score);
        MFREE(m);
        return OK;
ERROR:
        return FAIL;

}


int test_bq(char* filename)
{
        struct file_handler* f_hand = NULL;
        struct tl_seq_buffer* sb = NULL;
        struct bq_score* score = NULL;
        struct alphabet* a = NULL;
        struct qsubscore* m = NULL;


        float r;

        int i,j,c,l;

        RUN(create_alphabet(&a, 0, TLALPHABET_DEFAULT_DNA));
        RUN(init_bq_score(&score, 0.05f, 0.01f));

        calc_score_matrix(&m,0.05f,0.01f);
        RUN(open_fasta_fastq_file(&f_hand, filename, TLSEQIO_READ));

        RUN(read_fasta_fastq_file(f_hand, &sb, 100));
        for(i = 0; i < 1;i++){

                fprintf(stdout,"%s %s %s\n", sb->sequences[i]->name,sb->sequences[i]->seq,sb->sequences[i]->qual);
                for(j = 0; j < sb->sequences[i]->len;j++){
                        l = (int) sb->sequences[i]->qual[j] - sb->base_quality_offset;
                        /*if(sb->sequences[i]->seq[j] == 'N'){

                                sb->sequences[i]->seq[j] = 'A';
                                }*/


                        fprintf(stdout,"%d %c %c %f\t",i,sb->sequences[i]->seq[j],sb->sequences[i]->qual[j],scaledprob2prob(score->qual[l]));

                        for(c = 0; c < 4;c++){

                                calc_score_simple(score, c, tlalphabet_get_code(a, sb->sequences[i]->seq[j]),l,&r);
                                fprintf(stdout,"%c %f ","ACGT"[c],scaledprob2prob(r));
                        }

                        for(c = 0; c < 4;c++){
                                r = get_qsubscore(m, c, tlalphabet_get_code(a,sb->sequences[i]->seq[j]),l);
                                fprintf(stdout,"%c %f ","ACGT"[c],scaledprob2prob(r));
                        }

                        fprintf(stdout,"\n");

                }
        }
        LOG_MSG("Base Q: %d",sb->base_quality_offset);
        close_seq_file(&f_hand);

        free_tl_seq_buffer(sb);
        free_alphabet(a);
        free_bq_score(score);
        MFREE(m);
        return OK;
ERROR:
        return FAIL;
}


int calc_score_simple(struct bq_score*b, int t, int q , int l, float* ret)
{

        float sum = prob2scaledprob(0.0);

        float p = b->qual[l];
        float n = b->nqual[l];

        int i;
        if(q == 4){
                *ret = b->diag;
                return OK;
        }

        for(i = 0; i < 4; i++){
                if(i == q){
                        sum = logsum(sum, b->m[t][i]+ p);
                }else{
                        sum = logsum(sum, b->m[t][i]+ n);
                }
        }
        *ret = sum;
        return OK;
}


int init_bq_score(struct bq_score** b, float base_error, float indel_freq)
{
        struct bq_score* s = NULL;
        int i,j;
        MMALLOC(s, sizeof(struct bq_score));
        for(i = 0; i < 128;i++){
                s->qual[i] = prob2scaledprob(1.0 - powf(10.0f, (-1.0 * (float)i) / 10.0f));

                s->nqual[i] = prob2scaledprob((powf(10.0f, (-1.0 * (float)i) / 10.0f) )/ 3.0f);


                //fprintf(stdout,"%d %f %f \n",i, scaledprob2prob(qual[i]), scaledprob2prob(nqual[i]));
        }

        for(i = 0;i < 4;i++){
                for(j = 0;j < 4;j++){
                        if(i ==j){
                                s->m[i][j]  = prob2scaledprob((1.0  - base_error* (1.0- indel_freq)) / 4.0);
                        }else{
                                s->m[i][j] = prob2scaledprob( base_error* (1.0- indel_freq)/ 12.0);
                        }
                }
                s->m[i][4] = prob2scaledprob(0.0f);//  prob2scaledprob((1.0  - base_error* (1.0- indel_freq)) / 4.0);
                s->m[4][i] =prob2scaledprob(0.0f);// buffer for start and stop

        }

        s->m[4][4] = prob2scaledprob(0.0f);// buffer for start and stop
        s->diag =  prob2scaledprob((1.0  - base_error* (1.0- indel_freq)) / 4.0);
        s->off_d =  prob2scaledprob( base_error* (1.0- indel_freq)/ 12.0);
        for(i = 0;i < 5;i++){
                for(j = 0;j < 5;j++){
                        fprintf(stdout,"%f ", scaledprob2prob(s->m[i][j]));
                }
                fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
        LOG_MSG("On:%f", scaledprob2prob(s->diag));
        LOG_MSG("Off:%f",scaledprob2prob(s->off_d));
        *b = s;
        return OK;
ERROR:
        return FAIL;
}

void free_bq_score(struct bq_score*b)
{
        if(b){
                MFREE(b);
        }
}

#endif
