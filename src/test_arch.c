#include <math.h>

#include "tllogsum.h"
#include "tlseqio.h"
#include "tlalphabet.h"

#include "test_arch.h"

#include "poahmm.h"
#include "poahmm_structs.h"
#include "init_poahmm.h"


#include "lpst.h"

#include "hmm_model_bag.h"
//#include "io.h"
#include "core_hmm_functions.h"

#include <omp.h>

#define NUM_TEST_SEQ 100

static int test_arch(struct tl_seq_buffer** rb, struct arch_library* al, struct seq_stats* si,int i_file,int i_hmm);

int test_architectures(struct arch_library* al, struct seq_stats* si, struct parameters* param)
{
        struct arch_bag* ab = NULL;
        struct file_handler* f_hand = NULL;

        struct tl_seq_buffer** rb = NULL;
        //struct tl_seq** ri = NULL;

        struct alphabet* a = NULL;

        RUN(create_alphabet(&a, 0, TLALPHABET_DEFAULT_DNA));

        float sum;
        int best;
        int i,j;
        //int (*fp)(struct read_info_buffer* rb, struct file_handler* f_handle) = NULL;

        /* alloc matrix for posteriors  */
        al->arch_posteriors = NULL;

        MMALLOC(al->arch_posteriors, sizeof(float*) * al->num_arch);
        for(i = 0; i < al->num_arch;i++){
                al->arch_posteriors[i] = NULL;
                MMALLOC(al->arch_posteriors[i], sizeof(float) * param->num_infiles);
                for(j = 0; j < param->num_infiles;j++){
                        al->arch_posteriors[i][j] = prob2scaledprob(0.0);
                }
        }

        MMALLOC(al->arch_to_read_assignment, sizeof(int) * param->num_infiles);

        //omp_set_num_threads(param->num_threads);

        /* here I combine: architectures with sequence parameters of individual input files  */
        /* to: check which architecture belongs to which read */

        ASSERT(al != NULL, "No arch library");
        ASSERT(param != NULL, "No parameters");
        ASSERT(si != NULL, "no seqstats");

        MMALLOC(ab, sizeof(struct arch_bag));
        ab->archs = 0;
        ab->arch_posterior = 0;

        MMALLOC(ab->archs,sizeof(struct model_bag*) * al->num_arch);
        MMALLOC(ab->arch_posterior,sizeof(float) *  al->num_arch);
        ab->num_arch = al->num_arch; /* I'm not going to miss hmms in this process */

        MMALLOC(rb, sizeof(struct tl_seq_buffer*) * param->num_infiles);
        for(i = 0; i < param->num_infiles;i++){
                rb[i] = NULL;
                RUN(open_fasta_fastq_file(&f_hand, param->infile[i], TLSEQIO_READ));
                RUN(read_fasta_fastq_file(f_hand, &rb[i],NUM_TEST_SEQ));
                RUN(close_seq_file(&f_hand));
        }


#ifdef HAVE_OPENMP
#pragma omp parallel default(shared)
#pragma omp for collapse(2) private(i, j)
#endif
        for(i = 0; i < param->num_infiles;i++){
                for(j = 0; j < al->num_arch;j++){
                        test_arch(rb,al,si,i,j);
                }
        }



//LOG_MSG("Testing %d models,", al->num_arch);

/* convert posteriors into Phred scaled quality values */

        for(i = 0; i < param->num_infiles;i++){
                sum = prob2scaledprob(0.0f);
                for(j = 0; j < al->num_arch;j++){
                        //fprintf(stdout,"%f ",al->arch_posteriors[j][i]);
                        sum = logsum(sum,al->arch_posteriors[j][i]);
                }

                fprintf(stdout,"FILE: %s\t",param->infile[i]);
                for(j = 0; j < al->num_arch;j++){

                        al->arch_posteriors[j][i] = 1.0 - scaledprob2prob(al->arch_posteriors[j][i] - sum);
                        if(al->arch_posteriors[j][i] < 0.00001f){
                                al->arch_posteriors[j][i] = 0.00001f;
                        }
                        al->arch_posteriors[j][i] = -10.0f * log10f(al->arch_posteriors[j][i]);
                        fprintf(stdout,"%f ",al->arch_posteriors[j][i]);
                }
                fprintf(stdout,"\n");
        }


        for(i = 0; i < param->num_infiles;i++){
                al->arch_to_read_assignment[i] = -1;
                sum = -1.0;
                best = -1;
                for(j = 0; j < al->num_arch;j++){
                        //LOG_MSG("%s" , al->spec_line[j]);
                        if(al->arch_posteriors[j][i] >= sum){

                                best =j;
                                sum = al->arch_posteriors[j][i];

                        }
                }
                ASSERT(best != -1,"No best arch found");
                LOG_MSG("Best architecture for %s is (Q = %f):", param->infile[i], al->arch_posteriors[best][i]);
                LOG_MSG("%s" , al->spec_line[best]);
                if(al->arch_posteriors[best][i] < 10.0f){
                        WARNING_MSG("Two or more architectures have a high probability of matching %s", param->infile[i]);
                        WARNING_MSG("Attempting to select best arch using probabilistic suffix trees.");
                        sum = prob2scaledprob(0.0);
                        for(j = 0; j < al->num_arch;j++){
                                if(al->arch_posteriors[j][i] < 10.0f){
                                        lpst_score_read(al->read_structure[j], rb[i],si->ssi[i], &al->arch_posteriors[j][i]);
                                        //LOG_MSG("%s Re-scored as: %f", al->spec_line[j], al->arch_posteriors[j][i]);
                                        sum = logsum(sum, al->arch_posteriors[j][i]);
                                }else{
                                        al->arch_posteriors[j][i] = prob2scaledprob(0.0);
                                }
                        }
                        for(j = 0; j < al->num_arch;j++){
                                al->arch_posteriors[j][i] = 1.0 - scaledprob2prob(al->arch_posteriors[j][i] - sum);
                                if(al->arch_posteriors[j][i] < 0.00001f){
                                        al->arch_posteriors[j][i] = 0.00001f;
                                }
                                al->arch_posteriors[j][i] = -10.0f * log10f(al->arch_posteriors[j][i]);
                                //fprintf(stdout,"%f ",al->arch_posteriors[j][i]);
                        }
                        //fprintf(stdout,"\n");
                        sum = -1.0;
                        best = -1;
                        for(j = 0; j < al->num_arch;j++){
                                //LOG_MSG("%s" , al->spec_line[j]);
                                if(al->arch_posteriors[j][i] >= sum){

                                        best =j;
                                        sum = al->arch_posteriors[j][i];

                                }
                        }
                        ASSERT(best != -1,"No best arch found");
                        if(al->arch_posteriors[best][i] < 10.0f){
                                ERROR_MSG("Giving up!");
                        }
                        LOG_MSG("Best architecture for %s is (Q = %f):", param->infile[i], al->arch_posteriors[best][i]);
                        //score_pst(p, test_seq, test_len, &out);
                }

                al->arch_to_read_assignment[i] = best;
                al->arch_posteriors[0][i] = al->arch_posteriors[best][i];
        }

        al->num_file = param->num_infiles;

        for(i = 0; i < param->num_infiles;i++){
                free_tl_seq_buffer(rb[i]);
                //free_read_info_buffer(rb[i]);
        }
        MFREE(rb);

        MFREE(ab->arch_posterior);
        MFREE(ab->archs);
        MFREE(ab);
        free_alphabet(a);
        return OK;
ERROR:
        for(i = 0; i < param->num_infiles;i++){
                free_tl_seq_buffer(rb[i]);
                //free_read_info_buffer(rb[i]);
        }
        MFREE(rb);

        MFREE(ab->arch_posterior);
        MFREE(ab->archs);
        MFREE(ab);
        free_alphabet(a);
        return FAIL;
}

int test_arch(struct tl_seq_buffer** rb, struct arch_library* al, struct seq_stats* si,int i_file,int i_hmm)
{
        //fprintf(stdout,"Thread %d working on file %d; hmm %d\n",omp_get_thread_num(),i_file,i_hmm);
        //printf ("num_thds=%d, max_thds=%d\n",omp_get_thread_num(),omp_get_max_threads());
        struct poahmm* poahmm = NULL;

        struct global_poahmm_param* p = NULL;
        //struct model_bag* mb = NULL;
        struct tl_seq** ri = NULL;
        struct alphabet* a = NULL;

        uint8_t* tmp_seq = NULL;
        uint8_t* tmp_qual = NULL;
        int num_seq;
        int i,j;
        int base_q_offset = 0;
        float score = prob2scaledprob(1.0);

        base_q_offset = rb[i_file]->base_quality_offset;

        MMALLOC(p, sizeof(struct global_poahmm_param));
        p->min_seq_len = si->ssi[i_file]->average_length; /* FIXMEEEEEEE  */
        p->max_seq_len = si->ssi[i_file]->max_seq_len;
        p->base_error = 0.05f;
        p->indel_freq = 0.1f;
        for(i =0; i < 5;i++){
                p->back[i] = si->ssi[i_file]->background[i];
        }

        RUN(poahmm_from_read_structure(&poahmm, p,al->read_structure[i_hmm],  si->a));
        //LOG_MSG("Did I geta poa: %p", poahmm);

        //set_terminal_gap_prob(poahmm, si->ssi[i_file]->max_seq_len );
        //RUN(init_model_bag(&mb,al->read_structure[i_hmm], si->ssi[i_file], si->a, i_hmm));
        //if(!mb){                /* no model returned - probably because it is to long for reads on this file */


                //al->arch_posteriors[i_hmm][i_file] = prob2scaledprob(0.0f);
                //return OK;
        //}
        num_seq = rb[i_file]->num_seq;
        ri = rb[i_file]->sequences;

        a = si->a;

        MMALLOC(tmp_seq, sizeof(uint8_t) * (si->ssi[i_file]->max_seq_len+1));
        MMALLOC(tmp_qual, sizeof(uint8_t) * (si->ssi[i_file]->max_seq_len+1));

        for(i = 0; i < num_seq;i++){
                for(j = 0; j < ri[i]->len;j++){
                        tmp_seq[j] = tlalphabet_get_code(a, ri[i]->seq[j]);
                }
                tmp_seq[ri[i]->len] = 0;
                for(j = 0; j < ri[i]->len;j++){
                        tmp_qual[j] = ri[i]->qual[j] - base_q_offset;
                }
                //tmp_qualt[ri[i]->len] = 0;
                //l = (int) sb->sequences[i]->qual[j] - sb->base_quality_offset;
                RUN(viterbi_poahmm_banded(poahmm, tmp_seq, tmp_qual,  ri[i]->len, NULL, 0));
                score += poahmm->f_score;// - poahmm->random_scores[ ri[i]->len];
                //RUN(backward(mb, tmp_seq,ri[i]->len));
                //score +=  mb->b_score;
                //fprintf(stdout,"%d %d %d %f %d\n",i,i_file,i_hmm, mb->b_score, ri[i]->len);
        }
        //free_model_bag(mb);
        //fprintf(stdout,"Thread %d working on file %d; hmm %d\t%f score\n",omp_get_thread_num(),i_file,i_hmm,score);
        MFREE(tmp_seq);
        MFREE(tmp_qual);

        free_poahmm(poahmm);
        MFREE(p);
        //fprintf(stdout,"F:%d HMM:%d SCORE:%f\n",i_file,i_hmm,score);
        al->arch_posteriors[i_hmm][i_file] = score;
        return OK;
ERROR:
        return FAIL;
}


