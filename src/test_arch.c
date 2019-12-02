#include "test_arch.h"

#include "hmm_model_bag.h"
#include "io.h"
#include "core_hmm_functions.h"

#include <omp.h>

#define NUM_TEST_SEQ 1000


static int test_arch(struct read_info_buffer** rb, struct arch_library* al, struct seq_stats* si,int file,int hmm);

int test_architectures(struct arch_library* al, struct seq_stats* si, struct parameters* param)
{
        struct arch_bag* ab = NULL;
        struct file_handler* f_hand = NULL;
        struct read_info_buffer** rb = NULL;

        float sum;
        int best;
        int i,j;
        int (*fp)(struct read_info_buffer* rb, struct file_handler* f_handle) = NULL;

        /* alloc matrix for posteriors  */
        al->arch_posteriors = NULL;

        MMALLOC(al->arch_posteriors, sizeof(float*) * al->num_arch);
        for(i = 0; i < al->num_arch;i++){
                al->arch_posteriors[i] = NULL;
                MMALLOC(al->arch_posteriors[i], sizeof(float) * param->num_infiles);
                for(j = 0; j < param->num_infiles;j++){
                        al->arch_posteriors[i][j] = prob2scaledprob(1.0);
                }
        }

        MMALLOC(al->arch_to_read_assignment, sizeof(int) * param->num_infiles);

        omp_set_num_threads(param->num_threads);

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

        MMALLOC(rb, sizeof(struct read_info_buffer*) * param->num_infiles);
        for(i = 0; i < param->num_infiles;i++){

                rb[i] = NULL;
                RUN(alloc_read_info_buffer(&rb[i], NUM_TEST_SEQ));
                RUN(io_handler(&f_hand, param->infile[i]));
                if(f_hand->sam == 0){
                        fp = &read_fasta_fastq;
                }else {
                        fp = &read_sam_chunk;
                }
                /* read in sequences  */
                RUN(fp(rb[i],f_hand));//  param,file,&numseq));
                pclose(f_hand->f_ptr);
        }

#pragma omp parallel default(shared)
#pragma omp for collapse(2) private(i, j)
        for(i = 0; i < param->num_infiles;i++){
                for(j = 0; j < al->num_arch;j++){
                        test_arch(rb,al,si,i,j);
                }
        }


        /* convert posteriors into Phred scaled quality values */

        for(i = 0; i < param->num_infiles;i++){
                sum = prob2scaledprob(0.0f);
                for(j = 0; j < al->num_arch;j++){
                        sum = logsum(sum,al->arch_posteriors[j][i]);
                }
                //fprintf(stdout,"FILE: %s\t",param->infile[i]);
                for(j = 0; j < al->num_arch;j++){

                        al->arch_posteriors[j][i] = 1.0 - scaledprob2prob(al->arch_posteriors[j][i] - sum);
                        if(al->arch_posteriors[j][i] < 0.00001f){
                                al->arch_posteriors[j][i] = 0.00001f;
                        }
                        al->arch_posteriors[j][i] = -10.0f * log10f(al->arch_posteriors[j][i]);
                        //fprintf(stdout,"%f ",al->arch_posteriors[j][i]);
                }


        }

        for(i = 0; i < param->num_infiles;i++){
                al->arch_to_read_assignment[i] = -1;
                sum = -1.0;
                best = -1;
                for(j = 0; j < al->num_arch;j++){
                        if(al->arch_posteriors[j][i] > sum){
                                best =j;
                                sum = al->arch_posteriors[j][i];
                        }
                }
                ASSERT(best != -1,"No best arch found");
                LOG_MSG("Best architecture for %s is (Q = %f):", param->infile[i], al->arch_posteriors[best][i]);
                LOG_MSG("%s" , al->spec_line[best]);
                al->arch_to_read_assignment[i] = best;
                al->arch_posteriors[0][i] = al->arch_posteriors[best][i];
        }

        al->num_file = param->num_infiles;
        for(i = 0; i < param->num_infiles;i++){
                free_read_info_buffer(rb[i]);
        }
        MFREE(rb);
        MFREE(f_hand);

        MFREE(ab->arch_posterior);
        MFREE(ab->archs);
        MFREE(ab);
        return OK;
ERROR:
        return FAIL;
}

int test_arch(struct read_info_buffer** rb, struct arch_library* al, struct seq_stats* si,int i_file,int i_hmm)
{
        fprintf(stdout,"Thread %d working on file %d; hmm %d\n",omp_get_thread_num(),i_file,i_hmm);
        //printf ("num_thds=%d, max_thds=%d\n",omp_get_thread_num(),omp_get_max_threads());

        struct model_bag* mb = NULL;
        struct read_info** ri = NULL;
        int num_seq;
        int i;
        float score = prob2scaledprob(1.0);
        mb = init_model_bag(al->read_structure[i_hmm], si->ssi[i_file], i_hmm);
        num_seq = rb[i_file]->num_seq;
        ri = rb[i_file]->ri;
        for(i = 0; i < num_seq;i++){
                RUN(backward(mb, ri[i]->seq,ri[i]->len));
                score +=  mb->b_score;
                //fprintf(stdout,"%f\n", m->b_score);

        }
        free_model_bag(mb);

        al->arch_posteriors[i_hmm][i_file] = score;
        return OK;
ERROR:
        return FAIL;
}


