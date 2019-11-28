#include "test_arch.h"

#include "hmm_model_bag.h"
#include "io.h"
#include "core_hmm_functions.h"

#include <omp.h>


static int run_test_arch(const struct read_info_buffer* rb,struct model_bag* model_bag,float* res);

int test_architectures(struct arch_library* al, struct seq_stats* si, struct parameters* param)
{
        struct arch_bag* ab = NULL;

        struct file_handler* f_hand = NULL;
        struct read_info_buffer* rb = NULL;
        //struct read_info** ri = NULL;
        float sum;

        int i,j;
        int (*fp)(struct read_info_buffer* rb, struct file_handler* f_handle) = NULL;

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

        RUN(alloc_read_info_buffer(&rb,1000000));

        for(i = 0; i < param->num_infiles;i++){

                for(j = 0; j < al->num_arch;j++){
                        ab->archs[j] = NULL;
                        ab->archs[j] = init_model_bag( al->read_structure[j], si->ssi[i], j);
                        ab->arch_posterior[j] = prob2scaledprob(1.0);
                }
                RUN(io_handler(&f_hand, param->infile[i]));
                if(f_hand->sam == 0){
                        fp = &read_fasta_fastq;
                }else {
                        fp = &read_sam_chunk;
                }
                /* read in sequences  */
                RUN(fp(rb,f_hand));//  param,file,&numseq));
                        //if((status = fp(ri, param,file,&numseq)) != OK)  exit(status);
                /* run HMM code */
#pragma omp parallel for
                for(j = 0; j < al->num_arch;j++){
                        //ab->arch_posterior[j] = prob2scaledprob(1.0);
                        run_test_arch(rb, ab->archs[j], &ab->arch_posterior[j]);
                }
                /* this was the old call:  */
                //RUN(run_pHMM(ab,mb,ri,param,0,numseq,MODE_ARCH_COMP));

                pclose(f_hand->f_ptr);


                sum = prob2scaledprob(0.0f);
                for(j = 0; j < al->num_arch;j++){
                        //	fprintf(stderr,"%d %s %f	%s", file_num, param->infile[file_num] ab->arch_posterior[i],ab->command_line[i]);
                        sum = logsum(sum, ab->arch_posterior[j]);
                }
                //best_architecture = -1;
                //best_score = -1.0;

                for(j = 0; j < ab->num_arch;j++){
                        ab->arch_posterior[j] = scaledprob2prob(ab->arch_posterior[j] - sum);
                        fprintf(stderr,"model%d %f \n",j,  ab->arch_posterior[j]);
                        /*if(ab->arch_posterior[i]  > best_score){
                                best_score =ab->arch_posterior[i] ;
                                best_architecture = i;
                                }*/
                }
                for(j = 0; j < ab->num_arch;j++){
                        free_model_bag(ab->archs[j]);

                }


        }

        free_read_info_buffer(rb);
        MFREE(f_hand);

        MFREE(ab->arch_posterior);
        MFREE(ab->archs);
        MFREE(ab);
        return OK;
ERROR:
        return FAIL;
}



int run_test_arch(const struct read_info_buffer* rb,struct model_bag* model_bag,float* res)
{
        struct model_bag* mb = NULL;
        struct read_info* ri = NULL;
        int i;
        float score = prob2scaledprob(1.0);
        printf ("num_thds=%d, max_thds=%d\n",omp_get_thread_num(),omp_get_max_threads());

/* make local copy_ */
        mb = copy_model_bag(model_bag);
        for(i = 0; i < rb->num_seq;i++){
                ri = rb->ri[i];


                RUN(backward(mb, ri->seq,ri->len));


                score +=  mb->b_score;
                //fprintf(stdout,"%f\n", m->b_score);

        }
        free_model_bag(mb);
        *res = score;
        return OK;
ERROR:
        return FAIL;
}
