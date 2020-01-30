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

static int test_arch_file_order(struct arch_library* al,struct read_groups* rg);

static int test_files(struct arch_library* al,struct tl_seq_buffer** rb, struct read_ensembl* e, int* ambigious);
static int ambiguous_read_order_warning(struct arch_library* al, struct read_ensembl* e);
static int double_assignment_warning(struct arch_library* al,int x);
//static int test_files_poahmm(struct arch_library* al,struct tl_seq_buffer** rb, struct seq_stats* si,int n_files, int* ambigious);
//static int test_files_lpst(struct arch_library* al,struct tl_seq_buffer** rb, struct seq_stats* si,int n_files, int* ambigious);




static int test_arch(struct tl_seq_buffer** rb, struct arch_library* al, struct seq_stats* si,int i_file,int i_hmm,float* score);



int test_architectures(struct cookbook* cb, struct read_groups* rg)
{
        float sum;
        float max;
        int best;
        int n_best;
        int i;
        int g;
        int c;
        for(i = 0 ;i < cb->num_lib;i++){
                //LOG_MSG("------");
                //LOG_MSG("Testing: %s", cb->lib[i]->name);
                RUN(test_arch_file_order(cb->lib[i],rg));
                cb->scores[i] = cb->lib[i]->P;

                //LOG_MSG("score: %f", cb->lib[i]->P);
                //LOG_MSG("------");
        }
        sum = prob2scaledprob(0.0f);
        for(i = 0 ;i < cb->num_lib;i++){
                sum = logsum(sum, cb->scores[i]);
        }
        best = -1;
        max = -1.0;
        n_best = 0;
        for(i = 0 ;i < cb->num_lib;i++){
                //LOG_MSG("%s -> %f", cb->lib[i]->name, scaledprob2prob(cb->scores[i] - sum));
                if(max < scaledprob2prob(cb->scores[i] - sum)){
                        best = i;
                        max = scaledprob2prob(cb->scores[i] - sum);
                        n_best = 1;
                }else if(max == scaledprob2prob(cb->scores[i] - sum)){
                        n_best++;
                }
        }
        if( n_best > 1){

                WARNING_MSG("These %d recipes match the input:",n_best);
                WARNING_MSG("");
                for(i = 0 ;i < cb->num_lib;i++){
                        if(max == scaledprob2prob(cb->scores[i] - sum)){
                                WARNING_MSG("   %s",cb->lib[i]->name);
                        }
                }
                WARNING_MSG("");
                ERROR_MSG("Can't determine which one to use!");

        }
        if(n_best == 0){
                ERROR_MSG("No recipes matched the input");
        }
        if(!cb->lib[best]->read_order_check){
                RUN(ambiguous_read_order_warning(cb->lib[best], rg->e[0]));
                WARNING_MSG("will assume that the input file order matches the order in the recipe:");
                for(g = 0 ; g < rg->num_groups;g++){
                        //LOG_MSG("Read group %d:", g);
                        for(i = 0; i < cb->lib[best]->num_arch;i++){
                                rg->e[g]->arch_to_read_assignment[i] = i;
                                //LOG_MSG("%s -> %s",rg->e[g]->filenames[i],cb->lib[best]->spec_line[i]);
                        }
                }

                //for(i = 0; i < cb->lib[best]->num_arch;i++){
                //cb->lib[best]->arch_to_read_assignment[i] = i;
                //WARNING_MSG("%s -> %s", param->infile[i],cb->lib[best]->spec_line[i]);
                //}
        }

        LOG_MSG("Selecting: %s" , cb->lib[best]->name);
        for(g = 0 ; g < rg->num_groups;g++){
                LOG_MSG("Read group %d:", g);
                for(i = 0; i < cb->lib[best]->num_arch;i++){
                        c= rg->e[g]->arch_to_read_assignment[i];
                        LOG_MSG("%s -> %s",rg->e[g]->filenames[i],cb->lib[best]->spec_line[c]);
                }
        }

        /* selecting this library  */

        cb->best = best;
        return OK;
ERROR:
        return FAIL;
}

int test_arch_file_order(struct arch_library* al,struct read_groups* rg)
{
        //struct arch_bag* ab = NULL;
        struct file_handler* f_hand = NULL;
        struct tl_seq_buffer** rb = NULL;
        struct read_ensembl* e = NULL;
        //struct alphabet* a = NULL;


        int ambigious = 1;

        
        int i;
        int n_rg;
        if(al->num_arch != rg->e[0]->num_files){
                //WARNING_MSG("numarch != infiles %d %d", al->num_arch, rg->e[0]->num_files);
                al->P = prob2scaledprob(0.0);
                return OK;
        }
        al->P =  prob2scaledprob(1.0);
        //int (*fp)(struct read_info_buffer* rb, struct file_handler* f_handle) = NULL;

        //RUN(create_alphabet(&a, 0, TLALPHABET_DEFAULT_DNA));
        /* alloc matrix for posteriors  */
        al->arch_posteriors = NULL;

        /*MMALLOC(al->arch_posteriors, sizeof(float*) * al->num_arch);
        for(i = 0; i < al->num_arch;i++){
                al->arch_posteriors[i] = NULL;
                MMALLOC(al->arch_posteriors[i], sizeof(float) * param->num_infiles);
                for(j = 0; j < param->num_infiles;j++){
                        al->arch_posteriors[i][j] = 1.0f;// prob2scaledprob(0.0);
                }
        }

        MMALLOC(al->arch_to_read_assignment, sizeof(int) * param->num_infiles);*/

        //omp_set_num_threads(param->num_threads);

        /* here I combine: architectures with sequence parameters of individual input files  */
        /* to: check which architecture belongs to which read */

        ASSERT(al != NULL, "No arch library");
        //ASSERT(param != NULL, "No parameters");
        //ASSERT(si != NULL, "no seqstats");


        //MMALLOC(ab, sizeof(struct arch_bag));
        //ab->archs = 0;
        //ab->arch_posterior = 0;

        ///MMALLOC(ab->archs,sizeof(struct model_bag*) * al->num_arch);
        //MMALLOC(ab->arch_posterior,sizeof(float) *  al->num_arch);
        //ab->num_arch = al->num_arch; /* I'm not going to miss hmms in this process */

        MMALLOC(rb, sizeof(struct tl_seq_buffer*) * rg->e[0]->num_files);
        for(n_rg = 0; n_rg < rg->num_groups;n_rg++){
                LOG_MSG("Running readgroup %d", n_rg);
                e = rg->e[n_rg];
                for(i = 0; i < e->num_files;i++){
                        rb[i] = NULL;
                        LOG_MSG("Opening %s", e->filenames[i]);
                        RUN(open_fasta_fastq_file(&f_hand, e->filenames[i], TLSEQIO_READ));
                        RUN(read_fasta_fastq_file(f_hand, &rb[i],NUM_TEST_SEQ));
                        RUN(close_seq_file(&f_hand));
                }
                RUN(test_files(al, rb, e, &ambigious));
        }
        //LOG_MSG("%f poa ", al->P);
        //if( ambigious){
        //        RUN(test_files_lpst(al, rb, si, param->num_infiles, &ambigious));
        //LOG_MSG("%f lpst ", al->P);
        //}

        if( ambigious){
                al->read_order_check = 0;
                //ERROR_MSG("Giving up!");
        }



        al->num_file = rg->e[0]->num_files;

        for(i = 0; i < rg->e[0]->num_files;i++){
                free_tl_seq_buffer(rb[i]);
                //free_read_info_buffer(rb[i]);
        }
        MFREE(rb);

        //MFREE(ab->arch_posterior);
        //MFREE(ab->archs);
        //MFREE(ab);
        //free_alphabet(a);
        return OK;
ERROR:
        for(i = 0; i < rg->e[0]->num_files;i++){
                free_tl_seq_buffer(rb[i]);
                //free_read_info_buffer(rb[i]);
        }
        MFREE(rb);

        //MFREE(ab->arch_posterior);
        //MFREE(ab->archs);
        //MFREE(ab);
        //free_alphabet(a);
        return FAIL;
}

int test_files(struct arch_library* al,struct tl_seq_buffer** rb, struct read_ensembl* e, int* ambigious)
{
        float** post = NULL;
        float sum;
        float max;
        int best;
        int n_best;

        int i,j;
        int n_files;

        *ambigious = 0;
        n_files = e->num_files;
        RUN(galloc(&post,al->num_arch,n_files));
#ifdef HAVE_OPENMP
#pragma omp parallel default(shared)
#pragma omp for collapse(2) private(i, j)
#endif
        for(j = 0; j < al->num_arch;j++){
                for(i = 0; i < e->num_files;i++){
                        //if( al->arch_posteriors[j][i] == 1.0f){
                        test_arch(rb,al,e->si,i,j,&max);
                        post[j][i] = max;
                        lpst_score_read(al->read_structure[j], rb[i],e->si->ssi[i], &max);
                        post[j][i] += max;
                        //al->arch_posteriors[j][i] += max;
                                //}
                }
        }

        /* convert posteriors into Phred scaled quality values */
        //al->P = 1.0;
        for(i = 0; i < n_files ;i++){
                sum = prob2scaledprob(0.0f);
                for(j = 0; j < al->num_arch;j++){
                        //fprintf(stdout,"%f ",al->arch_posteriors[j][i]);
                        sum = logsum(sum,post[j][i]);
                }
                e->arch_to_read_assignment[i] = -1;
                //al->arch_to_read_assignment[i] = -1;
                max = -1.0;
                best = -1;
                n_best = 0;
                for(j = 0; j < al->num_arch;j++){
                        //fprintf(stdout,"%5f ",scaledprob2prob(post[j][i] - sum));
                        //LOG_MSG("%s" , al->spec_line[j]);
                        if(scaledprob2prob(post[j][i] - sum) > max){
                                n_best = 1;
                                best = j;
                                max = scaledprob2prob(post[j][i] - sum);

                        }else if(scaledprob2prob(post[j][i] - sum) == max){
                                n_best++;
                        }
                }
                //fprintf(stdout,"\n");
                /* score of matching the complete architecture */
                //al->P = al->P * sum;
                //ASSERT(best != -1,"No best arch found");
                if(best == -1){
                        al->P = prob2scaledprob(0.0f);
                        return OK;
                }
                al->P = al->P + post[best][i];

                e->arch_to_read_assignment[i] = best;
                //al->arch_to_read_assignment[i] = best;
                if(n_best == 1){

                        //al->arch_to_read_assignment[i] = best;
                        //al->arch_posteriors[0][i] = al->arch_posteriors[best][i];
                }else{
                        WARNING_MSG("%d or more architectures have a high probability of matching reads in file  %d.",n_best, i);

                        for(j = 0; j < al->num_arch;j++){
                                if(scaledprob2prob(post[j][i] - sum) == max){
                                        LOG_MSG("   %s" , al->spec_line[j]);
                                        //al->arch_posteriors[j][i] = 1.0f;
                                }else{
                                        //al->arch_posteriors[j][i] = 0.0f;
                                }
                        }
                        *ambigious = 1;
                }
        }
        fprintf(stdout,"\n");
        for(j = 0; j < al->num_arch;j++){
                post[j][0] = 0;
        }

        for(j = 0; j < n_files ;j++){
                post[ e->arch_to_read_assignment[j]][0]++;
        }
        for(j = 0; j < al->num_arch;j++){
                if(post[j][0] > 1.1f){
                        RUN(double_assignment_warning(al, j));
                        //al->P = prob2scaledprob(0.0);
                        break;
                }
                //LOG_MSG("%d %f",j , al->arch_posteriors[0][j]);
        }
        gfree(post);
        return OK;
ERROR:
        return FAIL;
}

int ambiguous_read_order_warning(struct arch_library* al, struct read_ensembl* e)
{
        int i;
        WARNING_MSG("");
        WARNING_MSG("Minor problem with this recipe:");
        for(i = 0;i < al->num_arch;i++){
                WARNING_MSG("%s", al->spec_line[i]);

        }
        WARNING_MSG("");
        WARNING_MSG("Can't determine which one of the above");
        WARNING_MSG("should be applied to which one of these input files:");
        for(i = 0;i < e->num_files;i++){
                WARNING_MSG("%s",e->filenames[i]);
        }
        return OK;
}

int double_assignment_warning(struct arch_library* al,int x)
{
        int i;
        WARNING_MSG("");
        WARNING_MSG("Problem with this: recipe:");
        for(i = 0;i < al->num_arch;i++){
                WARNING_MSG("%s", al->spec_line[i]);

        }
        WARNING_MSG("");
        WARNING_MSG("This read definition matches multiple reads:");
        WARNING_MSG("%s", al->spec_line[x]);

        return OK;
}

/*int test_files_lpst(struct arch_library* al,struct tl_seq_buffer** rb, struct seq_stats* si,int n_files, int* ambigious)
{
        float sum;
        float max;
        int best;
        int n_best;
        int i,j;

        *ambigious = 0;

#ifdef HAVE_OPENMP
#pragma omp parallel default(shared)
#pragma omp for collapse(2) private(i, j)
#endif
        for(j = 0; j < al->num_arch;j++){
                for(i = 0; i < n_files;i++){
                        //if( al->arch_posteriors[j][i] == 1.0f){
                        lpst_score_read(al->read_structure[j], rb[i],si->ssi[i], &al->arch_posteriors[j][i]);
                        //}
                }
        }

        for(i = 0; i < n_files ;i++){
                sum = prob2scaledprob(0.0f);
                for(j = 0; j < al->num_arch;j++){
                        //fprintf(stdout,"%f ",al->arch_posteriors[j][i]);
                        sum = logsum(sum,al->arch_posteriors[j][i]);
                }

                al->arch_to_read_assignment[i] = -1;
                max = -1.0;
                best = -1;
                n_best = 0;
                for(j = 0; j < al->num_arch;j++){
                        //LOG_MSG("%s" , al->spec_line[j]);
                        if(scaledprob2prob(al->arch_posteriors[j][i] - sum) > max){
                                n_best = 1;
                                best = j;
                                max = scaledprob2prob(al->arch_posteriors[j][i] - sum);

                        }else if(al->arch_posteriors[j][i] == max){
                                n_best++;
                        }
                }
                ASSERT(best != -1,"No best arch found");
                if(n_best == 1){
                        al->P = al->P + al->arch_posteriors[best][i];
                        al->arch_to_read_assignment[i] = best;
                        al->arch_posteriors[0][i] = al->arch_posteriors[best][i];
                }else{
                        WARNING_MSG("%d or more architectures have a high probability of matching reads in file  %d.",n_best, i);

                        for(j = 0; j < al->num_arch;j++){
                                if(al->arch_posteriors[j][i] == sum){
                                        LOG_MSG("   %f %s" ,al->arch_posteriors[j][i], al->spec_line[j]);
                                        al->arch_posteriors[j][i] = 1.0f;
                                }else{
                                        al->arch_posteriors[j][i] = 0.0f;
                                }

                        }
                        *ambigious = 1;
                }
        }
        return OK;
ERROR:
        return NULL;
}
*/




int test_arch(struct tl_seq_buffer** rb, struct arch_library* al, struct seq_stats* si,int i_file,int i_hmm,float* score)
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
        float l_score = prob2scaledprob(1.0);

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

        if(poahmm == NULL){
                *score = prob2scaledprob(0.0f);
                return OK;
        }
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
                l_score += poahmm->f_score;// - poahmm->random_scores[ ri[i]->len];
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
        *score = l_score;
        //al->arch_posteriors[i_hmm][i_file] = score;
        return OK;
ERROR:
        //*score = prob2scaledprob(0.0f);
        return FAIL;
}


