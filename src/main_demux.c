#include "interface.h"
#include "arch_lib.h"
#include "seq_stats.h"
#include "hmm_model_bag.h"
#include "test_arch.h"
#include "calibrate_hmm.h"
#include "extract_reads.h"

#include "tlrng.h"

#include "tllogsum.h"

int main (int argc,char * argv[]) {
        struct parameters* param = NULL;
        struct arch_library* al = NULL;
        struct seq_stats* si = NULL;
        struct rng_state* main_rng = NULL;
        int i,j;

        RUN(interface(&param,argc,argv));

        if(!param){
                return EXIT_SUCCESS;
        }
        ASSERT(param->num_infiles > 0, "Number of inputs has to be greater than 0");

        ASSERT(param->outfile != NULL, "No output file suffix");


        /* set top level rng generator  */
        RUNP(main_rng = init_rng(param->seed));

        /* create architecture library */
        RUN(alloc_arch_lib(&al));

        if(param->segments[0]){
                RUN(read_arch_into_lib(al, param->segments, param->num_segments));
        }

        if(param->arch_file){
                RUN(read_architecture_files(al, param->arch_file));
        }

        /* QC on architecture ?? */

        /* Start HMM stuff */
        init_logsum();


        /* get all sequence stats  */
        LOG_MSG("Got here");
        for( i = 0; i < param->num_infiles;i++){
                fprintf(stdout," File %d max: %s\n",i, param->infile[i]);
        }

        RUN(get_sequence_stats(&si,al, param->infile, param->num_infiles, main_rng));
        //si->ssi[0]->average_length;

        /* here there have to be sanity checks  */
        for(i = 0; i < param->num_infiles;i++){
                for(j = i+1; j < param->num_infiles;j++){
                        ASSERT(si->ssi[i]->total_num_seq == si->ssi[j]->total_num_seq,"File %s and %s contain different number of sequences", param->infile[i],param->infile[j]);
                }
        }
        //exit(0);
        for(i = 0; i < si->num;i++){
                LOG_MSG("%d %f", i,si->ssi[i]->average_length);
                for(j = 0; j < 5;j++){
                        fprintf(stdout,"%f ", scaledprob2prob(si->ssi[i]->background[j]));
                }
                fprintf(stdout,"\n");
        }
        //exit(0);
        /* allocate model for each architecture and each input file. */

        /* figure out which architectures belong to which read infiles */
        /* need to work on each read sequentially OR on all if enough memory ... */


        RUN(test_architectures(al,si,param));

        RUN(calibrate_architectures(al,si, main_rng));


        //int extract_reads(struct arch_library* al, struct seq_stats* si,struct parameters* param)
        RUN(extract_reads(al,si,param));
        exit(0);
//sprintf(param->buffer,"Start Run\n--------------------------------------------------\n");
        //param->messages = append_message(param->messages, param->buffer);
        //hmm_controller_multiple(param);
        free_arch_lib(al);
        free_sequence_stats(si);
        free_param(param);
        free_rng(main_rng);

        return EXIT_SUCCESS;
ERROR:
        if(param){
                //fprintf(stdout,"%s",param->errmsg);
                free_param(param);
        }
        return EXIT_SUCCESS;
}

