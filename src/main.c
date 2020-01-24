#ifdef HAVE_OPENMP
#include <omp.h>
#endif

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
        struct cookbook* cookbook = NULL;
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

        if(param->book_file && param->recipe){
                /* Format:
                   ./tagdust -r "@BC:S:ACAGTG,ACTTGA,TTAGGC;READ1:E:N{20-30};@READ2:E:N+;" ../dev/casava_read1_small.fastq.gz -o ~/tmp/gg
                */
                RUN(read_cookbook_command_line(&cookbook,param->recipe));
        }else if(param->book_file){
                RUN(read_cookbook_file (&cookbook,param->book_file));
        }else if(param->recipe){
                /* make book with one entry from command line...  */
                RUN(read_cookbook_command_line(&cookbook,param->recipe));
        }else{
                ERROR_MSG("No recipe specified. Either use --book <file> or --recipe <r>");
        }


        //exit(0);
        /* Sanity checks */
        /* are we dealing with multiple lanes? */
        /* are the sequences sorted? */
        /* are the segment specifications valid? */

        /* set top level rng generator */
        RUNP(main_rng = init_rng(param->seed));

        /* create architecture library */
        //RUN(alloc_arch_lib(&al));

        /*if(param->segments[0]){
                RUN(read_arch_into_lib(al, param->segments, param->num_segments));
        }

        if(param->arch_file){
                RUN(read_architecture_files(al, param->arch_file));
                }*/
        /* QC on architecture ?? */
#ifdef HAVE_OPENMP
        omp_set_num_threads(param->num_threads);
#else
        WARNING_MSG("Running without OpenMP!");
#endif
        /* Start HMM stuff */
        init_logsum();
        RUN(get_sequence_stats(&si, param->infile, param->num_infiles, main_rng));
        for(i = 0; i < param->num_infiles;i++){
                for(j = i+1; j < param->num_infiles;j++){
                        ASSERT(si->ssi[i]->total_num_seq == si->ssi[j]->total_num_seq,"File %s and %s contain different number of sequences", param->infile[i],param->infile[j]);
                }
        }


        RUN(test_architectures(cookbook,si,param));
        exit(0);
        RUN(calibrate_architectures(al,si, main_rng));
        //exit(0);
        //int extract_reads(struct arch_library* al, struct seq_stats* si,struct parameters* param)
        RUN(extract_reads(al,si,param,main_rng));

        //free_arch_lib(al);
        free_cookbook(&cookbook);
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

