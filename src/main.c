#include "interface.h"
#include "arch_lib.h"
#include "seq_stats.h"
#include "nuc_code.h"

int main (int argc,char * argv[]) {
        struct parameters* param = NULL;
        struct arch_library* al = NULL;
        struct seq_stats* si = NULL;
        int i,j;

        RUN(interface(&param,argc,argv));

        if(!param){
                return EXIT_SUCCESS;
        }
        ASSERT(param->num_infiles > 0, "Number of inputs has to be greater than 0");

        ASSERT(param->outfile != NULL, "No output file suffix");

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
        RUN(init_nuc_code());

        /* get all sequence stats  */

        RUN(get_sequence_stats(&si,al, param->infile, param->num_infiles));
        /* allocate model for each architecture and each input file. */

        /* figure out which architectures belong to which read infiles */
        /* need to work on each read sequentially OR on all if enough memory ... */



        //sprintf(param->buffer,"Start Run\n--------------------------------------------------\n");
        //param->messages = append_message(param->messages, param->buffer);
        //hmm_controller_multiple(param);
        free_arch_lib(al);
        free_sequence_stats(si);
        free_param(param);

        return EXIT_SUCCESS;
ERROR:
        if(param){
                //fprintf(stdout,"%s",param->errmsg);
                free_param(param);
        }
        return EXIT_SUCCESS;
}












