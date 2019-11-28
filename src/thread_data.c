#include "thread_data.h"



int alloc_thread_data(struct thread_data** td,int num_threads, int num_arch)
{
        struct thread_data* thread_data = NULL;
        int i;
        MMALLOC(thread_data,sizeof(struct thread_data)* param->num_threads);

        //interval =  (int)((double)numseq /(double)param->num_threads);

        for(i= 0;i< param->num_threads ;i++) {
                thread_data[i].fasta = NULL;//ce_fasta;
                thread_data[i].ri = NULL;
                thread_data[i].mb = NULL;copy_model_bag(mb);
                //thread_data[i].start = t*interval;
                //thread_data[i].end = t*interval + interval;
                thread_data[i].param = NULL;//param;
                thread_data[i].ab = NULL;
                MMALLOC(thread_data[i].ab,sizeof(struct arch_bag));
                thread_data[i].ab->num_arch = num_arch;
                thread_data[i].ab->arch_posterior = 0;
                thread_data[i].ab->archs = NULL;
                thread_data[i].ab->command_line = 0;
                MMALLOC(thread_data[i].ab->archs,sizeof(struct model_bag*) * num_arch );
                MMALLOC(thread_data[i].ab->arch_posterior,sizeof(float) * num_arch );
                for(i = 0; i < ab->num_arch;i++){
                        thread_data[i].ab->archs[i] = NULL;//copy_model_bag(ab->archs[i]);
                        thread_data[i].ab->arch_posterior[i] = prob2scaledprob(1.0);
                }

        }
        *td = thread_data;

        return OK;
ERROR:
        return FAIL;
}


void free_thread_data(struct thread_data* thread_data)
{
        int t;
        if(thread_data){

                for(t = 0;t < param->num_threads ;t++) {
                        for(i = 0; i < ab->num_arch;i++){
                                free_model_bag(thread_data[t].ab->archs[i]);// =copy_model_bag(ab->archs[i]);
                                //here I sum the posteriors from the different runs!!!
                                //ab->arch_posterior[i] = logsum(ab->arch_posterior[i] , thread_data[t].ab->arch_posterior[i]);
                                ab->arch_posterior[i]  += thread_data[t].ab->arch_posterior[i];
                        }

                        MFREE(thread_data[t].ab->archs);// = malloc(sizeof(struct model_bag*) * ab->num_arch );
                        MFREE(thread_data[t].ab->arch_posterior);//  = malloc(sizeof(float) * ab->num_arch );
                        MFREE(thread_data[t].ab);// = malloc(struct arch_bag);
                        free_model_bag(thread_data[t].mb);
                }
        }
        MFREE(thread_data);
}
