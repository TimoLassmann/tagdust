
#include "tldevel.h"
#include "tlseqio.h"



int main(int argc, char *argv[])
{
        struct file_handler* f = NULL;
        struct file_handler* f_out = NULL;
        struct tl_seq_buffer* sb = NULL;

        if(argc > 1){
                RUN(open_fasta_fastq_file(&f, argv[1], TLSEQIO_READ));


                RUN(open_fasta_fastq_file(&f_out, "Dummy", TLSEQIO_WRITE ));
                //int total_r = 0;
                //int total_w = 0;
                while(1){

                        RUN(read_fasta_fastq_file(f, &sb, 1000));
                        //total_r+= sb->num_seq;
                        //LOG_MSG("Finished reading chunk: found %d ",sb->num_seq);
                        /*for(i = 0; i < sb->num_seq;i++){

                                fprintf(stdout,"%s\n", sb->sequences[i]->name);
                                for(j = 0; j < sb->sequences[i]->len;j++){
                                        fprintf(stdout,"%c",(char)sb->sequences[i]->seq[j]);
                                }
                                fprintf(stdout,"\n");

                                for(j = 0; j < sb->sequences[i]->len;j++){
                                        fprintf(stdout,"%c",sb->sequences[i]->qual[j]);
                                }
                                fprintf(stdout,"\n");

                        }*/


                        if(sb->num_seq == 0){
                                break;
                        }
                        //total_w+= sb->num_seq;
                        //LOG_MSG("%d %d",total_r,total_w);
                        RUN(write_fasta_fastq(sb, f_out));
                }

                free_tl_seq_buffer(sb);
                RUN(close_fasta_fastq_file(&f));
                RUN(close_fasta_fastq_file(&f_out));
                //fprintf(stdout,"%p",f);
        }

        return EXIT_SUCCESS;
ERROR:
        return EXIT_FAILURE;
}
