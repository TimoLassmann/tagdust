
#include <string.h>
#include <zlib.h>

#include "tldevel.h"

#include "tlmisc.h"
#include "tlrng.h"
#include "tlseqio.h"

#define PST_IMPORT
#include "pst.h"


static int read_10x_white_list(struct tl_seq_buffer** b,char* filename);

int main(int argc, char *argv[])
{
        char alphabet[] = "ACGT";
        LOG_MSG("Hello World");
        char* filename = NULL;
        char* test_seq = NULL;
        struct tl_seq_buffer* sb = NULL;
        struct pst* p = NULL;
        struct rng_state* rng = NULL;
        int i,j;
        float out;
        LOG_MSG("%d",argc);
        if(argc == 2){
                filename = argv[1];
                LOG_MSG("%s",filename);
                if(!my_file_exists(filename)){
                        ERROR_MSG("File %s not found");
                }
                RUN(read_10x_white_list(&sb, filename));

                RUN(run_build_pst(&p, sb));

                int num = 0;
                //print_pst(p, p->pst_root, &num);
                LOG_MSG("Found %d leaves ",num);

                num = 0;
                //print_pst(p, p->ppt_root, &num);
                LOG_MSG("Found %d leaves ",num);

                //exit(0);
                for(i = 0; i < 10;i++){
                        fprintf(stdout,">%s\n%s\n",sb->sequences[i]->name,sb->sequences[i]->seq);
                }
                rng = init_rng(0);
                MMALLOC(test_seq, sizeof(char) * (sb->sequences[0]->len+1));
                for(i = 0; i < 10;i++){


                        scan_read_with_pst(p, sb->sequences[i]->seq, sb->sequences[i]->len,&out);
                        LOG_MSG("%f\t%s (overfit)",sb->sequences[i]->seq,out);
                        sb->sequences[i]->seq[6] = 'T';
                        scan_read_with_pst(p, sb->sequences[i]->seq, sb->sequences[i]->len,&out);
                        LOG_MSG("%f\t%s (Added T)",sb->sequences[i]->seq,out);


                        for(j = 0; j < sb->sequences[i]->len;j++){
                                test_seq[j] = alphabet[ tl_random_int(rng,4)];
                        }
                        test_seq[sb->sequences[i]->len]= 0;
                        scan_read_with_pst(p, test_seq, sb->sequences[i]->len,&out);
                        LOG_MSG("%f\t%s (random)",sb->sequences[i]->seq,out);
                }
                //free_error_correct_seq(&e);

                MFREE(test_seq);

                //if()
                free_pst(p);
                free_tl_seq_buffer(sb);
                free_rng(rng);
        }

        return EXIT_SUCCESS;
ERROR:
        return EXIT_FAILURE;
}


int read_10x_white_list(struct tl_seq_buffer** b,char* filename)
{
        struct tl_seq_buffer* sb = NULL;
        struct tl_seq* s = NULL;
        gzFile f_ptr;
        char* buffer = NULL;
        char* tmp = NULL;
        int buffer_len = 256;

        MMALLOC(buffer, sizeof(char) * buffer_len);

        RUN(alloc_tl_seq_buffer(&sb, 1000000));


        RUNP(f_ptr = gzopen(filename, "r"));
        while((tmp =  gzgets(f_ptr, buffer, buffer_len)) != NULL){
                //fprintf(stdout,"%s",buffer);
                s = sb->sequences[sb->num_seq];
                s->len = strnlen(buffer, buffer_len) -1;

                snprintf(s->name, TL_SEQ_MAX_NAME_LEN, "Bar%d", sb->num_seq+1);

                while(s->malloc_len < s->len){
                        RUN(resize_tl_seq(s));
                }
                strncpy(s->seq, buffer, s->len);
                s->seq[s->len] = 0;
                //snprintf(s->seq, s->malloc_len, "%s",buffer);
                sb->num_seq++;
                //if(sb->num_seq == 100000){
                //break;
                //}
                if(sb->num_seq == sb->malloc_num){
                        RUN(resize_tl_seq_buffer(sb));
                }
        }

        *b = sb;

        gzclose(f_ptr);


        MFREE(buffer);
        return OK;

ERROR:
        if(buffer){
                MFREE(buffer);
        }
        if(f_ptr){
                gzclose(f_ptr);
        }
        return FAIL;

}
