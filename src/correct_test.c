#include <string.h>
#include <zlib.h>

#include "tldevel.h"
#include "tlmisc.h"
#include "tlrng.h"
#include "tlseqio.h"


#define CORRECT_IMPORT
#include "correct.h"

int test_correct(char* infile);
static int read_10x_white_list(struct tl_seq_buffer** b,char* filename);


int main(int argc, char *argv[])
{
        LOG_MSG("%d",argc);
        if(argc == 2){
                test_correct(argv[1]);
        }else{

        }
        //GGTTTACT
        return EXIT_SUCCESS;
ERROR:
        return EXIT_FAILURE;
}

int test_correct(char* infile)
{
        struct tl_seq_buffer* sb = NULL;
        khash_t(exact)* exact = NULL;
        struct rng_state* rng = NULL;

        double s[5];
        char* test_seq = NULL;

        float out;
        int test_len = 0;
        int i_point;
        int i,j;
        uint32_t k;
        LOG_MSG("%s",infile);
        if(!my_file_exists(infile)){
                ERROR_MSG("File %s not found");
        }
        RUN(read_10x_white_list(&sb, infile));
        LOG_MSG("read: %d barcodes" , sb->num_seq);
        RUN(fill_exact_hash(&exact, sb));
        //k = seq_to_code(sb->sequences[0]->seq, sb->sequences[0]->len);
        //LOG_MSG("%s key: %x",  sb->sequences[0]->seq,k);
        //code_to_seq(k,  sb->sequences[0]->len);
        free_tl_seq_buffer(sb);
        return OK;
ERROR:
        return FAIL;
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
