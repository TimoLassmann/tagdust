#include <string.h>
#include <math.h>

#include "seq_stats.h"

#include "tldevel.h"
#include "tllogsum.h"
#include "tlseqio.h"
#include "tlalphabet.h"

static int alloc_sequence_stats_info(struct sequence_stats_info** si, int n);
static void free_sequence_stats_info(struct sequence_stats_info* si);

static int five_prime_exact_match(char* seq,char*p,int seq_len, double* res);
static int three_prime_exact_match(char* seq,char*p,int seq_len, double* res);

int get_sequence_stats(struct seq_stats** sequence_stats, struct arch_library* al,char** infiles,int numfiles,struct rng_state* main_rng)
{
        struct seq_stats* si = NULL;
        struct file_handler* f_hand = NULL;

        struct tl_seq_buffer* rb = NULL;
        struct tl_seq** ri = NULL;

        struct alphabet* a = NULL;


        int i,j,c;
        int total_read;
        int len;
        int last;

        char** five_test_sequence = NULL;
        char** three_test_sequence = NULL;

        double* five_s0 = NULL;
        double* five_s1 = NULL;
        double* five_s2 = NULL;
        double* three_s0 = NULL;
        double* three_s1 = NULL;
        double* three_s2 = NULL;
        double res;
        double sum;



        MMALLOC(five_s0, sizeof(double) * al->num_arch);
        MMALLOC(five_s1, sizeof(double) * al->num_arch);
        MMALLOC(five_s2, sizeof(double) * al->num_arch);
        MMALLOC(three_s0, sizeof(double) * al->num_arch);
        MMALLOC(three_s1, sizeof(double) * al->num_arch);
        MMALLOC(three_s2, sizeof(double) * al->num_arch);

        MMALLOC(si, sizeof(struct seq_stats));
        si->num = numfiles;
        si->ssi = NULL;
        MMALLOC(si->ssi, sizeof(struct sequence_stats_info*) * si->num);
        si->a = NULL;

        RUN(create_alphabet(&si->a, main_rng, TLALPHABET_DEFAULT_DNA));
        a = si->a;

        for(i = 0; i < si->num;i++){
                si->ssi[i] = NULL;
                RUN(alloc_sequence_stats_info(&si->ssi[i], al->num_arch));
        }


        //RUN(alloc_read_info_buffer(&rb,100000));

        /* copy5' and 3' sequences for matching in case of partial segments */
        MMALLOC(five_test_sequence, sizeof(char*) * al->num_arch);
        MMALLOC(three_test_sequence, sizeof(char*) * al->num_arch);
        for(i = 0;i < al->num_arch;i++){
                five_test_sequence[i] = NULL;
                three_test_sequence[i] = NULL;
                if(al->read_structure[i]->type[0] == 'P'){
                        len =   al->read_structure[i]->segment_length[0];//  strlen(al->read_structure[i]->sequence_matrix[0][0]);
                        MMALLOC(five_test_sequence[i], sizeof(char) * (len+1));
                        for(c = 0; c < numfiles;c++){
                                si->ssi[c]->expected_5_len[i] = (double)len;
                        }
                        for(j = 0; j < len;j++){
                                five_test_sequence[i][j] =  tlalphabet_get_code(a,al->read_structure[i]->sequence_matrix[0][0][j]);
                        }
                        five_test_sequence[i][len] = 0;
                }
                last = al->read_structure[i]->num_segments -1;
                if(al->read_structure[i]->type[last] == 'P'){
                        len =  al->read_structure[i]->segment_length[last];// strlen(al->read_structure[i]->sequence_matrix[last][0]);
                        MMALLOC(three_test_sequence[i], sizeof(char) * (len+1));
                        for(c = 0; c < numfiles;c++){
                                si->ssi[c]->expected_3_len[i] = (double)len;
                        }
                        for(j = 0; j < len;j++){
                                three_test_sequence[i][j] =  tlalphabet_get_code(a,al->read_structure[i]->sequence_matrix[0][0][j]);
                        }
                        three_test_sequence[i][len] = 0;
                }
        }

        /* Do stuff */

        for(i = 0; i < numfiles;i++){
                RUN(open_fasta_fastq_file(&f_hand, infiles[i], TLSEQIO_READ));
                //aopen_fasta_fastq_file(struct file_handler **fh, char *filename, int mode)
                //RUN(io_handler(&f_hand,infiles[i]));
                for(c = 0;c < al->num_arch;c++){
                        five_s0[c] = 0.0;
                        five_s1[c] = 0.0;
                        five_s2[c] = 0.0;
                        three_s0[c] = 0.0;
                        three_s1[c] = 0.0;
                        three_s2[c] = 0.0;
                }

                total_read = 0;
                //LOG_MSG("FILE:%s", infiles[i]);
                while(1){
                        RUN(read_fasta_fastq_file(f_hand, &rb,100000));
                        //read_fasta_fastq(struct read_info_buffer *rb, struct file_handler *f_handle)
                        //RUN(fp(rb,f_hand));//  param,file,&numseq));
                        //if((status = fp(ri, param,file,&numseq)) != OK)  exit(status);

                        if(!rb->num_seq ){
                                break;
                        }
                        ri = rb->sequences ;
                        for(j = 0; j < rb->num_seq;j++){
                                if(ri[j]->len > si->ssi[i]->max_seq_len){
                                        si->ssi[i]->max_seq_len = ri[j]->len;
                                }

                                //print_sequence(rb->ri[j], stdout);
                                si->ssi[i]->average_length += ri[j]->len;
                                for(c = 0;c < ri[i]->len;c++){

                                        si->ssi[i]->background[tlalphabet_get_code(a,ri[j]->seq[c])] += 1.0f;
                                }
                                for(c = 0;c < al->num_arch;c++){
                                        if(five_test_sequence[c]){
                                                five_prime_exact_match(ri[j]->seq, five_test_sequence[c], ri[j]->len, &res);
                                                if(res){
                                                        five_s0[c]++;
                                                        five_s1[c] += res;
                                                        five_s2[c] += res*res;
                                                }
                                        }

                                        if(three_test_sequence[c]){
                                                three_prime_exact_match(ri[j]->seq, three_test_sequence[c], ri[j]->len, &res);
                                                if(res){
                                                        three_s0[c]++;
                                                        three_s1[c] += res;
                                                        three_s2[c] += res*res;
                                                }
                                        }
                                }

                        }
                        total_read += rb->num_seq;
                        //LOG_MSG("total: %d", total_read);
#if DEBUG
                        //if(total_read > 10){
                        //break;
                        //}
#else
                        //if(total_read > 10){
                        //break;
                        //}
#endif


                }
                RUN(close_fasta_fastq_file(&f_hand));
                //pclose(f_hand->f_ptr);
                //LOG_MSG("total: %d", total_read);

                si->ssi[i]->total_num_seq = total_read;
                si->ssi[i]->average_length = (int) floor((double) si->ssi[i]->average_length / (double) total_read   + 0.5);

                sum = 0.0;
                for(j = 0; j < 5;j++){
                        sum += si->ssi[i]->background[j];
                }

                for(j = 0; j < 5;j++){
                        si->ssi[i]->background[j] = prob2scaledprob(si->ssi[i]->background[j]  / sum);
                }

                for(c = 0;c < al->num_arch;c++){
                        if(five_test_sequence[c]){
                                if(five_s0[c] <= 1){
                                        WARNING_MSG("there seems to e not a single read containing the 5' partial sequence.\n");
                                        si->ssi[i]->mean_5_len[c] = si->ssi[i]->expected_5_len[c];
                                        si->ssi[i]->stdev_5_len[c] = 1.0;
                                }else{
                                        si->ssi[i]->mean_5_len[c] = five_s1[c] / five_s0[c];
                                        //ssi->mean_5_len = five_s1 / five_s0;

                                        si->ssi[i]->stdev_5_len[c]  = sqrt((five_s0[c] * five_s2[c] - pow(five_s1[c],2.0)) / (five_s0[c] *( five_s0[c] - 1.0)));
                                        if(!si->ssi[i]->stdev_5_len[c]){
                                                si->ssi[i]->stdev_5_len[c] = 10000.0;
                                        }
                                        //fprintf(stderr,"5: %f %f	%f\n", ssi->mean_5_len,  ssi->stdev_5_len,five_s0);
                                        //if(ssi->stdev_5_len < 1){
                                        //	ssi->stdev_5_len = 1;
                                        //}

                                        //fprintf(stderr,"5: %f %f	%f\n", ssi->mean_5_len,  ssi->stdev_5_len,five_s0);
                                        if(si->ssi[i]->mean_5_len[c] <= 1){
                                                WARNING_MSG("5' partial segment seems not to be present in the data (length < 1).\n");
                                        }

                                }

                        }else{
                                si->ssi[i]->mean_5_len[c] = -1.0;
                                si->ssi[i]->stdev_5_len[c] = -1.0;
                        }
                        if(three_test_sequence[c]){
                                if(three_s0[c] <= 1){
                                        WARNING_MSG("there seems to e not a single read containing the 5' partial sequence.\n");

                                        si->ssi[i]->mean_3_len[c] = si->ssi[i]->expected_3_len[c];
                                        si->ssi[i]->stdev_3_len[c] = 1.0;
                                }else{
                                        si->ssi[i]->mean_3_len[c] = three_s1[c] / three_s0[c];
                                        //ssi->mean_5_len = three_s1 / three_s0;

                                        si->ssi[i]->stdev_3_len[c]  = sqrt((three_s0[c] * three_s2[c] - pow(three_s1[c],2.0)) / (three_s0[c] *( three_s0[c] - 1.0)));
                                        if(!si->ssi[i]->stdev_3_len[c]){
                                                si->ssi[i]->stdev_3_len[c] = 10000.0;
                                        }
                                        //fprintf(stderr,"5: %f %f	%f\n", ssi->mean_3_len,  ssi->stdev_3_len,three_s0);
                                        //if(ssi->stdev_3_len < 1){
                                        //	ssi->stdev_3_len = 1;
                                        //}

                                        //fprintf(stderr,"5: %f %f	%f\n", ssi->mean_3_len,  ssi->stdev_3_len,three_s0);
                                        if(si->ssi[i]->mean_3_len[c] <= 1){
                                                WARNING_MSG("5' partial segment seems not to be present in the data (length < 1).\n");
                                        }
                                }
                        }else{
                                si->ssi[i]->mean_3_len[c] = -1.0;
                                si->ssi[i]->stdev_3_len[c] = -1.0;
                        }
                }
        }


        for(i = 0;i < al->num_arch;i++){
                if(five_test_sequence[i]){
                        MFREE(five_test_sequence[i]);
                }

                if(three_test_sequence[i]){
                        MFREE(three_test_sequence[i]);
                }
        }
        MFREE(five_test_sequence);
        MFREE(three_test_sequence);

        MFREE(five_s0);
        MFREE(five_s1);
        MFREE(five_s2);
        MFREE(three_s0);
        MFREE(three_s1);
        MFREE(three_s2);
        free_tl_seq_buffer(rb);
        //free_read_info_buffer(rb);
        //free_alphabet(a);
        //LOG_MSG("Got here");
        *sequence_stats = si;
        return OK;
ERROR:
        free_sequence_stats(si);
        return FAIL;
}

int five_prime_exact_match(char* seq,char*p,int seq_len, double* res)
{
        int i,j;
        int five_len;

        five_len = strlen(p);
        *res = 0.0;
        if(five_len >= seq_len){
                *res = 0.0;
        }

        for(i = 0;i <= five_len ;i++){
                for(j = 0;j < five_len-i;j++){
                        if(seq[j] != p[i +j]){
                                break;
                        }
                }
                if(j == five_len-i && j > 3 ){
                        *res = (double) j;
                        break;
                }
        }

        return OK;
}

int three_prime_exact_match(char* seq,char*p,int seq_len, double* res)
{
        int three_len;
        int i,j;
        three_len = strlen(p);
        *res = 0.0;
        if(three_len >= seq_len){
                *res = 0.0;
        }

        for(i = 0;i <= three_len ;i++){
                for(j = 0;j < three_len-i;j++){
                        if(seq[seq_len - (three_len-i-j)] != p[j]){
                                break;
                        }
                }
                if(j == three_len-i  && j > 3){
                        *res = (double) j;
                        //three_s0++;
                        //three_s1 += three_len -i;
                        //three_s2 += (three_len-i) * (three_len-i);
                        break;
                }
        }
        return OK;
}

void free_sequence_stats(struct seq_stats* si)
{
        int i;
        if(si){
                for(i = 0; i < si->num;i++){
                        free_sequence_stats_info(si->ssi[i]);
                }
                free_alphabet(si->a);
                MFREE(si->ssi);
                MFREE(si);
        }
}

int alloc_sequence_stats_info(struct sequence_stats_info** si, int n)
{
        struct sequence_stats_info* ssi = NULL;
        int i;


        MMALLOC(ssi, sizeof(struct sequence_stats_info));

        //ssi->average_length
        ssi->average_length = 0;
        for(i = 0; i < 5;i++){
                ssi->background[i] = 1.0;
        }

        ssi->expected_5_len = NULL;
        ssi->expected_3_len = NULL;
        ssi->mean_5_len = NULL;
        ssi->stdev_5_len = NULL;
        ssi->mean_3_len = NULL;
        ssi->stdev_3_len = NULL;
        ssi->average_length = 0.0f;
        ssi->max_seq_len = 0;

        MMALLOC(ssi->expected_5_len, sizeof(double) * n);
        MMALLOC(ssi->expected_3_len, sizeof(double) * n);
        MMALLOC(ssi->mean_5_len, sizeof(double) * n);
        MMALLOC(ssi->mean_3_len, sizeof(double) * n);
        MMALLOC(ssi->stdev_5_len, sizeof(double) * n);
        MMALLOC(ssi->stdev_3_len, sizeof(double) * n);

        for(i = 0; i < n;i++){
                ssi->expected_5_len[i] =  0.0f;
                ssi->expected_3_len[i] =  0.0f;
                ssi->mean_5_len[i] =  0.0f;
                ssi->mean_3_len[i] =  0.0f;
                ssi->stdev_5_len[i] =  0.0f;
                ssi->stdev_3_len[i] = 0.0f;
        }
        *si = ssi;
        return OK;
ERROR:
        free_sequence_stats_info(ssi);
        return FAIL;
}

void free_sequence_stats_info(struct sequence_stats_info* si)
{
        if(si){
                MFREE(si->expected_5_len);
                MFREE(si->expected_3_len);
                MFREE(si->mean_5_len);
                MFREE(si->mean_3_len);
                MFREE(si->stdev_5_len);
                MFREE(si->stdev_3_len);
                MFREE(si);
        }
}

/*struct sequence_stats_info* get_sequence_stats(struct parameters* param, struct read_info** ri, int file_num )
{
        struct sequence_stats_info* ssi = 0;
        int status;
        MMALLOC(ssi, sizeof(struct sequence_stats_info));
        FILE* file = 0;

        int (*fp)(struct read_info** ,struct parameters*,FILE*,int* buffer_count  ) = 0;
        int i,j,c,numseq,total_read;
        int five_len = 0;
        int three_len = 0;

        double sum = 0.0;

        double five_s0 = 0.0;
        double five_s1 = 0.0;
        double five_s2 = 0.0;
        double three_s0 = 0.0;
        double three_s1 = 0.0;
        double three_s2 = 0.0;

        char* five_test_sequence = 0;
        char* three_test_sequence = 0;

        ssi->average_length = 0;
        for(i = 0; i < 5;i++){
                ssi->background[i] = 1.0;
        }

        ssi->expected_5_len = 0;
        ssi->expected_3_len = 0;
        ssi->mean_5_len = 0.0f;
        ssi->stdev_5_len = 0.0f;
        ssi->mean_3_len = 0.0f;
        ssi->stdev_3_len = 0.0f;
        ssi->average_length = 0.0f;
        ssi->max_seq_len = 0;


        file =  io_handler(file, file_num,param);

        if(param->sam == 0){
                fp = &read_fasta_fastq;
        }else {
                fp = &read_sam_chunk;
        }
        numseq = 0;
        total_read= 0;

        //Do I need to test for a 5' partial sequence?
        if(param->read_structure->type[0] == 'P'){
                five_len = (int) strlen(param->read_structure->sequence_matrix[0][0]);
                ssi->expected_5_len = five_len;
                MMALLOC(five_test_sequence, sizeof(char) * (five_len+1));

                for(i = 0; i < five_len;i++){
                        five_test_sequence[i] = nuc_code[(int) param->read_structure->sequence_matrix[0][0][i]];
                }
                five_test_sequence[five_len] = 0;
        }
        //Do I need to test for a 3' partial sequence?
        if(param->read_structure->type[param->read_structure->num_segments-1] == 'P'){
                three_len = (int) strlen(param->read_structure->sequence_matrix[ param->read_structure->num_segments-1][0]);
                ssi->expected_3_len = three_len;
                MMALLOC(three_test_sequence ,sizeof(char) * (three_len+1));
                for(i = 0; i < three_len;i++){
                        three_test_sequence[i] = nuc_code[(int) param->read_structure->sequence_matrix[ param->read_structure->num_segments-1][0][i]];
                }
                three_test_sequence[three_len] = 0;

        }
        while(1){
                RUN(fp(ri, param,file,&numseq));
                //if((status = fp(ri, param,file,&numseq)) != OK)  exit(status);
                if(!numseq){
                        break;
                }
//	while ((numseq = fp(ri, param,file)) != 0){
                for(i = 0; i < numseq;i++){
                        if(ri[i]->len > ssi->max_seq_len){
                                ssi->max_seq_len = ri[i]->len;
                        }
                        ssi->average_length += ri[i]->len;
                        for(j = 0;j < ri[i]->len;j++){
                                //fprintf(stderr,"%d ",(int)ri[i]->seq[j] );
                                ssi->background[(int)ri[i]->seq[j]] += 1.0f;
                        }
                        //check length of exact 5' matching sequence...
                        if(five_len){
                                for(j = 0;j <= five_len ;j++){
                                        for(c = 0;c < five_len-j;c++){
                                                if(ri[i]->seq[c] != five_test_sequence[j +c]){
                                                        break;
                                                }
                                        }
                                        if(c == five_len-j && c > 3 ){
                                                five_s0++;
                                                five_s1 += five_len -j;
                                                five_s2 += (five_len-j) * (five_len-j);
                                                break;
                                        }
                                }
                        }
                        //check length of exact 3' matching sequence...
                        if(three_len){
                                for(j = 0;j <= three_len ;j++){
                                        for(c = 0;c < three_len-j;c++){
                                                if(ri[i]->seq[ri[i]->len - (three_len-j -c)] != three_test_sequence[c]){
                                                        break;
                                                }
                                        }
                                        if(c == three_len-j  && c > 3){
                                                three_s0++;
                                                three_s1 += three_len -j;
                                                three_s2 += (three_len-j) * (three_len-j);
                                                break;
                                        }
                                }
                        }
                }

                total_read += numseq;
#if DEBUG
                if(total_read > 10001){
                        break;
                }
#else
                if(total_read > 1000000){
                        break;
                }
#endif
        }


        if(five_len){
                if(five_s0 <= 1){
                        sprintf(param->buffer,"WARNING: there seems to e not a single read containing the 5' partial sequence.\n");
                        param->messages = append_message(param->messages, param->buffer);
                        ssi->mean_5_len  = ssi->expected_5_len;
                        ssi->stdev_5_len  = 1.0;
                }else{
                        ssi->mean_5_len = five_s1 / five_s0;
                        ssi->stdev_5_len = sqrt(  (five_s0 * five_s2 - pow(five_s1,2.0))   /  (  five_s0 *(five_s0-1.0) )) ;
                        if(!ssi->stdev_5_len){
                                ssi->stdev_5_len = 10000.0;
                        }
                        //fprintf(stderr,"5: %f %f	%f\n", ssi->mean_5_len,  ssi->stdev_5_len,five_s0);
                        //if(ssi->stdev_5_len < 1){
                        //	ssi->stdev_5_len = 1;
                        //}

                        //fprintf(stderr,"5: %f %f	%f\n", ssi->mean_5_len,  ssi->stdev_5_len,five_s0);
                        if(ssi->mean_5_len <= 1){
                                sprintf(param->buffer,"WARNING: 5' partial segment seems not to be present in the data (length < 1).\n");
                                param->messages = append_message(param->messages, param->buffer);
                                //free_param(param);
                                //exit(EXIT_FAILURE);
                        }
                }
        }else{
                ssi->mean_5_len =  -1.0;
                ssi->stdev_5_len = -1.0;
        }

        if(three_len){
                if(three_s0 <= 1){
                        sprintf(param->buffer,"WARNING: 3' partial segment seems not to be present in the data.\n");
                        param->messages = append_message(param->messages, param->buffer);
                        ssi->mean_3_len  = ssi->expected_3_len;
                        ssi->stdev_3_len  = 1.0;

                }else{

                        ssi->mean_3_len = three_s1 / three_s0;
                        ssi->stdev_3_len = sqrt(  (three_s0 * three_s2 - pow(three_s1,2.0))   /  (  three_s0 *(three_s0-1.0) )) ;
                        if(!ssi->stdev_3_len){
                                ssi->stdev_3_len = 10000.0;
                        }
                        //fprintf(stderr,"3: %f %f	%f\n", ssi->mean_3_len,  ssi->stdev_3_len,three_s0);
                        //if(ssi->stdev_3_len < 1){
                        //	ssi->stdev_3_len = 1;
                        //}
                        //fprintf(stderr,"3: %f %f	%f\n", ssi->mean_3_len,  ssi->stdev_3_len,three_s0);
                        if(ssi->mean_3_len <= 1){
                                sprintf(param->buffer,"WARNING: 3' partial segment seems not to be present in the data (length < 1).\n");
                                //	fprintf(stderr,"%s",param->buffer);
                                param->messages = append_message(param->messages, param->buffer);
                                //free_param(param);
                                //exit(EXIT_FAILURE);
                        }
                }
        }else{
                ssi->mean_3_len =  -1.0;
                ssi->stdev_3_len = -1.0;
        }

        if(param->matchstart!= -1 || param->matchend !=-1){
                ssi->average_length = (param->matchend - param->matchstart )* total_read;
        }
        ssi->average_length =  (int) floor((double)  ssi->average_length / (double) total_read   + 0.5);

        sum = 0.0;
        for(i = 0; i < 5;i++){
                sum += ssi->background[i];
        }

        for(i = 0; i < 5;i++){
                ssi->background[i] = prob2scaledprob(ssi->background[i]  / sum);
        }

        if(five_test_sequence){
                MFREE(five_test_sequence);
        }
        if(three_test_sequence){
                MFREE(three_test_sequence);
        }
#ifdef DEBUG
        fprintf(stderr,"Backgound:\n");
        fprintf(stderr,"A:%f\n",scaledprob2prob( ssi->background[0]));
        fprintf(stderr,"C:%f\n",scaledprob2prob( ssi->background[1]));
        fprintf(stderr,"G:%f\n",scaledprob2prob( ssi->background[2]));
        fprintf(stderr,"T:%f:\n",scaledprob2prob( ssi->background[3]));
        fprintf(stderr,"N:%f\n",scaledprob2prob( ssi->background[4]));

        fprintf(stderr,"\nExpected Length:5':%f	3':%f\n",ssi->expected_5_len,ssi->expected_3_len );
        fprintf(stderr,"Observed Length:5':%f	3':%f\n",ssi->mean_5_len,ssi->mean_3_len );
        fprintf(stderr,"STDEV:5':%f	3':%f\n",ssi->stdev_5_len,ssi->stdev_3_len );


#endif

        pclose(file);
        return ssi;
ERROR:
        Return NULL;
}
*/
