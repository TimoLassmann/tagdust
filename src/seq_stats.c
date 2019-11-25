#include "seq_stats.h"

#include "tldevel.h"
#include "interface.h"
#include "nuc_code.h"

struct sequence_stats_info* get_sequence_stats(struct parameters* param, struct read_info** ri, int file_num )
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
        return NULL;
}
