#include "hmm_model_bag.h"

#include "hmm.h"
#include "init_hmm.h"
struct model_bag* init_model_bag(struct parameters* param,struct sequence_stats_info* ssi)
{
        int i,j,c,status;
        //int average_length = 12;
        int read_length = 1;
        int segment_length;

        struct model_bag* mb = 0;

        struct  model* model_p = 0;
        //struct hmm_column* col = 0;

        MMALLOC(mb,sizeof(struct model_bag));
        mb->path = 0;
        mb->dyn_prog_matrix = 0;
        mb->transition_matrix = 0;
        mb->label = 0;
        mb->average_raw_length = ssi->average_length;
        mb->current_dyn_length = ssi->max_seq_len + 10;
        mb->model = 0;
        MMALLOC(mb->model,sizeof(struct model* ) * param->read_structure->num_segments);



        mb->f_score = prob2scaledprob(0.0f);
        mb->b_score = prob2scaledprob(0.0f);
        mb->num_models = param->read_structure->num_segments;
        // get read length estimate...
        read_length = ssi->average_length;
        //fprintf(stderr,"READlength: %d\n",read_length);
        for(i = 0; i < mb->num_models;i++){
                //mb->model[i] = malloc_model_according_to_read_structure(param->read_structure,i);
                //fprintf(stderr," %d\n",read_length );
                if(param->read_structure->type[i] == 'G'){
                        read_length = read_length -2;
                }else if(param->read_structure->type[i] == 'R'){
                }else if(param->read_structure->type[i] == 'P'){
                        read_length = read_length - (int)strlen(param->read_structure->sequence_matrix[i][0])/2; // Initial guess - we don't know how much of the linker is present at this stage
                }else{
                        //	fprintf(stderr,"%s : %d \n", param->read_structure->sequence_matrix[i][0], (int)strlen(param->read_structure->sequence_matrix[i][0]));
                        read_length = read_length - (int)strlen(param->read_structure->sequence_matrix[i][0]);
                }

                //	fprintf(stderr,"READlength: %d\n",read_length);

        }
        //fprintf(stderr,"READlength: %d\n",read_length);

        if(read_length < 20){
                read_length = 20; // the expected read length should never be lower than 20!
        }
        mb->total_hmm_num = 0;


        mb->previous_silent = NULL;
        mb->next_silent = NULL;

        MMALLOC(mb->previous_silent,sizeof(float) * mb->current_dyn_length );
        MMALLOC(mb->next_silent,sizeof(float) * mb->current_dyn_length );


        for(i = 0; i < mb->num_models;i++){
                mb->model[i] = malloc_model_according_to_read_structure(param->read_structure->numseq_in_segment[i],(int)strlen(param->read_structure->sequence_matrix[i][0]),mb->current_dyn_length);
                segment_length = 0;
                if(param->read_structure->type[i] == 'G'){
                        segment_length = 2;
                }
                if(param->read_structure->type[i]  == 'R'){
                        segment_length = read_length;
                }
                mb->model[i] = init_model_according_to_read_structure(mb->model[i], param, i,ssi->background,segment_length);
                //print_model(mb->model[i] );
                mb->total_hmm_num += mb->model[i]->num_hmms;
        }


        //exit(0);

        double sum_prob = 0.0;
        double temp1;
        // 1) setting 5' parameters...
        if(ssi->expected_5_len){
                /*sum_prob = 0;

                  for(i = 1; i <=  ssi->expected_5_len;i++){
                  sum_prob +=gaussian_pdf(i , ssi->mean_5_len, ssi->stdev_5_len);
                  }*/
                model_p = mb->model[0];
                //model_p->skip = prob2scaledprob(  gaussian_pdf(0 ,ssi->mean_5_len, ssi->stdev_5_len));

                //sum_prob +=gaussian_pdf(0 ,ssi->mean_5_len, ssi->stdev_5_len);

                //temp1 = prob2scaledprob(0.0);
                //temp1 = logsum(temp1, model_p->skip);*/

                sum_prob = prob2scaledprob(0.0);

                for(i = 0 ; i < model_p->num_hmms;i++){
                        for(j = 0; j < ssi->expected_5_len;j++){
                                //			col = model_p->hmms[i]->hmm_column[j];
                                //fprintf(stderr,"%d MEAN:%f	STDEV:%f	g:%f\n",j,ssi->mean_5_len - ssi->expected_5_len, ssi->stdev_5_len,gaussian_pdf(j ,ssi->expected_5_len-ssi->mean_5_len , ssi->stdev_5_len));

                                model_p->silent_to_M[i][j]  = prob2scaledprob(1.0 / (float) model_p->num_hmms) + prob2scaledprob(   gaussian_pdf(j ,ssi->expected_5_len-ssi->mean_5_len , ssi->stdev_5_len));
                                //fprintf(stderr,"%f	%f	%f	%f\n",sum_prob, model_p->silent_to_M[i][j],scaledprob2prob(sum_prob), scaledprob2prob( model_p->silent_to_M[i][j]));
                                sum_prob = logsum(sum_prob, model_p->silent_to_M[i][j]);
                                //fprintf(stderr,"5': %d %f\n",j,gaussian_pdf(j ,ssi->mean_5_len - ssi->expected_5_len, ssi->stdev_5_len)  );

                                //	temp1 = logsum(temp1, model_p->silent_to_M[i][j]);
                        }
                        model_p->hmms[i] = set_hmm_transition_parameters(model_p->hmms[i],ssi->expected_5_len, param->sequencer_error_rate, param->indel_frequency, -1.0, -1.0);
                }
                model_p->skip = prob2scaledprob(  gaussian_pdf(ssi->expected_5_len,ssi->mean_5_len - ssi->expected_5_len, ssi->stdev_5_len));
                //fprintf(stderr,"5': skip: %f\n",gaussian_pdf(ssi->expected_5_len,ssi->mean_5_len - ssi->expected_5_len, ssi->stdev_5_len)  );
                sum_prob = logsum(sum_prob, model_p->skip);



                //fprintf(stderr,"Sanity: %f	%f\n",sum_prob, scaledprob2prob(sum_prob));

                temp1 = prob2scaledprob(0.0);

                for(i = 0 ; i < model_p->num_hmms;i++){
                        for(j = 0; j < ssi->expected_5_len;j++){
                                //			col = model_p->hmms[i]->hmm_column[j];
                                model_p->silent_to_M[i][j]  = model_p->silent_to_M[i][j]  - sum_prob;
                                temp1 = logsum(temp1, model_p->silent_to_M[i][j]);
                                //		fprintf(stderr,"5': %d %f\n",j,scaledprob2prob(model_p->silent_to_M[i][j] ) );

                                //	temp1 = logsum(temp1, model_p->silent_to_M[i][j]);
                        }

                }
                model_p->skip = model_p->skip- sum_prob;
                //fprintf(stderr,"5': skip: %f\n",scaledprob2prob(model_p->skip )  );
                temp1 = logsum(temp1, model_p->skip);

                //fprintf(stderr,"Sanity: %f	%f\n",temp1, scaledprob2prob(temp1));

                /*model_p->skip = model_p->skip - temp1;
                  for(i = 0 ; i < model_p->num_hmms;i++){
                  for(j = 0; j < ssi->expected_5_len;j++){
                  model_p->silent_to_M[i][j]  = model_p->silent_to_M[i][j]  - temp1;
                  }
                  }*/
        }

        // 2) setting 3' parameters...
        if(ssi->expected_3_len){
                sum_prob = 0;

                for(i = 0; i <  ssi->expected_3_len;i++){
                        sum_prob +=gaussian_pdf(i , ssi->mean_3_len ,ssi->stdev_3_len);
                }
                model_p = mb->model[mb->num_models-1];
                model_p->skip = prob2scaledprob(  gaussian_pdf(0 , ssi->mean_3_len ,ssi->stdev_3_len) / sum_prob    );
                temp1 = model_p->skip;
                for(i = 0 ; i < model_p->num_hmms;i++){
                        model_p->silent_to_M[i][0] = prob2scaledprob(1.0 / (float) model_p->num_hmms) + prob2scaledprob(1.0 -  gaussian_pdf(0 ,ssi->mean_3_len ,ssi->stdev_3_len) / sum_prob );
                        model_p->hmms[i] = set_hmm_transition_parameters(model_p->hmms[i],ssi->expected_3_len, param->sequencer_error_rate,  param->indel_frequency,ssi->mean_3_len ,ssi->stdev_3_len);
                }
        }
        // 3) sets parameters fot all internal P segments (note the 1 -> num_models -1 )
        for(c = 1; c < mb->num_models-1;c++){
                if(param->read_structure->type[c] == 'P'){

                        model_p = mb->model[c];
                        int len = model_p->hmms[0]->num_columns;
                        for(i = 0 ; i < model_p->num_hmms;i++){
                                j = 0;
                                model_p->hmms[i] = set_hmm_transition_parameters(model_p->hmms[i],len, param->sequencer_error_rate,  param->indel_frequency, 0.1, -1.0);
                        }
                }
        }
        //for(i = 0; i < mb->num_models;i++){
        //	print_model(mb->model[i]);
        //}

        //exit(0);

        MMALLOC(mb->path,sizeof(int*) * mb->current_dyn_length);
        MMALLOC(mb->dyn_prog_matrix,sizeof(float*) * mb->current_dyn_length );

        for (i = 0; i < mb->current_dyn_length;i++){
                mb->path[i] = 0;
                mb->dyn_prog_matrix[i] = 0;
                MMALLOC(mb->path[i],sizeof(int)* (mb->total_hmm_num +1) );
                MMALLOC(mb->dyn_prog_matrix[i],sizeof(float) * (mb->total_hmm_num +1) );
        }

        MMALLOC(mb->transition_matrix,sizeof(float*) * (mb->total_hmm_num +1));
        MMALLOC(mb->label,sizeof(int) *  (mb->total_hmm_num +1));

        mb->model_multiplier = 1.0f;

        c = 0;
        for(i = 0; i < mb->num_models ;i++){
                mb->model_multiplier  *= mb->model[i]->num_hmms;
                for(j = 0; j < mb->model[i]->num_hmms;j++){
                        mb->label[c] = (j << 16) | i ;
                        if(mb->model[i]->skip != prob2scaledprob(0.0)){
                                mb->label[c]  |= 0x80000000;
                        }
                        c++;

                }
        }

        mb->model_multiplier = prob2scaledprob(mb->model_multiplier);

        for(i = 0; i < mb->total_hmm_num+1 ;i++){
                mb->transition_matrix[i]  = 0;
                MMALLOC(mb->transition_matrix[i] ,sizeof(float) * (mb->total_hmm_num +1));
                for(j = 0; j <  mb->total_hmm_num+1 ;j++){
                        mb->transition_matrix[i][j] = 0;
                }
        }


        for(i = 0; i < mb->total_hmm_num ;i++){
                c = 1;
                for(j = i+1; j <  mb->total_hmm_num ;j++){
                        mb->transition_matrix[i][j] = 0;



                        if(i == j){
                                mb->transition_matrix[i][j] = 1;
                        }

                        if((mb->label[i] & 0xFFFF)+1 == ((mb->label[j] & 0xFFFF) ) ){
                                mb->transition_matrix[i][j] = 1;
                        }


                        if(((mb->label[i] & 0xFFFF) < ((mb->label[j] & 0xFFFF) ) )&& c ){
                                mb->transition_matrix[i][j] = 1;
                        }

                        if(!(mb->label[j] & 0x80000000)){
                                c =0;
                        }
                        //fprintf(stderr,"%d, %d, %d %d\n ", j,   mb->label[j],mb->label[j] & 0xFFFF, (mb->label[j] >> 16) & 0x7FFF);
                }

                // remain in the same state....
                mb->transition_matrix[i][i] = 1;
        }
        return mb;
ERROR:
        //KSLIB_MESSAGE(status,"Something wrong in init_model_bag.\n");
        return NULL;
}



/** \fn void free_model_bag(struct model_bag* mb)

    \brief Frees whole HMM model.


    \param mb @ref model_bag.

*/

void free_model_bag(struct model_bag* mb)
{
        int i;


        //mb->transition_matrix = malloc(sizeof(float*) * (mb->total_hmm_num +1));
        //mb->label = malloc(sizeof(int) *  (mb->total_hmm_num +1));

        for (i = 0; i < mb->current_dyn_length;i++){
                MFREE(mb->path[i]);// = malloc(sizeof(int)* (mb->total_hmm_num +1) );
                MFREE(mb->dyn_prog_matrix[i]);// = malloc(sizeof(float) * (mb->total_hmm_num +1) );
        }

        MFREE(mb->path);// = malloc(sizeof(int*) * MAX_SEQ_LEN);
        MFREE(mb->dyn_prog_matrix);// = malloc(sizeof(float*) * MAX_SEQ_LEN );


        for(i = 0; i < mb->total_hmm_num+1 ;i++){
                MFREE(mb->transition_matrix[i]);//  = malloc(sizeof(float) * (mb->total_hmm_num +1));

        }
        MFREE(mb->transition_matrix);
        MFREE(mb->label);
        MFREE(mb->previous_silent);//,sizeof(float) * mb->current_dyn_length );
        MFREE(mb->next_silent);//,sizeof(float) * mb->current_dyn_length );


        for(i = 0; i < mb->num_models;i++){
                free_model(mb->model[i]);
        }


        MFREE(mb->model);// = malloc(sizeof(struct model* ) * param->read_structure->num_segments);


        MFREE(mb);// = malloc(sizeof(struct model_bag));
}
