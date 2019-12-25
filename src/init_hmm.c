
#include "tllogsum.h"

#include "init_hmm.h"
#include "core_hmm_functions.h"



struct model* init_model_according_to_read_structure(struct model* model,struct read_structure* rs ,struct alphabet* a, int key,const  double* background,int assumed_length)
{

        //struct read_structure* rs = param->read_structure;
        float base_error = 0.05;
        float indel_freq = 0.1;
        struct hmm_column* col =0;
        int i,j,c,len;
        int current_nuc;
        char* tmp = 0;

        for(i= 0;i < 5;i++){
                model->background_nuc_frequency[i]= background[i];
                //fprintf(stderr,"%f\n",background[i]);
        }

        for(i = 0; i < model->num_hmms;i++){
                len = model->hmms[i]->num_columns;
                tmp = rs->sequence_matrix[key][i];

                //sets emission probabilities...
                for(j = 0; j < len;j++){
                        col = model->hmms[i]->hmm_column[j];
                        col->transition_e[MM] =  prob2scaledprob(0.0);
                        col->transition_e[MI] =  prob2scaledprob(0.0);
                        col->transition_e[MD] =  prob2scaledprob(0.0);
                        col->transition_e[MSKIP] =  prob2scaledprob(0.0);

                        col->transition_e[II] =  prob2scaledprob(0.0);
                        col->transition_e[IM] =  prob2scaledprob(0.0);
                        col->transition_e[ISKIP] =  prob2scaledprob(0.0);

                        col->transition_e[DD] =  prob2scaledprob(0.0);
                        col->transition_e[DM] =  prob2scaledprob(0.0);
                        //fprintf(stdout,"%d %d code:%d\n", j, tmp[j],nuc_code[(int) tmp[j]]);

                        current_nuc = tlalphabet_get_code(a, tmp[j]);
                        col->identifier = -1;
                        if(current_nuc < 4){
                                // before distributed the error probabilityequally to all remaining 4 letters. Now: distribute to 3 and set probability to emit 'N' to background...
                                for(c = 0; c < 4;c++){
                                        if(c == current_nuc){
                                                col->m_emit[c] = prob2scaledprob(1.0 -scaledprob2prob(background[4])  - base_error* (1.0- indel_freq));
                                        }else{
                                                col->m_emit[c] =  prob2scaledprob( base_error* (1.0- indel_freq)/ 3.0);
                                        }
                                        col->i_emit[c] = background[c];
                                        col->i_emit_e[c] =  prob2scaledprob(0.0f);
                                        col->m_emit_e[c] =  prob2scaledprob(0.0f);
                                }
                                col->m_emit[4] =  background[4];
                                col->i_emit[4] = background[4];
                                col->i_emit_e[4] =  prob2scaledprob(0.0f);
                                col->m_emit_e[4] =  prob2scaledprob(0.0f);
                        }else if(current_nuc == 4){
                                for(c = 0; c < 5;c++){
                                        col->m_emit[c] =  background[c];
                                        col->i_emit[c] =  background[c];
                                        col->i_emit_e[c] =  prob2scaledprob(0.0f);
                                        col->m_emit_e[c] =  prob2scaledprob(0.0f);
                                }
                        }else{
                                current_nuc = 4;
                                for(c = 0; c < 5;c++){
                                        if(c == current_nuc){
                                                col->m_emit[c] = prob2scaledprob(1.0);
                                        }else{
                                                col->m_emit[c] =  prob2scaledprob(0.0);
                                        }
                                        col->i_emit[c] = background[c];
                                        col->i_emit_e[c] =  prob2scaledprob(0.0f);
                                        col->m_emit_e[c] =  prob2scaledprob(0.0f);
                                }
                        }
                }
                //sets transition probabilities....
                model->hmms[i] = set_hmm_transition_parameters(model->hmms[i],len, base_error, indel_freq, -1.0, -1.0);
                /*if(len == 1){
                //single state - only silent to / from M everything else set to zero....
                col = model->hmms[i]->hmm_column[0];
                col->transition[MM] = prob2scaledprob(0.0f);
                col->transition[MI] = prob2scaledprob(0.0f);
                col->transition[MD] = prob2scaledprob(0.0f);
                col->transition[MSKIP] = prob2scaledprob(1.0f);

                col->transition[II] = prob2scaledprob(0.0f);
                col->transition[IM] = prob2scaledprob(0.0f);
                col->transition[ISKIP] = prob2scaledprob(0.0f);

                col->transition[DD] = prob2scaledprob(0.0);
                col->transition[DM] = prob2scaledprob(0.0);
                }else if(len == 2){

                //first column
                col = model->hmms[i]->hmm_column[0];
                col->transition[MM] = prob2scaledprob( 1.0 - base_error * indel_freq);
                col->transition[MI] = prob2scaledprob(base_error * indel_freq);

                col->transition[MD] =  prob2scaledprob(0.0);
                col->transition[MSKIP] = prob2scaledprob(0.0);

                col->transition[II] = prob2scaledprob(1.0 - 0.999);
                col->transition[IM] = prob2scaledprob(0.999);
                col->transition[ISKIP] = prob2scaledprob(0.0f);

                col->transition[DD] = prob2scaledprob(0.0);
                col->transition[DM] = prob2scaledprob(0.0);

                //second column
                col = model->hmms[i]->hmm_column[1];
                col->transition[MM] = prob2scaledprob(0.0f );// 1.0 - base_error * indel_freq);
                col->transition[MI] = prob2scaledprob(0.0f );//(base_error * indel_freq) +  prob2scaledprob(0.5);
                col->transition[MD] = prob2scaledprob(0.0f );//(base_error * indel_freq) +  prob2scaledprob(0.5);
                col->transition[MSKIP] = prob2scaledprob(1.0);

                col->transition[II] = prob2scaledprob(0.00);
                col->transition[IM] = prob2scaledprob(0.0);
                col->transition[ISKIP] = prob2scaledprob(0.0f);

                col->transition[DD] = prob2scaledprob(0.0f );//(1.0 - 0.999);
                col->transition[DM] = prob2scaledprob(0.0f );//0.999);

                }else{
                //first column....
                col = model->hmms[i]->hmm_column[0];
                col->transition[MM] = prob2scaledprob( 1.0 - base_error * indel_freq);
                col->transition[MI] = prob2scaledprob(base_error * indel_freq) +  prob2scaledprob(0.5);
                col->transition[MD] = prob2scaledprob(base_error * indel_freq) +  prob2scaledprob(0.5);
                col->transition[MSKIP] = prob2scaledprob(0.0);

                col->transition[II] = prob2scaledprob(1.0 - 0.999);
                col->transition[IM] = prob2scaledprob(0.999);
                col->transition[ISKIP] = prob2scaledprob(0.0f);

                col->transition[DD] = prob2scaledprob(0.0);
                col->transition[DM] = prob2scaledprob(0.0);

                //middle columns...
                for(j = 1; j < len-2;j++){
                col = model->hmms[i]->hmm_column[j];

                col->transition[MM] = prob2scaledprob( 1.0 - base_error * indel_freq);
                col->transition[MI] = prob2scaledprob(base_error * indel_freq) +  prob2scaledprob(0.5);
                col->transition[MD] = prob2scaledprob(base_error * indel_freq) +  prob2scaledprob(0.5);
                col->transition[MSKIP] = prob2scaledprob(0.0);

                col->transition[II] = prob2scaledprob(1.0 - 0.999);
                col->transition[IM] = prob2scaledprob(0.999);
                col->transition[ISKIP] = prob2scaledprob(0.0f);


                col->transition[DD] = prob2scaledprob(1.0 - 0.999);
                col->transition[DM] = prob2scaledprob(0.999);
                }

                //second last...
                col = model->hmms[i]->hmm_column[len -2];
                col->transition[MM] = prob2scaledprob( 1.0 - base_error * indel_freq);
                col->transition[MI] = prob2scaledprob(base_error * indel_freq);

                col->transition[MD] =  prob2scaledprob(0.0);

                col->transition[MSKIP] = prob2scaledprob(0.0);

                col->transition[II] = prob2scaledprob(1.0 - 0.999);
                col->transition[IM] = prob2scaledprob(0.999);
                col->transition[ISKIP] = prob2scaledprob(0.0f);

                col->transition[DD] = prob2scaledprob(0.0);
                col->transition[DM] = prob2scaledprob(1.0);
                //col->transition[DD] = prob2scaledprob(0.0);
                //col->transition[DM] = prob2scaledprob(0.0);

                col = model->hmms[i]->hmm_column[len -1];

                col->transition[MM] = prob2scaledprob(0.0f );// 1.0 - base_error * indel_freq);
                col->transition[MI] = prob2scaledprob(0.0f );//(base_error * indel_freq) +  prob2scaledprob(0.5);
                col->transition[MD] = prob2scaledprob(0.0f );//(base_error * indel_freq) +  prob2scaledprob(0.5);
                col->transition[MSKIP] = prob2scaledprob(1.0);

                col->transition[II] = prob2scaledprob(0.00);
                col->transition[IM] = prob2scaledprob(0.0);
                col->transition[ISKIP] = prob2scaledprob(0.0f);

                col->transition[DD] = prob2scaledprob(0.0f );//(1.0 - 0.999);
                col->transition[DM] = prob2scaledprob(0.0f );//0.999);





                }*/

        }

        // init all probs to 0

        for(i = 0 ; i < model->num_hmms;i++){
                len = model->hmms[i]->num_columns;
                for(j = 0; j < len;j++){

                        model->silent_to_M[i][j] = prob2scaledprob(0.0f);
                        model->silent_to_I[i][j] = prob2scaledprob(0.0f);
                        model->silent_to_M_e[i][j] = prob2scaledprob(0.0f);
                        model->silent_to_I_e[i][j] = prob2scaledprob(0.0f);
                }
        }
        model->skip = prob2scaledprob(0.0f);
        model->skip_e = prob2scaledprob(0.0f);

        if(rs->type[key] == 'B'){// barcodes all have same length & equal prior probability...
                for(i = 0 ; i < model->num_hmms;i++){
                        model->silent_to_M[i][0] = prob2scaledprob(1.0 / (float) model->num_hmms);// + prob2scaledprob(0.9);
                        model->silent_to_I[i][0] = prob2scaledprob(0.0f);
                }
                model->skip = prob2scaledprob(0.0);
        }

        if(rs->type[key] == 'F'){// fingerprint all have same length & equal prior probability... (of course we specify 1 with NNNNNNNN
                for(i = 0 ; i < model->num_hmms;i++){
                        model->silent_to_M[i][0] = prob2scaledprob(1.0 / (float) model->num_hmms);
                        //model->M_to_silent[i] = prob2scaledprob(1.0);
                }
                model->skip = prob2scaledprob(0.0);
        }


        if(rs->type[key] == 'S'){// fingerprint all have same length & equal prior probability... (of course we specify 1 with NNNNNNNN
                for(i = 0 ; i < model->num_hmms;i++){
                        model->silent_to_M[i][0] = prob2scaledprob(1.0 / (float) model->num_hmms);
                        model->silent_to_I[i][0] = prob2scaledprob(0.0f);
                        //model->M_to_silent[i] = prob2scaledprob(1.0);
                }
                model->skip = prob2scaledprob(0.0);
        }

        if(rs->type[key] == 'P'){// Partial - can skip and exit at every M / I state....
                len = model->hmms[0]->num_columns;
                for(i = 0 ; i < model->num_hmms;i++){
                        model->silent_to_M[i][0] = prob2scaledprob(1.0 / (float) model->num_hmms) + prob2scaledprob(1.0 - 0.01);
                        //model->M_to_silent[i] = prob2scaledprob(1.0);

                        for(j = 0; j < len;j++){
                                col = model->hmms[i]->hmm_column[j];
                                col->transition[MM] = prob2scaledprob( 1.0 - base_error * indel_freq ) + prob2scaledprob(0.99f);
                                col->transition[MI] = prob2scaledprob(base_error * indel_freq) +  prob2scaledprob(0.5)+ prob2scaledprob(0.99f);
                                col->transition[MD] = prob2scaledprob(base_error * indel_freq) +  prob2scaledprob(0.5)+ prob2scaledprob(0.99f);
                                col->transition[MSKIP] = prob2scaledprob(0.01f);

                                col->transition[II] = prob2scaledprob(1.0 - 0.999) + prob2scaledprob(0.99f);
                                col->transition[IM] = prob2scaledprob(0.999) + prob2scaledprob(0.99f);
                                col->transition[ISKIP] = prob2scaledprob(0.01f);
                        }
                }
                model->skip = prob2scaledprob(0.01);
        }



        if(rs->type[key] == 'O'){ // optional - like a G, GG or GGG priot probability set to 0.5  - assume length 2 for now,
                len = model->hmms[0]->num_columns;
                for(i = 0 ; i < model->num_hmms;i++){
                        //model->silent_to_M[i] = prob2scaledprob(1.0 / (float) model->num_hmms) + prob2scaledprob(0.5);
                        //model->M_to_silent[i] = prob2scaledprob(1.0);

                        model->silent_to_I[i][0] = prob2scaledprob(1.0 / (float) model->num_hmms) + prob2scaledprob(0.5);
                        //model->I_to_silent[i] = prob2scaledprob(1.0 / (float) (len+1));

                        //len = model->hmms[i]->num_columns;
                        //tmp = rs->sequence_matrix[key][i];
                        for(j = 0; j < len;j++){
                                col = model->hmms[i]->hmm_column[j];
                                for(c = 0; c < 5;c++){
                                        col->i_emit[c] = col->m_emit[c];
                                        col->m_emit[c] = prob2scaledprob(0.0);
                                }
                        }
                }
                model->skip = prob2scaledprob(0.5);
                col = model->hmms[0]->hmm_column[0];
                col->transition[MM] = prob2scaledprob( 0.0 );
                col->transition[MI] = prob2scaledprob(0.0);
                col->transition[MD] = prob2scaledprob(0.0);
                col->transition[MSKIP] = prob2scaledprob(0.0);

                //col->transition[MQUIT] = prob2scaledprob(1.0 / (float) 2);

                col->transition[II] = prob2scaledprob(1.0 - 1.0 / (float)(len+1) );
                col->transition[IM] = prob2scaledprob(0.0);
                col->transition[ISKIP] =  prob2scaledprob(1.0 / (float) (len+1));


                col->transition[DD] = prob2scaledprob(0.0);
                col->transition[DM] = prob2scaledprob(0.0);


                col->transition_e[MM] =  prob2scaledprob(0.0);
                col->transition_e[MI] =  prob2scaledprob(0.0);
                col->transition_e[MD] =  prob2scaledprob(0.0);

                col->transition_e[II] =  prob2scaledprob(0.0);
                col->transition_e[IM] =  prob2scaledprob(0.0);

                col->transition_e[DD] =  prob2scaledprob(0.0);
                col->transition_e[DM] =  prob2scaledprob(0.0);



        }

        if(rs->type[key] == 'G'){ // optional - like a G, GG or GGG priot probability set to 0.5  - assume length 2 for now,
                len = model->hmms[0]->num_columns;
                for(i = 0 ; i < model->num_hmms;i++){
                        //model->silent_to_M[i] = prob2scaledprob(1.0 / (float) model->num_hmms) + prob2scaledprob(0.5);
                        //model->M_to_silent[i] = prob2scaledprob(1.0);

                        model->silent_to_I[i][0] = prob2scaledprob(0.8935878);
                        //model->I_to_silent[i] = prob2scaledprob(1.0 - 0.195);

                        //len = model->hmms[i]->num_columns;
                        //tmp = rs->sequence_matrix[key][i];
                        for(j = 0; j < len;j++){
                                col = model->hmms[i]->hmm_column[j];
                                for(c = 0; c < 5;c++){
                                        col->i_emit[c] = col->m_emit[c];
                                        col->m_emit[c] = prob2scaledprob(0.0);
                                }
                        }
                }
                model->skip = prob2scaledprob(1.0 - 0.8935878);
                col = model->hmms[0]->hmm_column[0];
                col->transition[MM] = prob2scaledprob( 0.0 );
                col->transition[MI] = prob2scaledprob(0.0);
                col->transition[MD] = prob2scaledprob(0.0);

                //col->transition[MQUIT] = prob2scaledprob(1.0 / (float) 2);

                col->transition[II] = prob2scaledprob(0.195);
                col->transition[IM] = prob2scaledprob(0.0);

                col->transition[DD] = prob2scaledprob(0.0);
                col->transition[DM] = prob2scaledprob(0.0);


                col->transition_e[MM] =  prob2scaledprob(0.0);
                col->transition_e[MI] =  prob2scaledprob(0.0);
                col->transition_e[MD] =  prob2scaledprob(0.0);

                col->transition_e[II] =  prob2scaledprob(0.0);
                col->transition_e[IM] =  prob2scaledprob(0.0);

                col->transition_e[DD] =  prob2scaledprob(0.0);
                col->transition_e[DM] =  prob2scaledprob(0.0);
        }

        if(rs->type[key] == 'R'){// read - skip impossible;
                for(i = 0 ; i < model->num_hmms;i++){
                        model->silent_to_I[i][0] = prob2scaledprob(1.0 / (float) model->num_hmms);
                        //		model->I_to_silent[i] = prob2scaledprob(1.0 / (float) assumed_length);
                }
                col = model->hmms[0]->hmm_column[0];
                for(c = 0; c < 5;c++){

                        col->m_emit[c] =background[c];

                        col->i_emit[c] = background[c];
                        col->i_emit_e[c] =  prob2scaledprob(0.0f);
                        col->m_emit_e[c] =  prob2scaledprob(0.0f);
                }
                col->transition[MM] = prob2scaledprob( 0.0);
                col->transition[MI] = prob2scaledprob(0.0);
                col->transition[MD] = prob2scaledprob(0.0);
                col->transition[MSKIP] = prob2scaledprob(0.0);

                //col->transition[MQUIT] = prob2scaledprob(1.0 / (float) assumed_length);
                col->transition[II] = prob2scaledprob(1.0 - 1.0 / (float) assumed_length );
                col->transition[IM] = prob2scaledprob(0.0);
                col->transition[ISKIP] = prob2scaledprob(1.0 / (float) assumed_length);

                col->transition[DD] = prob2scaledprob(0.0);
                col->transition[DM] = prob2scaledprob(0.0);


                col->transition_e[MM] =  prob2scaledprob(0.0);
                col->transition_e[MI] =  prob2scaledprob(0.0);
                col->transition_e[MD] =  prob2scaledprob(0.0);

                col->transition_e[II] =  prob2scaledprob(0.0);
                col->transition_e[IM] =  prob2scaledprob(0.0);

                col->transition_e[DD] =  prob2scaledprob(0.0);
                col->transition_e[DM] =  prob2scaledprob(0.0);

                model->skip = prob2scaledprob(0.0);

        }
        return model;
}


struct model* malloc_model_according_to_read_structure(int num_hmm, int length,int dyn_length)
{
        struct model* model = NULL;
        //int status;
        int i = 0;
        int j = 0;
        int len = 0;

        MMALLOC(model,sizeof(struct model));

        model->average_length =0;
        model->hmms =0;
        model->silent_backward =0;
        model->silent_forward = 0;
        model->silent_to_I = 0;
        model->silent_to_I_e = 0;
        model->silent_to_M = 0;
        model->silent_to_M_e = 0;

        model->num_hmms = num_hmm;// (rs->numseq_in_segment[key]);
        MMALLOC(model->hmms,sizeof(struct hmm*) * model->num_hmms  );//(rs->numseq_in_segment[key]));

        for(i = 0; i < model->num_hmms;i++){
                model->hmms[i] = 0;
                MMALLOC(model->hmms[i],sizeof(struct hmm) );
        }

        MMALLOC(model->silent_to_M,sizeof(float*) * model->num_hmms);
        MMALLOC(model->silent_to_I,sizeof(float*) * model->num_hmms);
        MMALLOC(model->silent_to_M_e,sizeof(float*) * model->num_hmms);
        MMALLOC(model->silent_to_I_e,sizeof(float*) * model->num_hmms);
        MMALLOC(model->silent_forward,sizeof(float) * (dyn_length+1));
        MMALLOC(model->silent_backward,sizeof(float) * (dyn_length+1));

        len = length;// (int)strlen(rs->sequence_matrix[key][0]);

        for(i = 0 ;i  < model->num_hmms;i++){
                model->silent_to_M[i] = 0;
                model->silent_to_M_e[i] = 0;
                model->silent_to_I[i] = 0;
                model->silent_to_I_e[i] = 0;
                MMALLOC(model->silent_to_M[i],sizeof(float) * len);
                MMALLOC(model->silent_to_I[i],sizeof(float) * len);

                MMALLOC(model->silent_to_M_e[i],sizeof(float) * len);
                MMALLOC(model->silent_to_I_e[i],sizeof(float) * len);

                for(j = 0; j < len;j++){
                        model->silent_to_M[i][j] =  0.0f;
                        model->silent_to_I[i][j] = 0.0f;
                        model->silent_to_M_e[i][j] =  0.0f;
                        model->silent_to_I_e[i][j] = 0.0f;
                }
        }

        for(i = 0; i < model->num_hmms;i++){
                model->hmms[i]->num_columns = len;
                model->hmms[i]->hmm_column = 0;
                MMALLOC(model->hmms[i]->hmm_column,sizeof(struct hmm_column*) * len);
                for(j = 0; j < len;j++){
                        model->hmms[i]->hmm_column[j] = 0;


                        MMALLOC(model->hmms[i]->hmm_column[j] ,sizeof(struct hmm_column));
                        model->hmms[i]->hmm_column[j]->D_backward = 0;
                        model->hmms[i]->hmm_column[j]->D_foward = 0;
                        model->hmms[i]->hmm_column[j]->I_backward = 0;
                        model->hmms[i]->hmm_column[j]->I_foward = 0;
                        model->hmms[i]->hmm_column[j]->M_backward = 0;
                        model->hmms[i]->hmm_column[j]->M_foward = 0;

                        MMALLOC(model->hmms[i]->hmm_column[j]->M_foward,sizeof(float) * (dyn_length+1));
                        MMALLOC(model->hmms[i]->hmm_column[j]->M_backward ,sizeof(float) * (dyn_length+1));
                        MMALLOC(model->hmms[i]->hmm_column[j]->I_foward,sizeof(float) * (dyn_length+1));
                        MMALLOC(model->hmms[i]->hmm_column[j]->I_backward,sizeof(float) * (dyn_length+1));
                        MMALLOC(model->hmms[i]->hmm_column[j]->D_foward,sizeof(float) * (dyn_length+1));
                        MMALLOC(model->hmms[i]->hmm_column[j]->D_backward, sizeof(float) * (dyn_length+1));
                }
        }
        return model;
ERROR:
        //KSLIB_MESSAGE(status,"Something wrong in malloc_model_according_to_read_structure.\n");
        return NULL;
}




void free_model(struct model* model)
{
        int i = 0;
        int j = 0;


        for(i = 0; i < model->num_hmms;i++){

                for(j = 0; j < model->hmms[i]->num_columns;j++){
                        MFREE(model->hmms[i]->hmm_column[j]->M_foward);// = malloc(sizeof(float) * (dyn_length+1));
                        MFREE(model->hmms[i]->hmm_column[j]->M_backward);// =malloc(sizeof(float) * (dyn_length+1));
                        MFREE(model->hmms[i]->hmm_column[j]->I_foward);//  = malloc(sizeof(float) * (dyn_length+1));
                        MFREE(model->hmms[i]->hmm_column[j]->I_backward);// = malloc(sizeof(float) * (dyn_length+1));
                        MFREE(model->hmms[i]->hmm_column[j]->D_foward);// = malloc(sizeof(float) * (dyn_length+1));
                        MFREE(model->hmms[i]->hmm_column[j]->D_backward);// = malloc(sizeof(float) * (dyn_length+1));
                        MFREE(model->hmms[i]->hmm_column[j]);// = malloc(sizeof(struct hmm_column));
                }
                //model->hmms[i]->num_columns = sub_length;
                MFREE(model->hmms[i]->hmm_column);// = malloc(sizeof(struct hmm_column*) * sub_length);
        }

        for(i = 0; i < model->num_hmms;i++){
                MFREE(model->hmms[i]);// = malloc(sizeof(struct hmm) );
        }

        MFREE(model->hmms);// = malloc(sizeof(struct hmm*) * (1+ number_sub_models));

        for(i = 0; i < model->num_hmms;i++){
                MFREE(model->silent_to_M[i]);
                MFREE(model->silent_to_M_e[i]);
                MFREE(model->silent_to_I[i]);
                MFREE(model->silent_to_I_e[i]);
        }
        MFREE(model->silent_to_M);
        MFREE(model->silent_to_M_e);
        MFREE(model->silent_to_I);
        MFREE(model->silent_to_I_e);
        //}
        MFREE(model->silent_forward);// = malloc(sizeof(float) * (dyn_length+1));
        MFREE(model->silent_backward);// = malloc(sizeof(float) * (dyn_length+1));
        MFREE(model);// = malloc(sizeof(struct model));
}
