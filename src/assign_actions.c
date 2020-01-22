

#include "tllogsum.h"
#include "tldevel.h"

#include "assign_data.h"

static int qsort_seq_bit_vec(const void *a, const void *b);

int ref_correct(khash_t(exact) * h ,struct qsubscore* subm, struct seq_bit* sb, int q_offset);

int post_process_assign(struct assign_struct* as)
{
        struct seq_bit_vec* bv = NULL;
        struct seq_bit* bit = NULL;
        struct demux_struct* tmp_ptr = NULL;
        struct bit_annotation* ann = NULL;
        uint8_t etype;
        char* code =  NULL;
        //char* barcode;
        //char* umi;
        int i,j,c;//g;
        int correct_index;
        int len;
        //int umi_len;
        int tmp_len = 256;

        MMALLOC(code, sizeof(char) * tmp_len);

        for(i = 0; i < as->num_reads ;i++){
                bv = as->bit_vec[i];
                ///* qsort by file identifier  */
                //qsort(bv->bits, bv->num_bit,sizeof(struct seq_bit*), qsort_seq_bits_by_type_file);
                len = 0;
                //umi_len = 0;
                //fprintf(stdout,"%s ",bv->name);
                for(j = 0; j < bv->num_bit;j++){


                        etype = bv->bits[j]->type;
                        switch (etype) {
                        case ARCH_ETYPE_SPLIT: {
                                code[len] = bv->bits[j]->code;
                                //len += bv->bits[j]->len;// strnlen(bv->bits[j]->p,256);
                                len++;
                                break;
                        }
                        case ARCH_ETYPE_APPEND:{
                                bit = bv->bits[j];
                                ann = as->bit_ann[j];
                                //LOG_MSG("%s %s %s", ann->name.s, ann->c_name.s,ann->q_name.s);
                                kputs(ann->name.s, &bv->append);
                                kputs(":Z:", &bv->append);
                                kputs(bit->p.s , &bv->append);
                                kputc(' ', &bv->append);
                                break;
                        }
                        case ARCH_ETYPE_APPEND_CORRECT: {
                                        //LOG_MSG("correcyt");
                                bit = bv->bits[j];
                                ann = as->bit_ann[j];
                                //LOG_MSG("%s %s %s", ann->name.s, ann->c_name.s,ann->q_name.s);
                                correct_index =  bit->correct_index;
                                RUN(ref_correct(as->bit_ann[j]->bar_hash,as->subm, bit, 33));

                                kputs(ann->name.s, &bv->append);
                                kputs(":Z:", &bv->append);
                                kputs(bit->p.s , &bv->append);
                                kputc(' ', &bv->append);

                                kputs(ann->c_name.s  , &bv->append);
                                kputs(":Z:", &bv->append);
                                kputs(bit->p_corr.s , &bv->append);
                                kputc(' ', &bv->append);

                                kputs(ann->q_name.s  , &bv->append);
                                kputs(":Z:", &bv->append);
                                kputs(bit->q.s , &bv->append);
                                kputc(' ', &bv->append);
                                break;

                        }
                        default:
                                break;
                        }
                        bv->fail |= bv->bits[j]->fail;
                }
                //fprintf(stdout,"\n");
                code[len] =0;
                c = 0;
                for(j = 0; j < bv->num_bit;j++){
                        //LOG_MSG("Bit %d %d %d", j,bv->bits[j]->type,bv->sample_group);
                        if(bv->bits[j]->type == ARCH_ETYPE_EXTRACT){
                                code[len] = c+33;
                                //len += bv->bits[j]->len;// strnlen(bv->bits[j]->p,256);
                                len++;
                                code[len] =0;
                                //fprintf(stdout,"CODE::: %s\n",code);
                                RUNP(tmp_ptr = as->demux_names->tree_get_data(as->demux_names,code));
                                //LOG_MSG("Setting: %d to id: %d code: %s %s", c,tmp_ptr->id,code,tmp_ptr->out_filename);
                                bv->out_file_id[c] = tmp_ptr->id;
                                len--;
                                c++;
                        }
                }
                //exit(0);
                /*RUNP(tmp_ptr = as->demux_names->tree_get_data(as->demux_names,tmp));
                        //if(tmp_ptr){
                        tmp_ptr->count++;
                        //fprintf(stdout,"%s -> %d (%d)\n", tmp_ptr->name ,tmp_ptr->id,tmp_ptr->count);
                        bv->sample_group = tmp_ptr->id;
                if(len){
                        g = 0;
                        for(j = 0; j < bv->num_bit;j++){
                                if(bv->bits[j]->type == BAR_TYPE){
                                        barcode = bv->bits[j]->p;
                                        len = bv->bits[j]->len;
                                        //len = strnlen(barcode,256);
                                        for(c = 0; c < len;c++){
                                                tmp[g] = barcode[c];
                                                g++;
                                                if(g == tmp_len){
                                                        tmp_len = tmp_len + tmp_len /2;
                                                        MREALLOC(tmp,sizeof(char) * tmp_len);
                                                }
                                        }
                                        tmp[g] = '_';
                                        g++;
                                        if(g == tmp_len){
                                                tmp_len = tmp_len + tmp_len /2;
                                                MREALLOC(tmp,sizeof(char) * tmp_len);
                                        }
                                }
                        }
                        tmp[g-1] = 0;

                        //LOG_MSG("searchfor >%s<\n",tmp);
                        RUNP(tmp_ptr = as->demux_names->tree_get_data(as->demux_names,tmp));
                        //if(tmp_ptr){
                        tmp_ptr->count++;
                        //fprintf(stdout,"%s -> %d (%d)\n", tmp_ptr->name ,tmp_ptr->id,tmp_ptr->count);
                        bv->sample_group = tmp_ptr->id;
                        //bv->bc = tmp;
                        //tmp = NULL;
                }else{
                        bv->sample_group = 0;
                }
                if(umi_len){
                        tmp = NULL;
                        MMALLOC(tmp, sizeof(char) * umi_len);
                        g = 0;
                        for(j = 0; j < bv->num_bit;j++){
                                if(bv->bits[j]->type == UMI_TYPE){
                                        umi = bv->bits[j]->p;
                                        len = bv->bits[j]->len;
                                        for(c = 0; c < len;c++){
                                                tmp[g] = umi[c];
                                                g++;
                                        }
                                        tmp[g] = '_';
                                        g++;
                                }
                        }
                        tmp[g-1] = 0;
                        bv->umi = tmp;
                        tmp = NULL;
                }
                if(bv->sample_group == -1){
                        fprintf(stdout,"%d %d\n",i, bv->num_bit);
                        for(j = 0; j < bv->num_bit;j++){
                                fprintf(stdout,"bit%d: type: %d\n", j, bv->bits[j]->type);
                                if(bv->bits[j]->type == BAR_TYPE){
                                        barcode = bv->bits[j]->p;
                                        fprintf(stdout,"%s\n",barcode);
                                }
                        }
                        ERROR_MSG("asda");
                }
                */
        }
        //LOG_MSG("%d %d", e_fail,d_fail);
        MFREE(code);

        qsort(as->bit_vec,as->num_reads ,sizeof(struct seq_bit_vec*) , qsort_seq_bit_vec);
        return OK;
ERROR:
        return FAIL;
}

int ref_correct(khash_t(exact) * h ,struct qsubscore* subm, struct seq_bit* sb, int q_offset)
{
        float scores[48]; /* 16 (maxlen * 3 for every edit ) */
        uint8_t pos[48];
        khiter_t k;
        uint32_t i,j,l;
        uint32_t edit;
        uint32_t search;
        uint32_t key;
        uint32_t m_index;
        float s;
        float max;
        l = MACRO_MIN(16, sb->p.l);
        key = seq_to_code(sb->p.s , l);
        //LOG_MSG("Searching for: ");

        /* No error */
        k = kh_get(exact, h, key);
        if(k != kh_end(h)){
                kputsn( sb->p.s, sb->p.l , &sb->p_corr);
                return OK;
        }
        /* 1 error */
        m_index = 0;
        for(i = 0; i < l;i++){
                for(edit = 1; edit < 4;edit++){
                        search = key ^ ( edit << (i << 1));
                        k = kh_get(exact, h, search);
                        if(k != kh_end(h)){
                                //LOG_MSG("Found! with 1M at %d  %d -> %d",i, (key >> (i << 1)) & 3, (search >> (i << 1)) & 3);
                                //code_to_seq(key, 16);
                                //code_to_seq(search, 16);
                                //viterbi



                                //LOG_MSG("Scoring:");
                                s = prob2scaledprob(1.0);
                                for(j =0; j < l;j++){
                                        //fprintf(stdout,"%c", "ACGT"[(key >> (j << 1)) & 3]);
                                        s+= get_qsubscore(subm,  (search >> (j << 1)) & 3,  (key >> (j << 1)) & 3, sb->q.s[j]- q_offset);
                                }
                                /*fprintf(stdout,"\n");
                                for(j =0; j < l;j++){
                                        fprintf(stdout,"%c", "ACGT"[(search >> (j << 1)) & 3]);
                                }
                                fprintf(stdout,"\n");
                                for(j =0; j < l;j++){
                                        fprintf(stdout,"%c", sb->q.s[j]);
                                }


                                fprintf(stdout,"\t%f\n", s);*/
                                pos[m_index] =  (edit << 6u) | i;
                                scores[m_index] = s;
                                m_index++;
                        }
                }
        }
        if(m_index){
                //LOG_MSG("DONE");
                s= prob2scaledprob(0.0);
                for(i = 0; i < m_index;i++){
                        s = logsum(s, scores[i]);
                }
                max = 0.0;
                edit = 0;
                for(i = 0; i < m_index;i++){
                        scores[i] = scaledprob2prob(scores[i] / s);
                        if(scores[i] > max){
                                edit = i;
                                max = scores[i];
                        }

//      LOG_MSG("%d %f",i, scaledprob2prob(scores[i]  -s ));
                }

                if(max > 0.9){
                        i = pos[edit] & 0x3F;
                        edit = (pos[edit] >> 6u) & 3;
                        search = key ^ ( edit << (i << 1));
                        for(j =0; j < l;j++){
                                kputc("ACGT"[(search >> (j << 1)) & 3], &sb->p_corr);
                        }
                }else{
                        sb->fail |= READ_FAILC;
                        LOG_MSG("%d candidates found; best score is %f", m_index, max);
                }
        }else{
                kputsn("NA", 2,&sb->p_corr);
                sb->fail |= READ_FAILC;
        }
        return OK;
}




int qsort_seq_bit_vec(const void *a, const void *b)
{
        const struct seq_bit_vec **elem1 = (const struct seq_bit_vec**) a;
        const struct seq_bit_vec **elem2 = (const struct seq_bit_vec**) b;
        if ( (*elem1)->fail < (*elem2)->fail){
                return -1;
        }else if ((*elem1)->fail > (*elem2)->fail){
                return 1;
        }else{

                if((*elem1)->out_file_id[0] > (*elem2)->out_file_id[0]){
                //if ( (*elem1)->sample_group > (*elem2)->sample_group){
                        return 1;
                }else if((*elem1)->out_file_id[0] < (*elem2)->out_file_id[0]){
                        //}else if ((*elem1)->sample_group < (*elem2)->sample_group){
                        return -1;
                }else{
                        return 0;
                }
        }
}




