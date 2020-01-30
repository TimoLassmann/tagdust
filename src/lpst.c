#include "lpst.h"

#include "tldevel.h"
#include "tlseqio.h"
#include "tllogsum.h"

#include "pst_structs.h"
struct lpst_model{
        struct pst** pst;
        uint8_t* jmptbl;
        int len;

};


static int set_len_of_unknown(const struct read_structure*rs, int* min_plus_len,int* max_plus_len, int min_seq_len, int max_seq_len);

static inline int nuc_to_internal(const char c);

int lpst_score_read(struct tl_seq_buffer* sb, struct read_structure* rs, struct sequence_stats_info* si, float* score)
{
        struct lpst_model* lpst = NULL;
        struct segment_specs* spec;

        int i,j,c;
        int n_pst;
        int plus_min_len;
        int plus_max_len;
        float** s_prob;
        int** s;

        float P_S;
        int l,pos,n;

        init_logsum();
        RUN(set_len_of_unknown(rs,&plus_min_len,&plus_max_len, si->min_seq_len,si->max_seq_len));
        if(plus_max_len == -1 || plus_min_len == -1){
                WARNING_MSG("LPST: Model too long from sequences.");
                *score = prob2scaledprob(0.0);
                return OK;
        }
        n_pst = 0;
        for(i= 0; i < rs->num_segments;i++){
                spec = rs->seg_spec[i];
                if(spec->extract == ARCH_ETYPE_APPEND_CORRECT){
                        n_pst++;
                }
        }
        if(!n_pst){

                P_S = prob2scaledprob(1.0);
                for(i = 0; i < sb->num_seq;i++){
                        for(j = 0; j < sb->sequences[i]->len;j++){
                                l = nuc_to_internal(sb->sequences[i]->seq[j]);
                                P_S += si->background[l];
                        }
                }
                *score= P_S;
                return OK;
        }

        MMALLOC(lpst, sizeof(struct lpst_model));

        lpst->pst = NULL;
        lpst->len =  si->max_seq_len;
        lpst->jmptbl = NULL;
        MMALLOC(lpst->pst, sizeof(struct pst* ) * n_pst);
        MMALLOC(lpst->jmptbl, sizeof(uint8_t) * lpst->len);

        c = 0;
        n_pst = 0;
        for(i= 0; i < rs->num_segments;i++){
                spec = rs->seg_spec[i];
                if(spec->extract == ARCH_ETYPE_APPEND_CORRECT){

                        if(spec->max_len == INT32_MAX  ){
                                for(j = 0; j < plus_max_len;j++){
                                        lpst->jmptbl[c] = n_pst;
                                        c++;
                                }
                        }else{
                                for(j = 0; j < spec->max_len;j++){
                                        lpst->jmptbl[c] = n_pst;
                                        c++;
                                }
                        }
                        lpst->pst[n_pst] = spec->pst;
                        n_pst++;
                }else{
                        if(spec->max_len == INT32_MAX  ){
                                for(j = 0; j < plus_max_len;j++){
                                        lpst->jmptbl[c] = 255;
                                        c++;
                                }
                        }else{
                                for(j = 0; j < spec->max_len;j++){
                                        lpst->jmptbl[c] = 255;
                                        c++;
                                }
                        }
                }
        }
        P_S = prob2scaledprob(1.0);
        for(i = 0; i < sb->num_seq;i++){
                for(j = 0; j < sb->sequences[i]->len;j++){
                        l = nuc_to_internal(sb->sequences[i]->seq[j]);
                        if(lpst->jmptbl[j] == 255){
                                P_S += si->background[l];
                        }else{

                                s = lpst->pst[lpst->jmptbl[j]]->fpst_root->links;
                                s_prob = lpst->pst[lpst->jmptbl[j]]->fpst_root->prob;
                                pos = j;
                                n = 0;
                                while(pos){
                                        pos= pos-1;
                                        c = nuc_to_internal(sb->sequences[i]->seq[pos]);
                                        if(!s[n][c]){
                                                break;
                                        }
                                        n = s[n][c];
                                }
                                P_S += s_prob[n][l];
                        }
                }
        }
        *score = P_S;


        MFREE(lpst->pst);
        MFREE(lpst->jmptbl);

        MFREE(lpst);
        return OK;
ERROR:
        return FAIL;
}


int set_len_of_unknown(const struct read_structure*rs, int* min_plus_len,int* max_plus_len, int min_seq_len, int max_seq_len)
{
        int i;
        int known_len;
        int num_unknown;
        int min;

        known_len = 0;
        num_unknown = 0;

        for(i = 0; i < rs->num_segments  ;i++){
                //LOG_MSG("%d: %d->%d %s", i, rs->seg_spec[i]->min_len, rs->seg_spec[i]->max_len, rs->seg_spec[i]->seq[0]);
                if(rs->seg_spec[i]->max_len == INT32_MAX){
                        num_unknown++;
                }else{
                        known_len += rs->seg_spec[i]->max_len;
                }
        }
        *min_plus_len = 0;
        *max_plus_len = 0;
        //LOG_MSG("NumUnknown: %d len known:%d    %d %d ", num_unknown,known_len, min_seq_len, max_seq_len);


        if(num_unknown){
                if(known_len >= max_seq_len){
                        /* oh dear sequence too short to fit this models  */
                        *max_plus_len = -1;
                        *min_plus_len = -1;
                }else{
                        *max_plus_len = (int) ceilf (((float) max_seq_len - (float) known_len) / (float) num_unknown);
                        *min_plus_len = (int) floorf(((float) min_seq_len - (float) known_len) / (float) num_unknown);
                }
        }
        return OK;
ERROR:
        return FAIL;
}




static inline int nuc_to_internal(const char c)
{
        switch (c) {
        case 'A':
        case 'a':
                return 0;
                break;
        case 'C':
        case 'c':
                return 1;
                break;
        case 'G':
        case 'g':
                return 2;
                break;
        case 'T':
        case 't':
                return 3;
                break;
        case 'N':
        case 'n':
                return 0;
                break;
        default:
                return 0;
                break;
        }
        return -1;
}


