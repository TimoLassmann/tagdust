#include "math.h"


#include "tldevel.h"


#include "arch_lib_sim.h"




static int set_len_of_unknown_rs(const struct read_structure*rs, int* min_plus_len,int* max_plus_len, int min_seq_len, int max_seq_len);


int emit_from_rs(const struct read_structure* rs, struct rng_state* rng, struct alphabet* a,  uint8_t** seq, uint8_t** qual, int* len, int sim_len)
{
        struct segment_specs* spec = NULL;
        uint8_t* s = NULL;
        uint8_t* q = NULL;
        int l;
        int m;
        int i,j,b,g;
        int plus_min_len, plus_max_len;
        RUN(set_len_of_unknown_rs(rs,&plus_min_len,&plus_max_len, sim_len,sim_len));

        m = 128;
        l = 0;

        MMALLOC(s, sizeof(uint8_t) * m);
        MMALLOC(q, sizeof(uint8_t) * m);


        for(i = 0; i < rs->num_segments;i++){
                spec = rs->seg_spec[i];
                //print_segment_spec(spec);
                if(spec->max_len == INT32_MAX){
                        g = tl_random_int(rng, plus_max_len - plus_min_len) + plus_min_len;
                }else{
                        g = tl_random_int(rng, spec->max_len- spec->min_len) + spec->min_len;
                }
                //LOG_MSG("G: %d", g);
                b = 0;
                if(spec->num_seq >1){
                        b = tl_random_int(rng, spec->num_seq-1) + 1;
                }
                for(j = 0; j < g;j++){
                        //LOG_MSG("j:%d",j);
                        if(spec->max_len == INT32_MAX){
                                s[l] = tl_random_int(rng, 4);
                        }else{

                                if( spec->seq[b][j] == 'N'){
                                        s[l] = tl_random_int(rng, 4);

                                }else{
                                        s[l] = spec->seq[b][j];
                                }

                        }
                        q[l] = tl_random_gaussian(rng, 20.0, 2.0);
                        //s->label[l] = etype[spec->extract];
                        l++;
                        if(l == m){
                                m = m + m /2;
                                MREALLOC(s, sizeof(uint8_t) * m);
                                MREALLOC(q, sizeof(uint8_t) * m);

                        }
                }
        }
        *seq = s;
        *qual = q;
        *len = l;
        return OK;
ERROR:
        return FAIL;
}

int set_len_of_unknown_rs(const struct read_structure*rs, int* min_plus_len,int* max_plus_len, int min_seq_len, int max_seq_len)
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
