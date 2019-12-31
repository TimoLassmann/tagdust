
#include "tldevel.h"

#include "tlseqio.h"
#include "tllogsum.h"
#include <string.h>


#include "pst.h"
#define MAX_PST_LEN 12


struct pst_node{
        struct pst_node* next[4];
        float nuc_probability[4];
        char* label;
};

struct pst {
        struct pst_node* pst_root;
        struct pst_node* ppt_root;
        uint32_t** counts;

        float p_min;
        float gamma_min;
        float r;
        int L;
};




static int init_pst(struct pst** pst, struct tl_seq_buffer* sb);

static float get_pst_prob(struct pst_node* n, char* string,int target, int pos);
static float get_ppt_prob(struct pst_node* n, char* string,int target, int pos);


static struct pst_node* build_pst(struct pst* pst,struct pst_node* n );
static struct pst_node* build_ppt(struct pst* pst,struct pst_node* n );
//static void print_pst(struct pst* pst,struct pst_node* n,int* num);
static void free_pst_node(struct pst_node* n);

static struct pst_node* alloc_node(struct pst_node* n,char* string,int len);
//static void free_node(struct pst_node*n);

static inline int nuc_to_internal(const char c);



int scan_read_with_pst(struct pst* pst, char* seq, int len, float*r)
{
        int i;
        int  c;
        //float P_PT;
        float P_S;
        float P_P;
        float P_R;

        float* base_p = pst->pst_root->nuc_probability;

        float A;
        float B;

        P_S = prob2scaledprob(1.0);
        P_P = prob2scaledprob(1.0);
        P_R = prob2scaledprob(1.0);
        //fprintf(stdout,"%s - query sequence\n", seq);
        for(i = 0; i < len; i++ ){
                c = nuc_to_internal(seq[i]);
                //fprintf(stdout,"%f %f %d\n", P,P_R,c);
                P_R = P_R + prob2scaledprob(base_p[c]);

                A = get_pst_prob(pst->pst_root, seq,c, i);



                P_S = P_S + prob2scaledprob(A);

                B = get_ppt_prob(pst->ppt_root, seq,c, i);
                P_P = P_P + prob2scaledprob(B);
                //LOG_MSG("%c %f %f %f %f ",seq[i],A,B, P_S,P_P);
        }
        /*A = exp2f(P_S-P_R) / (1.0 + exp2f(P_S-P_R));
        fprintf(stdout,"%f\t%f,%s\n", P_S - P_R,A, seq);
        A = exp2f(P_P-P_R) / (1.0 + exp2f(P_P-P_R));
        fprintf(stdout,"%f\t%f,%s\n", P_P - P_R,A, seq);*/
        A = MACRO_MAX(P_S-P_R, P_P-P_R);
        A = exp2f(A) / (1.0 + exp2f(A));
        *r = A;

        return OK;
}




float get_pst_prob(struct pst_node* n, char* string,int target, int pos)
{
        if(pos == 0){
                return n->nuc_probability[target];
        }
        pos = pos -1;
        int c;
        c = nuc_to_internal(string[pos]);
        if(n->next[c]){
                return get_pst_prob(n->next[c], string, target,pos);
        }else{
                return n->nuc_probability[target];
        }
}

float get_ppt_prob(struct pst_node* n, char* string,int target, int pos)
{

        if(string[pos+1] == 0){
                return n->nuc_probability[target];
        }
        int c;
        pos = pos +1;
        c = nuc_to_internal(string[pos]);
        //c = nuc_code[(int)string[pos]];
        if(n->next[c]){
                return get_ppt_prob(n->next[c], string, target,pos);
        }else{

                return n->nuc_probability[target];
        }
}


 int run_build_pst(struct pst** pst, struct tl_seq_buffer* sb)
{
        struct pst* p = NULL;
        //char alphabet[] = "ACGT";
        //char tmp[MAX_PST_LEN+4];
        float sum;
        int i;
        int c;
        uint32_t x;
        RUN(init_pst(&p, sb));

        init_logsum();

        sum = 0.0;
        for(i = 0;i < 4;i++){
                x = i;

                c = p->counts[0][x];
                //fprintf(stdout," %c: %d\n", alphabet[i],c);
                //tmp[0] = alphabet[i];
                //tmp[1] = 0;//alphabet[i];
                //c = count_string(tmp,(const char**)p->suffix_array,p->suffix_len-1,1);
                //fprintf(stdout," %c: %d\n", alphabet[i],c);
                p->pst_root->nuc_probability[i] = c;
                //p->ppt_root->nuc_probability[i] = c;
                sum+= c;
        }
        //exit(0);
        for(i = 0;i < 4;i++){
                p->pst_root->nuc_probability[i] = p->pst_root->nuc_probability[i]/ sum;
                p->ppt_root->nuc_probability[i] = p->pst_root->nuc_probability[i];
                //fprintf(stdout,"%c %f\n", alphabet[i],p->pst_root->nuc_probability[i]);
        }


        p->pst_root = build_pst(p,p->pst_root );
        p->ppt_root = build_ppt(p,p->ppt_root );


        *pst = p;
        return OK;
ERROR:
        return FAIL;
}

struct pst_node* build_pst(struct pst* pst,struct pst_node* n )
{

        char alphabet[] = "ACGT";

        char tmp[MAX_PST_LEN+4];
        int i;
        int j;
        int c;
        int add;

        int len = (int) strlen(n->label);
        float sum = 0.0f;
        float tmp_counts_s[4];

        uint32_t x;

        //fprintf(stderr,"NODE: %s\n", n->label);

        //step 2 test expansion

        //loop though letters at present node
        if(len + 1 < MAX_PST_LEN ){
                for(i = 0; i < 4;i++){

                        if(n->nuc_probability[i] >= pst->p_min){ /// string (suffix  + letter is not rare compared to all other letters  ) // e.g. acgC = 0.2 and not acgC = 0.0000001. hence it makes sense to "extend the suffix - test [A,C,G,T] - acgC

                                //init longer suffix
                                x = i;
                                tmp[0] = alphabet[i];
                                for(j = 1; j < len+1;j++){
                                        x = x << 2 | nuc_to_internal(n->label[j-1]);
                                        tmp[j] = n->label[j-1];
                                }

                                sum = 0.0;
                                for(j = 0; j < 4;j++){
                                        x = x << 2 | j;
                                        c = pst->counts[len+1][x];
                                        tmp[len+1]  = alphabet[j];
                                        tmp[len+2] = 0;
                                        //c = count_string(tmp,(const char**)pst->suffix_array,pst->suffix_len-1,len+2);
                                        tmp_counts_s[j] = c;
                                        sum+= c;

                                        //LOG_MSG("Counts: %d %d %s", c, pst->counts[len+1][x], tmp);
                                        //ASSERT(c ==pst->counts[len+1][x],"Counts: %d %d %s", c, pst->counts[len+1][x], tmp);

                                        x = x >> 2;
                                }
                                //fprintf(stdout,"\n");

                                for(j = 0; j < 4;j++){
                                        tmp_counts_s[j] = tmp_counts_s[j]/sum;
                                }

                                add = 0;
                                for(j = 0; j < 4;j++){
                                        if(tmp_counts_s[j] / n->nuc_probability[j] >= pst->r){
                                                add = 1;
                                                //fprintf(stdout,"Adding because of %d %d\n", i ,add);
                                                //fprintf(stdout, "%f %f\n",tmp_counts_s[j] / n->nuc_probability[j],pst->r);

                                        }

                                        if(tmp_counts_s[j] / n->nuc_probability[j] <= (1.0f/ pst->r)){
                                                add = 1;
                                                //fprintf(stdout,"Adding because of %d %d\n", i ,add);
                                        }



                                }
                                if(add){
                                        n->next[i] = alloc_node(n->next[i] ,tmp,len+1);
                                        sum = 0;
                                        for(j = 0; j < 4;j++){

                                                n->next[i]->nuc_probability[j] = tmp_counts_s[j] *(1.0f  - 4.0f *  pst->gamma_min) + pst->gamma_min;
                                                //fprintf(stdout,"%f ", n->next[i]->nuc_probability[j]);
                                                sum += n->next[i]->nuc_probability[j];
                                        }
                                        //fprintf(stdout," sum: %f\n",sum);

                                        n->next[i] = build_pst(pst,n->next[i]  );
                                        //n->next[i]->in_T = 1;
                                }

                        }
                }
        }
        return n;

}

struct pst_node* build_ppt(struct pst* pst,struct pst_node* n )
{
        char alphabet[] = "ACGT";

        char tmp[MAX_PST_LEN+4];
        int i;
        int j;
        int c;
        int add;
        int len = (int) strlen(n->label);
        float sum = 0.0f;

        float tmp_counts_s[4];

        uint32_t x;
        uint32_t y;

        //fprintf(stderr,"NODE: %s\n", n->label);
        //for(i = 0;i < 5;i++){
        //	fprintf(stderr,"%c+%s\t%f\n",alphabet[i], n->label, n->nuc_probability[i]);
        //}

        //step 2 test expansion

        //loop though letters at present node
        if(len + 1 < MAX_PST_LEN ){
                /// search for all strings and record probabilities S+ACGT...
                /// don't search rare strings...
                /// - super mega simple ...


                for(i = 0; i < 4;i++){
                        if(n->nuc_probability[i] >= pst->p_min){ /// string (suffix  + letter is not rare compared to all other letters  ) // e.g. acgC = 0.2 and not acgC = 0.0000001. hence it makes sense to "extend the suffix - test [A,C,G,T] - acgC





                                //init longer prefix!!!!
                                x = 0;
                                for(j = 0; j < len;j++){
                                        tmp[j+1] = n->label[j];
                                        x = x << 2 |nuc_to_internal(n->label[j]);

                                        //x = x |  (nuc_to_internal(n->label[j]) << (2*(j+1)));
                                }
                                tmp[len+1] = alphabet[i];
                                x = x << 2 | i;

                                //x = x |  (i << (2*(len+1)));

                                sum = 0.0;
                                for(j = 0; j < 4;j++){
                                        y = j << (2 * (len+1));
//x = (x & 0xFFFFFFFC) | j;
                                        //x | j
                                        tmp[0]  = alphabet[j];
                                        tmp[len+2] = 0;
                                        c = pst->counts[len+1][x+y];
                                        //c = pst->counts[len+1][x];
                                        //c = count_string(tmp,(const char**)pst->suffix_array,pst->suffix_len-1,len+2);
                                        tmp_counts_s[j] = c;
                                        //if(c != pst->counts[len+1][x+ y]){
                                        //LOG_MSG("Counts: %d %d %s", c, pst->counts[len+1][x+ y], tmp);
                                        //}
                                        sum+= c;
                                }

                                for(j = 0; j < 4;j++){
                                        tmp_counts_s[j] = tmp_counts_s[j] /sum;
                                }

                                        // here I know that probablility of 'X' is non-neglible AND that there existsa string 'X' - [A,C,G,T,N] which is frequent in the data - hence I add...
                                        //n->next[i] = alloc_node(n->next[i],tmp+1,len+1);
                                add = 0;

                                for(j = 0; j < 4;j++){
                                        if(tmp_counts_s[j] / n->nuc_probability[j] >= pst->r){
                                                add++;
                                        }

                                        if(tmp_counts_s[j] / n->nuc_probability[j] <= 1.0/ pst->r){
                                                add++;
                                        }

                                }

                                if(add){
                                        n->next[i] = alloc_node(n->next[i],tmp+1,len+1);
                                        sum = 0;
                                        for(j = 0; j < 4;j++){

                                                n->next[i]->nuc_probability[j] = tmp_counts_s[j] *(1.0f  - 4.0f *  pst->gamma_min) + pst->gamma_min;
                                                //fprintf(stdout,"%f ", n->next[i]->nuc_probability[j]);
                                                sum += n->next[i]->nuc_probability[j];
                                        }
                                        //fprintf(stdout," sum: %f\n",sum);

                                        n->next[i] = build_ppt(pst,n->next[i]  );
                                        //n->next[i]->in_T = 1;
                                }
                        }
                }
        }
        return n;
}




int init_pst(struct pst** pst, struct tl_seq_buffer* sb)
{
        struct pst* p = NULL;
        int i;
        int j;
        int c;
        int l;
        int len;
        uint32_t x;
        uint32_t mask[MAX_PST_LEN];
        char* seq;
        MMALLOC(p, sizeof(struct pst));
        p->L = MAX_PST_LEN;
        p->gamma_min = 0.02f;

        p->p_min = 0.0001;
        //p->lamba = 0.001;
        p->r = 1.02f;

        p->pst_root = NULL;
        p->ppt_root = NULL;
        p->counts = NULL;
        RUNP(p->pst_root = alloc_node(p->pst_root,"",0));
        RUNP(p->ppt_root = alloc_node(p->ppt_root,"",0));


        DECLARE_TIMER(t);
        START_TIMER(t);


        MMALLOC(p->counts, sizeof(uint32_t*) * MAX_PST_LEN);
        c = 1;

        for(i = 0;i < MAX_PST_LEN;i++){

                p->counts[i] = NULL;
                c *= 4;
                MMALLOC(p->counts[i], sizeof(uint32_t) * c);
                //LOG_MSG("Allocing: %d",c);
                for(j = 0; j < c;j++){
                        p->counts[i][j] = 0;
                }
        }
        mask[0] = 0x3;
        //LOG_MSG("MASK: %x", mask[0]);
        for(i = 1;i < MAX_PST_LEN;i++){
                mask[i] = (mask[i-1] << 2) | 0x3;
                //LOG_MSG("MASK: %x", mask[i]);
        }
        //sb->num_seq = 1;
        int total_len = 0;
        for(i = 0; i < sb->num_seq;i++){
                x = 0;
                seq = sb->sequences[i]->seq;
                len = sb->sequences[i]->len;
                //fprintf(stdout,"%s\n", seq);
                l = 0;
                for(j = 0; j < len;j++){
                        x = (x << 2) | nuc_to_internal(seq[j]);
                        for(c = 0; c <= l;c++){
                                //LOG_MSG("Inserting len:%d seq: %d mask %x", c,x,x &mask[c]);
                                p->counts[c][x & mask[c]]++;
                        }
                        l++;
                        l = MACRO_MIN(l, MAX_PST_LEN-1);
                }
                total_len+= sb->sequences[i]->len;
        }

        STOP_TIMER(t);
        LOG_MSG("counting took %f (%d)", GET_TIMING(t), total_len);
        /*START_TIMER(t);
        p->suffix_len = 0;
        for(i =0; i < sb->num_seq;i++){
                p->suffix_len += sb->sequences[i]->len;
        }
        //p->suffix_len = sb->num_seq;
        p->suffix_array = NULL;
        MMALLOC(p->suffix_array, sizeof(char*) * p->suffix_len);

        c = 0;
        for(i = 0; i < sb->num_seq;i++){
                for(j = 0; j < sb->sequences[i]->len;j++){
                        p->suffix_array[c] = sb->sequences[i]->seq + j;
                        c++;
                }
                //p->mean_length += sb->sequences[i]->len;
        }
        //p->mean_length = p->mean_length / (float) sb->num_seq;

        p->numseq = sb->num_seq;
        ASSERT(c == p->suffix_len,"lengths differ");
        //p->suffix_len = sb->num_seq;
        ///exit(0);
        qsort(p->suffix_array, p->suffix_len, sizeof(char *), qsort_string_cmp);

        STOP_TIMER(t);
        LOG_MSG("counting took %f", GET_TIMING(t));
        //exit(0);*/
        *pst = p;
        return OK;
ERROR:
        return FAIL;
}

void free_pst(struct pst* p)
{
        if(p){
                int i;
                for(i = 0;i < MAX_PST_LEN;i++){
                        MFREE(p->counts[i]);
                }
                MFREE(p->counts);

                free_pst_node(p->pst_root);
                free_pst_node(p->ppt_root);
                //MFREE(p->suffix_array);
                MFREE(p);
        }
}

void free_pst_node(struct pst_node* n)
{
        int i;
        if(n){
                for(i = 0;i < 4;i++){
                        if(n->next[i]){
                                free_pst_node(n->next[i]);
                        }
                }
                MFREE(n->label);
                MFREE(n);
        }
}


struct pst_node* alloc_node(struct pst_node* n,char* string,int len)
{
        int i;
        MMALLOC(n, sizeof(struct pst_node));
        n->label = NULL;
        MMALLOC(n->label, sizeof(char) * (len+1));
        strncpy(n->label, string, len);
        n->label[len] = 0;
        //n->in_T 0;

        //n->last_seen = -1;
        //n->bit_occ = 0;
        for(i =0; i < 4;i++){
                n->next[i] = NULL;
                n->nuc_probability[i] = 0.25f;
        }
        return n;
ERROR:
        return NULL;
}






/*void print_pst(struct pst* pst,struct pst_node* n,int* num)
{
        int i;
        int internal;
        //char alphabet[] = "ACGTN";
        internal = 0;
        for(i = 0;i < 4;i++){
                if(n->next[i]){
                        internal++;
                }
        }
        if(!internal){
                *num = *num +1;
                //fprintf(stderr,"%s\n",n->label);
                //for(i = 0;i < 5;i++){
                //fprintf(stderr,"%c+%s\t%f\n",alphabet[i], n->label, n->nuc_probability[i]);
                //}

        }
        //}


        for(i = 0;i < 4;i++){
                if(n->next[i]){
                        //if(n->next[i]->in_T){
                                //fprintf(stderr,"Going:%d\n",i);
                        print_pst(pst,n->next[i],num);
                                //}
                }
        }
        }*/



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
                break;
        }
        return -1;
}


/* error  correct stuff  */
/*
struct error_correct_seq{
        char** seq;
        int* error;
        float* P;
        int alloc_num_seq;
        int num_seq;
        int len;
        int pos;
};


static int init_error_correct_seq(struct error_correct_seq** err_seq,int n, int len);
static int resize_error_correct_seq(struct error_correct_seq* e);
static int free_error_correct_seq(struct error_correct_seq** err_seq);

int error_correct(struct pst* pst, struct error_correct_seq* e);


static struct pst_node* get_pst_node(struct pst_node* n,char* s,int pos);

struct pst_node* get_pst_node(struct pst_node* n,char* s,int pos)
{
        int c;
        if(pos == 0){
                return n;
        }
        pos = pos -1;
        c = nuc_to_internal(s[pos]);
        if(n->next[c]){
                return get_pst_node(n->next[c], s,pos);
        }else{
                return n;
        }
}


int init_error_correct_seq(struct error_correct_seq** err_seq,int n, int len)
{
        struct error_correct_seq* e = NULL;
        int i;
        MMALLOC(e,sizeof(struct error_correct_seq));
        e->seq = NULL;
        e->error = NULL;
        e->P = NULL;
        e->alloc_num_seq = n;
        e->num_seq = 0;
        e->len = len;
        e->pos = 0;

        RUN(galloc(&e->seq, e->alloc_num_seq ,e->len+1));
        RUN(galloc(&e->P,e->alloc_num_seq ));
        RUN(galloc(&e->error,e->alloc_num_seq));

        *err_seq = e;
        return OK;
ERROR:
        return FAIL;

}

int resize_error_correct_seq(struct error_correct_seq* e)
{
        e->alloc_num_seq = e->alloc_num_seq + e->alloc_num_seq /2;
        RUN(galloc(&e->seq, e->alloc_num_seq ,e->len+1));
        RUN(galloc(&e->P,e->alloc_num_seq ));
        RUN(galloc(&e->error,e->alloc_num_seq));

        return OK;
ERROR:
        return FAIL;
}

int free_error_correct_seq(struct error_correct_seq** err_seq)
{
        struct error_correct_seq* e = NULL;
        e = *err_seq;
        if(e){
                gfree(e->P);
                gfree(e->seq);
                gfree(e->error);
                MFREE(e);
                e = NULL;
        }
        *err_seq = e;
        return OK;
}

int error_correct(struct pst* pst, struct error_correct_seq* e)
{
        struct pst_node* n;
        float p;
        int i,j,c;

        int cur_numseq;

        cur_numseq = e->num_seq;

        for(i = 0; i < cur_numseq;i++){
                fprintf(stdout,"%d\t",e->pos);
                for(j = 0; j < e->len;j++){
                        fprintf(stdout,"%c", e->seq[i][j]);
                }
                fprintf(stdout,"\n");


                n = get_pst_node(pst->pst_root, e->seq[i], e->pos);

                c = nuc_to_internal(e->seq[i][e->pos]);
                p = n->nuc_probability[c];

                if(e->error[i] < 3){
                        for(j = 0; j < 5;j++){
                                if(j!= c){
                                if(n->nuc_probability[j] * 0.05f > p){
                                        //LOG_MSG("Adding pos%d because %f > %f (%d not %d),", e->pos, n->nuc_probability[j],p,j,c);
                                        strncpy(e->seq[e->num_seq], e->seq[i], e->len);
                                        e->seq[e->num_seq][e->len] = 0;
                                        e->seq[e->num_seq][e->pos] = "ACGTN"[j];
                                        e->error[e->num_seq] = e->error[i] +1;
                                        e->P[e->num_seq] = e->P[i] + prob2scaledprob(n->nuc_probability[j] * 0.05f);
                                        e->num_seq++;
                                        if(e->num_seq == e->alloc_num_seq){
                                                resize_error_correct_seq(e);
                                        }
                                }
                                }
                        }
                }
                p = prob2scaledprob(n->nuc_probability[c]);
                e->P[i] += p;

        }
        e->pos++;
        if(e->pos == e->len){
                return OK;

        }
        return error_correct(pst, e);
}
*/