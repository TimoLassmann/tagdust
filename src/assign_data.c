#include "assign_data.h"
#include <string.h>
#include "tldevel.h"


#include "tlseqio.h"

//static int qsort_seq_bits_by_file(const void *a, const void *b);
static int qsort_seq_bits_by_type_file(const void *a, const void *b);
static int qsort_seq_bits_by_file(const void *a, const void *b);
static int qsort_seq_bit_vec(const void *a, const void *b);


static int setup_assign_structure(struct arch_library* al,struct assign_struct* as,int** plan);
static int setup_barcode_files(struct arch_library* al, struct assign_struct* as, char* prefix);
static int setup_output_files(struct assign_struct* as);

static int alloc_seq_bit(struct seq_bit** sb);

static int alloc_demux_struct(struct demux_struct** demux);

static int add_file_name_options(struct rbtree_root** r, struct arch_library* al);
static int add_options(struct rbtree_root** r, char** seq, int num_seq,int max_len,char sep);
//static int add_options(struct rbtree_root** r, char** seq, int num_seq,int max_len);
static int create_demux_tree(struct rbtree_root** r);

int sort_as_by_file_type(struct assign_struct* as)
{
        struct seq_bit_vec* bv = NULL;
        int i;
        for(i = 0; i < as->num_reads ;i++){
                bv = as->bit_vec[i];
                qsort(bv->bits, bv->num_bit,sizeof(struct seq_bit*), qsort_seq_bits_by_type_file);
        }

        return OK;
}


int post_process_assign(struct assign_struct* as)
{
        struct seq_bit_vec* bv = NULL;
        struct demux_struct* tmp_ptr = NULL;

        char* tmp =  NULL;
        //char* barcode;
        //char* umi;
        int i,j,c;//g;
        int len;
        //int umi_len;
        int tmp_len = 256;

        int e_fail;
        int d_fail;
        //char f_code;

        e_fail = 0;
        d_fail = 0;
        MMALLOC(tmp, sizeof(char) * tmp_len);

        for(i = 0; i < as->num_reads ;i++){
                bv = as->bit_vec[i];
                ///* qsort by file identifier  */
                //qsort(bv->bits, bv->num_bit,sizeof(struct seq_bit*), qsort_seq_bits_by_type_file);
                len = 0;
                //umi_len = 0;
                for(j = 0; j < bv->num_bit;j++){

                        if(bv->bits[j]->type == BAR_TYPE){
                                tmp[len] = bv->bits[j]->code;
                                //len += bv->bits[j]->len;// strnlen(bv->bits[j]->p,256);
                                len++;

                        }
                        //if(bv->bits[j]->type == UMI_TYPE){
                        //umi_len += bv->bits[j]->len;
                        //}
                        bv->fail |= bv->bits[j]->fail;
                }
                tmp[len] =0;
                c = 0;
                for(j = 0; j < bv->num_bit;j++){
                        if(bv->bits[j]->type == READ_TYPE){
                                tmp[len] = bv->bits[j]->code;
                                //len += bv->bits[j]->len;// strnlen(bv->bits[j]->p,256);
                                len++;
                                tmp[len] =0;
                                tmp_ptr = as->demux_names->tree_get_data(as->demux_names,tmp);
                                if(!tmp_ptr){
                                        WARNING_MSG("%s not found!", tmp);
                                        exit(0);
                                }else{
                                LOG_MSG("bit:%d in out file id: %d belongs to %s",j, bv->bits[j]->file, tmp_ptr->out_filename);
                                bv->out_file_id[c] = tmp_ptr->id;
                                }
                                len--;
                                c++;
                        }
                }

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
        LOG_MSG("%d %d", e_fail,d_fail);
        MFREE(tmp);
        exit(0);
        qsort(as->bit_vec,as->num_reads ,sizeof(struct seq_bit_vec*) , qsort_seq_bit_vec);
        return OK;
ERROR:
        return FAIL;
}

int init_assign_structure(struct assign_struct** assign,struct arch_library* al,char* prefix, int total)
{
        struct assign_struct* as = NULL;
        //ASSERT(num_files >= 1,"no infiles");

        int* plan;
        int i,j;
        MMALLOC(as, sizeof(struct assign_struct));
        as->num_files = al->num_file;
        as->demux_names = NULL;
        as->bit_vec = NULL;
        as->num_bits = 0;
        as->max_bar_len = 0;
        as->max_seq_len = 0;
        as->out_reads = 0;
        as->loc_out_reads = NULL;
        as->file_index = NULL;
        RUN(setup_assign_structure(al,as,&plan));
        RUN(setup_barcode_files(al,as,prefix));
        //exit(0);
        //RUN(setup_output_files(as));
        //LOG_MSG("%d", as->out_reads);
        //exit(0);
        as->alloc_total = total;
        as->num_reads = 0;

        as->bit_vec = NULL;
        MMALLOC(as->bit_vec, sizeof(struct seq_bit_vec*)* as->alloc_total );
        //MMALLOC(as->active_bits,sizeof(uint8_t) * as->total));
        for(i = 0; i < as->alloc_total;i++){

                as->bit_vec[i] = NULL;
                MMALLOC(as->bit_vec[i], sizeof(struct seq_bit_vec));
                as->bit_vec[i]->out_file_id = NULL;
                MMALLOC(as->bit_vec[i]->out_file_id, sizeof(int) * as->out_reads);
                for(j = 0; j < as->out_reads;j++){
                        as->bit_vec[i]->out_file_id[j] = -1;
                }
                as->bit_vec[i]->sample_group = -1;
                as->bit_vec[i]->num_bit = as->num_bits;
                as->bit_vec[i]->bits = NULL;
                //as->bit_vec[i]->Q = NULL;
                //as->bit_vec[i]->bc = NULL;
                as->bit_vec[i]->umi = NULL;
                //MMALLOC(as->bit_vec[i]->Q,sizeof(float) * as->num_files);
                as->bit_vec[i]->fail = 0;
                MMALLOC(as->bit_vec[i]->bits , sizeof(struct seq_bit) * as->num_bits);
                for(j = 0; j < as->num_bits;j++){
                        as->bit_vec[i]->bits[j] = NULL;

                        RUN(alloc_seq_bit(&as->bit_vec[i]->bits[j]));
                        //MMALLOC(as->bit_vec[i]->bits[j], sizeof(struct seq_bit));
                        as->bit_vec[i]->bits[j]->file = plan[j];
                }
        }
        MFREE(plan);
        *assign = as;
        return OK;
ERROR:
        free_assign_structure(as);
        return FAIL;
}

int alloc_seq_bit(struct seq_bit** sb)
{
        struct seq_bit* b = NULL;

        MMALLOC(b, sizeof(struct seq_bit));
        b->fail = 0;
        b->file = 0;
        b->type = 0;
        b->p = NULL;
        b->q = 0;
        b->len = 0;
        *sb = b;
        return OK;
ERROR:
        return FAIL;
}

int reset_assign_structute(struct assign_struct* as)
{
        int i;
        int j;
        for(i = 0; i < as->alloc_total;i++){
                //qsort(  as->bit_vec[i]->bits,  as->bit_vec[i]->num_bit,sizeof(struct seq_bit*), qsort_seq_bits_by_file);
                as->bit_vec[i]->fail = 0;
                as->bit_vec[i]->sample_group = -1;
                if(as->bit_vec[i]->umi){
                        MFREE(as->bit_vec[i]->umi);
                }
                for(j = 0; j < as->out_reads;j++){
                        as->bit_vec[i]->out_file_id[j] = -1;
                }
                //as->bits[i]->num_bit = 0;
        }
        return OK;
}

int setup_assign_structure(struct arch_library* al,struct assign_struct* as,int** plan)
{
        struct read_structure* read_structure = NULL;
        int* p = NULL;
        int i,j;
        //int f;
        uint8_t c;
        int p_size;
        int p_index;
        ASSERT(al != NULL,"No archlib");

        ASSERT(as != NULL,"No assign struct ");

        as->num_bits = 0;
        as->out_reads = 0;
        p_size = 256;
        p_index = 0;
        MMALLOC(p, sizeof(int) * p_size);
        MMALLOC(as->loc_out_reads, sizeof(int) * p_size);
        MMALLOC(as->file_index , sizeof(int) * p_size);

        for(i = 0; i < al->num_file;i++){
                as->file_index[i] = as->num_bits;
                read_structure = al->read_structure[al->arch_to_read_assignment[i]];
                //fprintf(stdout,"Read %d:\n",i);
                for(j = 0; j < read_structure->num_segments;j++){
                        //c = read_structure->type[j];

                        c = read_structure->seg_spec[j]->extract;
                        //print_segment_spec(read_structure->seg_spec[j]);

                        switch (c) {
                        case ARCH_ETYPE_EXTRACT:
                                as->loc_out_reads[as->out_reads] = as->num_bits;
                                p[p_index] = as->out_reads; /* might be wrong if 2 reads come from the same file; set to out_reads(?)  */
                                p_index++;
                                if(p_index == p_size){
                                        p_size = p_size + p_size /2;
                                        MREALLOC(p, sizeof(int)* p_size);
                                        MREALLOC(as->loc_out_reads , sizeof(int) * p_size);
                                }
                                as->out_reads++;

                                break;

                        default:
                                p[p_index] = i;
                                p_index++;
                                if(p_index == p_size){
                                        p_size = p_size + p_size /2;
                                        MREALLOC(p, sizeof(int)* p_size);
                                        MREALLOC(as->loc_out_reads , sizeof(int) * p_size);
                                }
                                break;

                        }
                        as->num_bits++;
                }


                //fprintf(stdout,"\n");
        }
        //exit(0);
        *plan = p;
        return OK;
ERROR:
        return FAIL;
}


void free_assign_structure(struct assign_struct* as)
{
        if(as){
                int i,j;
                if(as->bit_vec){
                        //MMALLOC(as->active_bits,sizeof(uint8_t) * as->total));
                        for(i = 0; i < as->alloc_total;i++){
                                for(j = 0; j < as->num_bits;j++){
                                        MFREE(as->bit_vec[i]->bits[j]);
                                }
                                MFREE(as->bit_vec[i]->out_file_id);
//MFREE(as->bits[i]->Q);
                                MFREE(as->bit_vec[i]->bits);
                                MFREE(as->bit_vec[i]);
                        }
                        MFREE(as->bit_vec);

                }
                if(as->demux_names){
                        as->demux_names->free_tree(as->demux_names);
                }
                if(as->file_names){
                        as->file_names->free_tree(as->file_names);
                }
                MFREE(as->file_index);
                MFREE(as->loc_out_reads);
                MFREE(as);
        }
}


int qsort_seq_bits_by_file(const void *a, const void *b)
{
        const struct seq_bit **elem1 = (const struct seq_bit**) a;
        const struct seq_bit **elem2 = (const struct seq_bit**) b;
        if ( (*elem1)->file > (*elem2)->file){
                return 1;
        }else if ((*elem1)->file < (*elem2)->file){
                return -1;
        }else{
                return 0;
        }
}

int qsort_seq_bits_by_type_file(const void *a, const void *b)
{
        const struct seq_bit **elem1 = (const struct seq_bit**) a;
        const struct seq_bit **elem2 = (const struct seq_bit**) b;
        if((*elem1)->type > (*elem2)->type){
                return 1;
        }else if((*elem1)->type < (*elem2)->type){
                return -1;
        }else{
                if ( (*elem1)->file > (*elem2)->file){
                        return 1;
                }else if ((*elem1)->file < (*elem2)->file){
                        return -1;
                }else{
                        return 0;
                }
        }
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
                if ( (*elem1)->sample_group > (*elem2)->sample_group){
                        return 1;
                }else if ((*elem1)->sample_group < (*elem2)->sample_group){
                        return -1;
                }else{
                        return 0;
                }
        }
}




static void* get_name(void* ptr);
static void* get_key(void* ptr);
static long int compare_name(void* keyA, void* keyB);
static void print_demux_struct(void* ptr,FILE* out_ptr);
static int resolve_default(void* ptr_a,void* ptr_b);
static void free_demux_struct(void* ptr);

int setup_output_files(struct assign_struct* as)
{
        struct rbtree_root* root = NULL;
        void*  (*fp_get)(void* ptr) = NULL;
        long int (*fp_cmp)(void* keyA, void* keyB)= NULL;
        int (*fp_cmp_same)(void* ptr_a,void* ptr_b);
        void (*fp_print)(void* ptr,FILE* out_ptr) = NULL;

        void (*fp_free)(void* ptr) = NULL;


        fp_get = &get_name;
        fp_cmp = &compare_name;
        fp_print = &print_demux_struct;
        fp_cmp_same = &resolve_default;
        fp_free = &free_demux_struct;
        RUNP(root = init_tree(fp_get,fp_cmp,fp_cmp_same,fp_print,fp_free));
        as->file_names = root;
        return OK;
ERROR:
        return FAIL;
}

int setup_barcode_files(struct arch_library* al, struct assign_struct* as, char* prefix)
{
        struct read_structure* read_structure = NULL;
        int i,j,f;
        int p_len;
        //uint8_t c;

        ASSERT(al != NULL,"No archlib");
        struct rbtree_root* root = NULL;
        //struct rbtree_root* root_new = NULL;
        struct demux_struct* tmp_ptr = NULL;




        RUN(create_demux_tree(&root));
        /* add prefix  */
        p_len = strlen(prefix);
        RUN(alloc_demux_struct(&tmp_ptr));


        MMALLOC(tmp_ptr->out_filename,sizeof(char) *(p_len+1));
        snprintf(tmp_ptr->out_filename,p_len+1,"%s",prefix);
        tmp_ptr->out_filename[p_len]=0;


        MMALLOC(tmp_ptr->key, sizeof(char)*1);

        tmp_ptr->key[0] = 0;
        tmp_ptr->count = 1;
        RUN(root->tree_insert(root,tmp_ptr));
        tmp_ptr = NULL;



        for(i = 0; i < al->num_file;i++){
                read_structure = al->read_structure[al->arch_to_read_assignment[i]];
                //fprintf(stdout,"Read %d: ",i);
                for(j = 0; j < read_structure->num_segments;j++){

                        //c = read_structure->type[j];
                        if(read_structure->seg_spec[j]->extract == ARCH_ETYPE_SPLIT){
                                RUN(add_options(&root,
                                                read_structure->seg_spec[j]->seq,
                                                read_structure->seg_spec[j]->num_seq,
                                                read_structure->seg_spec[j]->max_len,
                                                '_'
                                            ));

                        }
                }
                //fprintf(stdout,"\n");
        }
        RUN(add_file_name_options(&root,al));

        RUN(root->flatten_tree(root));

                //root_new = init_tree(fp_get,fp_cmp,fp_cmp_same,fp_print,fp_free);
        for(f = 0;f < root->num_entries;f++){
                tmp_ptr = root->data_nodes[f];
                tmp_ptr->id = f;
                fprintf(stdout,"%s %s %d %d\n",tmp_ptr->out_filename, tmp_ptr->key ,tmp_ptr->count, tmp_ptr->id);
        }

        //root->free_tree(root);
        as->demux_names = root;

        /*char* query = NULL;

        MMALLOC(query, sizeof(char) * 5);
        query[0] = '#';

        query[1] = '\"';
        //query[3] = '\"';
        query[2] = 0;


        RUNP(tmp_ptr = as->demux_names->tree_get_data(as->demux_names,query));
        LOG_MSG("Found:");
        fprintf(stdout,"%s %s %d %d\n",tmp_ptr->out_filename, tmp_ptr->key ,tmp_ptr->count, tmp_ptr->id);*/

        return OK;
ERROR:
        return FAIL;
}

/* This monstrosity:
1) looks for named segments that are extracted
2) makes sure there aren't any duplicate names - if there are append a number;
3) adds the .fastq.gz suffix
4) adds all the possibilities to the prefix output files (which is either the -o option ot -o + barcode sequences))
- ufff */

int add_file_name_options(struct rbtree_root** r, struct arch_library* al)
{
        struct rbtree_root* root = NULL;
        struct rbtree_root* root2 = NULL;
        struct demux_struct* tmp_ptr = NULL;
        struct demux_struct* new_ptr = NULL;
        struct read_structure* read_structure = NULL;
        char* suffix[] = {".fastq.gz", NULL};
        char** additional = NULL;

        int i,j,c;

        int len;
        RUN(create_demux_tree(&root));

        for(i = 0; i < al->num_file;i++){
                read_structure = al->read_structure[al->arch_to_read_assignment[i]];
                //fprintf(stdout,"Read %d: ",i);

                for(j = 0; j < read_structure->num_segments;j++){
                        if(read_structure->seg_spec[j]->extract == ARCH_ETYPE_EXTRACT){
                                RUN(alloc_demux_struct(&tmp_ptr));
                                if(read_structure->seg_spec[j]->name){
                                        len = strlen(read_structure->seg_spec[j]->name);
                                        MMALLOC(tmp_ptr->out_filename,sizeof(char) * (len+1));
                                        snprintf(tmp_ptr->out_filename, len+1,"%s",read_structure->seg_spec[j]->name);
                                        tmp_ptr->out_filename[len] = 0;
                                }else{
                                        WARNING_MSG("This should never happen -> archlib token code");
                                        MMALLOC(tmp_ptr->out_filename,sizeof(char) * 2);
                                        tmp_ptr->out_filename[0] = 'R';
                                        tmp_ptr->out_filename[1] = 0;

                                }

                                MMALLOC(tmp_ptr->key, sizeof(char)*2);
                                tmp_ptr->key[0] = (char)(j+33);
                                tmp_ptr->key[1] = 0;
                                tmp_ptr->count = 1;
                                RUN(root->tree_insert(root,tmp_ptr));
                                tmp_ptr = NULL;
                                //fprintf(stdout,"%s NAME\n",read_structure->seg_spec[j]->name);
                        }
                }


        }
        RUN(create_demux_tree(&root2));
        c = 0;
        //LOG_MSG("%d",root->num_entries);
        RUN(root->flatten_tree(root));
        for(i = 0;i < root->num_entries;i++){
                tmp_ptr = root->data_nodes[i];
                tmp_ptr->id = i;
                ///fprintf(stdout,"%s %s %d\n",tmp_ptr->out_filename, tmp_ptr->key ,tmp_ptr->count);
                if(tmp_ptr->count > 1){
                        for(j = 0; j < tmp_ptr->count;j++){
                                RUN(alloc_demux_struct(&new_ptr));

                                len = strlen(tmp_ptr->out_filename);
                                MMALLOC(new_ptr->out_filename,sizeof(char) * (len+10));
                                snprintf(new_ptr->out_filename, len+10,"%s%d",tmp_ptr->out_filename,j+1);
                                MMALLOC(new_ptr->key, sizeof(char)*2);
                                new_ptr->key[0] = (char)(c+33);
                                new_ptr->key[1] = 0;


                                RUN(root2->tree_insert(root2,new_ptr));
                                new_ptr = NULL;
                                c++;

                        }
                }else{
                        RUN(alloc_demux_struct(&new_ptr));

                        len = strlen(tmp_ptr->out_filename);
                        MMALLOC(new_ptr->out_filename,sizeof(char) * (len+1));
                        snprintf(new_ptr->out_filename, len+1,"%s",tmp_ptr->out_filename);
                        MMALLOC(new_ptr->key, sizeof(char)*2);
                        new_ptr->key[0] = (char)(c+33);
                        new_ptr->key[1] = 0;
                        RUN(root2->tree_insert(root2,new_ptr));
                        new_ptr = NULL;
                        c++;

                }
        }
        root->free_tree(root);
        root= NULL;


        len = 0;

        RUN(add_options(&root2, suffix, 1,  9,0));

        RUN(root2->flatten_tree(root2));
        MMALLOC(additional,sizeof(char*) * root2->num_entries);

        for(i = 0;i < root2->num_entries;i++){
                tmp_ptr = root2->data_nodes[i];
                tmp_ptr->id = i;
                //fprintf(stdout,"%s %s %d\n",tmp_ptr->out_filename, tmp_ptr->key ,tmp_ptr->count);
                additional[i] = tmp_ptr->out_filename;
                c = strlen(tmp_ptr->out_filename);
                if(c > len){
                        len = c;
                }
        }

        root= *r;
        RUN(add_options(&root, additional, root2->num_entries,  len,'_'));
        root2->free_tree(root2);

        root2 = NULL;

        MFREE(additional);

        *r = root;
        return OK;
ERROR:
        return FAIL;
}


int add_options(struct rbtree_root** r, char** seq, int num_seq,int max_len,char sep)
{
        struct rbtree_root* root = NULL;
        struct rbtree_root* root_new = NULL;
        struct demux_struct* tmp_ptr = NULL;
        struct demux_struct* new_ptr = NULL;

        int f;
        int i;
        int len;
        if(*r){
                root = *r;

                RUN(root->flatten_tree(root));

                RUN(create_demux_tree(&root_new));

                for(f = 0;f < root->num_entries;f++){
                        tmp_ptr = root->data_nodes[f];
                        for(i = 0; i < num_seq;i++){
                                //for(g = 0;g < read_structure->seg_spec[j]->num_seq;g++){
                                len = strlen(tmp_ptr->out_filename);
                                len++;
                                len += max_len;//  segment_length[j];
                                //len += strnlen(read_structure->sequence_matrix[j][g], 256);
                                len++;
                                RUN(alloc_demux_struct(&new_ptr));

                                MMALLOC(new_ptr->out_filename,sizeof(char) * len);
                                if(sep){
                                        snprintf(new_ptr->out_filename, len, "%s%c%s",tmp_ptr->out_filename,sep,seq[i]);
                                        new_ptr->out_filename[len-1] = 0;
                                }else{
                                        snprintf(new_ptr->out_filename, len, "%s%s",tmp_ptr->out_filename,seq[i]);
                                        new_ptr->out_filename[len-2] = 0;
                                }
                                len = strlen(tmp_ptr->key);
                                len+=2;

                                MMALLOC(new_ptr->key, sizeof(char) * len);

                                snprintf(new_ptr->key, len, "%s%c", tmp_ptr->key, (char) (i+33));
                                new_ptr->key[len-1] = 0;

                                RUN(root_new->tree_insert(root_new,new_ptr));
                                new_ptr = NULL;
                        }
//root->data_nodes[f]

                }
                root->free_tree(root);
                root = root_new;
                root_new = NULL;

        }else{
                ERROR_MSG("No red black tree!");
        }

        *r = root;
        return OK;
ERROR:
        return FAIL;
}




int create_demux_tree(struct rbtree_root** r)
{
        struct rbtree_root* root = NULL;

        void*  (*fp_get)(void* ptr) = NULL;
        long int (*fp_cmp)(void* keyA, void* keyB)= NULL;
        int (*fp_cmp_same)(void* ptr_a,void* ptr_b);
        void (*fp_print)(void* ptr,FILE* out_ptr) = NULL;

        void (*fp_free)(void* ptr) = NULL;



        if(*r){
                ERROR_MSG("Tree already exists!");
        }else{
                fp_get = &get_key;
                fp_cmp = &compare_name;
                fp_print = &print_demux_struct;
                fp_cmp_same = &resolve_default;
                fp_free = &free_demux_struct;
                RUNP(root = init_tree(fp_get,fp_cmp,fp_cmp_same,fp_print,fp_free));
        }
        *r = root;
        return OK;
ERROR:
        return FAIL;
}


int alloc_demux_struct(struct demux_struct** demux)
{
        struct demux_struct* d = NULL;

        MMALLOC(d, sizeof(struct demux_struct));

        d->f_hand = NULL;
        d->name = NULL;
        d->out_filename = NULL;
        d->key = NULL;
        d->id = 0;
        d->count = 0;

        *demux = d;
        return OK;
ERROR:
        return FAIL;
}

void* get_name(void* ptr)
{
        struct demux_struct* tmp = (struct demux_struct*)  ptr;
        return tmp->name;
}

void* get_key(void* ptr)
{
        struct demux_struct* tmp = (struct demux_struct*)  ptr;
        return tmp->key;
}

long int compare_name(void* keyA, void* keyB)
{
        return strcmp(keyA,keyB);
}


void print_demux_struct(void* ptr,FILE* out_ptr)
{
        struct demux_struct* tmp = (struct demux_struct*)  ptr;
        fprintf(out_ptr,"%s\t%d\t%d\n",tmp->name, tmp->id,tmp->count);
}

int resolve_default(void* ptr_a,void* ptr_b)
{
        struct demux_struct* a = (struct demux_struct*)  ptr_a;
        struct demux_struct* b = (struct demux_struct*)  ptr_b;
        a->count++;

        free_demux_struct(b);
        return 0;
}

void free_demux_struct(void* ptr)
{
        struct demux_struct* tmp = (struct demux_struct*)  ptr;
        if(tmp){
                if(tmp->f_hand){
                        close_seq_file(&tmp->f_hand);
                }
                if(tmp->name){
                        MFREE(tmp->name);
                }
                MFREE(tmp);
        }
}
