#include "assign_data.h"
#include <string.h>

#include "kstring.h"
#include "tldevel.h"

static int setup_assign_structure(struct arch_library* al,struct assign_struct* as);
//static int setup_barcode_files(struct arch_library* al, struct assign_struct* as, char* prefix);
static int setup_barcode_files(struct arch_library* al, struct assign_struct* as, char* prefix,int bam);

static int alloc_seq_bit(struct seq_bit** sb);
static int free_seq_bit(struct seq_bit** sb);

static int alloc_bit_ann(struct bit_annotation** ba);
static int free_bit_annotation(struct bit_annotation** sb);

static int alloc_demux_struct(struct demux_struct** demux);

static int add_file_name_options(struct rbtree_root** r, struct arch_library* al,int bam);
static int add_options(struct rbtree_root** r, char** seq, int num_seq,int max_len,char sep);
//static int add_options(struct rbtree_root** r, char** seq, int num_seq,int max_len);
static int create_demux_tree(struct rbtree_root** r);

/*int sort_as_by_file_type(struct assign_struct* as)
{
        struct seq_bit_vec* bv = NULL;
        int i;
        for(i = 0; i < as->num_reads ;i++){
                bv = as->bit_vec[i];
                //qsort(bv->bits, bv->num_bit,sizeof(struct seq_bit*), qsort_seq_bits_by_type_file);
        }
        return OK;
        }*/

int init_assign_structure(struct assign_struct** assign,struct arch_library* al,char* prefix, int total,int bam)
{
        struct assign_struct* as = NULL;
        struct seq_bit_vec* bv = NULL;
        //ASSERT(num_files >= 1,"no infiles");

        //int* plan;
        int i,j;
        MMALLOC(as, sizeof(struct assign_struct));
        as->num_files = al->num_file;
        as->demux_names = NULL;
        as->bit_vec = NULL;
        as->num_bits = 0;
        as->max_bar_len = 0;
        as->max_seq_len = 0;
        as->out_reads = 0;
        as->n_exact_hash = 0;
        //as->exact = NULL;
        as->loc_out_reads = NULL;
        as->file_index = NULL;
        as->subm = NULL;
        as->bit_ann = NULL;

        RUN(calc_score_matrix(&as->subm, 0.02f, 0.1f));
        RUN(setup_assign_structure(al,as));
        RUN(setup_barcode_files(al,as,prefix,bam));
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
                bv = as->bit_vec[i];
                bv->out_file_id = NULL;
                MMALLOC(bv->out_file_id, sizeof(int) * as->out_reads);
                for(j = 0; j < as->out_reads;j++){
                        bv->out_file_id[j] = -1;
                }
                bv->sample_group = -1;
                bv->num_bit = as->num_bits;
                bv->bits = NULL;
                //bv->Q = NULL;
                //bv->bc = NULL;
                bv->append.l = 0;
                bv->append.m = 0;
                bv->append.s = NULL;

                //bv->a_len = 0;

                //MMALLOC(bv->append, sizeof(char) * as->append_len);



                //MMALLOC(bv->Q,sizeof(float) * as->num_files);
                bv->fail = 0;
                MMALLOC(bv->bits , sizeof(struct seq_bit*) * as->num_bits);
                for(j = 0; j < as->num_bits;j++){
                        bv->bits[j] = NULL;

                        RUN(alloc_seq_bit(&bv->bits[j]));
                        //MMALLOC(bv->bits[j], sizeof(struct seq_bit));
                        //bv->bits[j]->file = plan[j];
                }
        }
        //MFREE(plan);
        *assign = as;
        return OK;
ERROR:
        free_assign_structure(as);
        return FAIL;
}

int alloc_bit_ann(struct bit_annotation** ba)
{
        struct bit_annotation* b = NULL;


        MMALLOC(b, sizeof(struct bit_annotation));
        b->name.l = 0;
        b->name.m = 0;
        b->name.s = NULL;

        b->c_name.l = 0;
        b->c_name.m = 0;
        b->c_name.s = NULL;

        b->q_name.l = 0;
        b->q_name.m = 0;
        b->q_name.s = NULL;

        *ba = b;
        return OK;
ERROR:
        return FAIL;
}

int free_bit_annotation(struct bit_annotation** sb)
{
        struct bit_annotation* b = NULL;

        b = *sb;
        if(b){
                if(b->name.m){
                        MFREE(b->name.s);
                }
                if(b->c_name.m){
                        MFREE(b->c_name.s);
                }
                if(b->q_name.m){
                        MFREE(b->q_name.s);
                }

                MFREE(b);
                b = NULL;
        }
        *sb = b;
        return OK;
//ERROR:
        //return FAIL;
}


int alloc_seq_bit(struct seq_bit** sb)
{
        struct seq_bit* b = NULL;

        MMALLOC(b, sizeof(struct seq_bit));
        b->p.l = 0;
        b->p.m = 0;
        b->p.s = NULL;

        b->p_corr.l = 0;
        b->p_corr.m = 0;
        b->p_corr.s = NULL;

        b->q.l = 0;
        b->q.m = 0;
        b->q.s = NULL;


        b->fail = 0;
        b->type = 0;
        b->correct_index = 0;
        //b->p = NULL;
        //b->q = 0;

        *sb = b;
        return OK;
ERROR:
        return FAIL;
}

int free_seq_bit(struct seq_bit** sb)
{
        struct seq_bit* b = NULL;

        b = *sb;
        if(b){
                if(b->p_corr.m){
                        MFREE(b->p_corr.s);
                }
                if(b->p.m){
                        MFREE(b->p.s);
                }
                if(b->q.m){
                        MFREE(b->q.s);
                }
                MFREE(b);
                b = NULL;
        }
        *sb = b;
        return OK;
//ERROR:
        //return FAIL;
}

int reset_assign_structute(struct assign_struct* as)
{
        int i;
        int j;
        struct seq_bit* s;
        for(i = 0; i < as->alloc_total;i++){
                //qsort(  as->bit_vec[i]->bits,  as->bit_vec[i]->num_bit,sizeof(struct seq_bit*), qsort_seq_bits_by_file);
                //ASSERT(as->bit_vec[i]->append.l == 0 ," Should be zero");
                as->bit_vec[i]->append.l = 0;
                as->bit_vec[i]->fail = 0;
                as->bit_vec[i]->sample_group = -1;
                //if(as->bit_vec[i]->append){
                //MFREE(as->bit_vec[i]->append);
                //a}

                for(j = 0; j < as->out_reads;j++){
                        as->bit_vec[i]->out_file_id[j] = -1;
                }
                for(j = 0; j < as->num_bits;j++){

                        s = as->bit_vec[i]->bits[j];
                        s->p.l = 0;
                        s->q.l = 0;
                        s->p_corr.l = 0;
                        //ASSERT(s->p.l == 0,"p is not null as %d  %d %d %s fail? %d",i, s->file,s->type, s->p.s, s->fail);
                        //ASSERT(s->q.l == 0,"p is not null");
                }
                //as->bits[i]->num_bit = 0;
        }
        return OK;
//ERROR:
        //      return FAIL;
}

int setup_assign_structure(struct arch_library* al,struct assign_struct* as)
{
        struct read_structure* read_structure = NULL;
        struct segment_specs* s = NULL;
        //int* p = NULL;
        int i,j;
        //int f;
        uint8_t c;
        int p_size;
        //int p_index;
        ASSERT(al != NULL,"No archlib");

        ASSERT(as != NULL,"No assign struct ");

        as->num_bits = 0;
        as->out_reads = 0;
        p_size = 0;
        for(i = 0; i < al->num_file;i++){
                read_structure = al->read_structure[al->arch_to_read_assignment[i]];
                //fprintf(stdout,"Read %d:\n",i);
                p_size += read_structure->num_segments;
        }
        //p_index = 0;
        //MMALLOC(p, sizeof(int) * p_size);
        MMALLOC(as->loc_out_reads, sizeof(int) * p_size);
        MMALLOC(as->file_index , sizeof(int) * p_size);

        MMALLOC(as->bit_ann, sizeof(struct bit_annotation*) * p_size);
        for(i = 0; i < p_size;i++){
                //as->bit_ann[i]
                as->bit_ann[i] = NULL;

                alloc_bit_ann(&as->bit_ann[i]);
                //MMALLOC(as->bit_ann[i], sizeof(struct bit_annotation));
                //as->bit_ann[i]->name;

        }
        //struct bit_annotation* a;
        //a->c_name
        for(i = 0; i < al->num_file;i++){
                as->file_index[i] = as->num_bits;
                read_structure = al->read_structure[al->arch_to_read_assignment[i]];
                //fprintf(stdout,"Read %d:\n",i);

                for(j = 0; j < read_structure->num_segments;j++){
                        //c = read_structure->type[j];
                        s = read_structure->seg_spec[j];
                        c = s->extract;

                        if(s->name){
                                kputs(s->name, &as->bit_ann[as->num_bits]->name);
                        }
                        if(s->name){
                                kputs(s->correct_name, &as->bit_ann[as->num_bits]->c_name);
                        }
                        if(s->name){
                                kputs(s->qual_name, &as->bit_ann[as->num_bits]->q_name);
                        }
                        if(s->bar_hash){
                                as->bit_ann[as->num_bits]->bar_hash = s->bar_hash;
                        }
                        //as->bit_ann[0]->name = 0;

                        //print_segment_spec(read_structure->seg_spec[j]);

                        /*p[p_index] = i;
                          p_index++;*/
                        /*if(p_index == p_size){
                                p_size = p_size + p_size /2;
                                MREALLOC(p, sizeof(int)* p_size);
                                MREALLOC(as->loc_out_reads , sizeof(int) * p_size);
                                }*/
                        if(c == ARCH_ETYPE_EXTRACT){
                                as->loc_out_reads[as->out_reads] = as->num_bits;
                                //p[p_index] = as->out_reads;
                                as->out_reads++;
                        }


                        as->num_bits++;
                }


                //fprintf(stdout,"\n");
        }
        /*as->n_exact_hash = 0;


        for(i = 0; i < al->num_file;i++){
                as->file_index[i] = as->num_bits;
                read_structure = al->read_structure[al->arch_to_read_assignment[i]];
                //fprintf(stdout,"Read %d:\n",i);
                for(j = 0; j < read_structure->num_segments;j++){
                        c = read_structure->seg_spec[j]->extract;
                        if(c == ARCH_ETYPE_APPEND_CORRECT){
                                as->n_exact_hash++;
                        }
                }
        }
        if(as->n_exact_hash){

                MMALLOC(as->exact, sizeof(khash_t(exact)*) * as->n_exact_hash);
                for(i = 0; i < as->n_exact_hash;i++){
                        as->exact[i] = NULL;
                }
                as->n_exact_hash = 0;
                for(i = 0; i < al->num_file;i++){
                        as->file_index[i] = as->num_bits;
                        read_structure = al->read_structure[al->arch_to_read_assignment[i]];
                        //fprintf(stdout,"Read %d:\n",i);
                        for(j = 0; j < read_structure->num_segments;j++){
                                c = read_structure->seg_spec[j]->extract;
                                if(c == ARCH_ETYPE_APPEND_CORRECT){
                                        LOG_MSG("reading in barcodes from: %s",read_structure->seg_spec[j]->seq[0]);
                                        RUN(read_white_list(&sb, read_structure->seg_spec[j]->seq[0]));

                                        LOG_MSG("read: %d barcodes" , sb->num_seq);
                                        //ASSERT(
                                        RUN(fill_exact_hash(& as->exact[as->n_exact_hash] , sb));
                                        as->n_exact_hash++;

                            }
                        }
                }
        }

        free_tl_seq_buffer(sb);*/
        //exit(0);
        //*plan = p;
        return OK;
ERROR:
        return FAIL;
}


void free_assign_structure(struct assign_struct* as)
{
        if(as){
                int i,j;
                if(as->subm){
                        MFREE(as->subm);
                }
                if(as->bit_vec){
                        //MMALLOC(as->active_bits,sizeof(uint8_t) * as->total));
                        for(i = 0; i < as->alloc_total;i++){
                                for(j = 0; j < as->num_bits;j++){

                                        free_seq_bit(&as->bit_vec[i]->bits[j]);
                                }
                                MFREE(as->bit_vec[i]->out_file_id);
//MFREE(as->bits[i]->Q);
                                if(as->bit_vec[i]->append.m){
                                        MFREE(as->bit_vec[i]->append.s);
                                }
                                MFREE(as->bit_vec[i]->bits);
                                MFREE(as->bit_vec[i]);
                        }
                        MFREE(as->bit_vec);

                }
                if(as->demux_names){
                        as->demux_names->free_tree(as->demux_names);
                }
                for(j = 0; j < as->num_bits;j++){
                        free_bit_annotation(&as->bit_ann[j]);
                }
                MFREE(as->bit_ann);
                //for(i = 0; i < as->n_exact_hash;i++){

                //}
                MFREE(as->file_index);
                MFREE(as->loc_out_reads);
                MFREE(as);
        }
}








static void* get_key(void* ptr);
static long int compare_name(void* keyA, void* keyB);
static void print_demux_struct(void* ptr,FILE* out_ptr);
static int resolve_default(void* ptr_a,void* ptr_b);
static void free_demux_struct(void* ptr);


int setup_barcode_files(struct arch_library* al, struct assign_struct* as, char* prefix,int bam)
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
        RUN(add_file_name_options(&root,al,bam));

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

int add_file_name_options(struct rbtree_root** r, struct arch_library* al,int bam)
{
        struct rbtree_root* root = NULL;
        struct rbtree_root* root2 = NULL;
        struct demux_struct* tmp_ptr = NULL;
        struct demux_struct* new_ptr = NULL;
        struct read_structure* read_structure = NULL;
        char* suffix[] = {".fastq.gz", NULL};

        char* suffix_bam[] = {".bam", NULL};
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

        if(bam){
                RUN(add_options(&root2, suffix_bam, 1,4,0));
        }else{

                RUN(add_options(&root2, suffix, 1,  9,0));
        }

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
        d->out_filename = NULL;
        d->key = NULL;
        d->id = 0;
        d->count = 0;

        *demux = d;
        return OK;
ERROR:
        return FAIL;
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
        fprintf(out_ptr,"%s\t%d\t%d\n",tmp->key, tmp->id,tmp->count);
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
                if(tmp->out_filename){
                        MFREE(tmp->out_filename);
                }
                if(tmp->key){
                        MFREE(tmp->key);
                }
                MFREE(tmp);
        }
}
