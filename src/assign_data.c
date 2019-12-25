#include "assign_data.h"
#include <string.h>

#include "tldevel.h"





//static int qsort_seq_bits_by_file(const void *a, const void *b);
static int qsort_seq_bits_by_type_file(const void *a, const void *b);
static int qsort_seq_bits_by_file(const void *a, const void *b);
static int qsort_seq_bit_vec(const void *a, const void *b);

static int setup_assign_structure(struct arch_library* al,struct assign_struct* as,int** plan);
static int setup_barcode_files(struct arch_library* al, struct assign_struct* as);
static int setup_output_files(struct assign_struct* as);

int sort_as_by_file_type(struct assign_struct* as)
{
        struct seq_bit_vec* bv = NULL;
        int i;
        for(i = 0; i < as->num_reads ;i++){
                bv = as->bits[i];
                qsort(bv->bits, bv->num_bit,sizeof(struct seq_bit*), qsort_seq_bits_by_type_file);
        }

        return OK;
}

int post_process_assign(struct assign_struct* as)
{
        struct seq_bit_vec* bv = NULL;
        struct demux_struct* tmp_ptr = NULL;

        char* tmp =  NULL;
        char* barcode;
        char* umi;
        int i,j,c,g;
        int len;
        int umi_len;
        int tmp_len = 256;

        int e_fail;
        int d_fail;


        e_fail = 0;
        d_fail = 0;
        MMALLOC(tmp, sizeof(char) * tmp_len);

        for(i = 0; i < as->num_reads ;i++){
                bv = as->bits[i];
                for(j = 0; j < as->num_files;j++){
                        if(bv->Q[j] == -1.0f){
                                e_fail++;
                        }

                        if(bv->Q[j] == -2.0f){
                                d_fail++;
                        }
                        if(bv->Q[j] < 0.0){
                                bv->fail = 1;
                                break;
                        }
                }

                ///* qsort by file identifier  */
                //qsort(bv->bits, bv->num_bit,sizeof(struct seq_bit*), qsort_seq_bits_by_type_file);
                len = 0;
                umi_len = 0;
                for(j = 0; j < bv->num_bit;j++){
                        if(bv->bits[j]->type == BAR_TYPE){
                                len += bv->bits[j]->len;// strnlen(bv->bits[j]->p,256);
                                len++;

                        }
                        if(bv->bits[j]->type == UMI_TYPE){
                                umi_len += bv->bits[j]->len;
                        }
                }
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

                        RUNP(tmp_ptr = as->demux_names->tree_get_data(as->demux_names,tmp));
                        //if(tmp_ptr){
                        tmp_ptr->count++;
                        //fprintf(stdout,"%s -> %d (%d)\n", tmp_ptr->name ,tmp_ptr->id,tmp_ptr->count);
                        bv->sample_group = tmp_ptr->id;
                        //bv->bc = tmp;
                        //tmp = NULL;
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
        }
        LOG_MSG("%d %d", e_fail,d_fail);
        MFREE(tmp);
        qsort(as->bits,as->num_reads ,sizeof(struct seq_bit_vec*) , qsort_seq_bit_vec);
        return OK;
ERROR:
        return FAIL;
}


int init_assign_structure(struct assign_struct** assign,struct arch_library* al, int total)
{
        struct assign_struct* as = NULL;
        //ASSERT(num_files >= 1,"no infiles");

        int* plan;
        int i,j;
        MMALLOC(as, sizeof(struct assign_struct));
        as->num_files = al->num_file;
        as->demux_names = NULL;
        as->bits = NULL;
        as->num_bits = 0;
        as->max_bar_len = 0;
        as->max_seq_len = 0;
        as->out_reads = 0;
        RUN(setup_assign_structure(al,as,&plan));
        RUN(setup_barcode_files(al,as));
        RUN(setup_output_files(as));
        //LOG_MSG("%d", as->out_reads);
        //exit(0);
        as->alloc_total = total;
        as->num_reads = 0;

        as->bits = NULL;
        MMALLOC(as->bits, sizeof(struct seq_bit_vec*)* as->alloc_total );
        //MMALLOC(as->active_bits,sizeof(uint8_t) * as->total));
        for(i = 0; i < as->alloc_total;i++){
                as->bits[i] = NULL;
                MMALLOC(as->bits[i], sizeof(struct seq_bit_vec));
                as->bits[i]->sample_group = -1;
                as->bits[i]->num_bit = as->num_bits;
                as->bits[i]->bits = NULL;
                as->bits[i]->Q = NULL;
                //as->bits[i]->bc = NULL;
                as->bits[i]->umi = NULL;
                MMALLOC(as->bits[i]->Q,sizeof(float) * as->num_files);
                as->bits[i]->fail = 0;
                MMALLOC(as->bits[i]->bits , sizeof(struct seq_bit) * as->num_bits);
                for(j = 0; j < as->num_bits;j++){
                        as->bits[i]->bits[j] = NULL;
                        MMALLOC(as->bits[i]->bits[j], sizeof(struct seq_bit));
                        as->bits[i]->bits[j]->file = plan[j];
                }
        }
        MFREE(plan);
        *assign = as;
        return OK;
ERROR:
        free_assign_structure(as);
        return FAIL;
}

int reset_assign_structute(struct assign_struct* as)
{
        int i;
        for(i = 0; i < as->alloc_total;i++){
                qsort(  as->bits[i]->bits,  as->bits[i]->num_bit,sizeof(struct seq_bit*), qsort_seq_bits_by_file);
                as->bits[i]->fail =0;
                as->bits[i]->sample_group = -1;
                if(as->bits[i]->umi){
                        MFREE(as->bits[i]->umi);
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
        char c;
        int p_size;
        int p_index;
        ASSERT(al != NULL,"No archlib");

        ASSERT(as != NULL,"No assign struct ");

        as->num_bits = 0;
        as->out_reads = 0;
        p_size = 256;
        p_index = 0;
        MMALLOC(p, sizeof(int) * p_size);
        for(i = 0; i < al->num_file;i++){
                read_structure = al->read_structure[al->arch_to_read_assignment[i]];
                //fprintf(stdout,"Read %d: ",i);
                for(j = 0; j < read_structure->num_segments;j++){
                        c = read_structure->type[j];
                        switch (c) {
                        case 'B':
                        case 'F':
                                p[p_index] = i;
                                p_index++;
                                if(p_index == p_size){
                                        p_size = p_size + p_size /2;
                                        MREALLOC(p, sizeof(int)* p_size);
                                }
                                as->num_bits++;
                                break;
                        case 'R':
                                p[p_index] = i;
                                p_index++;
                                if(p_index == p_size){
                                        p_size = p_size + p_size /2;
                                        MREALLOC(p, sizeof(int)* p_size);
                                }
                                as->out_reads++;
                                as->num_bits++;
                                break;

                        default:
                                break;
                        }
                }
                //fprintf(stdout,"\n");
        }
        /* create offsets  */
        *plan = p;
        return OK;
ERROR:
        return FAIL;
}


void free_assign_structure(struct assign_struct* as)
{
        if(as){
                int i,j;
                if(as->bits){
                        //MMALLOC(as->active_bits,sizeof(uint8_t) * as->total));
                        for(i = 0; i < as->alloc_total;i++){
                                for(j = 0; j < as->num_bits;j++){
                                        MFREE(as->bits[i]->bits[j]);
                                }
                                MFREE(as->bits[i]->Q);
                                MFREE(as->bits[i]->bits);
                                MFREE(as->bits[i]);
                        }
                        MFREE(as->bits);

                }
                if(as->demux_names){
                        as->demux_names->free_tree(as->demux_names);
                }
                if(as->file_names){
                        as->file_names->free_tree(as->file_names);
                }

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

int setup_barcode_files(struct arch_library* al, struct assign_struct* as)
{
        struct read_structure* read_structure = NULL;
        int i,j,g,f,len;
        char c;
        int num_barcodes = 0;
        ASSERT(al != NULL,"No archlib");
        struct rbtree_root* root = NULL;
        struct rbtree_root* root_new = NULL;
        struct demux_struct* tmp_ptr = NULL;
        struct demux_struct* new_ptr = NULL;
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

        //int test;

        //for(test = 0; test < 3;test++){
        for(i = 0; i < al->num_file;i++){
                read_structure = al->read_structure[al->arch_to_read_assignment[i]];
                //fprintf(stdout,"Read %d: ",i);
                for(j = 0; j < read_structure->num_segments;j++){
                        c = read_structure->type[j];
                        switch (c) {
                        case 'B':
                                if(!num_barcodes){
                                        for(g = 0;g < read_structure->numseq_in_segment[j];g++){
                                                //len = strnlen(read_structure->sequence_matrix[j][g], 256);
                                                len = read_structure->segment_length[j];
                                                if(len > as->max_bar_len){
                                                        as->max_bar_len = len;
                                                }
                                                MMALLOC(tmp_ptr, sizeof(struct demux_struct));
                                                tmp_ptr->name = NULL;
                                                MMALLOC(tmp_ptr->name,sizeof(char) * (len+1));
                                                //for(f = 0;f < len;f++){
                                                //tmp_ptr->name =
                                                //}
                                                strncpy(tmp_ptr->name, (char*)read_structure->sequence_matrix[j][g], len+1);

                                                tmp_ptr->name[len] =0;
                                                //snprintf(sample1->name , 10,"ABBBB");
                                                tmp_ptr->id = 0;
                                                tmp_ptr->count = 0;



                                                RUN(root->tree_insert(root,tmp_ptr));
                                                tmp_ptr = NULL;
                                        }

                                }else{
                                        RUN(root->flatten_tree(root));
                                        RUNP(root_new = init_tree(fp_get,fp_cmp,fp_cmp_same,fp_print,fp_free));
                                        for(f = 0;f < root->num_entries;f++){
                                                tmp_ptr = root->data_nodes[f];
                                                for(g = 0;g < read_structure->numseq_in_segment[j];g++){
                                                        len = strnlen(tmp_ptr->name, 256);
                                                        len++;
                                                        len += read_structure->segment_length[j];
                                                        //len += strnlen(read_structure->sequence_matrix[j][g], 256);
                                                        len++;
                                                        if(len > as->max_bar_len){
                                                                as->max_bar_len = len;
                                                        }

                                                        MMALLOC(new_ptr, sizeof(struct demux_struct));
                                                        new_ptr->name = NULL;

                                                        MMALLOC(new_ptr->name,sizeof(char) * (len+1));

                                                        snprintf(new_ptr->name, len, "%s_%s",tmp_ptr->name,read_structure->sequence_matrix[j][g]);
                                                        tmp_ptr->name[len] =0;
                                                        //snprintf(sample1->name , 10,"ABBBB");
                                                        new_ptr->id =0;
                                                        new_ptr->count = 0;
                                                        RUN(root_new->tree_insert(root_new,new_ptr));
                                                        new_ptr = NULL;
                                                }
//root->data_nodes[f]

                                        }
                                        root->free_tree(root);
                                        root = root_new;
                                        root_new = NULL;

                                }
                                num_barcodes++;
                                //as->num_bits++;
                                break;
                        default:
                                break;
                        }
                }
                //fprintf(stdout,"\n");
        }
        //}
        //root->print_tree(root,NULL);


        //tmp_ptr =  root->tree_get_data(root, "ACTTGA_ACAGTG_TTAGGC");


        //tmp_ptr->count++;
        RUN(root->flatten_tree(root));

        //root_new = init_tree(fp_get,fp_cmp,fp_cmp_same,fp_print,fp_free);
        for(f = 0;f < root->num_entries;f++){

                tmp_ptr = root->data_nodes[f];
                tmp_ptr->id = f;
                //fprintf(stdout,"%s %d\n",tmp_ptr->name,tmp_ptr->count);
        }
        //root->free_tree(root);
        as->demux_names = root;
        /* create offsets  */
        return OK;
ERROR:
        return FAIL;
}

void* get_name(void* ptr)
{
        struct demux_struct* tmp = (struct demux_struct*)  ptr;
        return tmp->name;
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
        return 1;
}

void free_demux_struct(void* ptr)
{
        struct demux_struct* tmp = (struct demux_struct*)  ptr;
        if(tmp){
                if(tmp->name){
                        MFREE(tmp->name);
                }
                MFREE(tmp);
        }
}
