#include "assign_data.h"


#include <string.h>

#include "tldevel.h"

#include "tlrbtree.h"


static int qsort_seq_bits_by_file(const void *a, const void *b);
static int qsort_seq_bit_vec(const void *a, const void *b);


static int set_up_barcode_files(struct arch_library* al);
int init_assign_structure(struct assign_struct** assign,struct arch_library* al, int total)
{
        struct assign_struct* as = NULL;
        //ASSERT(num_files >= 1,"no infiles");
        int i,j;
        MMALLOC(as, sizeof(struct assign_struct));
        as->num_files = al->num_file;
        as->bits = NULL;
        as->num_bits = 0;


        RUN(set_up_assign_structure(al,as));
        RUN(set_up_barcode_files(al));
        exit(0);
        as->total = total;

        as->bits = NULL;
        MMALLOC(as->bits, sizeof(struct seq_bit_vec*)* as->total);
        //MMALLOC(as->active_bits,sizeof(uint8_t) * as->total));
        for(i = 0; i < as->total;i++){
                as->bits[i] = NULL;
                MMALLOC(as->bits[i], sizeof(struct seq_bit_vec));
                as->bits[i]->num_bit = 0;
                as->bits[i]->bits = NULL;
                as->bits[i]->Q = NULL;
                as->bits[i]->bc = NULL;
                as->bits[i]->umi = NULL;
                MMALLOC(as->bits[i]->Q,sizeof(float) * as->num_files);
                as->bits[i]->pass = 1;
                MMALLOC(as->bits[i]->bits , sizeof(struct seq_bit) * as->num_bits);
                for(j = 0; j < as->num_bits;j++){
                        as->bits[i]->bits[j] = NULL;
                        MMALLOC(as->bits[i]->bits[j], sizeof(struct seq_bit));
                }
        }
        *assign = as;
        return OK;
ERROR:
        free_assign_structure(as);
        return FAIL;
}



int set_up_assign_structure(struct arch_library* al,struct assign_struct* as)
{
        struct read_structure* read_structure = NULL;
        int i,j;
        char c;
        ASSERT(al != NULL,"No archlib");

        ASSERT(as != NULL,"No assign struct ");

        as->num_bits = 0;
        for(i = 0; i < al->num_file;i++){
                read_structure = al->read_structure[al->arch_to_read_assignment[i]];
                //fprintf(stdout,"Read %d: ",i);
                for(j = 0; j < read_structure->num_segments;j++){
                        c = read_structure->type[j];
                        switch (c) {
                        case 'B':
                        case 'R':
                        case 'F':
                                as->num_bits++;
                                break;
                        default:
                                break;
                        }
                }
                //fprintf(stdout,"\n");
        }
        /* create offsets  */
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
                        for(i = 0; i < as->total;i++){
                                for(j = 0; j < as->num_bits;j++){
                                        MFREE(as->bits[i]->bits[j]);
                                }
                                MFREE(as->bits[i]->Q);
                                MFREE(as->bits[i]->bits);
                                MFREE(as->bits[i]);
                        }
                        MFREE(as->bits);
                }

                MFREE(as);
        }
}


int post_process_assign(struct assign_struct* as)
{
        struct seq_bit_vec* bv = NULL;
        char* tmp =  NULL;
        char* barcode;
        char* umi;
        int i,j,c,g;
        int len;
        int umi_len;
        for(i = 0; i < as->total;i++){
                bv = as->bits[i];

                /* qsort by file identifier  */

                qsort(bv->bits, bv->num_bit,sizeof(struct seq_bit*), qsort_seq_bits_by_file);
                len = 0;
                umi_len = 0;
                for(j = 0; j < bv->num_bit;j++){
                        if(i < 10){
                                fprintf(stdout,"%d %d %d\n",i,bv->bits[j]->type, bv->bits[j]->file);
                        }
                        if(bv->bits[j]->type == BAR_TYPE){
                                len += strnlen(bv->bits[j]->p,256);
                                len++;

                        }
                        if(bv->bits[j]->type == UMI_TYPE){

                                umi_len += bv->bits[j]->len;
                        }

                }
                if(len){
                        tmp = NULL;
                        MMALLOC(tmp, sizeof(char) * len);
                        g = 0;
                        for(j = 0; j < bv->num_bit;j++){
                                if(bv->bits[j]->type == BAR_TYPE){
                                        barcode = bv->bits[j]->p;
                                        len = strnlen(barcode,256);
                                        for(c = 0; c < len;c++){
                                                tmp[g] = barcode[c];
                                                g++;
                                        }
                                        tmp[g] = '_';
                                        g++;
                                }
                        }
                        tmp[g-1] = 0;
                        bv->bc = tmp;
                        tmp = NULL;
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

        }
        qsort(as->bits,as->total,sizeof(struct seq_bit_vec*) , qsort_seq_bit_vec);
        return OK;
ERROR:
        return FAIL;
}

int qsort_seq_bits_by_file(const void *a, const void *b)
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
        if ( (*elem1)->pass > (*elem2)->pass){
                return -1;
        }else if ((*elem1)->pass < (*elem2)->pass){
                return 1;
        }else{
                int c;
                c = strcmp((*elem1)->bc, (*elem2)->bc);
                if(c < 0){
                        return -1;

                }else if (c >= 0){
                        return 1;
                }
                return 0;
        }
}


struct test_struct{
        char* name;
        int number;
        float flt_num;
};



static void* get_name(void* ptr);
static long int compare_name(void* keyA, void* keyB);
static void print_test_struct(void* ptr,FILE* out_ptr);
static int resolve_default(void* ptr_a,void* ptr_b);
static void free_test_struct(void* ptr);

int set_up_barcode_files(struct arch_library* al)
{
        struct read_structure* read_structure = NULL;
        int i,j,g,f,len;
        char c;
        int num_barcodes = 0;
        ASSERT(al != NULL,"No archlib");
        struct rbtree_root* root = NULL;
        struct rbtree_root* root_new = NULL;
        struct test_struct* tmp_ptr = NULL;
        void*  (*fp_get)(void* ptr) = NULL;
        long int (*fp_cmp)(void* keyA, void* keyB)= NULL;
        int (*fp_cmp_same)(void* ptr_a,void* ptr_b);
        void (*fp_print)(void* ptr,FILE* out_ptr) = NULL;

        void (*fp_free)(void* ptr) = NULL;


        fp_get = &get_name;
        fp_cmp = &compare_name;
        fp_print = &print_test_struct;
        fp_cmp_same = &resolve_default;
        fp_free = &free_test_struct;

        root = init_tree(fp_get,fp_cmp,fp_cmp_same,fp_print,fp_free);


        for(i = 0; i < al->num_file;i++){
                read_structure = al->read_structure[al->arch_to_read_assignment[i]];
                //fprintf(stdout,"Read %d: ",i);
                for(j = 0; j < read_structure->num_segments;j++){
                        c = read_structure->type[j];
                        switch (c) {
                        case 'B':
                                if(!num_barcodes){
                                        for(g = 0;g < read_structure->numseq_in_segment[j];g++){
                                                len = strnlen(read_structure->sequence_matrix[j][g], 256);
                                                MMALLOC(tmp_ptr, sizeof(struct test_struct));
                                                tmp_ptr->name = NULL;
                                                MMALLOC(tmp_ptr->name,sizeof(char) * (len+1));
                                                strncpy(tmp_ptr->name, read_structure->sequence_matrix[j][g], len+1);
                                                //snprintf(sample1->name , 10,"ABBBB");
                                                tmp_ptr->number = 0;
                                                tmp_ptr->flt_num = 0.0f;
                                                root->tree_insert(root,tmp_ptr);
                                                tmp_ptr = NULL;
                                        }

                                }else{
                                        root->flatten_tree(root);
                                        root_new = init_tree(fp_get,fp_cmp,fp_cmp_same,fp_print,fp_free);
                                        for(f = 0;f < root->num_entries;f++){
                                                //root->data_nodes[f]
                                        }
                                }
                                //as->num_bits++;
                                break;
                        default:
                                break;
                        }
                }
                //fprintf(stdout,"\n");
        }
        root->print_tree(root,NULL);
        root->flatten_tree(root);
        //root_new = init_tree(fp_get,fp_cmp,fp_cmp_same,fp_print,fp_free);
        for(f = 0;f < root->num_entries;f++){
                tmp_ptr = root->data_nodes[f];
                fprintf(stdout,"%s\n",tmp_ptr->name);
        }
        root->free_tree(root);
        /* create offsets  */
        return OK;
ERROR:
        return FAIL;
}

void* get_name(void* ptr)
{
        struct test_struct* tmp = (struct test_struct*)  ptr;
        return tmp->name;
}

long int compare_name(void* keyA, void* keyB)
{
        return strcmp(keyA,keyB);
}


void print_test_struct(void* ptr,FILE* out_ptr)
{
        struct test_struct* tmp = (struct test_struct*)  ptr;
        fprintf(out_ptr,"%s\t%d\t%f\n",tmp->name, tmp->number,tmp->flt_num );
}

int resolve_default(void* ptr_a,void* ptr_b)
{
        return 1;
}

void free_test_struct(void* ptr)
{
        struct test_struct* tmp = (struct test_struct*)  ptr;
        if(tmp){
                if(tmp->name){
                        MFREE(tmp->name);
                }
                MFREE(tmp);
        }
}
