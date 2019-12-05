
#include <string.h>

#include "tldevel.h"

#include "tlrbtree.h"

struct test_struct{
        char* name;
        int number;
        float flt_num;
};

long int compare_name(void* keyA, void* keyB);
long int compare_number(void* keyA, void* keyB);
long int compare_flt(void* keyA, void* keyB);

static int resolve_same_flt(void* a,void* b);


static int resolve_default(void* ptr_a,void* ptr_b);

static void* get_name(void* ptr);
static void* get_number(void* ptr);
static void* get_flt(void* ptr);

void print_test_struct(void* ptr,FILE* out_ptr);
void free_test_struct(void* ptr);


static void* get_name(void* ptr)
{
        struct test_struct* tmp = (struct test_struct*)  ptr;
        return tmp->name;
}

static void* get_number(void* ptr)
{
        struct test_struct* tmp = (struct test_struct*)  ptr;
        return &tmp->number;
}

static void* get_flt(void* ptr)
{
        struct test_struct* tmp = (struct test_struct*)  ptr;
        return &tmp->flt_num;
}


static int resolve_same_flt(void* ptr_a,void* ptr_b)
{
        struct test_struct* one = (struct test_struct*)ptr_a;
        struct test_struct* two =(struct test_struct*)ptr_b;
        //DPRINTF3("resolving: %f +  %f \n",one->flt_num,two->flt_num );
        one->number  = one->number+two->number;
        MFREE(two->name);
        MFREE(two);
        return 0;
}

static int resolve_default(void* ptr_a,void* ptr_b)
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


long int compare_name(void* keyA, void* keyB)
{
        return strcmp(keyA,keyB);
}


long int compare_number(void* keyA, void* keyB)
{
        int* num1 = (int*)keyA;
        int* num2 = (int*)keyB;

        return *num1 - *num2;
}

long int compare_flt(void* keyA, void* keyB)
{
        float* num1 = (float*)keyA;
        float* num2 = (float*)keyB;

        if(*num1 == *num2){
                return 0;
        }
        if(*num1 > *num2){
                return 1;
        }else{
                return -1;
        }
}



void print_test_struct(void* ptr,FILE* out_ptr)
{
        struct test_struct* tmp = (struct test_struct*)  ptr;
        fprintf(out_ptr,"%s\t%d\t%f\n",tmp->name, tmp->number,tmp->flt_num );
}


int main (int argc,char * argv[])
{
        struct test_struct* sample1 = NULL;
        struct test_struct* sample2 = NULL;
        struct test_struct* sample3 = NULL;

        struct test_struct* tmp_ptr = NULL;
        //int status;
        int x,i;

        float y;


        fprintf(stdout,"%d debug", DEBUGLEVEL);

        //DPRINTF3("LEVEL 3 !!!!!!");

        MMALLOC(sample1, sizeof(struct test_struct));
        sample1->name = NULL;
        MMALLOC(sample1->name,sizeof(char) * 10);
        snprintf(sample1->name , 10,"ABBBB");
        sample1->number = 6032;
        sample1->flt_num = 3.1;

        MMALLOC(sample2, sizeof(struct test_struct));
        sample2->name = NULL;
        MMALLOC(sample2->name,sizeof(char) * 10);
        snprintf(sample2->name , 10,"CCC");
        sample2->number = 23;
        sample2->flt_num = 2.1;

        MMALLOC(sample3, sizeof(struct test_struct));
        sample3->name = NULL;
        MMALLOC(sample3->name,sizeof(char) * 10);
        snprintf(sample3->name , 10,"ZAAA");
        sample3->number = 367;
        sample3->flt_num = 2.1;



        struct rbtree_root* root = NULL;
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
        root->tree_insert(root,sample1);
//	tree_insert(root, sample1);
        sample1 = NULL;

        root->tree_insert(root, sample2);
        sample2 = NULL;

        root->tree_insert(root, sample3);
        sample3 = NULL;





        //fprintf(stdout,"Sorted by names: (%d entries)\n", root->node->num);

        root->print_tree(root,NULL);

        /*MMALLOC(sample1, sizeof(struct test_struct));
          sample1->name = NULL;
          MMALLOC(sample1->name,sizeof(char) * 10);
          snprintf(sample1->name , 10,"ZAAA");
          sample1->number = 23;
          sample1->flt_num = 2.1;


          root->tree_delete(root, sample1);
          fprintf(stdout,"DELETED CCC (%d entries)\n", root->node->num);

          root->print_tree(root,NULL);
        */

        fprintf(stdout,"Search for CCC:\n");
        tmp_ptr = root->tree_get_data(root,"CCC");
        //tmp_ptr = tree_get_data(root,"CCC");
        if(tmp_ptr){
                print_test_struct(tmp_ptr,stdout);
        }

        fprintf(stdout,"\n");

        root->free_tree(root);


        sample1 = NULL;
        sample2 = NULL;
        sample3 = NULL;



        MMALLOC(sample1, sizeof(struct test_struct));
        sample1->name = NULL;
        MMALLOC(sample1->name,sizeof(char) * 10);
        snprintf(sample1->name , 10,"ABBBB");
        sample1->number = 6032;
        sample1->flt_num = 3.1;

        MMALLOC(sample2, sizeof(struct test_struct));
        sample2->name = NULL;
        MMALLOC(sample2->name,sizeof(char) * 10);
        snprintf(sample2->name , 10,"CCC");
        sample2->number = 23;
        sample2->flt_num = 2.1;

        MMALLOC(sample3, sizeof(struct test_struct));
        sample3->name = NULL;
        MMALLOC(sample3->name,sizeof(char) * 10);
        snprintf(sample3->name , 10,"ZAAA");
        sample3->number = 367;
        sample3->flt_num = 2.1;



        fp_get = &get_number;
        fp_cmp_same = &resolve_default;
        fp_cmp = &compare_number;
        fp_print = &print_test_struct;
        fp_free = &free_test_struct;

        root = init_tree(fp_get,fp_cmp,fp_cmp_same,fp_print,fp_free);
        root->tree_insert(root, sample1);
        root->tree_insert(root, sample2);
        root->tree_insert(root, sample3);
        fprintf(stdout,"Sorted by numbers: (%d entries)\n",root->num_entries);
        root->print_tree(root,NULL);


        fprintf(stdout,"Search for 6032:\n");
        x = 6032;
        tmp_ptr= root->tree_get_data(root,&x);
        if(tmp_ptr){
                print_test_struct(tmp_ptr,stdout);
        }
        fprintf(stdout,"\n");

        root->free_tree(root);
        sample1 = NULL;
        sample2 = NULL;
        sample3 = NULL;
        MMALLOC(sample1, sizeof(struct test_struct));
        sample1->name = NULL;
        MMALLOC(sample1->name,sizeof(char) * 10);
        snprintf(sample1->name , 10,"ABBBB");
        sample1->number = 6032;
        sample1->flt_num = 3.1;

        MMALLOC(sample2, sizeof(struct test_struct));
        sample2->name = NULL;
        MMALLOC(sample2->name,sizeof(char) * 10);
        snprintf(sample2->name , 10,"CCC");
        sample2->number = 23;
        sample2->flt_num = 2.1;

        MMALLOC(sample3, sizeof(struct test_struct));
        sample3->name = NULL;
        MMALLOC(sample3->name,sizeof(char) * 10);
        snprintf(sample3->name , 10,"ZAAA");
        sample3->number = 367;
        sample3->flt_num = 2.1;


        fp_get = &get_flt;
        fp_cmp = &compare_flt;
        fp_cmp_same = &resolve_same_flt;
        fp_print = &print_test_struct;
        fp_free = &free_test_struct;

        root = init_tree(fp_get,fp_cmp,fp_cmp_same,fp_print,fp_free);
        root->tree_insert(root, sample1);
        root->tree_insert(root, sample2);
        root->tree_insert(root, sample3);
        fprintf(stdout,"Sorted by floats:\n");
        root->print_tree(root,NULL);
        fprintf(stdout,"Search for 2.1:\n");
        y = 2.1f;
        tmp_ptr= root->tree_get_data(root,&y);
        if(tmp_ptr){
                print_test_struct(tmp_ptr,stdout);
        }
        root->flatten_tree(root);
        fprintf(stdout,"After Flatten:\n");
        ASSERT((root->cur_data_nodes == root->num_entries),"fail");

        for(i = 0; i < root->cur_data_nodes;i++){
                print_test_struct(root->data_nodes[i],stdout);
        }

        root->free_tree(root);

        //MFREE(sample3);

        return EXIT_SUCCESS;
ERROR:
        return EXIT_FAILURE;
}
