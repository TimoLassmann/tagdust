

#include <stdio.h>
#include <stdlib.h>

#include "tldevel.h"

#define TLRBTREE_IMPORT
#include "tlrbtree.h"



#define RED 1
#define BLACK 0

#define ID_MODE_SUM 0
#define ID_MODE_MAX 1
#define ID_MODE_MIN 2
#define ID_MODE_NOTMERGE 3

struct rbtree_node{
        struct rbtree_node* right;
        struct rbtree_node* left;
        void* data_node;
        int color;
        unsigned int num;
};



static struct rbtree_root* flatten_tree_worker(struct rbtree_root* root,struct rbtree_node* n);
static struct rbtree_node* search_node(struct rbtree_root* root, struct rbtree_node* n,void* key);
static struct rbtree_node* insert_val_start(struct rbtree_root* root, struct rbtree_node* n,void* datanode);
static struct rbtree_node* delete_val_start(struct rbtree_root* root, struct rbtree_node* n,void* datanode);
static struct rbtree_node* insert_val(struct rbtree_root* root, struct rbtree_node* n,void* datanode);

static void free_tree( struct rbtree_root* root);
static int tree_insert(struct rbtree_root* root, void* datanode);
static int tree_delete(struct rbtree_root* root, void* datanode);
static void print_tree(struct rbtree_root* root,FILE* out_ptr);
static int flatten_tree(struct rbtree_root* root);
static void* tree_get_data(struct rbtree_root* root, void* key);
static struct rbtree_node* tree_get_node(struct rbtree_root* root, void* key);



int rank(struct rbtree_node* n);
static struct rbtree_node* set_num(struct rbtree_node* n);
static int size(struct rbtree_node* n);
static int isred(struct rbtree_node* n);
static struct rbtree_node* rotateLeft(struct rbtree_node* n);
static struct rbtree_node* rotateRight(struct rbtree_node* n);
static struct rbtree_node* flipColors(struct rbtree_node* n);


static struct rbtree_node* moveredleft(struct rbtree_node* n);
static struct rbtree_node* moveredright(struct rbtree_node* n);
static void* get_tree_min(struct rbtree_root* root);
static void* min_rb_data_node(struct rbtree_node* n);
static struct rbtree_node*  deleteMin(struct rbtree_node* n);
static struct rbtree_node* fixUp(struct rbtree_node* n);
static struct rbtree_node* delete(struct rbtree_root* root, struct rbtree_node* n, void* datanode);

static void free_rbtree(struct rbtree_node* n,void (*free_function_pointer)(void* ptr));
static int print_tree2(struct rbtree_node* node,void (*print_function_pointer)(void* ptr,FILE* out_ptr),FILE* out_ptr, int d);


void print_tree(struct rbtree_root* root,FILE* out_ptr)
{
        int d = 0;
        if(out_ptr){
                d = print_tree2(root->node,root->fp_print,out_ptr,d);
        }else{
                d = print_tree2(root->node,root->fp_print,stdout,d);
        }
}

int print_tree2(struct rbtree_node* n,void (*print_function_pointer)(void* ptr,FILE* out_ptr),FILE* out_ptr,int d)
{
        int i;
        if(n){
                if(n->left){
                        d++;
                        d = print_tree2(n->left,print_function_pointer,out_ptr,d);
                        d--;
                }
                if(print_function_pointer!= NULL){
                        for(i = 0; i < d;i++){
                                fprintf(out_ptr, "-");
                        }
                        print_function_pointer(n->data_node,out_ptr);
                }

                if(n->right){
                        d++;
                        d = print_tree2(n->right,print_function_pointer,out_ptr,d);
                        d--;
                }

        }
        return d;
}

int flatten_tree(struct rbtree_root* root)
{

        ASSERT(root->cur_data_nodes == 0,"root->cur_data_nodes is not 0");
        MMALLOC(root->data_nodes,root->num_entries * sizeof(void*));

        root->cur_data_nodes = 0;

        root = flatten_tree_worker(root, root->node);


        return OK;
ERROR:
        MFREE(root->data_nodes);
        return FAIL;
}

struct rbtree_root* flatten_tree_worker(struct rbtree_root* root,struct rbtree_node* n)
{
        if(n){
                if(n->left){
                        root = flatten_tree_worker(root,n->left);
                }

                root->data_nodes[root->cur_data_nodes] = n->data_node;
                root->cur_data_nodes++;

                if(n->right){
                        root = flatten_tree_worker(root,n->right);
                }

        }
        return root;
}

struct rbtree_root* init_tree(void* (*key_function_pointer)(void* ptr), long int (*compare_function_pointer)(void* keyA, void* keyB),int (*resolve_same_pointer)(void* ptr_a,void* ptr_b),void (*fp_print)(void* ptr,FILE* out_ptr),void (*fp_free)(void* ptr))
{
        struct rbtree_root* root = NULL;


        MMALLOC(root, sizeof(struct rbtree_root));

        root->node = NULL;
        root->data_nodes = NULL;
        root->num_entries = 0;
        root->cur_data_nodes = 0;
        root->key = key_function_pointer;
        root->compare = compare_function_pointer;
        root->resolve_same = resolve_same_pointer;
        root->fp_print = fp_print;
        root->fp_free = fp_free;
        root->tree_insert = tree_insert;
        root->tree_delete = tree_delete;
        root->free_tree = free_tree;
        root->print_tree = print_tree;
        root->flatten_tree = flatten_tree;
        root->tree_get_data = tree_get_data;
        root->tree_get_node = tree_get_node;
        return root;

ERROR:
        MFREE(root);
        return NULL;
}

int tree_insert(struct rbtree_root* root, void* datanode)
{
        RUNP(root->node = insert_val_start(root, root->node, datanode));
        root->num_entries = root->node->num;
        return OK;
ERROR:
        return FAIL;
}

int tree_delete(struct rbtree_root* root, void* datanode)
{
        RUNP(root->node = delete_val_start(root, root->node, datanode));
        root->num_entries = root->node->num;
        return OK;
ERROR:
        return FAIL;
}


void* tree_get_data(struct rbtree_root* root, void* key)
{
        struct rbtree_node* n = root->node;
        n = search_node(root,n,key);
        if(n){
                return n->data_node;
        }else{
                return NULL;
        }
}

struct rbtree_node*  tree_get_node(struct rbtree_root* root, void* key)
{
        struct rbtree_node* n = root->node;
        n = search_node(root,n,key);
        if(n){
                return n;
        }else{
                return NULL;
        }
}


struct rbtree_node* search_node(struct rbtree_root* root, struct rbtree_node* n,void* key)
{
        //DPRINTF3("searching: %p.",n);
        long int c;
        if(n){
                c = root->compare(root->key(n->data_node),key);
                if(c == 0){
                        return n;
                }else if(c > 0){
                        //DPRINTF3("going right.");
                        if(n->right){

                                return search_node(root, n->right , key);
                        }else{
                                return NULL;
                        }
                }else{
                        //DPRINTF3("going left.");

                        if(n->left){
                                return search_node(root, n->left , key);
                        }else{
                                return NULL;
                        }
                }
        }
        return n;
}


struct rbtree_node* insert_val_start(struct rbtree_root* root, struct rbtree_node* n,void* datanode)
{
        RUNP(n = insert_val(root,n,datanode));
        n->color = BLACK;
        return n;
ERROR:
        return NULL;

}

struct rbtree_node* delete_val_start(struct rbtree_root* root, struct rbtree_node* n,void* datanode)
{
        RUNP(n = delete(root,n, datanode));
        n->color = BLACK;
        return n;
ERROR:
        return NULL;
}

struct rbtree_node* delete(struct rbtree_root* root, struct rbtree_node* n, void* datanode)
{
        long int c;
        c = root->compare(root->key(n->data_node),root->key(datanode));
        //int cmp = key.compareTo(h.key);
        if (c < 0){
                if(!isred(n->left) && !isred(n->left->left)){
                        n  = moveredleft(n);
                }
                n->left = delete(root, n->left, datanode);

                //if (!isRed(h.left) && !isRed(h.left.left))
                //	h = moveRedLeft(h);
                //h.left =  delete(h.left, key);
        } else {
                if(isred(n->left)){
                        n = rotateRight(n);
                        n->color = n->right->color;
                        n->right->color = RED;
                }

                //if (isRed(h.left)) h = leanRight(h);

                if(c == 0 && (n->right == NULL)){
                        return NULL;
                }
                //if (cmp == 0 && (h.right == null))
                //	return null;

                if(!isred(n->right) && !isred(n->right->left)){
                        n = moveredright(n);
                }
                //if (!isRed(h.right) && !isRed(h.right.left))
                //	h = moveRedRight(h);
                if (c == 0) {

                        n->data_node = get_tree_min(root);
                        //n->data_node = min_rb_data_node(root,n->right, n->data_node);

                        n->right = deleteMin(n->right);
                        //h .key = min(h.right);
                        //h.value = get(h.right, h.key);
                        //	h.right = deleteMin(h.right);
                        //	else h.right = delete(h.right, key);
                }else{
                        n->right = delete(root, n->right, datanode);
                }

        }
        return fixUp(n);
}

struct rbtree_node*  deleteMin(struct rbtree_node* n)
{
        if(n->left == NULL){
                return NULL;
        }
        if(!isred(n->left) && !isred(n->left->left)){
                n = moveredleft(n);
        }
        n->left = deleteMin(n->left);

        return  fixUp(n);
}

void* get_tree_min(struct rbtree_root* root)
{
        return min_rb_data_node(root->node);
}


void* min_rb_data_node(struct rbtree_node* n)
{
        void* node_ptr = NULL;
        for(; n != NULL;n = n->left){
                node_ptr = n->data_node;
        }
        return node_ptr;
}


struct rbtree_node* moveredleft(struct rbtree_node* n)
{
        n = flipColors(n);
        if(isred(n->right->left)){
                n->right = rotateRight(n->right);
                n = rotateLeft(n);
                n = flipColors(n);
        }
        return n;

}
struct rbtree_node* moveredright(struct rbtree_node* n)
{
        n->color = BLACK;
        n->right->color = RED;

        if(!isred(n->left->left)){
                n->left->color = RED;
        }else{
                n = rotateRight(n);
                n->color = RED;
                n->left->color = BLACK;
        }
        return n;
}

struct rbtree_node* fixUp(struct rbtree_node* n)
{
        if(isred(n->right)){
                n = rotateLeft(n);
        }
        if(isred(n->left) && isred(n->left->left)){
                n = rotateRight(n);
        }
        if(isred(n->left) &&  isred(n->right )){
                n = flipColors(n);
        }

        return n;
}
/*	flipColors(<#struct rbtree_node *n#>)
    if (isRed(h.right))
    h = rotateLeft(h);
    if (isRed(h.left) && isRed(h.left.left))
    h = rotateRight(h);
    if (isRed(h.left) && isRed(h.right))
    colorFlip(h);
    return h; }
*/

struct rbtree_node* insert_val (struct rbtree_root* root, struct rbtree_node* n,void* datanode)
{
        long int c;
        if(n == 0){
                MMALLOC(n,sizeof(struct rbtree_node));

                n->left = 0;
                n->right = 0;
                n->color = RED;
                n->data_node = datanode;
                n->num = 1;
                return n;
        }

        if(isred(n->left)&&isred(n->right)){
                n = flipColors(n);
        }

        c = root->compare(root->key(n->data_node),root->key(datanode));
        if(c == 0){
                c = root->resolve_same(n->data_node,datanode);
                if(c){
                        n->right = insert_val(root, n->right , datanode);
                }
        }else if(c > 0){
                n->right = insert_val(root, n->right , datanode);
        }else{
                n->left = insert_val(root, n->left , datanode);
        }

        if(isred(n->right) && !isred(n->left)){
                n = rotateLeft(n);
        }

        if(isred(n->left)){
                if(isred(n->left->left)){
                        n = rotateRight(n);
                }
        }
        return set_num(n);
ERROR:
        MFREE(n);
        return NULL;
}

int rank(struct rbtree_node* n)
{
        if(n->left == 0){
                return 0;
        }
        return n->left->num;
}

struct rbtree_node* set_num(struct rbtree_node* n)
{
        n->num = size(n->left) + size(n->right) + 1;
        return n;
}

int size(struct rbtree_node* n)
{
        if (n == 0){
                return 0;
        }
        return n->num;
}


int isred(struct rbtree_node* n)
{
        if(!n){
                return 0;
        }
        return n->color == RED;
}

struct rbtree_node* rotateLeft(struct rbtree_node* n)
{
        struct rbtree_node* x;
        x = n->right;
        n->right = x->left;
        x->left = set_num(n);
        x->color = n->color;
        n->color = RED;
        return set_num(x);
}

struct rbtree_node* rotateRight(struct rbtree_node* n)
{
        struct rbtree_node* x;
        x = n->left;
        n->left = x->right;
        x->right = set_num(n);
        x->color = n->color;
        n->color = RED;
        return set_num(x);
}

struct rbtree_node* flipColors(struct rbtree_node* n)
{
        n->color = !(n->color);
        n->left->color = !(n->left->color);
        n->right->color = !(n->right->color);
        return n;
}


void free_tree( struct rbtree_root* root)
{

        free_rbtree(root->node,root->fp_free);

        if(root->data_nodes){
                MFREE(root->data_nodes);
        }

        MFREE(root);

}

void free_rbtree(struct rbtree_node* n,void (*free_function_pointer)(void* ptr))
{

        if(n){
                if(n->left){

                        free_rbtree(n->left,free_function_pointer);
                }
                if(n->right){

                        free_rbtree(n->right,free_function_pointer);
                }
                if(free_function_pointer!= NULL){

                        free_function_pointer(n->data_node);
                }
                //MFREE(n->data_node);
                MFREE(n);
        }
}



#ifdef RB_TEST

#include <string.h>

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
        DPRINTF3("resolving: %f +  %f \n",one->flt_num,two->flt_num );
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

        DPRINTF3("LEVEL 3 !!!!!!");

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





        fprintf(stdout,"Sorted by names: (%d entries)\n", root->node->num);

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
        tmp_ptr = tree_get_data(root,"CCC");
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
        tmp_ptr= tree_get_data(root,&x);
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
        tmp_ptr= tree_get_data(root,&y);
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


#endif

