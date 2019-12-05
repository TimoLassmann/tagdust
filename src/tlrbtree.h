
#ifndef TLRBTREE_H
#define TLRBTREE_H

#ifdef TLRBTREE_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif

#include <stdio.h>

typedef struct rbtree_node rbtree_node;

struct rbtree_root{
        struct rbtree_node* node;
        void** data_nodes;
        int cur_data_nodes;
        int num_entries;
        int identical_mode;
        /* User metqhods to deal with their data  */
        void* (*key)(void* ptr);
        long int (*compare)(void* keyA, void* keyB);
        int (*resolve_same)(void* ptr_a,void* ptr_b);
        void (*fp_print)(void* ptr,FILE* out_ptr);
        void (*fp_free)(void* ptr);
        /* Methods to manipulate tree */

        int(*tree_insert)(struct rbtree_root* root,void* data);
        int(*tree_delete) (struct rbtree_root* root, void* datanode);
        void(*free_tree)(struct rbtree_root* root);
        int(*flatten_tree) (struct rbtree_root* root);
        void (*print_tree) (struct rbtree_root* root,FILE* out_ptr);
        void* (*tree_get_data)(struct rbtree_root* root, void* key);
        struct rbtree_node*(*tree_get_node)(struct rbtree_root* root, void* key);
};


EXTERN struct rbtree_root* init_tree(void* (*key_function_pointer)(void* ptr), long int (*compare_function_pointer)(void* keyA, void* keyB),int (*resolve_same_pointer)(void* ptr_a,void* ptr_b),void (*fp_print)(void* ptr,FILE* out_ptr),void (*fp_free)(void* ptr));

//void tree_insert(struct rbtree_root* root, void* datanode);
//extern int tree_insert(struct rbtree_root* root, void* datanode);


#undef TLRBTREE_IMPORT
#undef EXTERN
#endif

