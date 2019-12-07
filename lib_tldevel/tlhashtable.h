#ifndef TLHASH_TABLE_H
#define TLHASH_TABLE_H

#ifdef TLHASHTABLE_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif


#include <stdint.h>
#include <inttypes.h>
#include <stdio.h>

#define HT_GLOBAL_INIT(name, type)                                      \
                                                                        \
        typedef struct hash_table_node_##name{                          \
                type key;                                               \
                int count;                                              \
                void* data;                                             \
                struct hash_table_node_##name *next;                    \
        }hash_table_node_##name##_t;                                    \
                                                                        \
                                                                        \
        typedef struct {                                                \
                hash_table_node_##name##_t* head;                       \
                hash_table_node_##name##_t* tail;                       \
        }hash_table_item_##name##_t;                                    \
                                                                        \
        typedef struct{                                                 \
                hash_table_item_##name##_t** table;                     \
                struct hash_table_node_##name** flat;                   \
                int num_items;                                          \
                int total_count;                                        \
                int table_size;                                         \
        } hash_table_##name##_t;                                        \
                                                                        \
                                                                        \
                                                                        \
        hash_table_##name##_t* init_hash_table_##name(int size);        \
        void free_hash_table_##name(hash_table_##name##_t* ht);         \
        hash_table_node_##name##_t* alloc_hash_table_node_##name (type key,void* data); \
                                                                        \
        hash_table_##name##_t* init_hash_table_##name(int size)         \
        {                                                               \
                hash_table_##name##_t* ht = NULL;                       \
                int i;                                                  \
                MMALLOC(ht, sizeof(hash_table_##name##_t));             \
                ht->num_items = 0;                                      \
                ht->total_count = 0;                                    \
                ht->table = NULL;                                       \
                ht->table_size = size;                                  \
                ht->flat = NULL;                                        \
                MMALLOC(ht->table, sizeof(hash_table_item_##name##_t*) * ht->table_size); \
                for(i =0; i < ht->table_size;i++){                      \
                        ht->table[i] = NULL;                            \
                        MMALLOC(ht->table[i], sizeof(hash_table_item_##name##_t)); \
                        ht->table[i]->head = NULL;                      \
                        ht->table[i]->tail = NULL;                      \
                }                                                       \
                return ht;                                              \
        ERROR:                                                          \
                free_hash_table_##name(ht);                             \
                return NULL;                                            \
        }                                                               \
                                                                        \
        void free_hash_table_##name(hash_table_##name##_t* ht)          \
        {                                                               \
                hash_table_node_##name##_t* n = NULL;                   \
                hash_table_node_##name##_t* tmp = NULL;                 \
                int i;                                                  \
                if(ht){                                                 \
                        if(ht->table){                                  \
                                for(i =0; i < ht->table_size;i++){      \
                                        n = ht->table[i]->head;         \
                                        while(n){                       \
                                                tmp = n;                \
                                                n = n->next;            \
                                                gfree(tmp->key);        \
                                                MFREE(tmp);             \
                                        }                               \
                                        MFREE(ht->table[i]);            \
                                }                                       \
                                MFREE(ht->table);                       \
                        }                                               \
                        if(ht->flat){                                   \
                                MFREE(ht->flat);                        \
                        }                                               \
                        MFREE(ht);                                      \
                }                                                       \
        }                                                               \
                                                                        \
        int hash_table_flat__##name(hash_table_##name##_t* ht)          \
        {                                                               \
                int i,j;                                                \
                hash_table_node_##name##_t* n = NULL;                   \
                MMALLOC(ht->flat, sizeof(hash_table_node_##name##_t*) * ht->num_items); \
                j = 0;                                                  \
                for(i = 0; i < ht->table_size;i++){                     \
                        n = ht->table[i]->head;                         \
                        while(n){                                       \
                                ht->flat[j] = n;                        \
                                j++;                                    \
                                n = n->next;                            \
                        }                                               \
                }                                                       \
                return OK;                                              \
        ERROR:                                                          \
                return FAIL;                                            \
        }                                                               \
                                                                        \
        hash_table_node_##name##_t* get_entry_hash_table_linked_list_##name(hash_table_node_##name##_t* n, int pos) \
        {                                                               \
                int i = 0;                                              \
                hash_table_node_##name##_t* tmp = n;                    \
                while (i != pos){                                       \
                        tmp = tmp->next;                                \
                        i++;                                            \
                }                                                       \
                return tmp;                                             \
        }                                                               \
                                                                        \
        int search_hash_table_linked_list_##name(hash_table_node_##name##_t* n, type key) \
        {                                                               \
                int ret = 0;                                            \
                hash_table_node_##name##_t* tmp = n;                    \
                while (tmp != NULL){                                    \
                        if (!ht_compare_key(key,tmp->key)){             \
                                return ret;                             \
                        }                                               \
                        tmp = tmp->next;                                \
                        ret++;                                          \
                }                                                       \
                return -1;                                              \
        }                                                               \
                                                                        \
                                                                        \
        hash_table_node_##name##_t* search_##name(hash_table_##name##_t* ht, type key) \
        {                                                               \
        hash_table_node_##name##_t* n = NULL;                           \
        hash_table_node_##name##_t* ret = NULL;                         \
        uint32_t index;                                                 \
        index = get_hash_value(key,ht->table_size);                     \
        n = ht->table[index]->head;                                     \
        if(n){                                                          \
                int pos = search_hash_table_linked_list_##name(n, key); \
                if(pos != -1){                                          \
                        ret = get_entry_hash_table_linked_list_##name(n, pos); \
                        return ret;                                     \
                }                                                       \
        }                                                               \
        return NULL;                                                    \
        }                                                               \
         int insert_##name(hash_table_##name##_t* ht, type key, void* data) \
         {                                                              \
                 hash_table_node_##name##_t* n = NULL;                  \
                 hash_table_node_##name##_t* new = NULL;                \
                 uint32_t index;                                        \
                 index = get_hash_value(key,ht->table_size);            \
                 n = ht->table[index]->head;                            \
                 RUNP(new = alloc_hash_table_node_##name(key,data));    \
                 if(n == NULL){                                         \
                         ht->table[index]->head = new;                  \
                         ht->table[index]->tail = new;                  \
                         ht->num_items++;                               \
                 }else{                                                 \
                         int pos = search_hash_table_linked_list_##name(n, key); \
                         if (pos == -1){                                \
                                 ht->table[index]->tail->next = new;    \
                                 ht->table[index]->tail = new;          \
                                 ht->num_items++;                       \
                         }else{                                         \
                                 gfree(new->key);                       \
                                 MFREE(new);                            \
                                 new = get_entry_hash_table_linked_list_##name(n, pos); \
                                 new->count++;                          \
                         }                                              \
                 }                                                      \
                 ht->total_count++;                                     \
                 return OK;                                             \
         ERROR:                                                         \
                 return FAIL;                                           \
         }                                                              \
                                                                        \
         int print_hash_table_##name(hash_table_##name##_t* ht)         \
         {                                                              \
                 hash_table_node_##name##_t* n = NULL;                  \
                 int i;                                                 \
                 for (i = 0; i < ht->table_size;i++){                   \
                         n = ht->table[i]->head;                        \
                         if(n == NULL){                                 \
                                 fprintf(stdout,"%d\tno entry\n",i);    \
                         }else{                                         \
                                 fprintf(stdout,"%d\t",i);              \
                                 while(n){                              \
                                         ht_print_key(n->key);          \
                                         n = n->next;                   \
                                 }                                      \
                                 fprintf(stdout,"\n");                  \
                         }                                              \
                 }                                                      \
                 return OK;                                             \
         }                                                              \
                                                                        \
         hash_table_node_##name##_t* alloc_hash_table_node_##name (type key, void* data) \
         {                                                              \
                 hash_table_node_##name##_t* n = NULL;                  \
                 MMALLOC(n, sizeof(hash_table_node_##name##_t));        \
                 n->key = key;                                          \
                 n->data = data;                                        \
                 n->count = 1;                                          \
                 n->next = NULL;                                        \
                 return n;                                              \
         ERROR:                                                         \
                 return NULL;                                           \
         }                                                              \



uint32_t get_hash_value_int(const int x, const int table_size);
uint32_t get_hash_value_double(const double x,const int table_size);
uint32_t get_hash_value_string(const char* s, const int table_size);
uint32_t get_hash_value_int_array(int* x,const int table_size);

#define get_hash_value(type, size) _Generic ((type),                \
                                             char*: get_hash_value_string, \
                                             int*: get_hash_value_int_array, \
                                             int: get_hash_value_int,   \
                                             double: get_hash_value_double \
                )(type, size)

int ht_compare_key_int(const int a, const int b);
int ht_compare_key_int_star(int* a, int* b);
int ht_compare_key_double(const double a, const double b);
int ht_compare_key_strings(const char* a, const char* b);

#define ht_compare_key(a, b) _Generic ((a),                           \
                                       double: ht_compare_key_double, \
                                       int: ht_compare_key_int,       \
                                       char*: ht_compare_key_strings, \
                                       int*: ht_compare_key_int_star       \
                )(a, b)

int print_int(const int a);
int print_int_star(int* a);
int print_double(const double a);
int print_string(const char* a);

#define ht_print_key(a) _Generic ((a),                  \
                                  int: print_int,       \
                                  double: print_double, \
                                  char*: print_string,  \
                                  int*: print_int_star  \
                )(a)






#define HT_TYPE(name) hash_table_##name##_t
#define HT_INIT(name,size) init_hash_table_##name(size)
#define HT_INSERT(name,ht,key,data) insert_##name(ht,key,data)

#define HT_SEARCH(name,ht,key)  search_##name(ht,key);

#define HT_FLATTEN(name,ht)  hash_table_flat__##name(ht);
#define HT_PRINT(name,ht)  print_hash_table_##name(ht)
#define HT_FREE(name,ht) free_hash_table_##name(ht)

#undef TLHASHTABLE_IMPORT
#undef EXTERN
#endif
