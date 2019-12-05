
#include "tldevel.h"
#include "tlhashtable.h"

HT_GLOBAL_INIT(TEST, char*)
HT_GLOBAL_INIT(TESTINT, int*)
HT_GLOBAL_INIT(TEST_DOUBLE, double)
int double_cmp(const void *a, const void *b);

struct hash_table_node_TEST_DOUBLE_t;
int double_cmp(const void *a, const void *b)
{
        hash_table_node_TEST_DOUBLE_t* const *ia = a;
        hash_table_node_TEST_DOUBLE_t* const *ib = b;
        if((*ia)->key  > (*ib)->key){
                return 1;
        }
        if((*ia)->key  > (*ib)->key){
                return 0;
        }
        return -1;

/* integer comparison: returns negative if b > a
   and positive if a > b */
}


int main (int argc,char * const argv[])
{

        int i,j;
        fprintf(stdout,"Hello world\n");
        char* f = NULL;


        RUN(galloc(&f,10));

        f[0] = 'a';
        f[1] = 'c';
        f[2] = 'c';
        f[3] = 'g';
        f[4] = 't';
        f[5] = 0;

        int test_int = 32;

        get_hash_value(test_int,10);
        fprintf(stdout,"HASH = %d\n",get_hash_value(test_int,100));
        fprintf(stdout,"HASH = %d\n",get_hash_value(f,100));
        // fprintf(stdout,"%d\n",get_hash_value(test_int));
        gfree(f);


        /*struct hash_table* ht = NULL;

        RUNP(ht = init_hash_table(345));

        for(i = 0; i < 2000;i++){
                insert(ht,i,NULL);
        }

        print_hash_table(ht);
        free_hash_table(ht);*/

        HT_TYPE(TEST)* my_ht = NULL;

        my_ht = HT_INIT(TEST,123);
        for(i = 0; i < 2000;i++){
                char* test = NULL;
                //MMALLOC(test, sizeof(char)*100);
                RUN(galloc(&test,100));
                int j,c;
                c = i;
                for(j = 0; j < 10;j++){
                        test[j] = (c & 0x3) +65;
                        c = c >> 2;
                }
                test[10] = 0;

                RUN(HT_INSERT(TEST,my_ht,test,NULL));

        }

        for(i = 0; i < 2000;i++){
                char* test = NULL;
                //MMALLOC(test, sizeof(char)*100);
                RUN(galloc(&test,100));
                int j,c;
                c = i;
                for(j = 0; j < 10;j++){
                        test[j] = (c & 0x3) +65;
                        c = c >> 2;
                }
                test[10] = 0;

                RUN(HT_INSERT(TEST,my_ht,test,NULL));

        }
        HT_PRINT(TEST,my_ht);
        HT_FREE(TEST,my_ht);


        HT_TYPE(TEST_DOUBLE )* my_htt = NULL;

        my_htt = HT_INIT(TEST_DOUBLE,123);
        for(i = 0; i < 2000;i++){
                RUN(HT_INSERT(TEST_DOUBLE,my_htt,(double)i,NULL));

        }
        for(i = 0; i < 2000;i++){
                RUN(HT_INSERT(TEST_DOUBLE,my_htt,(double)i,NULL));

        }
        hash_table_node_TEST_DOUBLE_t* hashnode = NULL;



        HT_PRINT(TEST_DOUBLE,my_htt);
        HT_FLATTEN(TEST_DOUBLE,my_htt);

        qsort(my_htt->flat,my_htt->num_items, sizeof(hash_table_node_TEST_DOUBLE_t*), double_cmp);
        for(i = 0; i < my_htt->num_items;i++){
                fprintf(stdout,"%d %f %d \n", i, my_htt->flat[i]->key,my_htt->flat[i]->count);
        }
        fprintf(stdout,"Hash table has %d entries, %d nodes and %d item counts\n", my_htt->table_size, my_htt->num_items, my_htt->total_count);

        LOG_MSG("searching for %f", (double) 42);

        hashnode = HT_SEARCH(TEST_DOUBLE, my_htt, 42);
        if(hashnode){
                LOG_MSG("Found:%f %d %p\n", hashnode->key,hashnode->count, hashnode->data);
        }
        HT_FREE(TEST_DOUBLE,my_htt);

        HT_TYPE(TESTINT )* int_array_table = NULL;
        int_array_table = HT_INIT(TESTINT,123);
        int* tmp = NULL;
        for(i = 0; i < 200;i++){
                tmp = NULL;
                RUN(galloc(&tmp,10));
                for(j = 0 ; j < 10;j++){
                        tmp[j] = rand() % 151;
                }

                RUN(HT_INSERT(TESTINT,int_array_table,tmp,NULL));

        }
        HT_PRINT(TESTINT,int_array_table);

        HT_FREE(TESTINT,int_array_table);
        return EXIT_SUCCESS;
ERROR:
        return EXIT_FAILURE;

}
