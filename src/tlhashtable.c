
#include <string.h>
#include <stdio.h>

#include "tldevel.h"

#define TLHASHTABLE_IMPORT
#include "tlhashtable.h"




static uint32_t hash_uint32( uint32_t a);

/* print commands */
int print_int(const int a)
{
        fprintf(stdout,"%d ",a);
        return OK;
}

int print_int_star(int* a)
{

        int len;
        RUN(get_dim1(a,&len));
        int i;
        for(i = 0; i < len;i++){
                fprintf(stdout,"%d,",a[i]);
        }
        fprintf(stdout,"\n");
        return OK;
ERROR:
        return FAIL;
}

int print_double(const double a)
{
        fprintf(stdout,"%f ",a);
        return OK;
}

int print_string(const char* a)
{
        fprintf(stdout,"%s ",a);
        return OK;
}

int ht_compare_key_int(const int a, const int b)
{
        if(a > b){
                return -1;
        }
        if(a == b){
                return 0;
        }
        return 1;
}

int ht_compare_key_int_star(int* a, int* b)
{
        int len_a;
        int len_b;

        get_dim1(a, &len_a);
        get_dim1(b, &len_b);
        int min_l = MACRO_MIN(len_a, len_b);
        int i;
        for(i = 0.; i < min_l;i++){
                if(a[i] > b[i]){
                        return -1;
                }else if(a[i] < b[i]){
                        return 1;
                }
        }
        if(len_a > len_b){
                return -1;
        }

        return 0;
}

int ht_compare_key_strings(const char* a, const char* b)
{
        return strcmp(a,b);
}


int ht_compare_key_double(const double a, const double b)
{
        if(a > b){
                return -1;
        }

        if(a == b){
                return 0;
        }
        return 1;

}



/* Hash variable commands  */
uint32_t get_hash_value_double(const double x,const int table_size)
{
        ASSERT(table_size != 0, "Table size cannot be 0!");

        return (unsigned long long) x %table_size;
ERROR:
        return 0;
}

uint32_t get_hash_value_int(const int x,const int table_size)
{

        ASSERT(table_size != 0, "Table size cannot be 0!");
        return x % table_size;
ERROR:
        return 0;
}


uint32_t get_hash_value_string(const char* s,const int table_size)
{
        ASSERT(table_size != 0, "Table size cannot be 0!");
        uint32_t hash = 0;

        for(; *s; ++s)
        {
                // fprintf(stdout,"%c\n",*s);
                hash += *s;
                hash += (hash << 10);
                hash ^= (hash >> 6);
        }

        hash += (hash << 3);
        hash ^= (hash >> 11);
        hash += (hash << 15);

        return hash % table_size;
ERROR:
        return 0;
}

uint32_t hash_uint32( uint32_t a)
{
        a = (a+0x7ed55d16) + (a<<12);
        a = (a^0xc761c23c) ^ (a>>19);
        a = (a+0x165667b1) + (a<<5);
        a = (a+0xd3a2646c) ^ (a<<9);
        a = (a+0xfd7046c5) + (a<<3);
        a = (a^0xb55a4f09) ^ (a>>16);
        return a;
}

uint32_t get_hash_value_int_array(int* x,const int table_size)
{
        ASSERT(table_size != 0, "Table size cannot be 0!");
        uint32_t hash = 0;
        int len;
        int i;
        get_dim1(x,&len);
        hash = hash ^  hash_uint32(x[0]);
        for(i = 1; i < len;i++){
                hash = (hash << 5) ^ ( hash >> (27));
                hash ^= hash_uint32(x[i]);
        }
        return hash % table_size;
ERROR:
        return 0;
}
