
#include "tldevel.h"

#include <stdint.h>
#include <string.h>

#define TLBITVEC_IMPORT
#include "tlbitvec.h"


#define SHIFT 5
#define MASK 0x1F

struct bitvec{
        uint32_t* arr;
        int size;
};


int make_bitvector( bitvec** bv, int num_elem)
{
        struct bitvec* bitvec = NULL;

        ASSERT(num_elem > 0, "No elements");
        MMALLOC(bitvec, sizeof(struct bitvec));
        bitvec->size = num_elem;
        bitvec->arr = NULL;
        MMALLOC(bitvec->arr, sizeof(int) *((num_elem / 32) + 1));
        RUN(clear_bitvector(bitvec));
        //RUN(clear_bitvector(x, num_elem));
        *bv = bitvec;
        return OK;
ERROR:
        free_bitvector(&bitvec);
        return FAIL;
}

int clear_bitvector(bitvec* bv)
{
        ASSERT(bv!= NULL, "No bit vector allocated");
        memset( bv->arr, 0, sizeof(int) *((bv->size / 32) + 1));
        return OK;
ERROR:
        return FAIL;
}

int bit_set(bitvec* bv, int i)
{
        ASSERT(bv!= NULL, "No bit vector allocated");
        ASSERT(i >= 0, "i needs to be positive (is: %d)", i);
        ASSERT(i < bv->size, "i is larger than bit vector (is: %d > %d)", i,bv->size);
        bv->arr[i >> SHIFT] |= (1 << (i & MASK));
        return OK;
ERROR:
        return FAIL;
}



int bit_clr(bitvec* bv, int i)
{
        ASSERT(bv!= NULL, "No bit vector allocated");
        ASSERT(i >= 0, "i needs to be positive (is: %d)", i);
        ASSERT(i < bv->size, "i is larger than bit vector (is: %d > %d)", i,bv->size);

        bv->arr[i >> SHIFT] &= ~(1 << (i & MASK));
        return OK;
ERROR:
        return FAIL;
}

int bit_test(bitvec* bv, int i, int* ret)
{
        ASSERT(bv!= NULL, "No bit vector allocated");
        DASSERT(i >= 0, "i needs to be positive (is: %d)", i);
        DASSERT(i < bv->size, "i is larger than bit vector (is: %d > %d)", i,bv->size);
        *ret = (bv->arr[i >> SHIFT] & (1 << (i & MASK))) != 0 ;
        return OK;
ERROR:
        return FAIL;
}

int free_bitvector(bitvec** bv)
{
        if(*bv){
                MFREE((*bv)->arr);
                MFREE(*bv);
                *bv = NULL;
        }
        return OK;
}
