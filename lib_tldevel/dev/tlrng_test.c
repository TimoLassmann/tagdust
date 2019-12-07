
#include <stdio.h>

#include "tldevel.h"
#include "tlrng.h"

int main(int argc, char *argv[])
{
        struct rng_state* rng = NULL;
        struct rng_state* rng_second = NULL;
        int i;
        RUNP(rng = init_rng(0));
        RUNP(rng_second = init_rng_from_rng(rng));
        for(i = 0; i < 10;i++){
                fprintf(stdout,"%f\t%d\n", tl_random_double(rng), tl_random_int(rng,10));
                fprintf(stdout,"%f\t%d\n", tl_random_double(rng_second), tl_random_int(rng_second,10));
        }
        free_rng(rng);
        return EXIT_SUCCESS;
ERROR:
        return EXIT_FAILURE;
}
