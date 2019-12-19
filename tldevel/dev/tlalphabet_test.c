
#include "tldevel.h"

#include "tlalphabet.h"

int print_alphabet(struct alphabet* a);


int main(int argc, char *argv[])
{
        struct alphabet* a = NULL;

        RUN(create_alphabet(&a,TLALPHABET_DEFAULT_PROTEIN));

        print_alphabet(a);
        MFREE(a);
        a = NULL;
        RUN(create_alphabet(&a,TLALPHABET_REDUCED_PROTEIN));

        print_alphabet(a);
        MFREE(a);

        RUN(create_alphabet(&a,TLALPHABET_DEFAULT_DNA));

        print_alphabet(a);
        MFREE(a);

        return EXIT_SUCCESS;
ERROR:
        return EXIT_FAILURE;
}



int print_alphabet(struct alphabet* a)
{
        LOG_MSG("Type: %d", a->type);
        fprintf(stdout,"LEN: %d\n",a->L);
        int i;
        for(i = 64;i < 96;i++){
                if(a->to_internal[i] != -1){
                        fprintf(stdout,"%c\t%d\n",  (char)i, a->to_internal[i]);
                }
        }
        return OK;
}
