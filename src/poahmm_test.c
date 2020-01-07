
#include <stdint.h>
#include "tldevel.h"

#include "poahmm.h"

#include "sim_seq_lib.h"
int single_seq_test(void);

int main(int argc, char *argv[])
{
        single_seq_test();
        return EXIT_SUCCESS;
ERROR:
        return EXIT_FAILURE;
}

int single_seq_test(void)
{
        struct poahmm* poahmm = NULL;
        char* seq = NULL;
        uint8_t* i_seq = NULL;
        int len = 16;
        int nuc_count[4];

        int i;

        for(i = 0; i < 4;i++){
                nuc_count[i] = 1000;
        }
        RUNP(poahmm = init_poahmm(128,nuc_count,500.0));

        //RU

        //RUN(init_nodes_from_single_sequence(poahmm, uint8_t* seq, int len)

        free_poahmm(poahmm);
        return OK;
ERROR:
        return FAIL;

}
