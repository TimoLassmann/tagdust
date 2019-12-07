

#include "tlcheckpoint.h"

#include <stdio.h>
int dummy_func(int i);

int main(int argc, char *argv[])
{
        char buffer[200];


        DECLARE_CHK(MAIN,".");
        snprintf(buffer,200,"dummy func");
        RUN_CHECKPOINT(MAIN,dummy_func(42),buffer);


        DESTROY_CHK(MAIN);
        //MFREE(cmd);
        return EXIT_SUCCESS;
ERROR:
        return EXIT_FAILURE;
}


int dummy_func(int i)
{
        return OK;
}
