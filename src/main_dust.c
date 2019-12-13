#include "tllogsum.h"

#include "interface.h"
#include "nuc_code.h"

#include "io.h"
#include "filter.h"


int main (int argc,char * argv[])
{
        struct parameters* param = NULL;
        int i,j;

        RUN(interface(&param,argc,argv));

        if(!param){
                return EXIT_SUCCESS;
        }

        init_logsum();
        RUN(init_nuc_code());


        free_param(param);

        return EXIT_SUCCESS;
ERROR:
        if(param){
                //fprintf(stdout,"%s",param->errmsg);
                free_param(param);
        }
        return EXIT_SUCCESS;
}
