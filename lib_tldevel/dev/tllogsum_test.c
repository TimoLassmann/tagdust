

#include "tldevel.h"
#include "tllogsum.h"

int main(int argc, char *argv[])
{
        init_logsum();

        float a,b;
        ASSERT(1 != 0,"How odd");

        a = prob2scaledprob(0.4f);
        b = prob2scaledprob(0.4f);

        LOG_MSG("SUM: %f",scaledprob2prob(logsum(a, b)));
        return EXIT_SUCCESS;
ERROR:
        return EXIT_FAILURE;

}
