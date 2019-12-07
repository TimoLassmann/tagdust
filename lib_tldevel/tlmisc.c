#include <sys/stat.h>

#define TLMISC_IMPORT
#include "tlmisc.h"


int my_file_exists(char* name)
{
        struct stat buf;
        int ret,local_ret;
        ret = 0;
        local_ret= stat ( name, &buf );
        /* File found */
        if ( local_ret == 0 )
        {
                ret++;
        }
        return ret;
}
