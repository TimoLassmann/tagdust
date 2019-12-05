#include <stdio.h>
#include <time.h>
#include <string.h>

#include "tldevel.h"
#include "tlmisc.h"

#define TLCHECKPOINT_IMPORT
#include "tlcheckpoint.h"


#define BUFFER_LEN 128


int set_checkpoint_file(struct checkpoint* chk,char* function,char* location,char* cmd)
{
        char buffer[BUFFER_LEN];
        FILE* f_ptr = NULL;

        struct tm *ptr;

        char time_string[200];

        time_t current = time(NULL);
        ptr = localtime(&current);

        if(!strftime(time_string, 200, "%F %H:%M:%S", ptr)){
                ERROR_MSG("Write failed");
        }
        snprintf(buffer,BUFFER_LEN ,"%s/%s_%d.chk", chk->base_dir,chk->base_name,chk->test_num );
        RUNP(f_ptr = fopen(buffer , "w" ));
        fprintf(f_ptr,"%*s: %s\n",MESSAGE_MARGIN, "command", cmd);
        fprintf(f_ptr,"%*s: %d\n",MESSAGE_MARGIN, "checkpoint ID", chk->test_num);
        fprintf(f_ptr,"%*s: %s\n",MESSAGE_MARGIN, "function", function);
        fprintf(f_ptr,"%*s: %s\n",MESSAGE_MARGIN, "called in", location);
        fprintf(f_ptr,"%*s: %s\n",MESSAGE_MARGIN, "at time", time_string);

        fclose(f_ptr);

        return OK;
ERROR:
        if(f_ptr){
                fclose(f_ptr);
        }
        return FAIL;
}

int test_for_checkpoint_file(struct checkpoint* chk,char* function,char* location, char* cmd)
{
        FILE* f_ptr = NULL;
        char buffer[BUFFER_LEN];
        static int8_t found = 0;

        snprintf(buffer,BUFFER_LEN ,"%s/%s_%d.chk", chk->base_dir,chk->base_name,chk->test_num );
        if(my_file_exists(buffer) && !found){
                RUNP(f_ptr = fopen(buffer , "r" ));
                /* get first line and compare to  */
                buffer[0]= 0;
                if(fscanf(f_ptr,"%*s %99[^\n]s",buffer) != 1){
                        ERROR_MSG("fscanf failed.");
                }
                fclose(f_ptr);
                //fprintf(stdout,"%s\n%s\n",cmd,buffer);
                if(!strncmp(cmd,buffer,99)){
                        return 1;
                }

                LOG_MSG("   Re-running: %s (%s)",function,location);
                LOG_MSG("   arguments have changed from:");
                LOG_MSG("     %s",cmd);
                LOG_MSG("   to:");
                LOG_MSG("     %s",buffer);
                LOG_MSG("   will re-run everything from this point.");

                found = 1;
        }else{
                found = 1;
        }
        return 0;
ERROR:
        LOG_MSG("test_for_checkpoint file has failed.");
        return 0;
}

struct checkpoint* init_checkpoint(char* base_name,char* target_dir)
{
        struct checkpoint* chk = NULL;
        size_t i = 0;
        int j;
        MMALLOC(chk, sizeof(struct checkpoint));

        chk->test_num = 0;
        chk->base_dir = NULL;
        chk->base_name = NULL;

        i = strlen(target_dir);
        MMALLOC(chk->base_dir, sizeof(char) * (i+1));

        for(j = 0;j < i;j++){
                chk->base_dir[j] = target_dir[j];
        }
        chk->base_dir[i] = 0;

        i = strlen(base_name);
        MMALLOC(chk->base_name, sizeof(char) * (i+1));
        for(j = 0;j < i;j++){
                chk->base_name[j] = base_name[j];
        }
        chk->base_name[i] = 0;

        return chk;
ERROR:
        return NULL;
}

void free_checkpoint(struct checkpoint* chk)
{
        if(chk){
                MFREE(chk->base_dir);
                MFREE(chk->base_name);
                MFREE(chk);
                chk = NULL;
        }
}
