#ifndef TLCHECKPOINT_H
#define TLCHECKPOINT_H

#include "tldevel.h"

#ifdef TLCHECKPOINT_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif



struct checkpoint{
        char* base_dir;
        char* base_name;
        int test_num;
};


EXTERN struct checkpoint* init_checkpoint(char* base_name,char* target_dir);
EXTERN int set_checkpoint_file(struct checkpoint* chk,char* function,char* location,char* cmd);
EXTERN int test_for_checkpoint_file(struct checkpoint* chk,char* function,char* location, char* cmd);
EXTERN void free_checkpoint(struct checkpoint* chk);

#define DECLARE_CHK(n,dir) struct checkpoint* chk_##n = NULL;  RUNP( chk_##n =  init_checkpoint(TOSTRING(n),dir));

#define RUN_CHECKPOINT(n,EXP,CMD) do {                                  \
                if(test_for_checkpoint_file(chk_##n,TOSTRING(EXP),AT,CMD) ==0 ){ \
                        RUN(EXP);                                       \
                        RUN(set_checkpoint_file(chk_##n,TOSTRING(EXP),AT,CMD)); \
                }else{                                                  \
                        log_message("Skipping over: %s (%s)",TOSTRING(EXP),AT); \
                }                                                       \
                chk_##n->test_num += 1;                                 \
        }while (0)

#define DESTROY_CHK(n) if(chk_##n){free_checkpoint( chk_##n);};

#undef TLCHECKPOINT_IMPORT
#undef EXTERN

#endif
