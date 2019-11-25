#ifndef NUC_CODE_H
#define NUC_CODE_H


#if HAVE_CONFIG_H
#include "config.h"
#endif



#include "tldevel.h"


unsigned int nuc_code[256];

unsigned int rev_nuc_code[5];

extern int init_nuc_code(void);

#endif
