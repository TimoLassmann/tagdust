#ifndef TLALPHABET_H
#define TLALPHABET_H

#ifdef TLALPHABET_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif

#define TLALPHABET_DEFAULT_PROTEIN 1
#define TLALPHABET_DEFAULT_DNA 2

#define TLALPHABET_REDUCED_PROTEIN 3

#include <stdint.h>

struct alphabet{
        int8_t to_internal[128];
        int8_t to_external[32];
        uint8_t type;
        uint8_t L;
};

EXTERN int create_alphabet(struct alphabet** alphabet, int type);

#undef TLALPHABET_IMPORT
#undef EXTERN
#endif
