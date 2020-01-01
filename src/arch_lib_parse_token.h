#ifndef ARCH_LIB_PARSE_TOKEN_H
#define ARCH_LIB_PARSE_TOKEN_H


#include "arch_lib.h"

extern int alloc_segment_spec(struct segment_specs** s);
extern void free_segment_spec(struct segment_specs*s);

extern int parse_rs_token_message(char* token, struct segment_specs** s_spec);
#endif
