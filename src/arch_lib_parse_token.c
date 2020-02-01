
#include <string.h>

#include "tldevel.h"

#include "tlseqio.h"



#include "arch_lib_parse_token.h"



#define BUFFER_LEN 256

static int parse_rs_token(char* token, struct segment_specs** s_spec);
static int detect_extract_type(char c, uint8_t* t);



int parse_rs_token_message(char* token, struct segment_specs** s_spec)
{
        struct segment_specs* s = NULL;
        int status;

        status = parse_rs_token(token, &s);

        if(status == OK){
                *s_spec = s;
        }else{
                ASSERT(s == NULL,"s should be NULL");
                WARNING_MSG("%s - unable to parse", token);
                fprintf(stdout,"\n");
                fprintf(stdout,"Have a look at the manual to see how to specify segments.\n");
                fprintf(stdout,"\n");
                fprintf(stdout,"TL;DR\n");
                fprintf(stdout,"\n");
                fprintf(stdout,"Valid formats:\n");
                fprintf(stdout,"\n");
                fprintf(stdout,"   <NAME>:<TYPE>:<SEQ1>,<SEQ2>,...\n");
                fprintf(stdout,"   <TYPE>:<SEQ1>,<SEQ2>,...\n");
                fprintf(stdout,"\n");
                fprintf(stdout,"   NAME - will be used to identify the extracted sequences.\n");
                fprintf(stdout,"   TYPE - valid options are: \n");
                fprintf(stdout,"      \"E\" - extract this read\n");
                fprintf(stdout,"      \"A\" - attach these sequences to the read name.\n");
                fprintf(stdout,"      \"S\" - split output based on matched to these sequences.\n");
                fprintf(stdout,"      \"I\" - ignore these sequences.\n");
                fprintf(stdout,"      \"P\" - partial sequence.\n");
                fprintf(stdout,"      \"L\" - partial sequence - optional residues 5'.\n");
                fprintf(stdout,"      \"R\" - partial sequence - optional residues 3'.\n");
                fprintf(stdout,"   SEQ - sequence consisting of A,C,G,T\n");
                fprintf(stdout,"      or:\n");
                fprintf(stdout,"      N+ - match as many N's as possible.\n");
                fprintf(stdout,"      N{n} - match N exactly n times.\n");
                fprintf(stdout,"      N{n,m} - match N exactly n .. m times.\n");
                fprintf(stdout,"\n");
                fprintf(stdout,"Or:\n");
                fprintf(stdout,"   <NAME>:<CORRECTED>:<QUAL>:file.txt\n");
                fprintf(stdout,"   e.g. CR:CB:CZ:whitelist.txt \n");
                fprintf(stdout,"\n");
                fprintf(stdout,"Some examples of valid options:\n");
                fprintf(stdout,"   READ1:E:N+\n");
                fprintf(stdout,"     extract a read into file \"XXX_READ1_XXX\"\n");
                fprintf(stdout,"\n");
                fprintf(stdout,"   READ1:E:N{4,10}\n");
                fprintf(stdout,"     extract a read of length 4 .. 10 into file \"XXX_READ1_XXX\"\n");
                fprintf(stdout,"\n");
                fprintf(stdout,"   BARCODE:S:AAA,TTT\n");
                fprintf(stdout,"     split output based on barcode options \"AAA\", \"TTT\"\n");
                fprintf(stdout,"     the name \"BARCODE\" is ignored.\n");
                fprintf(stdout,"     this is identical to S:AAA,TTT\n");
                fprintf(stdout,"\n");

                return FAIL;
        }
        return OK;
ERROR:
        return FAIL;
}

int parse_rs_token(char* token, struct segment_specs** s_spec)
{
        struct segment_specs* spec = NULL;
        int i;
        int j;
        int len;
        int f,g;
        int state;
        int pos;

        int o_bracket;
        int c_bracket;
        int comma;
        int plus;

        int max_seq_len;

        RUN(alloc_segment_spec(&spec));

        ASSERT(token != NULL, "No token");

        len = strlen(token);
        ASSERT(len > 2,"token %s is too short to contain anything meaningful", token);

        state = 0;
        pos= 0;
        for(i = 0; i < len;i++){
                //fprintf(stdout,"%d %c\n", i, token[i]);
                if(token[i] == ':'){
                        if(state == 0){ /* name of segment */
                                for(j = pos; j < i;j++){
                                        spec->name[j-pos] = token[j];
                                        if(j - pos == BUFFER_LEN-1){
                                                spec->name[BUFFER_LEN-1] = 0;
                                                WARNING_MSG("Name %s is long; will only use first %d letters",spec->name,BUFFER_LEN);
                                                break;
                                        }
                                }
                                spec->name[j-pos] = 0;
                                if(j == 0){
                                        ERROR_MSG("token starts with ':' -> no name can be detected");
                                }
                                pos = i+1;
                                state++;
                        }else if(state == 1){ /* extract type  */

                                if(token[pos+1] == ':'){
                                        RUN(detect_extract_type(token[pos], &spec->extract));
                                        pos = i+1;
                                        state++;
                                        i = len+1; /* quit loop */

                                }else{
                                        spec->extract = ARCH_ETYPE_APPEND_CORRECT;
                                        for(j = pos; j < i;j++){
                                                spec->correct_name[j-pos] = token[j];
                                        }
                                        spec->correct_name[j-pos] = 0;
                                        pos = i+1;
                                        state++;
                                }

                        }else if(state == 2){

                                for(j = pos; j < i;j++){
                                        spec->qual_name[j-pos] = token[j];
                                }
                                spec->qual_name[j-pos] = 0;
                                pos = i+1;
                                state++;
                                i = len+1;
                        }

                        //state++;
                }
        }
        ASSERT(state != 0, "Could not reach sequence part = state is 0");
        ASSERT(state <= 3, "Could not reach sequence part = state is %d",state);

        //LOG_MSG("State: %d: %s %s", state,spec->name,spec->correct_name,spec->qual_name );
        //exit(0);
        if(state == 1){
                ASSERT(pos < len, "No sequence?");
                //LOG_MSG("name: %s %d %d",spec->name,pos,len);
                j = strlen(spec->name);
                ASSERT(j == 1,"Expecting length %d",j);

                RUN(detect_extract_type(spec->name[0], &spec->extract));
                if(spec->extract == ARCH_ETYPE_EXTRACT){
                        snprintf(spec->name, BUFFER_LEN, "R");
                }else{
                        snprintf(spec->name, BUFFER_LEN, "NA");
                }
                //exit(0);
        }else if(state == 2){
                ASSERT(pos < len, "No sequence?");
                if(!strncmp(spec->name, "NA", BUFFER_LEN)){
                        ERROR_MSG("Name: \"NA\" has a special meaning in tagdust");
                }
        }else if(state == 3){
                spec->extract = ARCH_ETYPE_APPEND_CORRECT;
                spec->num_seq =1;
                MMALLOC(spec->seq, sizeof(char*) * spec->num_seq);

                for(i = 0; i < spec->num_seq;i++){
                        spec->seq[i] = NULL;
                        MMALLOC(spec->seq[0], sizeof(char) * (len));

                }
                //MMALLOC(spec->seq[0], sizeof(char) * (j+1));
                //MMALLOC(spec->seq[1], sizeof(char) * (len));
                memcpy(spec->seq[0], token +pos, len- pos);
                spec->seq[0][len-pos] = 0;

//LOG_MSG("read: %d barcodes" , sb->num_seq);
                //ASSERT(
                //RUN(fill_exact_hash(&spec->bar_hash , sb));
                ///as->n_exact_hash++;

                //RUN(alloc_kmer_counts(&k, 16));
                //RUN(add_counts(k, sb));




                //RUN(run_build_pst(&p, k));
                //RUN(rm_counts(k,sb));

#ifdef DEBUG
                RUN(print_segment_spec(spec));
#endif

                *s_spec = spec;
                return OK;

        }
        ASSERT(pos != 0, "Could not reach sequence part - pos is 0");
        o_bracket = 0;
        c_bracket = 0;
        comma = 0;
        plus = 0;
        max_seq_len = 0;
        j = 0;
        for(i = pos; i < len;i++){
                if(token[i] == '{'){
                        o_bracket++;
                }
                if(token[i] == '}'){
                        c_bracket++;
                }
                if(token[i] == ','){
                        if(j+1 > max_seq_len){
                                max_seq_len = j+1;
                        }
                        comma++;
                        j = 0;
                }
                if(token[i] == '+'){
                        plus++;
                }
                j++;

        }
        ASSERT(o_bracket == c_bracket, "Brackets don't match");

        ASSERT(o_bracket <= 1, "To many brackets");


        if(o_bracket == 1){
                ASSERT(comma <= 1, "Only zero or one comma is allowed when specifying length of segment");
                ASSERT(plus == 0, "+ and brackets not allowed");
        }
        if(comma){
                ASSERT(plus == 0, "+ and commas not allowed");
                ASSERT(max_seq_len != 0,"sequence length is too short");
        }else{
                max_seq_len = len -pos+1 ;
        }

        /* start reading in sequences;  */
        if(o_bracket == 1){     /* a single variable length sequence */
                spec->num_seq = 1;
        }else if(o_bracket == 0 && comma != 0){
                spec->num_seq = comma +2;
        }else{
                spec->num_seq = 1;
        }
        MMALLOC(spec->seq, sizeof(char*) * spec->num_seq);

        for(i = 0; i < spec->num_seq;i++){
                spec->seq[i] = NULL;
                MMALLOC(spec->seq[i], sizeof(char) * max_seq_len);
        }

        if(o_bracket == 0 && comma ){
                f = 1;
                g = 0;
                spec->min_len = INT32_MAX;
                spec->max_len = INT32_MIN;
//fprintf(stderr,"%d\n",len);
                for(i = pos;i <  len;i++){
                        //fprintf(stdout,"%d %d %c\n",f,g,token[i]);
                        if(token[i] != ','){
                                spec->seq[f][g] = token[i];

                                g++;
                        }else{
                                spec->seq[f][g] = 0;
                                f++;
                                if(g < spec->min_len){
                                        spec->min_len =g;
                                }
                                if(g > spec->max_len){
                                        spec->max_len = g;
                                }
                                g = 0;
                        }
                }
                spec->seq[f][g] = 0;
                if(g < spec->min_len){
                        spec->min_len =g;
                }
                if(g > spec->max_len){
                        spec->max_len = g;
                }

                f = 0;
                g = 0;

                for(i = 0; i <  spec->max_len;i++){
                        spec->seq[f][g] = 'N';
                        g++;
                }
                spec->seq[f][g] = 0;
        }else if(o_bracket == 1){
                f = 0;
                g = 0;
                for(i = pos;i <  len;i++){
                        //fprintf(stdout,"%d %d %c\n",f,g,tmp[i]);
                        if(token[i] != '{'){
                                spec->seq[f][g] = token[i];
                                g++;
                        }else{
                                spec->seq[f][g] = 0;

                                g = sscanf(token + i, "{%d,%d}",&spec->min_len,&spec->max_len);
                                if(g == 0){
                                        ERROR_MSG("Could not read min/max len");
                                        LOG_MSG("unable to read %d %d", spec->min_len,spec->max_len);
                                }else if(g == 1){
                                        if(comma){
                                                ERROR_MSG("could not read maxlen");
                                        }
                                        spec->max_len = spec->min_len;
                                        //LOG_MSG("read1 %d %d", min_len,max_len);
                                }
                                break;
                        }
                }
                g = strnlen(spec->seq[f], max_seq_len);
                //spec->seq[f][g] = 0;
                //LOG_MSG("%s",spec->seq[f]);
                if(g != 1){
                        ERROR_MSG("Tagdust only accepts a single character (not %d) with options + {n} {n,m}",g);
                }

                MREALLOC(spec->seq[f], sizeof(char) *(spec->max_len+1));
                for(g = 1;g < spec->max_len;g++){
                        spec->seq[f][g] = spec->seq[f][0];
                }
                spec->seq[f][spec->max_len] = 0;
                //LOG_MSG("%s",spec->seq[f]);
        }else if(plus){
                f = 0;
                g = 0;
                for(i = pos;i <  len;i++){
                        //fprintf(stdout,"%d %d %c\n",f,g,tmp[i]);
                        if(token[i] != '+'){
                                spec->seq[f][g] = token[i];

                                g++;
                        }else{
                                spec->seq[f][g] = 0;
                                spec->min_len = g;
                                if(g != 1){
                                        ERROR_MSG("Tagdust only accepts a single character (not %d) with options + {n} {n,m}",g);
                                }
                        }
                }

                spec->max_len = INT32_MAX;

        }else{
                f = 0;
                g = 0;
                for(i = pos;i <  len;i++){
                        spec->seq[f][g] = token[i];
                        g++;
                }
                spec->seq[f][g] = 0;
                spec->max_len = g;
                spec->min_len = g;

        }
        ASSERT(spec->max_len >= spec->min_len,"max len is greater than min_len");

        if(spec->num_seq > 1){
                ASSERT(spec->min_len == spec->max_len, "sequences have different lengths");
        }
        /* set alloc_len */
        /* {n,m} -> longer one  */
        /* + - set to one  */
        spec->alloc_len = spec->max_len;
        if(spec->alloc_len == INT32_MAX){
                spec->alloc_len = 1;
        }
        /* set min_len to 1 for L and R segments  */


        if(spec->extract == ARCH_ETYPE_WOBBLE_LEFT){
                spec->min_len = 1;
        }

        if(spec->extract == ARCH_ETYPE_WOBBLE_RIGHT){
                spec->min_len = 1;
        }

        ASSERT(spec->alloc_len != 0, "Alloc len can't be 0");
#ifdef DEBUG
        fprintf(stdout,"%s input\n", token);

        RUN(print_segment_spec(spec));
#endif

        *s_spec = spec;

        return OK;
ERROR:
        free_segment_spec(spec);
        *s_spec = NULL;
        WARNING_MSG("Parsing error occurred in this token: %s",token);
        return FAIL;
}

int print_segment_spec(const struct segment_specs* spec)
{
        int i;
        ASSERT(spec != NULL, "No spec");
        fprintf(stdout,"   %s name\n", spec->name);
        fprintf(stdout,"   %s name\n", spec->correct_name);
        fprintf(stdout,"   %s name\n", spec->qual_name);
        fprintf(stdout,"   %d type\n", spec->extract);
        fprintf(stdout,"   %d %d (%d)  min,max len\n", spec->min_len,spec->max_len, spec->alloc_len);
        ASSERT(spec->num_seq != 0, "no sequence!");
        for(i = 0; i < spec->num_seq;i++){
                fprintf(stdout,"   %s\n", spec->seq[i]);

        }
        return OK;
ERROR:
        return FAIL;
}

int detect_extract_type(char c, uint8_t* t)
{
        *t = 0;
        switch (c) {
        case 'E': {
                *t =  ARCH_ETYPE_EXTRACT;
                break;
        }
        case 'A': {
                *t =  ARCH_ETYPE_APPEND;
                break;
        }
        case 'S': {
                *t =  ARCH_ETYPE_SPLIT;
                break;
        }
        case 'I': {
                *t =  ARCH_ETYPE_IGNORE;
                break;
        }
        case 'P': {
                *t =  ARCH_ETYPE_PARTIAL;
                break;
        }

        case 'L': {
                *t =  ARCH_ETYPE_WOBBLE_LEFT;
                break;
        }

        case 'R': {
                *t =  ARCH_ETYPE_WOBBLE_RIGHT;
                break;
        }

        default:

                ERROR_MSG("Extract type %c not recognised (allowed is: \"EASIPLR\")", c);
                break;
        }

        return OK;
ERROR:
        return FAIL;
}


int alloc_segment_spec(struct segment_specs** s)
{
        struct segment_specs* spec = NULL;

        MMALLOC(spec, sizeof(struct segment_specs));
        spec->max_len = 0;
        spec->min_len = 0;
        spec->extract = 0;
        spec->num_seq = 0;
        spec->alloc_len = 0;
        spec->seq = NULL;
        spec->name = NULL;
        spec->correct_name = NULL;
        spec->qual_name = NULL;
        spec->bar_hash = NULL;
        spec->pst = NULL;
        MMALLOC(spec->name, sizeof(char) * BUFFER_LEN);
        MMALLOC(spec->correct_name , sizeof(char) * BUFFER_LEN);
        MMALLOC(spec->qual_name , sizeof(char) * BUFFER_LEN);

        spec->name[0] = 0;
        spec->correct_name[0] =0;
        spec->qual_name[0] =0;
        *s = spec;
        return OK;
ERROR:

        return FAIL;
}

void free_segment_spec(struct segment_specs*s)
{
        int i;
        if(s){
                if(s->pst){
                        //free_pst(s->pst);
                }
                if(s->bar_hash){

                        //kh_destroy(exact,  s->bar_hash);
                }
                if(s->seq){
                        for(i = 0; i < s->num_seq;i++){
                                MFREE(s->seq[i]);
                        }
                        MFREE(s->seq);
                }
                if(s->name){
                        MFREE(s->name);
                }
                if(s->correct_name){
                        MFREE(s->correct_name);
                }
                if(s->qual_name){
                        MFREE(s->qual_name);
                }
                MFREE(s);
        }
}



