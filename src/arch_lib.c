
#include <string.h>
#include <ctype.h>

#include "arch_lib.h"

#include "tldevel.h"

#include "tlmisc.h"

#define BUFFER_LEN 256

static int malloc_read_structure(struct read_structure** read_structure);
static int resize_read_structure(struct read_structure* rs);
static void free_read_structure(struct read_structure* read_structure);

static int alloc_segment_spec(struct segment_specs** s);
static void free_segment_spec(struct segment_specs*s);
static int detect_extract_type(char c, uint8_t* t);

static int assign_segment_sequences(struct read_structure* read_structure, char* tmp, int segment);
static int resize_arch_lib(struct arch_library* al);

static int QC_read_structure(struct read_structure* read_structure );
static int sanity_check_arch_lib(struct arch_library* al);

static int parse_rs_token(char* token, struct segment_specs** s_spec);
static int parse_rs_token_message(char* token, struct segment_specs** s_spec);

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
                fprintf(stdout,"   SEQ - sequence consisting of A,C,G,T\n");
                fprintf(stdout,"      or:\n");
                fprintf(stdout,"      N+ - match as many N's as possible.\n");
                fprintf(stdout,"      N{n} - match N exactly n times.\n");
                fprintf(stdout,"      N{n,m} - match N exactly n .. m times.\n");
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

                                RUN(detect_extract_type(token[pos], &spec->extract));
                                pos = i+1;
                                i = len+1; /* quit loop */
                                state++;
                        }
                        //state++;
                }
        }
        ASSERT(state != 0, "Could not reach sequence part = state is 0");
        ASSERT(state <= 2, "Could not reach sequence part = state is %d",state);


        if(state == 1){
                ASSERT(pos < len, "No sequence?");
                //LOG_MSG("name: %s %d %d",spec->name,pos,len);
                j = strlen(spec->name);
                ASSERT(j == 1,"Expecting length %d",j);

                RUN(detect_extract_type(spec->name[0], &spec->extract));
                snprintf(spec->name, BUFFER_LEN, "NA");
                //exit(0);
        }else if(state == 2){
                ASSERT(pos < len, "No sequence?");
                if(!strncmp(spec->name, "NA", BUFFER_LEN)){
                        ERROR_MSG("Name: \"NA\" has a special meaning in tagdust");
                }
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
                spec->seq[f][g] = 0;

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
//#ifdef DEBUG
        fprintf(stdout,"%s input\n", token);

        fprintf(stdout,"   %s name\n", spec->name);
        fprintf(stdout,"   %d type\n", spec->extract);
        fprintf(stdout,"   %d %d  min,max len\n", spec->min_len,spec->max_len);
        for(i = 0; i < spec->num_seq;i++){
                fprintf(stdout,"   %s\n", spec->seq[i]);

        }

        *s_spec = spec;

//#endif
        return OK;
ERROR:
        free_segment_spec(spec);
        *s_spec = NULL;
        WARNING_MSG("Parsing error occurred in this token: %s",token);
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

        default:

                ERROR_MSG("Extract type %c not recognised (allowed is: \"EASIP\")", c);
                break;
        }

        return OK;
ERROR:
        return FAIL;
}


int read_architecture_files(struct arch_library* al, char* filename)
{
        FILE* f_ptr = NULL;
        char* line_buffer = NULL;
        size_t line_len = 0;
        ssize_t nread;
        char* s_ptr = NULL;
        char* token;
        int c;
        int target;
        //int read;
        int n_segment;
        int i,index;

        char* tmp = NULL;

        ASSERT(my_file_exists(filename),"Could not find file %s.",filename);

        RUNP(f_ptr = fopen(filename, "r"));
        while ((nread = getline(&line_buffer, &line_len, f_ptr)) != -1) {
                s_ptr = strstr( line_buffer, "tagdust");
                if(s_ptr){
                        line_buffer[strcspn(line_buffer, "\r\n")] = 0;
                        tmp = NULL;
                        MMALLOC(tmp,sizeof(char) * line_len);
                        index= 0;
                        //LOG_MSG("%s", line_buffer);
                        if(al->num_arch == al->alloc_num_arch){
                                resize_arch_lib(al);
                        }
                        target = al->num_arch;
                        n_segment = -1;
                        token = strtok(line_buffer, " ");
                        while(token != NULL){
                                c = strnlen(token,  line_len);
                                //fprintf(stdout,"t:%s\n",token);
                                if(token[0] == '-' && isdigit((int) token[1])){
                                        n_segment = atoi(token+1);
                                        n_segment--;
                                        //LOG_MSG("setsegment to %d", n_segment);
                                }else{
                                        if(n_segment != -1){
                                                RUN(assign_segment_sequences(al->read_structure[target], token, n_segment));
                                                //LOG_MSG("readinto %d", n_segment);
                                                n_segment = -1;
                                                for(i = 0;i < c;i++){
                                                        tmp[index] = token[i];
                                                        index++;
                                                }
                                                tmp[index] = ' ';
                                                index++;

                                        }
                                }
                                token = strtok(NULL, " ");
                        }
                        tmp[index-1] = 0;
                        al->spec_line[target] = tmp;
                        al->num_arch++;
                }
        }
        fclose(f_ptr);
        MFREE(line_buffer);
        RUN(sanity_check_arch_lib(al));
        return OK;
ERROR:
        if(line_buffer){
                MFREE(line_buffer);
        }
        return FAIL;
}

int read_arch_into_lib(struct arch_library* al, char** list, int len)
{
        int i,j,c;
        int target;
        char* buffer = NULL;
        char* tmp = NULL;
        int line_len;
        MMALLOC(buffer, sizeof(char) * BUFFER_LEN);

        ASSERT(al != NULL, "No architecture_library");

        if(al->num_arch == al->alloc_num_arch){
                resize_arch_lib(al);
        }
        target = al->num_arch;
        /* create spec line */
        line_len = 0;
        for(i = 0; i < len;i++){
                line_len += strlen(list[i]);
        }
        line_len += len + 1;
        MMALLOC(tmp,sizeof(char) * line_len);
        c = 0;
        for(i = 0; i < len;i++){
                line_len = strlen(list[i]);
                for(j = 0;j< line_len;j++){
                        tmp[c] = list[i][j];
                        c++;
                }
                tmp[c] = ' ';
                c++;
        }
        tmp[c-1] = 0;
        al->spec_line[target] = tmp;
        for(i = 0; i < len;i++){
                if(list[i]){
                        snprintf(buffer, BUFFER_LEN, "%s", list[i]);
                        RUN(assign_segment_sequences(al->read_structure[target], buffer, i));
                }
        }
        RUN(QC_read_structure(al->read_structure[target]));
        al->num_arch++;
        MFREE(buffer);
        RUN(sanity_check_arch_lib(al));
        return OK;
ERROR:
        if(buffer){
                MFREE(buffer);
        }
        return FAIL;
}

int assign_segment_sequences(struct read_structure* read_structure, char* tmp, int segment)
{

        int i,f,g;
        int count;
        int len;

        count = 0;

        switch (tmp[0]) {
        case 'R':
        case 'G':
        case 'O':
        case 'P':
        case 'S':
        case 'F':
        case 'B':
        case 'L':               /* DO I NEED THIS???  */
                break;

        default:
                ERROR_MSG("Segment type :%c not recognized.\n",tmp[0]);
                break;
        }

        len = strlen(tmp);
        count = 0;
        for(i = 0; i < len;i++){
                if(tmp[i] == ','){
                        count++;
                }
        }
        //count = byg_count(",", tmp);

        if(tmp[0] == 'B'){ // add extra space for all N barcode....
                count++;
        }

        if(tmp[0] == 'S'){ // add extra space for all N barcode....
                count++;
        }
        f = 0;

        //fprintf(stderr,"Segment %d: %d	X%sX sequences\n",segment,count+1,tmp);
        read_structure->numseq_in_segment[segment] = count+1;

        MMALLOC(read_structure->sequence_matrix[segment],sizeof(char*) * (count+1));
        for(i = 0; i < count+1;i++){
                read_structure->sequence_matrix[segment][i] = 0;
                MMALLOC(read_structure->sequence_matrix[segment][i],sizeof(char)* strlen(tmp));
        }

        read_structure->type[segment] = tmp[0];
        //read_structure->segment_length[segment] = strlen(tmp);

        if(tmp[0] == 'R'){
                read_structure->sequence_matrix[segment][0][0] = 'N';
                read_structure->sequence_matrix[segment][0][1] = 0;
                read_structure->segment_length[segment] = 1;
        }else if(tmp[0] == 'B'){
                f = 1;
                g = 0;
                len = (int)strlen(tmp) ;
//fprintf(stderr,"%d\n",len);
                for(i = 2;i <  len;i++){
                        //fprintf(stdout,"%d %d %c\n",f,g,tmp[i]);
                        if(tmp[i] != ','){
                                read_structure->sequence_matrix[segment][f][g] = tmp[i];
                                g++;
                        }else{
                                read_structure->sequence_matrix[segment][f][g] = 0;
                                f++;
                                g = 0;
                        }
                }
                read_structure->segment_length[segment] = g;
                read_structure->sequence_matrix[segment][f][g] = 0;
                f = 0;
                g = 0;

                for(i = 0; i <  read_structure->segment_length[segment];i++){
                        read_structure->sequence_matrix[segment][f][g] = 'N';
                        g++;
                }
                read_structure->segment_length[segment] = g;
                read_structure->sequence_matrix[segment][f][g] = 0;
                //f++;

        }else{
                f = 0;
                g = 0;
                len = (int)strlen(tmp) ;
                //fprintf(stderr,"%d\n",len);
                for(i = 2;i <  len;i++){
                        //fprintf(stdout,"%d %d %c\n",f,g,tmp[i]);
                        if(tmp[i] != ','){
                                read_structure->sequence_matrix[segment][f][g] = tmp[i];
                                g++;
                        }else{
                                read_structure->sequence_matrix[segment][f][g] = 0;
                                f++;
                                g = 0;
                        }
                }
                read_structure->segment_length[segment] = g;
                read_structure->sequence_matrix[segment][f][g] = 0;
        }
        //fprintf(stdout,"%s LEn:%d\n", tmp, read_structure->segment_length[segment]);
        /* if(tmp[0] == 'B'){ */
        /*         f = f + 1; */
        /*         g = 0; */
        /*         for(i = 0; i <  read_structure->segment_length[segment];i++){ */
        /*                 read_structure->sequence_matrix[segment][f][g] = 'N'; */
        /*                 g++; */
        /*         } */
        /*         read_structure->sequence_matrix[segment][f][g] = 0; */
        /* } */

        if(tmp[0] == 'S'){
                f = f + 1;
                g = 0;
                for(i = 0; i < read_structure->segment_length[segment];i++){
                        read_structure->sequence_matrix[segment][f][g] = 'N';
                        g++;
                }
                read_structure->sequence_matrix[segment][f][g] = 0;
        }


        if(segment+1 >read_structure->num_segments  ){
                read_structure->num_segments = segment+1;
        }
        return OK;
ERROR:
        if(read_structure->sequence_matrix[segment]){

                for(i = 0; i < count+1;i++){
                        MFREE(read_structure->sequence_matrix[segment][i]);//,sizeof(char)* strlen(tmp));
                }
                MFREE(read_structure->sequence_matrix[segment]);//,sizeof(char*) * (count+1));
        }

        return FAIL;
}

int alloc_arch_lib(struct arch_library** arch)
{
        struct arch_library* al = NULL;
        char* in[] = {
                "R:N"
        };

        int i;
        //ASSERT(num != 0, "Too few architectures requested");
        MMALLOC(al, sizeof(struct arch_library));
        al->arch_to_read_assignment = NULL;
        al->arch_posteriors = NULL;
        al->num_arch = 0;
        al->alloc_num_arch = 8;
        al->confidence_thresholds = NULL;
        al->read_structure = NULL;
        al->spec_line = NULL;
        MMALLOC(al->confidence_thresholds, sizeof(float) * al->alloc_num_arch);
        MMALLOC(al->spec_line, sizeof(char*) * al->alloc_num_arch);
        MMALLOC(al->read_structure, sizeof(struct read_structure*) * al->alloc_num_arch);
        for(i = 0; i < al->alloc_num_arch;i++){
                al->spec_line[i] = NULL;
                al->read_structure[i] = NULL;
                RUN(malloc_read_structure(&al->read_structure[i]));
        }

        /* load default model */
        RUN(read_arch_into_lib(al, in, 1));
        *arch = al;
        return OK;
ERROR:
        return FAIL;
}

int resize_arch_lib(struct arch_library* al)
{
        int old;
        int i;

        ASSERT(al != NULL, "No al lib");

        old = al->alloc_num_arch;

        al->alloc_num_arch = al->alloc_num_arch  << 1;
        MREALLOC(al->spec_line, sizeof(char*) * al->alloc_num_arch);
        MREALLOC(al->read_structure, sizeof(struct read_structure*) * al->alloc_num_arch);
        MREALLOC(al->confidence_thresholds, sizeof(float) * al->alloc_num_arch);
        for(i = old;i < al->alloc_num_arch;i++){
                al->spec_line[i] = NULL;
                al->read_structure[i] = NULL;
                RUN(malloc_read_structure(&al->read_structure[i]));
        }
        return OK;
ERROR:
        return FAIL;
}

void free_arch_lib(struct arch_library* al)
{
        int i;
        if(al){
                for(i = 0; i < al->alloc_num_arch;i++){
                        if(al->spec_line[i]){
                                MFREE(al->spec_line[i]);
                        }
                        free_read_structure(al->read_structure[i]);
                }
                if(al->arch_posteriors){
                        for(i = 0; i < al->num_arch ;i++){
                                MFREE(al->arch_posteriors[i]);
                        }
                        MFREE(al->arch_posteriors);
                }
                if(al->confidence_thresholds){
                        MFREE(al->confidence_thresholds);
                }
                if(al->arch_to_read_assignment){
                        MFREE(al->arch_to_read_assignment);
                }
                MFREE(al->spec_line);
                MFREE(al->read_structure);
                MFREE(al);
        }
}


int malloc_read_structure(struct read_structure** read_structure)
{
        struct read_structure* rs = 0;
        int i;
        MMALLOC(rs, sizeof(struct read_structure));
        rs->num_segments = 0;
        rs->alloc_num_seqments = 8;
        rs->segment_name = NULL;
        rs->sequence_matrix = NULL;
        rs->numseq_in_segment = NULL;
        rs->type = NULL;
        rs->segment_length = NULL;
        rs->segment_name = NULL;
        rs->extract = NULL;
        rs->max_len = NULL;
        rs->min_len = NULL;
        //rs->assignment_to_read = 0;
        MMALLOC(rs->sequence_matrix ,sizeof(char**) * rs->alloc_num_seqments );
        MMALLOC(rs->segment_name ,sizeof(char*) * rs->alloc_num_seqments );
        MMALLOC(rs->numseq_in_segment, sizeof(int) * rs->alloc_num_seqments);
        MMALLOC(rs->segment_length, sizeof(int) * rs->alloc_num_seqments);
        MMALLOC(rs->type ,sizeof(char) * rs->alloc_num_seqments );
        MMALLOC(rs->extract, sizeof(uint8_t) * rs->alloc_num_seqments);
        MMALLOC(rs->max_len, sizeof(int) * rs->alloc_num_seqments);
        MMALLOC(rs->min_len, sizeof(int) * rs->alloc_num_seqments);


        for(i = 0;i < rs->alloc_num_seqments;i++){
                rs->sequence_matrix[i] = NULL;
                rs->segment_name[i] = NULL;
                rs->numseq_in_segment[i] = 0;
                rs->segment_length[i] = 0;
                rs->type[i] = 0;
                rs->extract[i] = 0;
                rs->max_len[i] = 0;
                rs->min_len[i] = 0;
        }


        *read_structure = rs;
        return OK;
ERROR:
        free_read_structure(rs);
        return FAIL;
}

int resize_read_structure(struct read_structure* rs)
{
        int i;
        int old;
        ASSERT(rs!= NULL, "No read structure");
        old = rs->alloc_num_seqments;
        rs->alloc_num_seqments = rs->alloc_num_seqments + rs->alloc_num_seqments /2;
        MREALLOC(rs->sequence_matrix ,sizeof(char**) * rs->alloc_num_seqments );
        MREALLOC(rs->segment_name ,sizeof(char*) * rs->alloc_num_seqments );
        MREALLOC(rs->numseq_in_segment, sizeof(int) * rs->alloc_num_seqments);
        MREALLOC(rs->segment_length, sizeof(int) * rs->alloc_num_seqments);
        MREALLOC(rs->type ,sizeof(char) * rs->alloc_num_seqments );
        MREALLOC(rs->extract, sizeof(uint8_t) * rs->alloc_num_seqments);
        MREALLOC(rs->max_len, sizeof(int) * rs->alloc_num_seqments);
        MREALLOC(rs->min_len, sizeof(int) * rs->alloc_num_seqments);
        for(i = old;i < rs->alloc_num_seqments;i++){
                rs->sequence_matrix[i] = NULL;
                rs->segment_name[i] = NULL;
                rs->numseq_in_segment[i] = 0;
                rs->segment_length[i] = 0;
                rs->type[i] = 0;
                rs->extract[i] = 0;
                rs->max_len[i] = 0;
                rs->min_len[i] = 0;
        }
        return OK;
ERROR:
        free_read_structure(rs);
        return FAIL;
}

void free_read_structure(struct read_structure* read_structure)
{
        int i,j;
        if(read_structure){
                for(i = 0; i < read_structure->alloc_num_seqments;i++){
                        if(read_structure->sequence_matrix[i]){
                                for(j = 0; j < read_structure->numseq_in_segment[i];j++){
                                        MFREE(read_structure->sequence_matrix[i][j]);
                                }

                                MFREE(read_structure->sequence_matrix[i]);
                        }
                        if(read_structure->segment_name[i]){
                                MFREE(read_structure->segment_name[i]);
                        }
                }

                MFREE(read_structure->sequence_matrix);
                MFREE(read_structure->segment_name);
                MFREE(read_structure->segment_length);
                MFREE(read_structure->numseq_in_segment );
                MFREE(read_structure->type);
                MFREE(read_structure->extract);
                MFREE(read_structure->max_len);
                MFREE(read_structure->min_len);
                MFREE(read_structure);
        }
}


int QC_read_structure(struct read_structure* read_structure )
{
        int i,g,f;//min_error;//errors;
        int last = -1;
        //int num_pairs = 0;
        for(i = 0; i < read_structure->num_segments;i++){


                if(read_structure->sequence_matrix[i]){
                        if(last +1 != i){
                                ERROR_MSG("ERROR: a hmm building lock was skipped??\n");
                                //sprintf(param->buffer,"ERROR: a hmm building lock was skipped??\n");
                                //param->messages = append_message(param->messages, param->buffer);

                                //return FAIL;
                        }


                        //serious checking...
                        for(g = 0;g < read_structure->numseq_in_segment[i];g++){
                                for(f = g+1;f < read_structure->numseq_in_segment[i];f++){
                                        //if(read_structure->segment_length[i] != read_structure->segment_length[i]){
                                        //      ERROR_MSG("ERROR: the sequences in the same segment have to have the same length.\n");
                                        //}
                                }
                        }
                        last = i;
                }



                /*if(read_structure->type[i] == 'B'){
                        min_error = 1000;
                        for(g = 0;g <  read_structure->numseq_in_segment[i];g++){
p                                for(f = g+1;f < read_structure->numseq_in_segment[i];f++){
                                        errors = bpm(read_structure->sequence_matrix[i][g], read_structure->sequence_matrix[i][f], (int)strlen(read_structure->sequence_matrix[i][0]),(int)strlen(read_structure->sequence_matrix[i][0]));

                                        if(errors < min_error){
                                                min_error = errors;
                                                num_pairs = 1;
                                        }else if(errors == min_error ){
                                                num_pairs++;
                                        }
                                }
                        }
                        }*/
        }
        return OK;
ERROR:
        return FAIL;
}

int sanity_check_arch_lib(struct arch_library* al)
{
        int i,j,c;
        ASSERT(al != NULL, "No arch library");
        //LOG_MSG("Check sanity");
        for(i = 0 ;i < al->num_arch;i++){
                if(al->spec_line[i] != NULL){
                        for(j = i+1;j < al->num_arch;j++){
                                if(al->spec_line[j] != NULL){
                                        //fprintf(stdout,"%d X%sX\n%d X%sX\n",i,al->spec_line[i],j, al->spec_line[j]);
                                        if(!strcmp(al->spec_line[i], al->spec_line[j])){
                                                WARNING_MSG("two architectures are the same:\n%s\n%s\n", al->spec_line[i],al->spec_line[j]);
                                                MFREE(al->spec_line[j]);
                                                al->spec_line[j] = NULL;
                                                free_read_structure(al->read_structure[j]);
                                                al->read_structure[j] = NULL;
                                                RUN(malloc_read_structure(&al->read_structure[j]));
                                        }
                                }
                        }
                }
        }
        c = 0;
        for(i = 0 ;i < al->num_arch;i++){
                if(al->spec_line[i]){
                        //fprintf(stdout,"%d %d %s\n",i,c, al->spec_line[i]);
                        al->spec_line[c] = al->spec_line[i];
                        al->read_structure[c] = al->read_structure[i];
                        c++;
                }
        }
        al->num_arch = c;

        //LOG_MSG("Check sanity DONE");

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
        spec->seq = NULL;
        spec->name = NULL;

        MMALLOC(spec->name, sizeof(char) * BUFFER_LEN);

        *s = spec;
        return OK;
ERROR:

        return FAIL;
}

void free_segment_spec(struct segment_specs*s)
{
        int i;
        if(s){
                if(s->seq){
                        for(i = 0; i < s->num_seq;i++){
                                MFREE(s->seq[i]);
                        }
                        MFREE(s->seq);
                }
                if(s->name){
                        MFREE(s->name);
                }
                MFREE(s);
        }
}



#ifdef ARCH_TEST

int test_token_parsing(void);
int main(int argc, char *argv[])
{
        struct arch_library* al = NULL;
        char* in[] = {
                "O:N",
                "B:GTA,AAC",
                "R:N",
                "L:CCTTAA",
                "B:ACAGTG,ACTTGA,TTAGGC"
        };
        char* buffer = NULL;
        int size;
        int i;



        MMALLOC(buffer, sizeof(char) * BUFFER_LEN);
        RUN(alloc_arch_lib(&al));

        size= sizeof(in) / sizeof(char*);
        //LOG_MSG("%d",size);
        RUN(read_architecture_files(al, "../dev/casava_arch.txt"));
        read_arch_into_lib(al, in, size);
        read_arch_into_lib(al, in, size);
        read_arch_into_lib(al, in, size);
        read_arch_into_lib(al, in, size);

        RUN(read_architecture_files(al, "../dev/casava_arch.txt"));
        LOG_MSG("Read in %d architectures.",al->num_arch);
        for(i = 0; i < al->num_arch;i++){
                fprintf(stdout,"%d %s\n",i, al->spec_line[i]);
        }
        free_arch_lib(al);
        al = NULL;

        MFREE(buffer);

        RUN(test_token_parsing());
        return EXIT_SUCCESS;
ERROR:
        if(al){
                free_arch_lib(al);
        }
        if(buffer){
                MFREE(buffer);
        }
        return EXIT_FAILURE;
}

int test_token_parsing(void)
{
        char buffer[256];
        int status;
        int i;
        char types[4] = "EASI";
        struct segment_specs* s;
        LOG_MSG("Test format with name:");
        for(i = 0; i < 4;i++){
                snprintf(buffer, 256, "OTTO:%c:AA",types[i]);

                status = parse_rs_token_message(buffer,&s);
                if(status != OK){
                        ERROR_MSG("Parsing of %s failed", buffer);
                }
                free_segment_spec(s);
        }
        LOG_MSG("Test format without name:");
        for(i = 0; i < 4;i++){
                snprintf(buffer, 256, "%c:AA",types[i]);

                status = parse_rs_token_message(buffer,&s);
                if(status != OK){
                        ERROR_MSG("Parsing of %s failed", buffer);
                }
                free_segment_spec(s);
        }

        LOG_MSG("Test other valid formats:");
        snprintf(buffer, 256, "OTTO:A:AA,CC");
        status = parse_rs_token_message(buffer,&s);
        if(status != OK){
                ERROR_MSG("Parsing of %s failed", buffer);
        }
        free_segment_spec(s);

        snprintf(buffer, 256, "OTTO:A:AA,CC,CG");
        status = parse_rs_token_message(buffer,&s);
        if(status != OK){
                ERROR_MSG("Parsing of %s failed", buffer);
        }
        free_segment_spec(s);

        snprintf(buffer, 256, "OTTO:A:AA,CC,CG,TT");
        status = parse_rs_token_message(buffer,&s);
        if(status != OK){
                ERROR_MSG("Parsing of %s failed", buffer);
        }
        free_segment_spec(s);

        snprintf(buffer, 256, "OTTO:A:N{4}");
        status = parse_rs_token_message(buffer,&s);
        if(status != OK){
                ERROR_MSG("Parsing of %s failed", buffer);
        }
        free_segment_spec(s);

        snprintf(buffer, 256, "O:A:N{4,8}");
        status = parse_rs_token_message(buffer,&s);
        if(status != OK){
                ERROR_MSG("Parsing of %s failed", buffer);
        }
        free_segment_spec(s);


        /* negative */
        LOG_MSG("Test invalid formats:");
        snprintf(buffer, 256, "O:AN{4,8}");
        status = parse_rs_token_message(buffer,&s);
        if(status != FAIL){
                ERROR_MSG("Parsing of %s succeeded (but it should not have)", buffer);
                free_segment_spec(s);
        }
        //free_segment_spec(s);

        snprintf(buffer, 256, "O:A:N{4,8},AAA");
        status = parse_rs_token_message(buffer,&s);
        if(status != FAIL){
                ERROR_MSG("Parsing of %s succeeded (but it should not have)", buffer);
                free_segment_spec(s);
        }
        //free_segment_spec(s);

        snprintf(buffer, 256, "O:A:");
        status = parse_rs_token_message(buffer,&s);
        if(status != FAIL){
                ERROR_MSG("Parsing of %s succeeded (but it should not have)", buffer);
                free_segment_spec(s);
        }
        //free_segment_spec(s);

        snprintf(buffer, 256, "O:A:A,AA,AAA");
        status = parse_rs_token_message(buffer,&s);
        if(status != FAIL){
                ERROR_MSG("Parsing of %s succeeded (but it should not have)", buffer);
                free_segment_spec(s);
        }
        //free_segment_spec(s);

        snprintf(buffer, 256, ":A:A,AA,AAA");
        status = parse_rs_token_message(buffer,&s);
        if(status != FAIL){
                ERROR_MSG("Parsing of %s succeeded (but it should not have)", buffer);
                free_segment_spec(s);
        }
        //free_segment_spec(s);

        snprintf(buffer, 256, "GGGG:AA,AA,AAA");
        status = parse_rs_token_message(buffer,&s);
        if(status != FAIL){
                ERROR_MSG("Parsing of %s succeeded (but it should not have)", buffer);
                free_segment_spec(s);
        }


        snprintf(buffer, 256, "NA:AA,AA,AAA");
        status = parse_rs_token_message(buffer,&s);
        if(status != FAIL){
                ERROR_MSG("Parsing of %s succeeded (but it should not have)", buffer);
                free_segment_spec(s);
        }


        snprintf(buffer, 256, "NA:E:AA");
        status = parse_rs_token_message(buffer,&s);
        if(status != FAIL){
                ERROR_MSG("Parsing of %s succeeded (but it should not have)", buffer);
                free_segment_spec(s);
        }

        //free_segment_spec(s);

        return OK;
ERROR:
        if(s){
                free_segment_spec(s);
        }
        return FAIL;
}


#endif
