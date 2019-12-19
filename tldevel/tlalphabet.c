#include "tldevel.h"

#define TLALPHABET_IMPORT
#include "tlalphabet.h"

static int create_default_protein(struct alphabet* a);
static int create_default_DNA(struct alphabet* a);
static int create_reduced_protein(struct alphabet* a);

static int clean_and_set_to_extern(struct alphabet* a);
static int merge_codes(struct alphabet*a,const int X, const int Y);


int create_alphabet(struct alphabet** alphabet, int type)
{
        struct alphabet* a = NULL;
        int i;
        MMALLOC(a, sizeof(struct alphabet));

        for(i = 0; i < 128;i++){
                a->to_internal[i] = -1;
        }

        for(i = 0; i < 32;i++){
                a->to_external[i] = -1;
        }
        a->type = type;

        switch (type) {
        case TLALPHABET_DEFAULT_PROTEIN : {
                create_default_protein(a);
                break;
        }
        case TLALPHABET_DEFAULT_DNA : {
                create_default_DNA(a);
                break;
        }
        case TLALPHABET_REDUCED_PROTEIN : {
                create_reduced_protein(a);
                break;
        }
        default:
                ERROR_MSG("Alphabet type %d not recognised", type);
                break;
        }

        RUN(clean_and_set_to_extern(a));
        *alphabet = a;
        return OK;
ERROR:
        if(a){
                MFREE(a);
        }
        return FAIL;
}

int create_default_protein(struct alphabet* a)
{
        char aacode[20] = "ACDEFGHIKLMNPQRSTVWY";

        int code;
        int i;
        code = 0;
        for(i = 0; i < 20;i++){
                //fprintf(stdout,"%c %d CODE: %d\n", aacode[i], (int) aacode[i], code);
                a->to_internal[(int) aacode[i]] = code;

                code++;
        }
        /* ambiguity codes  */
        /* BZX  */

        a->to_internal[(int) 'B'] = code;

        a->to_internal[(int) 'Z'] = code;
        a->to_internal[(int) 'X'] = code;

        /* Some protein sequences contain 'U' - a non-IUPAC code
           I will treat these as an ambiguous aa
           e.g:
           >Q74EN2_GEOSL/108-206
           TRELEALVAKGTEEGGYLLIDSRPAGKYNEAHIPTAVSIPFAELEKNPALLTASKDRLLVFYCGGVTUVLSPKSAGLAKKSGYEKVRVYLDGEPEWKKA

        */
        a->to_internal[(int) 'U'] = code;
        code++;
        return OK;
}

int create_default_DNA(struct alphabet* a)
{

        char dnacode[16] = "ACGTUNRYSWKMBDHV";

        int code;
        int i;
        code = 0;
        for(i = 0; i < 16;i++){
                //fprintf(stdout,"%c %d CODE: %d\n", aacode[i], (int) aacode[i], code);
                a->to_internal[(int) dnacode[i]] = code;

                code++;
        }

        merge_codes(a,'U','T');

        /* R.................A or G */
        /* Y.................C or T */
        /* S.................G or C */
        /* W.................A or T */
        /* K.................G or T */
        /* M.................A or C */
        /* B.................C or G or T */
        /* D.................A or G or T */
        /* H.................A or C or T */
        /* V.................A or C or G */
        merge_codes(a,'N','R');
        merge_codes(a,'N','Y');
        merge_codes(a,'N','S');
        merge_codes(a,'N','W');
        merge_codes(a,'N','K');
        merge_codes(a,'N','M');
        merge_codes(a,'N','B');
        merge_codes(a,'N','D');
        merge_codes(a,'N','H');
        merge_codes(a,'N','V');

        return OK;
}


int create_reduced_protein(struct alphabet* a)
{
        char aacode[20] = "ACDEFGHIKLMNPQRSTVWY";

        int code;
        int i;
        code = 0;
        for(i = 0; i < 20;i++){
                a->to_internal[(int) aacode[i]] = code;
                code++;
        }
        /* ambiguity codes  */
        /* BZX  */

        a->to_internal[(int) 'B'] = code;
        code++;

        a->to_internal[(int) 'Z'] = code;
        code++;

        a->to_internal[(int) 'X'] = code;
        code++;
        /* From  Clustering huge protein sequence sets in linear time
Martin Steinegger 1, 2, 3 and Johannes Söding 1 */
        /* The default alphabet with A = 13 merges (L,M), (I,V), (K,R), (E, Q), (A,S,T), (N, D) and (F,Y).*/

        /* reduced codes */
        merge_codes(a,'L','M');
        merge_codes(a,'I','V');

        merge_codes(a,'K','R');
        merge_codes(a,'E','Q');

        merge_codes(a,'A','S');
        merge_codes(a,'A','T');
        merge_codes(a,'S','T');
        merge_codes(a,'N','D');

        merge_codes(a,'F','Y');

        /* merge ambiguity codes */
        merge_codes(a,'B','N');
        merge_codes(a,'B','D');

        merge_codes(a,'Z','E');
        merge_codes(a,'Z','Q');


        return OK;

}

int merge_codes(struct alphabet*a,const int X, const int Y)
{
        int min;

        min = MACRO_MIN(a->to_internal[X],a->to_internal[Y]);

        ASSERT(min != -1, "code not set!");

        a->to_internal[X] = min;
        a->to_internal[Y] = min;
        return OK;
ERROR:
        return FAIL;
}

int clean_and_set_to_extern(struct alphabet* a)
{
        int i;
        uint8_t code = 0;
        int8_t trans[32];
        for(i = 0; i < 32;i++){
                trans[i] = -1;

        }

        for(i = 64; i < 96;i++){
                if(a->to_internal[i] != -1){
                        trans[a->to_internal[i]] = 1;
                }
        }
        code = 0;
        for(i = 0; i < 32;i++){
                if(trans[i] == 1){
                        trans[i] = code;
                        code++;
                }
        }
        a->L = code;
        for(i = 64; i < 96;i++){
                if(a->to_internal[i] != -1){
                        a->to_internal[i] = trans[a->to_internal[i]];//a->to_internal[i]];
                        a->to_internal[i+32] = a->to_internal[i];

                }

        }

        for(i = 64;i < 96;i++){
                if(a->to_internal[i] != -1){
                        a->to_external[a->to_internal[i]] = i;
                }
        }
        return OK;
}
