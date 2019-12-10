#include "assign_data.h"



#include "tldevel.h"


int set_up_assign_structure(struct arch_library* al,struct assign_struct* as)
{
        struct read_structure* read_structure = NULL;
        int i,j;
        char c;
        ASSERT(al != NULL,"No archlib");

        ASSERT(as != NULL,"No assign struct ");

        as->num_bits = 0;
        for(i = 0; i < al->num_file;i++){
                read_structure = al->read_structure[al->arch_to_read_assignment[i]];
                //fprintf(stdout,"Read %d: ",i);
                for(j = 0; j < read_structure->num_segments;j++){
                        c = read_structure->type[j];
                        switch (c) {
                        case 'B':
                        case 'R':
                        case 'F':
                                as->num_bits++;
                                break;
                        default:
                                break;
                        }
                }
                //fprintf(stdout,"\n");
        }
        /* create offsets  */
        return OK;
ERROR:
        return FAIL;
}


int init_assign_structure(struct assign_struct** assign,struct arch_library* al, int total)
{
        struct assign_struct* as = NULL;
        //ASSERT(num_files >= 1,"no infiles");
        int i,j;
        MMALLOC(as, sizeof(struct assign_struct));
        as->num_files = al->num_file;
        as->bits = NULL;
        as->num_bits = 0;


        RUN(set_up_assign_structure(al,as));

        as->total = total;

        as->bits = NULL;
        MMALLOC(as->bits, sizeof(struct seq_bit_vec*)* as->total);
        //MMALLOC(as->active_bits,sizeof(uint8_t) * as->total));
        for(i = 0; i < as->total;i++){
                as->bits[i] = NULL;
                MMALLOC(as->bits[i], sizeof(struct seq_bit_vec));
                as->bits[i]->num_bit = 0;
                as->bits[i]->bits = NULL;
                as->bits[i]->Q = NULL;
                as->bits[i]->bc = NULL;
                as->bits[i]->umi = NULL;
                MMALLOC(as->bits[i]->Q,sizeof(float) * as->num_files);
                as->bits[i]->pass = 1;
                MMALLOC(as->bits[i]->bits , sizeof(struct seq_bit) * as->num_bits);
                for(j = 0; j < as->num_bits;j++){
                        as->bits[i]->bits[j] = NULL;
                        MMALLOC(as->bits[i]->bits[j], sizeof(struct seq_bit));
                }
        }
        *assign = as;
        return OK;
ERROR:
        free_assign_structure(as);
        return FAIL;
}

void free_assign_structure(struct assign_struct* as)
{
        if(as){
                int i,j;
                if(as->bits){
                        //MMALLOC(as->active_bits,sizeof(uint8_t) * as->total));
                        for(i = 0; i < as->total;i++){
                                for(j = 0; j < as->num_bits;j++){
                                        MFREE(as->bits[i]->bits[j]);
                                }
                                MFREE(as->bits[i]->Q);
                                MFREE(as->bits[i]->bits);
                                MFREE(as->bits[i]);
                        }
                        MFREE(as->bits);
                }

                MFREE(as);
        }
}
