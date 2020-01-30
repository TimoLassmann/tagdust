#include "read_groups.h"

#include <string.h>
#include <ctype.h>
#include "tldevel.h"

#include "tlseqio.h"
/* These are functions to figure out which input files belong together...  */


#define TOP_READS 1000000


struct rg_f_sort{
        struct sequence_stats_info* si;
        char* name;
        int index;
};

static int alloc_read_groups(struct read_groups** read_group, int infiles);
static int compare_read_names(char* name1, char* name2);
static int print_read_groups(struct read_groups* read_groups);

static int qsort_rg_f(const void *a, const void *b);




int generate_read_groups(struct read_groups** rg, char** infile, int num_infiles)
{

        /* Logic:
           check read names of input -> group together;
           if impossible; check via hashing whether read names are shared between files -> unsorted
           else assume all reads are different lanes and should be put together.
        */
        struct read_groups* read_groups = NULL;
        struct read_ensembl* e = NULL;

        struct tl_seq_buffer**  rb = NULL;
        struct file_handler** f_hand = NULL;

        int i,j,c;
        int r;
        int** dist_mat = NULL;

        uint8_t* active_seq = NULL;
        MMALLOC(rb, sizeof(struct tl_seq_buffer*) * num_infiles);
        MMALLOC(f_hand, sizeof(struct file_handler*) * num_infiles);

        RUN(alloc_read_groups(&read_groups, num_infiles));


        RUN(galloc(&dist_mat,num_infiles,num_infiles));

        RUN(galloc(&active_seq,num_infiles));
        for(i = 0; i < num_infiles;i++){
                active_seq[i] =1;
                for(j = 0; j < num_infiles;j++){
                        dist_mat[i][j] = 0;
                }
        }
        for(i = 0; i < num_infiles;i++){
                f_hand[i] = NULL;
                RUN(open_fasta_fastq_file(&f_hand[i], infile[i], TLSEQIO_READ));
                //RUN(io_handler(&f_hand[i], infile[i]));
        }
        for(i = 0; i < num_infiles;i++){
                rb[i] = NULL;
                //LOG_MSG("Reading : %s", infile[i]);
                RUN(read_fasta_fastq_file(f_hand[i], &rb[i], TOP_READS));


        }
        for(i = 0; i < num_infiles;i++){
                for(j = i+1;j < num_infiles;j++){
                        if(rb[i]->num_seq == rb[j]->num_seq){
                                for(c = 0; c < rb[i]->num_seq;c++){
                                        r = compare_read_names(rb[i]->sequences[c]->name, rb[j]->sequences[c]->name);
                                        if(r == OK){
                                                dist_mat[i][j]++;
                                                dist_mat[j][i]++;
                                        }
                                }
                        }
                }
        }

        for(i = 0; i < num_infiles;i++){
                if(active_seq[i]){
                        e = read_groups->e[read_groups->num_groups];
                        e->filenames[e->num_files] = infile[i];
                        e->num_files++;
                        //fprintf(stdout,"Start clu: %d (%s)\n", i, infile[i]);
                        active_seq[i] = 0;
                        for(j = i+1;j < num_infiles;j++){
                                if(active_seq[j]){
                                        if(dist_mat[i][j] == rb[i]->num_seq){
                                                e->filenames[e->num_files] = infile[j];
                                                e->num_files++;

                                                //fprintf(stdout,"Adding  %d (%s)\n", j, infile[j]);
                                                active_seq[j] = 0;
                                        }
                                }
                        }
                        //fprintf(stdout,"DONE\n\n");
                        read_groups->num_groups++;
                }
        }

        RUN(print_read_groups(read_groups));

        for(i = 0; i < read_groups->num_groups;i++){
                for(j = i+1; j < read_groups->num_groups;j++){
                        if(read_groups->e[i]->num_files != read_groups->e[j]->num_files){
                                WARNING_MSG("Problem with the input files:");
                                RUN(print_read_groups(read_groups));
                                fprintf(stdout,"\n");
                                fprintf(stdout,"Two (or more) groups of reads have a different number of associated files.\n");
                                fprintf(stdout,"You probably want to run tagcook separately on each group.\n\n");
                                break;
                        }
                }
        }

        *rg = read_groups;

        for(i = 0; i < num_infiles;i++){
                if(rb[i]){
                        free_tl_seq_buffer(rb[i]);
                }
        }
        MFREE(rb);

        for(i = 0; i < num_infiles;i++){
                RUN(close_seq_file(&f_hand[i]));
        }
        MFREE(f_hand);
        gfree(dist_mat);
        gfree(active_seq);


        return OK;
ERROR:
        free_read_groups(read_groups);
        return FAIL;
}

int print_read_groups(struct read_groups* read_groups)
{
        int i,j;
        struct read_ensembl* e = NULL;
        for(i = 0; i < read_groups->num_groups;i++){
                LOG_MSG("Group %d:", i);
                e = read_groups->e[i];

                for(j = 0; j < e->num_files;j++){
                        LOG_MSG("   %s", e->filenames[j]);
                }
                LOG_MSG("");
        }
        return OK;
}

int alloc_read_groups(struct read_groups** read_group, int infiles)
{
        struct read_groups* rg = NULL;
        int i,j;
        MMALLOC(rg, sizeof(struct read_groups));
        rg->num_groups = 0;
        rg->alloc_num_groups = infiles;
        rg->e = NULL;

        MMALLOC(rg->e, sizeof(struct read_ensembl*) * rg->alloc_num_groups);
        for(i = 0; i < rg->alloc_num_groups;i++){
                rg->e[i] = NULL;
                MMALLOC(rg->e[i], sizeof(struct read_ensembl));
                rg->e[i]->num_files = 0;
                rg->e[i]->filenames = NULL;
                rg->e[i]->arch_to_read_assignment = NULL;
                rg->e[i]->si = NULL;
                rg->e[i]->ssi = NULL;
                MMALLOC(rg->e[i]->ssi, sizeof(struct sequence_stats_info*) * infiles);
                MMALLOC(rg->e[i]->filenames, sizeof(char*) * infiles);
                MMALLOC(rg->e[i]->arch_to_read_assignment, sizeof(int) * infiles);

                for(j = 0; j < infiles;j++){
                        rg->e[i]->filenames[j] = NULL;
                        rg->e[i]->arch_to_read_assignment[j] = -1;
                }
        }
        *read_group = rg;
        return OK;
ERROR:
        free_read_groups(rg);
        *read_group = NULL;
        return FAIL;
}

void free_read_groups(struct read_groups* rg)
{
        if(rg){
                int i;
                for(i = 0; i < rg->alloc_num_groups;i++){
                        if(rg->e[i]->ssi){
                                MFREE(rg->e[i]->ssi);
                        }
                        //free_sequence_stats(rg->e[i]->si);
                        MFREE(rg->e[i]->filenames);
                        MFREE(rg->e[i]->arch_to_read_assignment);
                        MFREE(rg->e[i]);
                }
                MFREE(rg->e);
                MFREE(rg);
        }
}


int compare_read_names(char* name1, char* name2)
{

#ifdef READ_GROUP_ITEST
        int detected = -1;
#else
        static int detected = -1;
#endif
        char instrument_R1[100];
        int run_id_R1 = 0;
        char flowcell_R1[100];
        int flowcell_lane_R1= 0;
        int tile_number_R1= 0;
        int x_coordinate_R1= 0;
        int y_coordinate_R1= 0;

        char instrument_R2[100];
        int run_id_R2= 0;
        char flowcell_R2[100];
        int flowcell_lane_R2= 0;
        int tile_number_R2= 0;
        int x_coordinate_R2= 0;
        int y_coordinate_R2= 0;

        int i;

        int number_of_values_found = 0;

        instrument_R1[0] = 0;
        instrument_R2[0] = 0;

        flowcell_R1[0] = 0;
        flowcell_R2[0] = 0;



        if(detected == -1){
                //option 1: casava 1.8
                // name should look like this:@EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG
                number_of_values_found =sscanf(name1,"%[^:]:%d:%[^:]:%d:%d:%d:%d ", instrument_R1,&run_id_R1,flowcell_R1,&flowcell_lane_R1,&tile_number_R1,&x_coordinate_R1,&y_coordinate_R1 );
                if(number_of_values_found == 7){
                        detected = 1;
#ifdef DEBUG
//                        LOG_MSG("Detected casava 1.8 format.\n");
#endif
                        //param->messages = append_message(param->messages, param->buffer);
                }
                //fprintf(stderr,"casava 1.8?:%d %s\n",number_of_values_found, name1);
                //fprintf(stderr,"%s\n%d\n%s\n%d\n%d\n%d\n%d\n", instrument_R1,run_id_R1,flowcell_R1,flowcell_lane_R1,tile_number_R1,x_coordinate_R1,y_coordinate_R1);

        }

        if(detected == -1){
                //option 2: casava 1.7
                // name should look like this:@HWUSI-EAS100R:6:73:941:1973#0/1
                //HWUSI-EAS747_0040_FC64GRTAAXX:8:1:3268:1065#0/1
                number_of_values_found =sscanf(name1,"%[^:]:%d:%d:%d:%d", instrument_R1,&flowcell_lane_R1,&tile_number_R1,&x_coordinate_R1,&y_coordinate_R1);
                //fprintf(stderr,"casava 1.7?:%d %s\n", number_of_values_found,name1);
                //fprintf(stderr,"%s\n%d\n%d\n%d\n%d\n", instrument_R1,flowcell_lane_R1,tile_number_R1,x_coordinate_R1,y_coordinate_R1);

                if(number_of_values_found == 5){
                        detected = 2;
                        #ifdef DEBUG
                        //LOG_MSG("Detected casava <1.7 format.\n");
                        #endif
                        //param->messages = append_message(param->messages, param->buffer);
                }
        }

        if(detected == -1){
                detected = 1000;
        }

        if(detected == 1){
                number_of_values_found =sscanf(name1,"%[^:]:%d:%[^:]:%d:%d:%d:%d ", instrument_R1,&run_id_R1,flowcell_R1,&flowcell_lane_R1,&tile_number_R1,&x_coordinate_R1,&y_coordinate_R1 );
                if(number_of_values_found != 7){
                        #ifdef DEBUG
                        //WARNING_MSG("File name %s\n does not match detected casava 1.8 format.\n",name1);
                        #endif
                        return FAIL;
                }

                number_of_values_found =sscanf(name2,"%[^:]:%d:%[^:]:%d:%d:%d:%d ", instrument_R2,&run_id_R2,flowcell_R2,&flowcell_lane_R2,&tile_number_R2,&x_coordinate_R2,&y_coordinate_R2 );
                if(number_of_values_found != 7){
                        #ifdef DEBUG
                        //WARNING_MSG("File name %s\n does not match detected casava 1.8 format.\n",name2);
                        #endif
#ifdef READ_GROUP_ITEST
                        fprintf(stderr,"@%s:%d:%d:%d:%d\n", instrument_R1,flowcell_lane_R1,tile_number_R1,x_coordinate_R1,y_coordinate_R1);
                        fprintf(stderr,"@%s:%d:%d:%d:%d\n", instrument_R2,flowcell_lane_R2,tile_number_R2,x_coordinate_R2,y_coordinate_R2);
#endif
                        return FAIL;
                }

                if(y_coordinate_R1 != y_coordinate_R2){
                        #ifdef DEBUG
                        //WARNING_MSG("Files seem to contain reads in different order:\n%s\n%s\n",name1, name2);
                        #endif
#ifdef READ_GROUP_ITEST
                        fprintf(stderr,"@%s:%d:%d:%d:%d\n", instrument_R1,flowcell_lane_R1,tile_number_R1,x_coordinate_R1,y_coordinate_R1);
                        fprintf(stderr,"@%s:%d:%d:%d:%d\n", instrument_R2,flowcell_lane_R2,tile_number_R2,x_coordinate_R2,y_coordinate_R2);
#endif
                        return FAIL;
                }

                if(x_coordinate_R1 != x_coordinate_R2){
                        #ifdef DEBUG
                        //WARNING_MSG("Files seem to contain reads in different order:\n%s\n%s\n",name1, name2);
                        #endif
#ifdef READ_GROUP_ITEST
                        fprintf(stderr,"@%s:%d:%d:%d:%d\n", instrument_R1,flowcell_lane_R1,tile_number_R1,x_coordinate_R1,y_coordinate_R1);
                        fprintf(stderr,"@%s:%d:%d:%d:%d\n", instrument_R2,flowcell_lane_R2,tile_number_R2,x_coordinate_R2,y_coordinate_R2);
#endif
                        return FAIL;
                }

                if(tile_number_R1 !=  tile_number_R2){
                        #ifdef DEBUG
                        //WARNING_MSG("Files seem to contain reads in different order:\n%s\n%s\n",name1, name2);
                        #endif
#ifdef READ_GROUP_ITEST
                        fprintf(stderr,"@%s:%d:%d:%d:%d\n", instrument_R1,flowcell_lane_R1,tile_number_R1,x_coordinate_R1,y_coordinate_R1);
                        fprintf(stderr,"@%s:%d:%d:%d:%d\n", instrument_R2,flowcell_lane_R2,tile_number_R2,x_coordinate_R2,y_coordinate_R2);
#endif
                        return FAIL;
                }

                if(flowcell_lane_R1 !=  flowcell_lane_R2){
                        #ifdef DEBUG
                        //WARNING_MSG("Files seem to contain reads in different order:\n%s\n%s\n",name1, name2);
                        #endif
#ifdef READ_GROUP_ITEST
                        fprintf(stderr,"@%s:%d:%d:%d:%d\n", instrument_R1,flowcell_lane_R1,tile_number_R1,x_coordinate_R1,y_coordinate_R1);
                        fprintf(stderr,"@%s:%d:%d:%d:%d\n", instrument_R2,flowcell_lane_R2,tile_number_R2,x_coordinate_R2,y_coordinate_R2);
#endif
                        return FAIL;
                }

                if(strcmp(flowcell_R1,flowcell_R2)){
                        #ifdef DEBUG
                        //WARNING_MSG("Files seem to contain reads in different order:\n%s\n%s\n",name1, name2);
                        #endif
#ifdef READ_GROUP_ITEST
                        fprintf(stderr,"@%s:%d:%d:%d:%d\n", instrument_R1,flowcell_lane_R1,tile_number_R1,x_coordinate_R1,y_coordinate_R1);
                        fprintf(stderr,"@%s:%d:%d:%d:%d\n", instrument_R2,flowcell_lane_R2,tile_number_R2,x_coordinate_R2,y_coordinate_R2);
#endif
                        return FAIL;
                }

                if(run_id_R1 !=  run_id_R2){
                        #ifdef DEBUG
                        //WARNING_MSG("Files seem to contain reads in different order:\n%s\n%s\n",name1, name2);
                        #endif
#ifdef READ_GROUP_ITEST
                        fprintf(stderr,"@%s:%d:%d:%d:%d\n", instrument_R1,flowcell_lane_R1,tile_number_R1,x_coordinate_R1,y_coordinate_R1);
                        fprintf(stderr,"@%s:%d:%d:%d:%d\n", instrument_R2,flowcell_lane_R2,tile_number_R2,x_coordinate_R2,y_coordinate_R2);
#endif
                        return FAIL;
                }

                if(strcmp(instrument_R1,instrument_R2)){

                        #ifdef DEBUG
                        //WARNING_MSG("Files seem to contain reads in different order:\n%s\n%s\n",name1, name2);
                        #endif
#ifdef READ_GROUP_ITEST
                        fprintf(stderr,"@%s:%d:%d:%d:%d\n", instrument_R1,flowcell_lane_R1,tile_number_R1,x_coordinate_R1,y_coordinate_R1);
                        fprintf(stderr,"@%s:%d:%d:%d:%d\n", instrument_R2,flowcell_lane_R2,tile_number_R2,x_coordinate_R2,y_coordinate_R2);
#endif
                        return FAIL;
                }
        }

        if(detected == 2){
                number_of_values_found =sscanf(name1,"%[^:]:%d:%d:%d:%d", instrument_R1,&flowcell_lane_R1,&tile_number_R1,&x_coordinate_R1,&y_coordinate_R1);

                if(number_of_values_found != 5){
                        #ifdef DEBUG
                        //WARNING_MSG("File name %s\n does not match detected casava <1.8 format.\n",name1);
                        #endif
#ifdef READ_GROUP_ITEST
                        fprintf(stderr,"@%s:%d:%d:%d:%d\n", instrument_R1,flowcell_lane_R1,tile_number_R1,x_coordinate_R1,y_coordinate_R1);
                        fprintf(stderr,"@%s:%d:%d:%d:%d\n", instrument_R2,flowcell_lane_R2,tile_number_R2,x_coordinate_R2,y_coordinate_R2);
#endif
                        return FAIL;
                }
                number_of_values_found =sscanf(name2,"%[^:]:%d:%d:%d:%d", instrument_R2,&flowcell_lane_R2,&tile_number_R2,&x_coordinate_R2,&y_coordinate_R2);

                if(number_of_values_found != 5){
                        #ifdef DEBUG
                        //WARNING_MSG("File name %s\n does not match detected casava <1.8 format.\n",name2);
                        #endif
#ifdef READ_GROUP_ITEST
                        fprintf(stderr,"@%s:%d:%d:%d:%d\n", instrument_R1,flowcell_lane_R1,tile_number_R1,x_coordinate_R1,y_coordinate_R1);
                        fprintf(stderr,"@%s:%d:%d:%d:%d\n", instrument_R2,flowcell_lane_R2,tile_number_R2,x_coordinate_R2,y_coordinate_R2);
#endif
                        return FAIL;
                }

                if(y_coordinate_R1 != y_coordinate_R2){
                        #ifdef DEBUG
                        //WARNING_MSG("Files seem to contain reads in different order:\n%s\n%s\n",name1, name2);
                        #endif
#ifdef READ_GROUP_ITEST
                        fprintf(stderr,"@%s:%d:%d:%d:%d\n", instrument_R1,flowcell_lane_R1,tile_number_R1,x_coordinate_R1,y_coordinate_R1);
                        fprintf(stderr,"@%s:%d:%d:%d:%d\n", instrument_R2,flowcell_lane_R2,tile_number_R2,x_coordinate_R2,y_coordinate_R2);
#endif
                        return FAIL;
                }

                if(x_coordinate_R1 != x_coordinate_R2){
                        #ifdef DEBUG
                        //WARNING_MSG("Files seem to contain reads in different order:\n%s\n%s\n",name1, name2);
                        #endif

#ifdef READ_GROUP_ITEST
                        fprintf(stderr,"@%s:%d:%d:%d:%d\n", instrument_R1,flowcell_lane_R1,tile_number_R1,x_coordinate_R1,y_coordinate_R1);
                        fprintf(stderr,"@%s:%d:%d:%d:%d\n", instrument_R2,flowcell_lane_R2,tile_number_R2,x_coordinate_R2,y_coordinate_R2);
#endif
                        return FAIL;
                }

                if(tile_number_R1 !=  tile_number_R2){
                        #ifdef DEBUG
                        //WARNING_MSG("Files seem to contain reads in different order:\n%s\n%s\n",name1, name2);
                        #endif
#ifdef READ_GROUP_ITEST
                        fprintf(stderr,"@%s:%d:%d:%d:%d\n", instrument_R1,flowcell_lane_R1,tile_number_R1,x_coordinate_R1,y_coordinate_R1);
                        fprintf(stderr,"@%s:%d:%d:%d:%d\n", instrument_R2,flowcell_lane_R2,tile_number_R2,x_coordinate_R2,y_coordinate_R2);
#endif
                        return FAIL;
                }

                if(flowcell_lane_R1 !=  flowcell_lane_R2){
                        #ifdef DEBUG
                        //WARNING_MSG("Files seem to contain reads in different order:\n%s\n%s\n",name1, name2);
                        #endif
#ifdef READ_GROUP_ITEST
                        fprintf(stderr,"@%s:%d:%d:%d:%d\n", instrument_R1,flowcell_lane_R1,tile_number_R1,x_coordinate_R1,y_coordinate_R1);
                        fprintf(stderr,"@%s:%d:%d:%d:%d\n", instrument_R2,flowcell_lane_R2,tile_number_R2,x_coordinate_R2,y_coordinate_R2);
#endif
                        return FAIL;
                }

                if(strcmp(instrument_R1,instrument_R2)){
                        #ifdef DEBUG
                        ///WARNING_MSG("Files seem to contain reads in different order:\n%s\n%s\n",name1, name2);
                        #endif
#ifdef READ_GROUP_ITEST
                        fprintf(stderr,"@%s:%d:%d:%d:%d\n", instrument_R1,flowcell_lane_R1,tile_number_R1,x_coordinate_R1,y_coordinate_R1);
                        fprintf(stderr,"@%s:%d:%d:%d:%d\n", instrument_R2,flowcell_lane_R2,tile_number_R2,x_coordinate_R2,y_coordinate_R2);
#endif
                        return FAIL;
                }
        }
        if(detected == 1000){
                number_of_values_found = 0;
                for(i = 0; i < (int) strlen(name1);i++){
                        if(isspace(name1[i]) || name1[i] == ';'){
                                break;
                        }
                        if(name1[i] != name2[i]){
                                number_of_values_found = 1;
                                break;
                        }
                        /*if(isspace(name1[i])){
                          name1[i] = 0;
                          name2[i] = 0;
                          }

                          if(name1[i] == ';'){
                          name1[i] = 0;
                          name2[i] = 0;
                          }*/
                }
                if(number_of_values_found){

#ifdef DEBUG
                        WARNING_MSG("Files seem to contain reads in different order:\n%s\n%s\n",name1,name2);
                        #endif
#ifdef READ_GROUP_ITEST
                        fprintf(stderr,"@%s:%d:%d:%d:%d\n", instrument_R1,flowcell_lane_R1,tile_number_R1,x_coordinate_R1,y_coordinate_R1);
                        fprintf(stderr,"@%s:%d:%d:%d:%d\n", instrument_R2,flowcell_lane_R2,tile_number_R2,x_coordinate_R2,y_coordinate_R2);
                        return FAIL;
#endif
                }
        }


        return OK;
ERROR:
        return FAIL;
}


int sort_read_groups_based_on_arch_assign(struct read_groups* rg)
{
        struct rg_f_sort** list = NULL;
        int i;
        int g;


        MMALLOC(list, sizeof(struct rg_f_sort*) * rg->e[0]->num_files);

        for(i = 0;i < rg->e[0]->num_files;i++){
                list[i] = NULL;
                MMALLOC(list[i], sizeof(struct rg_f_sort));

        }
        for(g = 0; g < rg->num_groups;g++){
                for(i = 0; i < rg->e[g]->num_files;i++){
                        list[i]->name = rg->e[g]->filenames[i];
                        list[i]->index = rg->e[g]->arch_to_read_assignment[i];
                        list[i]->si = rg->e[g]->ssi[i];
                }
                qsort(list,rg->e[g]->num_files ,sizeof(struct rg_f_sort*) , qsort_rg_f);
                for(i = 0; i < rg->e[g]->num_files;i++){
                        rg->e[g]->filenames[i] = list[i]->name;
                        rg->e[g]->arch_to_read_assignment[i] = list[i]->index;
                        rg->e[g]->ssi[i] = list[i]->si;
                }
        }

        for(i = 0;i < rg->e[0]->num_files;i++){
                MFREE(list[i]);

        }
        MFREE(list);
        return OK;
ERROR:
        return FAIL;
}

int qsort_rg_f(const void *a, const void *b)
{
        const struct rg_f_sort **elem1 = (const struct rg_f_sort**) a;
        const struct rg_f_sort **elem2 = (const struct rg_f_sort**) b;
        if ( (*elem1)->index < (*elem2)->index){
                return -1;
        }else if ((*elem1)->index > (*elem2)->index){
                return 1;
        }else{

                return 0;
        }
}



#ifdef READ_GROUP_ITEST
int test_read_name_recognition(void);

int test_read_name_recognition(void)
{

        char* name1 = NULL;

        char* name2 = NULL;
        MMALLOC(name1, sizeof(char)* 1000);
        MMALLOC(name2, sizeof(char)* 1000);

        sprintf (name1, "EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG");
        sprintf (name2, "EAS139:136:FC706VJ:2:2104:15343:197393 2:Y:18:ATCACG");

        int i = 0;

        fprintf(stdout,"Running I/O Unit tests.\n");
        i = compare_read_names( name1, name2);
        fprintf(stdout,"	cmp1:\n	%s\n	%s\n	Result:%d\n",name1,name2,i);
        if(i == 1){
                ERROR_MSG("ERROR - names should match");
        }
        name1[0] = 0;
        name2[0] = 0;

        sprintf (name1, "EAS139:136:FC706VJ:2:2104:15343:197393 2:Y:18:ATCACG");
        sprintf (name2, "EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG");

        i = compare_read_names( name1, name2);
        fprintf(stdout,"	cmp2:\n	%s\n	%s\n	Result:%d\n",name1,name2,i);
        if(i == 1){
                ERROR_MSG("ERROR - names should match");
        }


        name1[0] = 0;
        name2[0] = 0;

        sprintf (name1, "HWUSI-EAS100R:6:73:941:1973#0/1");
        sprintf (name2, "HWUSI-EAS100R:6:73:941:1973#0/2");

        fprintf(stdout,"	cmp3:\n	%s\n	%s\n	Result:%d\n",name1,name2,i);
        i = compare_read_names( name1, name2);

        if(i == 1){
                ERROR_MSG("ERROR - names should match");
        }

        name1[0] = 0;
        name2[0] = 0;

        sprintf (name1, "HWUSI-EAS100R:6:73:941:1973#0/2");
        sprintf (name2, "HWUSI-EAS100R:6:73:941:1973#0/1");

        i = compare_read_names( name1, name2);
        fprintf(stdout,"	cmp4:\n	%s\n	%s\n	Result:%d\n",name1,name2,i);
        if(i == 1){
                ERROR_MSG("ERROR - names should match");
        }

        name1[0] = 0;
        name2[0] = 0;

        sprintf (name1, "EAS139:136:FC706VJ:2:2104:15343:197393 2:Y:18:ATCACG");
        sprintf (name2, "HWUSI-EAS100R:6:73:941:1973#0/1");

        i = compare_read_names( name1, name2);
        fprintf(stdout,"	cmp5:\n	%s\n	%s\n	Result:%d\n",name1,name2,i);
        if(i == 0){
                ERROR_MSG("ERROR - names should match");
        }

        name1[0] = 0;
        name2[0] = 0;

        sprintf (name1, "EAS139:136:FC706VJ:2:2104:15343:197393 2:Y:18:ATCACG");
        sprintf (name2, "EAS139:136:FC706VJ:2:2104:15343:197393 1:N:18:GGGACG");

        i = compare_read_names( name1, name2);
        fprintf(stdout,"	cmp6:\n	%s\n	%s\n	Result:%d\n",name1,name2,i);
        if(i == 1){
                ERROR_MSG("ERROR - names should match");
        }

        name1[0] = 0;
        name2[0] = 0;

        sprintf (name1, "EAS139:136:FC706VJ:2:2104:15343:197393 2:Y:18:ATCACG");
        sprintf (name2, "EAS139:136:FC706VJ:2:2104:15344:197393 1:N:18:GGGACG");

        i = compare_read_names(name1, name2);
        fprintf(stdout,"	cmp7:\n	%s\n	%s\n	Result:%d\n",name1,name2,i);
        if(i == 0){
                ERROR_MSG("ERROR - names should match");
        }

        name1[0] = 0;
        name2[0] = 0;

        sprintf (name1, "HWUSI-EAS747_0040_FC64GRTAAXX:8:1:3268:1065#0/1");
        sprintf (name2, "HWUSI-EAS747_0040_FC64GRTAAXX:8:1:3268:1065#0/2");
        i = compare_read_names( name1, name2);
        fprintf(stdout,"	cmp8:\n	%s\n	%s\n	Result:%d\n",name1,name2,i);
        if(i == 1){
                ERROR_MSG("ERROR - names should match");

        }

        MFREE(name1);
        MFREE(name2);

        return OK;
ERROR:


        if(name1){
                MFREE(name1);
        }
        if(name2){
                MFREE(name2);
        }

        return FAIL;
}

int main (int argc,char * argv[]) {
        //struct parameters* param = 0;
        test_read_name_recognition();

        if(argc > 1){
                int i;
                char** files = NULL;
                int num_files;
                struct read_groups* rg = NULL;
                MMALLOC(files, sizeof(char*) * argc);
                num_files = 0;
                for(i = 0; i < argc-1;i++){
                        files[num_files] = argv[i+1];
                        num_files++;
                }

                LOG_MSG("options: %s", argv[1]);
                RUN(generate_read_groups(&rg, files, num_files));

                free_read_groups(rg);
                MFREE(files);
        }
        //free_param(param);
        return EXIT_SUCCESS;
ERROR:
        return EXIT_FAILURE;
}


#endif
