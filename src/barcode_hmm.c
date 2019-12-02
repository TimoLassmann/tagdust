#include "barcode_hmm.h"

#include "core_hmm_functions.h"
#include "hmm.h"
#include "hmm_model_bag.h"
#include "seq_stats.h"
#include "misc.h"

int hmm_controller_multiple(struct parameters* param)
{
        struct read_info*** read_info_container = NULL;
        struct sequence_stats_info** sequence_stats_info_container = NULL;
        struct model_bag** model_bag_container = NULL;

        int* numseqs = NULL;
        int status;

        struct fasta* reference_fasta = NULL;
        struct log_information* li = NULL;

        int i,j,c;//,f;
        int num_out_reads = 0;
        long long int barcode_present = 0;

        char* read_present = NULL;

        MMALLOC(read_present,sizeof(char)* param->infiles);

        for(i = 0; i < param->infiles;i++){
                read_present[i] = 0;
        }


        //long long int read_present = 0;

        FILE** file_container = NULL;
        int (*fp)(struct read_info** ,struct parameters*,FILE* ,int* buffer_count) = NULL;

        MMALLOC(read_info_container, sizeof(struct read_info**) * param->infiles);
        MMALLOC(sequence_stats_info_container,sizeof(struct sequence_stats_info*)* param->infiles);
        MMALLOC(model_bag_container, sizeof(struct model_bag*) * param->infiles);
        MMALLOC(file_container,sizeof(FILE*) * param->infiles);

        MMALLOC(param->read_structures, sizeof(struct read_info*) * param->infiles);
        MMALLOC(param->confidence_thresholds, sizeof(float) * param->infiles);

        MMALLOC(numseqs, sizeof(int) * param->infiles);

        for(i= 0; i < param->infiles;i++){
                read_info_container[i] = NULL;
                sequence_stats_info_container[i] = NULL;
                model_bag_container[i] = NULL;
                file_container[i] = NULL;
                param->read_structures[i] = NULL;
                param->confidence_thresholds[i] = 0.0f;
        }
        /* Look through input files and find best architecture for each input file... */

        for(i = 0; i < param->infiles;i++){
                if( !i && param->read_structure->num_segments){
                        param->read_structures[0] = param->read_structure;
                        param->read_structure = NULL;
                }else	 if(param->arch_file){
                        RUN(test_architectures(param, i));

//param = test_architectures(param, i);
                        param->read_structures[i] = param->read_structure;
                        param->read_structure = NULL;
                        //param->read_structure = malloc_read_structure();
                }else{
                        //Resort top default R:N
                        RUN(malloc_read_structure(&param->read_structure));

                        RUN(assign_segment_sequences(param, "R:N" , 0 ));
                        //param->read_structure = assign_segment_sequences(param, "R:N" , 0 );
                        if(QC_read_structure(param)){
                                sprintf(param->buffer,"Something wrong with architecture....\n");
                                //param->messages = append_message(param->messages, param->buffer);
                                free_param(param);
                                exit(EXIT_FAILURE);
                        }
                        param->read_structures[i] = param->read_structure;
                        param->read_structure = NULL;
                }
                for(j = 0; j < param->read_structures[i]->num_segments;j++){
                        if(param->read_structures[i]->type[j] == 'B'){
                                barcode_present |= (1 << i);
                        }
                        if(param->read_structures[i]->type[j] == 'R'){
                                read_present[i]++;
                        }
                }
        }

        // sanity check - barcode present in multiple reads? - Can't handle at the moment
        if(bitcount64(barcode_present) > 1){
                sprintf(param->buffer,"Barcodes seem to be in both architectures... \n");
                //param->messages = append_message(param->messages, param->buffer);
                free_param(param);
                exit(EXIT_FAILURE);
        }

        for(i = 0; i < param->infiles;i++){
                num_out_reads += (int) read_present[i];
        }
        /* check if outputs already exists  */
        for(i = 0; i < param->infiles;i++){
                if(barcode_present &  (1 << i)){
                        param->read_structure = param->read_structures[i];
                        param->read_structures[i] = NULL;
                        RUN(check_for_existing_demultiplexed_files_multiple(param, num_out_reads));

                        //if(check_for_existing_demultiplexed_files_multiple(param, num_out_reads)) KSLIB_XFAIL(kslFAIL, param->errmsg,"Error: some output files already exists.\n");
                        param->read_structures[i] =param->read_structure;
                        param->read_structure  = NULL;
                }
        }
        // Get ready for log space...
        init_logsum();
#if DEBUG
        param->num_query = 1001;
#else
#if  RTEST
        param->num_query = 1000;
#else
        param->num_query = 1000001;
#endif

#endif

        //malloc read_info structure (buffer for reads) AND collect basic sequence statistics from each input file...

        for(i =0 ; i < param->infiles;i++){
                read_info_container[i] = malloc_read_info(read_info_container[i],param->num_query);
                param->read_structure = param->read_structures[i];
                sequence_stats_info_container[i] =get_sequence_stats(param,read_info_container[i], i);
        }


        // determine confidence thresholds
        // this should happen all the time...


        if(!param->confidence_threshold ){
                for(i = 0; i < param->infiles;i++){
                        sprintf(param->buffer,"Determining threshold for read%d.\n",i);
                        //param->messages = append_message(param->messages, param->buffer);

                        param->read_structure = param->read_structures[i];
                        RUN(estimateQthreshold(param, sequence_stats_info_container[i]));
                        //if((status =estimateQthreshold(param, sequence_stats_info_container[i])) != kslOK) KSLIB_XFAIL(kslFAIL,param->errmsg,"estimateQthreshold failed.\n");
                        //param = estimateQthreshold(param, sequence_stats_info_container[i]);
                        param->confidence_thresholds[i] = param->confidence_threshold;
                }
        }

        // Now malloc models...
        for(i = 0; i < param->infiles;i++){
                param->read_structure = param->read_structures[i];
                model_bag_container[i] =  init_model_bag(param, sequence_stats_info_container[i]);
        }

        //Read in known contaminants.
        if(param->reference_fasta){
                reference_fasta = get_fasta(reference_fasta,param->reference_fasta);
                MMALLOC(reference_fasta->mer_hash ,sizeof(int)* reference_fasta->numseq);
                for(i = 0; i < reference_fasta->numseq;i++){
                        reference_fasta->mer_hash[i] = 0;
                }
        }

        //start reading files and computation...
        for(i = 0;i < param->infiles;i++){
                file_container[i] =  io_handler(file_container[i] , i,param);
                numseqs[i] = 0;
        }
        if(param->sam == 0){
                fp = &read_fasta_fastq;
        }else {
                fp = &read_sam_chunk;
        }



        MMALLOC(li,sizeof(struct log_information));


        li->total_read = 0;

        li->num_EXTRACT_SUCCESS = 0;
        li->num_EXTRACT_FAIL_BAR_FINGER_NOT_FOUND = 0;
        li->num_EXTRACT_FAIL_READ_TOO_SHORT = 0;
        li->num_EXTRACT_FAIL_AMBIGIOUS_BARCODE = 0;
        li->num_EXTRACT_FAIL_ARCHITECTURE_MISMATCH = 0;
        li->num_EXTRACT_FAIL_MATCHES_ARTIFACTS = 0;
        li->num_EXTRACT_FAIL_LOW_COMPLEXITY= 0;

        // infinite loop...
        while(1){
                // read batches of sequences from input file(s)
                c = 0;
                for(i = 0; i < param->infiles;i++){

                        RUN(fp( read_info_container[i], param,file_container[i],&numseqs[i]));
                        //if((status = fp( read_info_container[i], param,file_container[i],&numseqs[i])) != kslOK) KSLIB_XEXCEPTION_SYS(kslEWRT,"Failed to read data chunk from file: %s", param->infile[i]);
                        //numseqs[i] = fp( read_info_container[i], param,file_container[i]);
                        c+=numseqs[i];
                }
                if(!c){
                        break;
                }

                // stop if files are not sorted...

                //check if number of reads read from each file are identical.
                for(i = 0; i < param->infiles-1;i++){
                        for(j = i +1; j < param->infiles;j++){
                                if(numseqs[i] != numseqs[j]){
                                        sprintf(param->buffer,"Input File:%s and %s differ in number of entries.\n", param->infile[i],param->infile[j]);
                                        //param->messages = append_message(param->messages, param->buffer);
                                        free_param(param);
                                        exit(EXIT_FAILURE);

                                }
                        }
                }
                //fprintf(stderr,"%d	%d\n",numseqs[0], li->total_read);

                // On first iteration compare top 1000 read names in all files...
                if(!li->total_read){
                        for(i = 0; i < param->infiles-1;i++){
                                for(j = i +1; j < param->infiles;j++){
                                        for(c = 0;c < MACRO_MIN(1000,numseqs[0] );c++){
                                                //fprintf(stderr,"%s\n%s\n", read_info_container[i][c]->name,read_info_container[j][c]->name);
                                                if(compare_read_names(param,read_info_container[i][c]->name,read_info_container[j][c]->name) ){
                                                        sprintf(param->buffer,"Files seem to contain reads in different order:\n%s\n%s\n", read_info_container[i][c]->name,read_info_container[j][c]->name);
                                                        //oparam->messages = append_message(param->messages, param->buffer);
                                                        free_param(param);
                                                        exit(EXIT_FAILURE);
                                                }
                                                //fprintf(stderr,"%s\n%s\n", read_info_container[i][c]->name,read_info_container[j][c]->name);
                                        }
                                }
                        }
                }


                // Check if any read in any file was longer than what was allocated based on the top X reads.
                for(i = 0; i < param->infiles;i++){

                        for(j = 0;j < numseqs[0];j++){
                                c = 0;
                                if(read_info_container[i][j]->len >=  sequence_stats_info_container[i]->max_seq_len){
                                        sequence_stats_info_container[i]->max_seq_len=read_info_container[i][j]->len ;
                                        c = 1;
                                }
                                if(c){
                                        fprintf(stderr," %d %d\n", read_info_container[i][j]->len ,sequence_stats_info_container[i]->max_seq_len);
                                        sprintf(param->buffer,"Long sequence found. Need to realloc model...\n");
                                        //param->messages = append_message(param->messages, param->buffer);

                                        free_model_bag(model_bag_container[i] );
                                        model_bag_container[i] = init_model_bag(param, sequence_stats_info_container[i]);

                                }
                        }
                }
                // Run HMM models over sequences
                for(i = 0;i < param->infiles;i++){
                        param->read_structure = param->read_structures[i];
                        param->confidence_threshold = param->confidence_thresholds[i];
                        if(param->read_structure->num_segments == 1 && param->read_structure->type[0] == 'R'){
                                //for(i = 0; i < numseq1;i++){
                                //	r1[i]->read_type = EXTRACT_SUCCESS;
                                //}
                                RUN(run_rna_dust(read_info_container[i],param,reference_fasta,numseqs[i]));

                                //if((status =run_rna_dust(read_info_container[i],param,reference_fasta,numseqs[i] )) != kslOK) KSLIB_XFAIL(kslFAIL,param->errmsg,"run_rna_dust failed\n");
                                //read_info_container[i] =  run_rna_dust(read_info_container[i],param,reference_fasta,numseqs[i]);
                        }else{

                                RUN(run_pHMM(0, model_bag_container[i] ,read_info_container[i], param, reference_fasta, numseqs[i], MODE_GET_LABEL));
                                //if((status =run_pHMM(0, model_bag_container[i] ,read_info_container[i], param, reference_fasta, numseqs[i], MODE_GET_LABEL)) != kslOK) KSLIB_XFAIL(kslFAIL,param->errmsg,"run_pHMM failed\n");
                                //model_bag_container[i] = run_pHMM(0, model_bag_container[i] ,read_info_container[i], param, reference_fasta, numseqs[i], MODE_GET_LABEL);
                                //mb_R1 =  run_pHMM(0,mb_R1,r1,param,reference_fasta,numseq1,MODE_GET_LABEL);
                        }
                }

                //make sure the read_arch is set to the one holding the  barcode AND copy which barcode was found to readinfocontainer 0
                for(i = 0; i< param->infiles;i++){
                        if(barcode_present & (1 << i)){
                                param->read_structure = param->read_structures[i];
                                //param->read_structure_R2 = 0;
                                if(i){
                                        for(j= 0;j < numseqs[0];j++){
                                                read_info_container[0][j]->barcode = read_info_container[i][j]->barcode ;
                                        }
                                }
                                break;
                        }
                }



                for(i = 0; i < numseqs[0] ;i++){
                        c = -100000;
                        for(j = 0; j < param->infiles;j++){
                                c = MACRO_MAX(read_info_container[j][i]->read_type , c);
                        }
                        read_info_container[0][i]->read_type = c;

                }
                print_all(read_info_container,param, numseqs[0], read_present);

                li->total_read += numseqs[0];

                for(i = 0; i <  numseqs[0];i++){
                        ///fprintf(stdout,"%d	%d	%d\n",i,ri[i]->read_type,ri[i]->barcode);
                        switch ((int) read_info_container[0][i]->read_type) {

                        case EXTRACT_SUCCESS:
                                //	print_seq(ri[i],outfile);
                                li->num_EXTRACT_SUCCESS++;
                                break;
                        case EXTRACT_FAIL_BAR_FINGER_NOT_FOUND:
                                li->num_EXTRACT_FAIL_BAR_FINGER_NOT_FOUND++;
                                break;
                        case  EXTRACT_FAIL_READ_TOO_SHORT:
                                li->num_EXTRACT_FAIL_READ_TOO_SHORT++;
                                break;
                        case  EXTRACT_FAIL_ARCHITECTURE_MISMATCH:
                                li->num_EXTRACT_FAIL_ARCHITECTURE_MISMATCH++;
                                break;
                        case  EXTRACT_FAIL_MATCHES_ARTIFACTS:
                                li->num_EXTRACT_FAIL_MATCHES_ARTIFACTS++;
                        case  EXTRACT_FAIL_LOW_COMPLEXITY:
                                li->num_EXTRACT_FAIL_LOW_COMPLEXITY++;
                                break;
                        default:
                                li->num_EXTRACT_FAIL_MATCHES_ARTIFACTS++;
                                //fprintf(stderr,"WTF: no reference but ended up here: %d\n",read_info_container[0][i]->read_type);
                                reference_fasta->mer_hash[ ((int)(read_info_container[0][i]->read_type) >> 8 ) -1] ++;
                                break;
                        }
                }
        }

        sprintf(param->buffer,"Done.\n\n");
        //param->messages = append_message(param->messages, param->buffer);

        for(i =0; i < param->infiles;i++){
                sprintf(param->buffer,"%s	Input file %d.\n",param->infile[i],i);
                //param->messages = append_message(param->messages, param->buffer);
        }


        sprintf(param->buffer,"%d	total input reads\n", li->total_read);
        //param->messages = append_message(param->messages, param->buffer);

        sprintf(param->buffer,"%0.2f	selected threshold\n", param->confidence_threshold);
        //param->messages = append_message(param->messages, param->buffer);

        sprintf(param->buffer,"%d	successfully extracted\n" ,li->num_EXTRACT_SUCCESS);
        ///param->messages = append_message(param->messages, param->buffer);

        sprintf(param->buffer,"%0.1f%%	extracted\n",  (float) li->num_EXTRACT_SUCCESS / (float) li->total_read  *100.0f);
        //param->messages = append_message(param->messages, param->buffer);

        sprintf(param->buffer,"%d	problems with architecture\n" , li->num_EXTRACT_FAIL_ARCHITECTURE_MISMATCH);
        //param->messages = append_message(param->messages, param->buffer);

        sprintf(param->buffer,"%d	barcode / UMI not found\n" ,li->num_EXTRACT_FAIL_BAR_FINGER_NOT_FOUND);
        //param->messages = append_message(param->messages, param->buffer);

        sprintf(param->buffer,"%d	too short\n" , li->num_EXTRACT_FAIL_READ_TOO_SHORT);
        //param->messages = append_message(param->messages, param->buffer);

        sprintf(param->buffer,"%d	low complexity\n" , li->num_EXTRACT_FAIL_LOW_COMPLEXITY);
        //param->messages = append_message(param->messages, param->buffer);

        sprintf(param->buffer,"%d	match artifacts:\n" , li->num_EXTRACT_FAIL_MATCHES_ARTIFACTS);
        //param->messages = append_message(param->messages, param->buffer);

        if(reference_fasta){
                for(i = 0; i < reference_fasta->numseq;i++){
                        if(reference_fasta->mer_hash[i]){
                                sprintf(param->buffer,"%d	%s\n" , reference_fasta->mer_hash[i], reference_fasta->sn[i]);
                                //param->messages = append_message(param->messages, param->buffer);
                        }

                }
                free_fasta(reference_fasta);

        }
        for(i = 0; i < param->infiles;i++){
                free_model_bag(model_bag_container[i]);
                free_read_info(read_info_container[i], param->num_query);
                MFREE(sequence_stats_info_container[i]);
                pclose(file_container[i]);
        }

        MFREE(li);

        param->read_structure = 0;

        MFREE(file_container);
        MFREE(model_bag_container);
        MFREE(sequence_stats_info_container);
        MFREE(read_info_container);
        MFREE(numseqs);
        MFREE(read_present);
        return OK;
ERROR:

        //fprintf(stderr,"%s\n",param->errmsg);
        return FAIL;
}

/*
  void hmm_controller_pe(struct parameters* param)
  {
	struct read_info** r1 = NULL;
	struct read_info** r2 = NULL;

	struct sequence_stats_info* ssi_R1 = NULL;
	struct sequence_stats_info* ssi_R2 = NULL;

	struct model_bag* mb_R1 = NULL;
	struct model_bag* mb_R2 = NULL;

	FILE* file1 = 0;
	FILE* file2 = 0;


	int (*fp)(struct read_info** ,struct parameters*,FILE* ) = NULL;

	int i,j;
	//int total_read = 0;
	//double sum = 0.0;
	// Check for output files....
	j = 0;
	for(i = 0; i < param->read_structure_R1->num_segments;i++){
  if(param->read_structure_R1->type[i] == 'B'){
  j |=1;
  }
  fprintf(stderr,"%c",  param->read_structure_R1->type[i]);
	}

	fprintf(stderr,"ARCH2: %d segments\n", param->read_structure_R2->num_segments);
	for(i = 0; i < param->read_structure_R2->num_segments;i++){
  if(param->read_structure_R2->type[i] == 'B'){
  j |=2;
  }
  fprintf(stderr,"%c",  param->read_structure_R2->type[i]);
	}


	if(j == 3){
  //param->read_structure = 0;
  sprintf(param->buffer,"Barcodes seem to be in both architectures... \n");
  param->messages = append_message(param->messages, param->buffer);
  free_param(param);
  exit(EXIT_FAILURE);
	}else if (j ==2){
  param->read_structure = param->read_structure_R2;
	}else{ // 0 or 1 - keep first arch
  param->read_structure = param->read_structure_R1;
  //param->read_structure_R1 = 0;
	}

	j = check_for_existing_demultiplexed_files(param);

	if(j){
  sprintf(param->buffer , "ERROR: %d output files with prefix %s already exist.\n", j,param->outfile);
  param->messages = append_message(param->messages, param->buffer  );
  free_param(param);
  exit(EXIT_FAILURE);
	}

	init_logsum();

	//double* back = 0;
	//int average_length = 0;


  #if DEBUG
	//printf("Debug\n");
  #if RTEST
	param->num_query = 1000;
  #else
	param->num_query = 1001;
  #endif
  #else
  #if RTEST
	param->num_query = 1000;
  #else
	param->num_query = 1000001;
  #endif

  #endif

	//param->num_query = 500000;



	r1 = malloc_read_info(r1, param->num_query );
	r2 = malloc_read_info(r2, param->num_query );



	param->read_structure = param->read_structure_R1;
	ssi_R1 = get_sequence_stats(param, r1, 0 );
	param->read_structure = param->read_structure_R2;
	ssi_R2 = get_sequence_stats(param, r2, 1 );

	if(!param->confidence_threshold ){
  sprintf(param->buffer,"Determining threshold for read1.\n");
  param->messages = append_message(param->messages, param->buffer);

  param->read_structure = param->read_structure_R1;
  param = estimateQthreshold(param,ssi_R1);
  param->confidence_threshold_R1 = param->confidence_threshold;

  sprintf(param->buffer,"Determining threshold for read2.\n");
  param->messages = append_message(param->messages, param->buffer);


  param->read_structure = param->read_structure_R2;
  param = estimateQthreshold(param,ssi_R2);
  param->confidence_threshold_R2 = param->confidence_threshold;

	}

	// Inits model.

	param->read_structure = param->read_structure_R1;
	mb_R1 = init_model_bag(param, ssi_R1);

	param->read_structure = param->read_structure_R2;
	mb_R2  = init_model_bag(param, ssi_R2);


	struct fasta* reference_fasta = 0;

	if(param->reference_fasta){
  reference_fasta = get_fasta(reference_fasta,param->reference_fasta);


  MMALLOC(reference_fasta->mer_hash ,sizeof(int)* reference_fasta->numseq);
  for(i = 0; i < reference_fasta->numseq;i++){
  reference_fasta->mer_hash[i] = 0;
  }

	}

	//file =  io_handler(file, file_num,param);
	if(param->sam == 0){
  fp = &read_fasta_fastq;
	}else {
  fp = &read_sam_chunk;
	}

	file1 =  io_handler(file1, 0,param);

	file2 =  io_handler(file2, 1,param);

	int numseq1,numseq2;


	struct log_information* li = 0;

	void* tmp = 0;

	MMALLOC(li,sizeof(struct log_information));


	li->total_read = 0;

	li->num_EXTRACT_SUCCESS = 0;
	li->num_EXTRACT_FAIL_BAR_FINGER_NOT_FOUND = 0;
	li->num_EXTRACT_FAIL_READ_TOO_SHORT = 0;
	li->num_EXTRACT_FAIL_AMBIGIOUS_BARCODE = 0;
	li->num_EXTRACT_FAIL_ARCHITECTURE_MISMATCH = 0;
	li->num_EXTRACT_FAIL_MATCHES_ARTIFACTS = 0;
	li->num_EXTRACT_FAIL_LOW_COMPLEXITY= 0;


	param->multiread = 2;

	while ((numseq1 = fp(r1, param,file1)) != 0){
  numseq2 = fp(r2, param,file2);
  if(numseq1 != numseq2){
  sprintf(param->buffer,"Two files seem to be of different length.\n");
  param->messages = append_message(param->messages, param->buffer);
  free_param(param);
  exit(EXIT_FAILURE);

  }
  j = 0;
  for(i = 0; i < numseq1;i++){
  if(r1[i]->len >=ssi_R1->max_seq_len){
  ssi_R1->max_seq_len = r1[i]->len;
  ///mb->current_dyn_length = ri[i]->len + 10;
  j = 1;
  }
  if(r2[i]->len >=ssi_R2->max_seq_len){
  ssi_R2->max_seq_len = r2[i]->len;
  ///mb->current_dyn_length = ri[i]->len + 10;
  j = 1;
  }

  }

  if(j){
  sprintf(param->buffer,"Long sequence found. Need to realloc model...\n");
  param->messages = append_message(param->messages, param->buffer);

  free_model_bag(mb_R1);

  mb_R1 = init_model_bag(param, ssi_R1);

  free_model_bag(mb_R2);

  mb_R2 = init_model_bag(param, ssi_R2);


  }
  if(!param->sim_numseq){
  for(i = 0; i < numseq1;i++){
  if(compare_read_names(param,r1[i]->name,r2[i]->name) ){
  sprintf(param->buffer,"Files seem to contain reads in different order:\n%s\n%s\n",r1[i]->name,r2[i]->name );
  param->messages = append_message(param->messages, param->buffer);
  free_param(param);
  exit(EXIT_FAILURE);
  }



  }
  }

  param->read_structure = param->read_structure_R1;
  param->confidence_threshold = param->confidence_threshold_R1;
  if(param->read_structure->num_segments == 1 && param->read_structure->type[0] == 'R'){
  //for(i = 0; i < numseq1;i++){
  //	r1[i]->read_type = EXTRACT_SUCCESS;
  //}

  r1 =  run_rna_dust(r1,param,reference_fasta,numseq1);
  }else{
  mb_R1 =  run_pHMM(0,mb_R1,r1,param,reference_fasta,numseq1,MODE_GET_LABEL);
  }


  param->read_structure = param->read_structure_R2;
  param->confidence_threshold = param->confidence_threshold_R2;
  if(param->read_structure->num_segments == 1 && param->read_structure->type[0] == 'R'){
  for(i = 0; i < numseq1;i++){
  r2[i]->read_type = EXTRACT_SUCCESS;
  }
  r2 =  run_rna_dust(r2,param,reference_fasta,numseq1);
  }else{
  mb_R2 =  run_pHMM(0,mb_R2,r2,param,reference_fasta,numseq1,MODE_GET_LABEL);
  }

  for(i = 0; i < numseq1;i++){
  if(r1[i]->read_type != EXTRACT_SUCCESS || r2[i]->read_type != EXTRACT_SUCCESS){
  if(r1[i]->read_type != r2[i]->read_type){
  if(r1[i]->read_type  == EXTRACT_FAIL_MATCHES_ARTIFACTS ||  r2[i]->read_type ==EXTRACT_FAIL_MATCHES_ARTIFACTS){
  r1[i]->read_type =EXTRACT_FAIL_MATCHES_ARTIFACTS;
  }else if(r1[i]->read_type  == EXTRACT_FAIL_LOW_COMPLEXITY ||  r2[i]->read_type ==EXTRACT_FAIL_LOW_COMPLEXITY){
  r1[i]->read_type =EXTRACT_FAIL_LOW_COMPLEXITY;
  }else if(r1[i]->read_type  == EXTRACT_FAIL_ARCHITECTURE_MISMATCH ||  r2[i]->read_type ==EXTRACT_FAIL_ARCHITECTURE_MISMATCH){
  r1[i]->read_type =EXTRACT_FAIL_ARCHITECTURE_MISMATCH;
  }
  }
  //r1[i]->read_type = EXTRACT_FAIL_ARCHITECTURE_MISMATCH;
  }



  MREALLOC(r1[i]->seq,tmp, sizeof(char) * (r1[i]->len +  r2[i]->len +2));
  MREALLOC(r1[i]->qual,tmp,sizeof(char) * (r1[i]->len +  r2[i]->len +2));
  r1[i]->seq[r1[i]->len] = 65;
  r1[i]->qual[r1[i]->len] = 65;
  r1[i]->len++;

  for(j = 0; j < r2[i]->len;j++){
  r1[i]->seq[r1[i]->len] = r2[i]->seq[j];
  r1[i]->qual[r1[i]->len] = r2[i]->qual[j];
  r1[i]->len++;
  }
  r1[i]->seq[r1[i]->len] = 65;
  r1[i]->qual[r1[i]->len] = 65;
  }

  //param->read_structure = 0;
  //fprintf(stderr,"ARCH1: %d segments\n", param->read_structure_R1->num_segments);
  j = 0;
  for(i = 0; i < param->read_structure_R1->num_segments;i++){
  if(param->read_structure_R1->type[i] == 'B'){
  j |=1;
  }
  //fprintf(stderr,"%c",  param->read_structure_R1->type[i]);
  }

  //fprintf(stderr,"ARCH2: %d segments\n", param->read_structure_R2->num_segments);
  for(i = 0; i < param->read_structure_R2->num_segments;i++){
  if(param->read_structure_R2->type[i] == 'B'){
  j |=2;
  }
  //fprintf(stderr,"%c",  param->read_structure_R2->type[i]);
  }


  if(j == 3){
  //param->read_structure = 0;
  sprintf(param->buffer,"Barcodes seem to be in both architectures... \n");
  param->messages = append_message(param->messages, param->buffer);
  free_param(param);
  exit(EXIT_FAILURE);
  }else if (j ==2){
  param->read_structure = param->read_structure_R2;
  //param->read_structure_R2 = 0;
  for(i= 0;i< numseq1;i++){
  r1[i]->barcode =r2[i]->barcode ;
  }

  }else{ // 0 or 1 - keep first arch
  param->read_structure = param->read_structure_R1;
  //param->read_structure_R1 = 0;
  }


  //param->read_structure = param->read_structure_R1;

  print_split_files(param, r1, numseq1);

  li->total_read += numseq1;

  for(i = 0; i < numseq1;i++){
  ///fprintf(stdout,"%d	%d	%d\n",i,ri[i]->read_type,ri[i]->barcode);
  switch ((int) r1[i]->read_type) {

  case EXTRACT_SUCCESS:
  //	print_seq(ri[i],outfile);
  li->num_EXTRACT_SUCCESS++;
  break;
  case EXTRACT_FAIL_BAR_FINGER_NOT_FOUND:
  li->num_EXTRACT_FAIL_BAR_FINGER_NOT_FOUND++;
  break;
  case  EXTRACT_FAIL_READ_TOO_SHORT:
  li->num_EXTRACT_FAIL_READ_TOO_SHORT++;
  break;
  case  EXTRACT_FAIL_ARCHITECTURE_MISMATCH:
  li->num_EXTRACT_FAIL_ARCHITECTURE_MISMATCH++;
  break;
  case  EXTRACT_FAIL_MATCHES_ARTIFACTS:
  li->num_EXTRACT_FAIL_MATCHES_ARTIFACTS++;
  case  EXTRACT_FAIL_LOW_COMPLEXITY:
  li->num_EXTRACT_FAIL_LOW_COMPLEXITY++;
  break;
  default:
  li->num_EXTRACT_FAIL_MATCHES_ARTIFACTS++;
  reference_fasta->mer_hash[ ((int)(r1[i]->read_type) >> 8 ) -1] ++;
  break;
  }
  }
	}


	sprintf(param->buffer,"Done.\n\n");
	param->messages = append_message(param->messages, param->buffer);

	sprintf(param->buffer,"%s	Input file name 1.\n",param->infile[0]);
	param->messages = append_message(param->messages, param->buffer);

	sprintf(param->buffer,"%s	Input file name 2.\n",param->infile[1]);
	param->messages = append_message(param->messages, param->buffer);

	sprintf(param->buffer,"%d	total input reads\n", li->total_read);
	param->messages = append_message(param->messages, param->buffer);

	sprintf(param->buffer,"%0.2f	selected threshold\n", param->confidence_threshold);
	param->messages = append_message(param->messages, param->buffer);

	sprintf(param->buffer,"%d	successfully extracted\n" ,li->num_EXTRACT_SUCCESS);
	param->messages = append_message(param->messages, param->buffer);

	sprintf(param->buffer,"%0.1f%%	extracted\n",  (float) li->num_EXTRACT_SUCCESS / (float) li->total_read  *100.0f);
	param->messages = append_message(param->messages, param->buffer);

	sprintf(param->buffer,"%d	problems with architecture\n" , li->num_EXTRACT_FAIL_ARCHITECTURE_MISMATCH);
	param->messages = append_message(param->messages, param->buffer);

	sprintf(param->buffer,"%d	barcode / UMI not found\n" ,li->num_EXTRACT_FAIL_BAR_FINGER_NOT_FOUND);
	param->messages = append_message(param->messages, param->buffer);

	sprintf(param->buffer,"%d	too short\n" , li->num_EXTRACT_FAIL_READ_TOO_SHORT);
	param->messages = append_message(param->messages, param->buffer);

	sprintf(param->buffer,"%d	low complexity\n" , li->num_EXTRACT_FAIL_LOW_COMPLEXITY);
	param->messages = append_message(param->messages, param->buffer);

	sprintf(param->buffer,"%d	match artifacts:\n" , li->num_EXTRACT_FAIL_MATCHES_ARTIFACTS);
	param->messages = append_message(param->messages, param->buffer);

	if(reference_fasta){
  for(i = 0; i < reference_fasta->numseq;i++){
  if(reference_fasta->mer_hash[i]){
  sprintf(param->buffer,"%d	%s\n" , reference_fasta->mer_hash[i], reference_fasta->sn[i]);
  param->messages = append_message(param->messages, param->buffer);
  }

  }
  free_fasta(reference_fasta);

	}

	free_model_bag(mb_R1);
	free_model_bag(mb_R2);
	MFREE(li);

	free_read_info(r1,param->num_query);
	free_read_info(r2,param->num_query);

	pclose(file1);
	pclose(file2);

	MFREE(ssi_R1);
	MFREE(ssi_R2);

	param->read_structure = 0;

  }

*/

/** \fn void filter_controller(struct parameters* param,int (*fp)(struct read_info** ,struct parameters*,FILE* ),int file_num)
    \brief Runs all functions.

    Constructs HMM, runs it on sequences.

    \param param @ref parameters.
    \param fp Pointer to function used to read sequences (either SAM or fastq).
    \param filenum Number of input files.

*/

/*
  void filter_controller(struct parameters* param, int file_num)
  {
	struct read_info** ri = 0;
	int (*fp)(struct read_info** ,struct parameters*,FILE* ) = 0;
	FILE* outfile = 0;
	FILE* artifact_file = 0;
	int i;
	int numseq;
	//int total_read = 0;
	struct tm *ptr;
	int hour;
	char am_or_pm;
	char logfile[1000];



  #if DEBUG
	//printf("Debug\n");
  #if RTEST
	param->num_query = 1000;
  #else
	param->num_query = 501;
  #endif
  #else
  #if RTEST
	param->num_query = 1000;
  #else
	param->num_query = 1000001;
  #endif
  #endif

	//param->num_query = 500000;

	FILE* file = 0;


	ri = malloc_read_info(ri, param->num_query );




	struct fasta* reference_fasta = 0;

	if(param->reference_fasta){
  reference_fasta = get_fasta(reference_fasta,param->reference_fasta);


  MMALLOC(reference_fasta->mer_hash,sizeof(int)* reference_fasta->numseq);
  for(i = 0; i < reference_fasta->numseq;i++){
  reference_fasta->mer_hash[i] = 0;
  }

	}

	file =  io_handler(file, file_num,param);
	if(param->sam == 0){
  fp = &read_fasta_fastq;
	}else {
  fp = &read_sam_chunk;
	}

	if(param->outfile){
  if ((outfile = fopen( param->outfile, "w")) == NULL){
  sprintf(param->buffer,"can't open output file: %s\n",  param->outfile);
  param->messages = append_message(param->messages, param->buffer);
  free_param(param);
  exit(EXIT_FAILURE);
  }
	}else{
  outfile= stdout;
	}


	if(param->print_artifact){

  if ((artifact_file = fopen( param->print_artifact, "w")) == NULL){
  sprintf(param->buffer,"can't open artifact file: %s\n",  param->print_artifact);
  param->messages = append_message(param->messages, param->buffer);
  free_param(param);
  exit(EXIT_FAILURE);
  }

	}

	struct log_information* li = 0;

	MMALLOC(li,sizeof(struct log_information));


	li->total_read = 0;

	li->num_EXTRACT_SUCCESS = 0;
	li->num_EXTRACT_FAIL_BAR_FINGER_NOT_FOUND = 0;
	li->num_EXTRACT_FAIL_READ_TOO_SHORT = 0;
	li->num_EXTRACT_FAIL_AMBIGIOUS_BARCODE = 0;
	li->num_EXTRACT_FAIL_ARCHITECTURE_MISMATCH = 0;
	li->num_EXTRACT_FAIL_MATCHES_ARTIFACTS = 0;
	li->num_EXTRACT_FAIL_LOW_COMPLEXITY= 0;

	//total_read = 0;
	while ((numseq = fp(ri, param,file)) != 0){
  //	numseq = fp(ri, param,file);
  ri =  run_rna_dust(ri,param,reference_fasta,numseq);
  li->total_read += numseq;

  for(i = 0; i < numseq;i++){
  switch ((int) ri[i]->read_type) {

  case EXTRACT_SUCCESS:
  print_seq(ri[i],outfile);
  li->num_EXTRACT_SUCCESS++;
  break;
  case EXTRACT_FAIL_BAR_FINGER_NOT_FOUND:
  li->num_EXTRACT_FAIL_BAR_FINGER_NOT_FOUND++;
  break;
  case  EXTRACT_FAIL_READ_TOO_SHORT:
  li->num_EXTRACT_FAIL_READ_TOO_SHORT++;
  break;
  case  EXTRACT_FAIL_ARCHITECTURE_MISMATCH:
  li->num_EXTRACT_FAIL_ARCHITECTURE_MISMATCH++;
  break;
  case  EXTRACT_FAIL_MATCHES_ARTIFACTS:
  li->num_EXTRACT_FAIL_MATCHES_ARTIFACTS++;
  case  EXTRACT_FAIL_LOW_COMPLEXITY:
  li->num_EXTRACT_FAIL_LOW_COMPLEXITY++;
  break;
  default:
  li->num_EXTRACT_FAIL_MATCHES_ARTIFACTS++;

  reference_fasta->mer_hash[ ((int)(ri[i]->read_type) >> 8 ) -1] ++;

  //fprintf(stderr,"Guessing it matches sequence %d (%s)\n",  ((int)(ri[i]->prob) >> 8 ) -1,reference_fasta->sn[((int)(ri[i]->prob) >> 8)-1]   ) ;

  break;
  }
  }
  if(param->print_artifact){
  for(i = 0; i < numseq;i++){
  if ((int) ri[i]->read_type != EXTRACT_SUCCESS) {
  print_seq(ri[i],artifact_file);
  }
  }
  }
	}

	if(param->print_artifact){
  fclose(artifact_file);
	}

	if(param->outfile){
  fclose(outfile);
	}


	time_t current = time(NULL);
	ptr = localtime(&current);
	hour = ptr->tm_hour;
	if (hour <= 11)
  am_or_pm = 'a';
	else {
  hour -= 12;
  am_or_pm = 'p';
	}
	if (hour == 0){
  hour = 12;
	}

	if(param->log){
  sprintf (logfile, "%s/%s_tagdust_log.txt",param->log,shorten_pathname(param->infile[file_num]));
  fprintf(stderr,"LOGFILE::::%s\n",logfile);
  if ((outfile = fopen( logfile, "w")) == NULL){
  fprintf(stderr,"can't open logfile\n");
  exit(-1);
  }

  fprintf(outfile,"%s	Input file name.\n",param->infile[file_num]);
  fprintf(outfile,"%.2d-%.2d-%d;%2d:%.2d%cm	Date and Time\n",ptr->tm_mon + 1,ptr->tm_mday, ptr->tm_year + 1900,hour,ptr->tm_min, am_or_pm );
  fprintf(outfile,"%f	selected threshold\n", param->confidence_threshold);

  fprintf(outfile,"%d	total input reads\n", li->total_read);

  fprintf(outfile,"%d	successfully extracted\n" ,li->num_EXTRACT_SUCCESS);
  fprintf(outfile,"%0.1f%%	extracted\n",  (float) li->num_EXTRACT_SUCCESS / (float) li->total_read  *100.0f);
  fprintf(outfile,"%d	barcode / UMI not found\n" ,li->num_EXTRACT_FAIL_BAR_FINGER_NOT_FOUND);
  fprintf(outfile,"%d	too short\n" , li->num_EXTRACT_FAIL_READ_TOO_SHORT);
  fprintf(outfile,"%d	ambiguous barcode\n" , li->num_EXTRACT_FAIL_AMBIGIOUS_BARCODE);
  fprintf(outfile,"%d	problems with architecture\n" , li->num_EXTRACT_FAIL_ARCHITECTURE_MISMATCH);
  fprintf(outfile,"%d	low complexity\n" , li->num_EXTRACT_FAIL_LOW_COMPLEXITY);
  fprintf(outfile,"%d	match artifacts:\n" , li->num_EXTRACT_FAIL_MATCHES_ARTIFACTS);
  if(reference_fasta){
  for(i = 0; i < reference_fasta->numseq;i++){
  if(reference_fasta->mer_hash[i]){
  fprintf(outfile,"%d	%s\n" , reference_fasta->mer_hash[i], reference_fasta->sn[i]);
  }

  }
  }
  fclose(outfile);
	}

	if(reference_fasta){
  free_fasta(reference_fasta);
	}


	MFREE(li);

	free_read_info(ri,param->num_query);
	pclose(file);
  }

*/
/** \fn void hmm_controller(struct parameters* param,int (*fp)(struct read_info** ,struct parameters*,FILE* ),int file_num)
    \brief Runs all functions.

    Constructs HMM, runs it on sequences.

    \param param @ref parameters.
    \param fp Pointer to function used to read sequences (either SAM or fastq).
    \param filenum Number of input files.

*/

/*
  void hmm_controller(struct parameters* param,int file_num)
  {
	struct read_info** ri = 0;
	int (*fp)(struct read_info** ,struct parameters*,FILE* ) = 0;

	int i,j;
	int numseq;
	//int total_read = 0;
	//double sum = 0.0;

	j = check_for_existing_demultiplexed_files(param);

	if(j){
  sprintf(param->buffer , "ERROR: %d output files with prefix %s already exist.\n", j,param->outfile);
  param->messages = append_message(param->messages, param->buffer  );
  free_param(param);


  exit(EXIT_FAILURE);
	}


	init_logsum();

	//double* back = 0;
	//int average_length = 0;

	//back = malloc(sizeof(double)*5);

  #if DEBUG
	//printf("Debug\n");
  #if RTEST
	param->num_query = 1000;
  #else
	param->num_query = 1001;
  #endif

  #else
  #if RTEST
	param->num_query = 1000;
  #else
	param->num_query = 1000001;
  #endif

  #endif

	//param->num_query = 500000;

	FILE* file = 0;

	ri = malloc_read_info(ri, param->num_query );





	struct sequence_stats_info* ssi = get_sequence_stats(param, ri, file_num );

	if(param->read_structure->num_segments == 1 && param->read_structure->type[0] == 'R'){
  sprintf(param->buffer,"When using \" -1 R:N\" TagDust will echo reads and apply post fileting steps.\n");
  param->messages = append_message(param->messages, param->buffer);
	}else{

  if(!param->confidence_threshold ){
  param = estimateQthreshold(param, ssi);
  }
	}


	// Inits model.

	struct model_bag* mb = init_model_bag(param, ssi);

	struct fasta* reference_fasta = 0;

	if(param->reference_fasta){
  reference_fasta = get_fasta(reference_fasta,param->reference_fasta);


  MMALLOC(reference_fasta->mer_hash,sizeof(int)* reference_fasta->numseq);
  for(i = 0; i < reference_fasta->numseq;i++){
  reference_fasta->mer_hash[i] = 0;
  }

	}

	file =  io_handler(file, file_num,param);
	if(param->sam == 0){
  fp = &read_fasta_fastq;
	}else {
  fp = &read_sam_chunk;
	}



	if(!param->train ){

	}else if( !strcmp( param->train , "full")){
  for(i = 0; i < 10;i++){
  fprintf(stderr,"Iteration %d\n",i);
  mb = set_model_e_to_laplace(mb);
  while ((numseq = fp(ri, param,file)) != 0){
  //	numseq = fp(ri, param,file);
  mb =  run_pHMM(0,mb,ri,param,reference_fasta,numseq,MODE_TRAIN);
  }
  pclose(file);
  file =  io_handler(file, file_num,param);
  //rewind(file);
  for(j = 0; j < mb->num_models;j++){
  print_model(mb->model[j]);
  mb->model[j] = reestimate(mb->model[j], 0);
  print_model(mb->model[j]);
  }
  }

	}else if (!strcmp(param->train, "half" )){
  for(i = 0; i < 10;i++){
  fprintf(stderr,"Iteration %d\n",i);
  while ((numseq = fp(ri, param,file)) != 0){
  //	numseq = fp(ri, param,file);
  mb =  run_pHMM(0,mb,ri,param,reference_fasta,numseq,MODE_TRAIN);
  }
  pclose(file);
  file =  io_handler(file, file_num,param);
  //rewind(file);
  for(j = 0; j < mb->num_models;j++){
  print_model(mb->model[i]);
  mb->model[j] = reestimate(mb->model[j], 2);
  print_model(mb->model[i]);
  }
  }
	}

	pclose(file);




	file =  io_handler(file, file_num,param);


	struct log_information* li = 0;

	MMALLOC(li,sizeof(struct log_information));


	li->total_read = 0;

	li->num_EXTRACT_SUCCESS = 0;
	li->num_EXTRACT_FAIL_BAR_FINGER_NOT_FOUND = 0;
	li->num_EXTRACT_FAIL_READ_TOO_SHORT = 0;
	li->num_EXTRACT_FAIL_AMBIGIOUS_BARCODE = 0;
	li->num_EXTRACT_FAIL_ARCHITECTURE_MISMATCH = 0;
	li->num_EXTRACT_FAIL_MATCHES_ARTIFACTS = 0;
	li->num_EXTRACT_FAIL_LOW_COMPLEXITY= 0;

	while ((numseq = fp(ri, param,file)) != 0){
  j = 0;
  for(i = 0; i < numseq;i++){
  if(ri[i]->len >=ssi->max_seq_len){
  ssi->max_seq_len = ri[i]->len;
  ///mb->current_dyn_length = ri[i]->len + 10;
  j = 1;
  }
  }
  if(j){
  sprintf(param->buffer,"Long sequence found. Need to realloc model...\n");
  param->messages = append_message(param->messages, param->buffer);

  free_model_bag(mb);

  mb = init_model_bag(param, ssi);

  }

  //	numseq = fp(ri, param,file);

  if(param->read_structure->num_segments == 1 && param->read_structure->type[0] == 'R'){
  ri =  run_rna_dust(ri,param,reference_fasta,numseq);
  }else{
  mb =  run_pHMM(0,mb,ri,param,reference_fasta,numseq,MODE_GET_LABEL);
  }

  print_split_files(param, ri, numseq);

  li->total_read += numseq;

  for(i = 0; i < numseq;i++){
  ///fprintf(stdout,"%d	%d	%d\n",i,ri[i]->read_type,ri[i]->barcode);
  switch ((int) ri[i]->read_type) {

  case EXTRACT_SUCCESS:
  //	print_seq(ri[i],outfile);
  li->num_EXTRACT_SUCCESS++;
  break;
  case EXTRACT_FAIL_BAR_FINGER_NOT_FOUND:
  li->num_EXTRACT_FAIL_BAR_FINGER_NOT_FOUND++;
  break;
  case  EXTRACT_FAIL_READ_TOO_SHORT:
  li->num_EXTRACT_FAIL_READ_TOO_SHORT++;
  break;
  case  EXTRACT_FAIL_ARCHITECTURE_MISMATCH:
  li->num_EXTRACT_FAIL_ARCHITECTURE_MISMATCH++;
  break;
  case  EXTRACT_FAIL_MATCHES_ARTIFACTS:
  li->num_EXTRACT_FAIL_MATCHES_ARTIFACTS++;
  case  EXTRACT_FAIL_LOW_COMPLEXITY:
  li->num_EXTRACT_FAIL_LOW_COMPLEXITY++;
  break;
  default:
  li->num_EXTRACT_FAIL_MATCHES_ARTIFACTS++;
  reference_fasta->mer_hash[ ((int)(ri[i]->read_type) >> 8 ) -1] ++;
  break;
  }
  }
	}

	sprintf(param->buffer,"Done.\n\n");
	param->messages = append_message(param->messages, param->buffer);

	sprintf(param->buffer,"%s	Input file name.\n",param->infile[file_num]);
	param->messages = append_message(param->messages, param->buffer);

	sprintf(param->buffer,"%d	total input reads\n", li->total_read);
	param->messages = append_message(param->messages, param->buffer);

	sprintf(param->buffer,"%0.2f	selected threshold\n", param->confidence_threshold);
	param->messages = append_message(param->messages, param->buffer);

	sprintf(param->buffer,"%d	successfully extracted\n" ,li->num_EXTRACT_SUCCESS);
	param->messages = append_message(param->messages, param->buffer);

	sprintf(param->buffer,"%0.1f%%	extracted\n",  (float) li->num_EXTRACT_SUCCESS / (float) li->total_read  *100.0f);
	param->messages = append_message(param->messages, param->buffer);

	sprintf(param->buffer,"%d	problems with architecture\n" , li->num_EXTRACT_FAIL_ARCHITECTURE_MISMATCH);
	param->messages = append_message(param->messages, param->buffer);

	sprintf(param->buffer,"%d	barcode / UMI not found\n" ,li->num_EXTRACT_FAIL_BAR_FINGER_NOT_FOUND);
	param->messages = append_message(param->messages, param->buffer);

	sprintf(param->buffer,"%d	too short\n" , li->num_EXTRACT_FAIL_READ_TOO_SHORT);
	param->messages = append_message(param->messages, param->buffer);

	sprintf(param->buffer,"%d	low complexity\n" , li->num_EXTRACT_FAIL_LOW_COMPLEXITY);
	param->messages = append_message(param->messages, param->buffer);

	sprintf(param->buffer,"%d	match artifacts:\n" , li->num_EXTRACT_FAIL_MATCHES_ARTIFACTS);
	param->messages = append_message(param->messages, param->buffer);

	if(reference_fasta){
  for(i = 0; i < reference_fasta->numseq;i++){
  if(reference_fasta->mer_hash[i]){
  sprintf(param->buffer,"%d	%s\n" , reference_fasta->mer_hash[i], reference_fasta->sn[i]);
  param->messages = append_message(param->messages, param->buffer);
  }

  }
  free_fasta(reference_fasta);

	}

	free_model_bag(mb);
	MFREE(li);

	free_read_info(ri,param->num_query);
	pclose(file);
	MFREE(ssi);
  }


*/





/** \fn float model_information_content(struct model_bag*mb)

    \brief Calculated the information content of a model based on the match states.

    \f[

    IC =  - \sum_{i,j} p_{i,j} \log \left( p_{i,j} b_j \right)
    \f]



    \param mb  @ref model_bag - contains the HMM model.

*/


float model_information_content(struct model_bag*mb)
{
        float IC =0.0;
        int i,j,c,f;
        struct hmm_column* col = 0;
        float test_prob = 0;
        for(i = 0; i < mb->num_models;i++){

                for(j = 0;j < mb->model[i]->hmms[0]->num_columns ;j++){
                        for(c= 0; c < 5;c++){
                                test_prob = 0.0;
                                for(f = 0;f < mb->model[i]->num_hmms;f++){
                                        col = mb->model[i]->hmms[f]->hmm_column[j];
                                        test_prob += scaledprob2prob( col->m_emit[c])  * 1.0 / (float)mb->model[i]->num_hmms;

                                }
                                IC +=test_prob* log2 (   test_prob  /  scaledprob2prob( mb->model[0]->background_nuc_frequency[c]));

                        }
                }
        }
        return  IC;
}




/** \fn struct model_bag* estimate_length_distribution_of_partial_segments(struct model_bag*mb,struct read_info** ri,struct parameters* param, int numseq)

    \brief Initialized length related probabilities for partial segments.

    This function matches partial sequences to reads exactly. Based on the exact hits HMM parameters are derived using a normal distribution.
    \param mb  @ref model_bag - contains the HMM model.
    \param ri @ref read_info - contains the sequences.
    \param param @ref parameters - contain the parameters.
    \param numseq number of sequences.

    \bug Seems to set some probabilities to dead states at the end of HMMs...
*/

struct model_bag* estimate_length_distribution_of_partial_segments(struct model_bag*mb,struct read_info** ri,struct parameters* param, int numseq)
{
        int i,j,c;
        char* test_sequence = 0;

        struct model* model = 0;
        //struct hmm_column* col = 0;

        double mean;
        double stdev;

        double sum_prob = 0;
        double s0,s1,s2;
        int len = 0;

        float base_error = param->sequencer_error_rate;
        float indel_freq = param->indel_frequency;

        //5'
        if(param->read_structure->type[0] == 'P'){
                test_sequence = param->read_structure->sequence_matrix[0][0];
                len = (int) strlen(test_sequence);

                for(i = 0; i < len;i++){
                        test_sequence[i] = nuc_code[(int) test_sequence[i]];
                }
                //for(c = 0;c < len;c++){
                //	fprintf(stderr,"%c",*(test_sequence +c) + 65);
                //
                //}
                //fprintf(stderr,"\n");

                mean = 0;
                s0 = 0;
                s1 = 0;
                s2 = 0;


                for(i = 0; i < numseq;i++){
                        for(j = 0;j <= len ;j++){
                                for(c = 0;c < len-j;c++){
                                        if(ri[i]->seq[c] != test_sequence[j +c]){
                                                break;
                                        }
                                }
                                if(c == len-j && c > 3 ){
                                        s0++;
                                        s1 += len -j;
                                        s2 += (len-j) * (len-j);
                                        break;
                                }

                        }
                }
                if(!s0){
                        //fprintf(stderr,"ERROR: there seems to e not a single read containing the 5' partial sequence.\n");
                        sprintf(param->buffer,"ERROR: there seems to e not a single read containing the 5' partial sequence.\n");
                        //param->messages = append_message(param->messages, param->buffer);
                        free_param(param);
                        exit(EXIT_FAILURE);

                }

                mean = s1 / s0;
                stdev = sqrt(  (s0 * s2 - pow(s1,2.0))   /  (  s0 *(s0-1.0) )) ;
                if(stdev < 1){
                        stdev = 1;
                }
                //fprintf(stderr,"5: %f %f	%f\n", mean,  stdev,s0);

                if(mean <= 1){
                        //fprintf(stderr,"
                        sprintf(param->buffer,"WARNING: 5' partial segment seems not to be present in the data (length < 1).\n");
                        //param->messages = append_message(param->messages, param->buffer);
                        //free_param(param);

                }

                sum_prob = 0;

                for(i = 0; i <  len;i++){
                        //	fprintf(stderr,"%f ",sum_prob );
                        sum_prob +=gaussian_pdf(i , mean ,stdev);
                }

                //fprintf(stderr,"\n");

                //Init model ....
                model = mb->model[0];
                model->skip = prob2scaledprob(  gaussian_pdf(0 , mean ,stdev) / sum_prob    );

                //fprintf(stderr,"5SKIP:%f\t%f\n", model->skip, scaledprob2prob(model->skip));
                s1 = prob2scaledprob(0.0);
                //fprintf(stderr,"SUMPROB:%f\n", s1);
                s1 = logsum(s1, model->skip);
                //fprintf(stderr,"SUMPROB:%f\n", s1);
                //if(rs->type[key] == 'P'){// Partial - can skip and exit at every M / I state....
                len = model->hmms[0]->num_columns;
                for(i = 0 ; i < model->num_hmms;i++){

                        //model->M_to_silent[i] = prob2scaledprob(1.0);

                        for(j = 0; j < len;j++){
                                //			col = model->hmms[i]->hmm_column[j];
                                model->silent_to_M[i][j]  = prob2scaledprob(1.0 / (float) model->num_hmms) + prob2scaledprob(   gaussian_pdf(len-j , mean ,stdev) / sum_prob);

                                s1 = logsum(s1, model->silent_to_M[i][j]);

                        }

                        model->hmms[i] = set_hmm_transition_parameters(model->hmms[i],len, base_error, indel_freq, -1.0, -1.0);
                }

                model->skip = model->skip - s1;
                s2 = model->skip;
                for(i = 0 ; i < model->num_hmms;i++){

                        //model->M_to_silent[i] = prob2scaledprob(1.0);

                        for(j = 0; j < len;j++){
                                model->silent_to_M[i][j]  = model->silent_to_M[i][j]  - s1;
                        }
                }



                /*fprintf(stderr,"%f	skip\n", scaledprob2prob(model->skip));
                  for(j = 0; j < len;j++){
                  fprintf(stderr,"%f	%d len \n", scaledprob2prob(model->silent_to_M[0][j] ),j);
                  }*/

                //}




//#endif

        }

        for(c = 1; c < mb->num_models-1;c++){
                if(param->read_structure->type[c] == 'P'){

                        model = mb->model[c];
                        len = model->hmms[0]->num_columns;
                        for(i = 0 ; i < model->num_hmms;i++){

                                j = 0;
                                model->hmms[i] = set_hmm_transition_parameters(model->hmms[i],len, base_error, indel_freq, 0.1, -1.0); // 0.1 is exit parameter for partial internal sequence (I don't know how to set this.....)
                        }
                }
        }


        //3'
        if(param->read_structure->type[mb->num_models-1] == 'P'){
                test_sequence = param->read_structure->sequence_matrix[ mb->num_models-1][0];
                len = (int) strlen(test_sequence);
                for(i = 0; i < len;i++){
                        test_sequence[i] = nuc_code[(int) test_sequence[i]];
                }

                mean = 0;
                s0 = 0;
                s1 = 0;
                s2 = 0;

                for(i = 0; i < numseq;i++){
                        for(j = 0;j <= len ;j++){

                                for(c = 0;c < len-j;c++){
                                        if(ri[i]->seq[ri[i]->len - (len-j -c)] != test_sequence[c]){
                                                break;
                                        }
                                }
                                if(c == len-j  && c > 3){
                                        s0++;
                                        s1 += len -j;
                                        s2 += (len-j) * (len-j);
                                        break;
                                }
                        }
                }
                if(!s0){
                        sprintf(param->buffer,"ERROR: there seems to e not a single read containing the 3' partial sequence.\n");
                        //param->messages = append_message(param->messages, param->buffer);
                        free_param(param);
                        exit(EXIT_FAILURE);
                }
                mean = s1 / s0;
                stdev = sqrt(  (s0 * s2 - pow(s1,2.0))   /  (  s0 *(s0-1.0) )) ;
                if(stdev < 1){
                        stdev = 1;
                }

                //fprintf(stderr,"3: %f %f\n", mean,  stdev);
                if(mean <= 1){
                        sprintf(param->buffer,"WARNING: 3' partial segment seems not to be present in the data (length < 1).\n");
                        //param->messages = append_message(param->messages, param->buffer);

                }

                sum_prob = 0;

                for(i = 0; i <  len;i++){
                        sum_prob +=gaussian_pdf(i , mean ,stdev);
                }

                //Init model ....
                model = mb->model[mb->num_models-1];
                model->skip = prob2scaledprob(  gaussian_pdf(0 , mean ,stdev) / sum_prob    );


                s1 = model->skip;

                //fprintf(stderr,"3SKIP:%f\t%f\n", model->skip, scaledprob2prob(model->skip));
                //if(rs->type[key] == 'P'){// Partial - can skip and exit at every M / I state....
                //len = model->hmms[mb->num_models-1]->num_columns;
                for(i = 0 ; i < model->num_hmms;i++){

                        //model->M_to_silent[i] = prob2scaledprob(1.0);
                        model->silent_to_M[i][0] = prob2scaledprob(1.0 / (float) model->num_hmms) + prob2scaledprob(1.0 -  gaussian_pdf(0 , mean ,stdev) / sum_prob );
                        //fprintf(stderr,"move:%d:%f\t%f\n", i,model->silent_to_M[i][0], scaledprob2prob(model->silent_to_M[i][0]));

                        model->hmms[i] = set_hmm_transition_parameters(model->hmms[i],len, base_error, indel_freq, mean, stdev);
                }

        }
        return mb;
}



/** \fnstruct hmm* set_hmm_transition_parameters(struct hmm* hmm, int len,double base_error, double indel_freq,  double mean, double stdev)

    \brief Initialized transition probabilities of HMM segments.

    This function sets the transition probabilities within HMMs.

    If mean and stdev = -1.0 MSKIP transitions are set to 0 - the sequence has to match the entire HMM.
    If mean is > -1.0 and stdev = -1.0 MSKIP is set to a constant throughout the model.
    If mean and stdev are not -1.0 MSKIP is set according to a gaussian distribiution - i.e. leaving the model towards the end is more probable than at the beginning...




    \param hmm  @ref hmm - contains the HMM model.
    \param len - length of the HMM.
    \param base_error - seuqencer error rate.
    \param indel_freq - InDel frequency.

    \param mean - mean length of modelled 3' sequence
    \param stdev - standard deviation of read length.

*/






/** \fn struct model_bag* run_pHMM(struct model_bag* mb,struct read_info** ri,struct parameters* param,int numseq, int mode)
    \brief Starts threads to run HMM functions on subsets of sequences.
    Depending on the specifies mode, this function startes threads with different function pointers.


    \param ri @ref read_info.
    \param param @ref parameters.
    \param numseq Number of sequences.
    \param mode Determines which function to run.

*/
int run_pHMM(struct arch_bag* ab,struct model_bag* mb,struct read_info** ri,struct parameters* param,struct fasta* reference_fasta ,int numseq, int mode)
{
        struct thread_data* thread_data = 0;
        int status;


        pthread_t threads[param->num_threads];
        pthread_attr_t attr;
        int i,t;

        int interval = 0;
        int rc;
        /*

        MMALLOC(thread_data,sizeof(struct thread_data)* param->num_threads);

        interval =  (int)((double)numseq /(double)param->num_threads);

        for(t = 0;t < param->num_threads ;t++) {
                thread_data[t].fasta = reference_fasta;
                thread_data[t].ri = ri;
                thread_data[t].mb = copy_model_bag(mb);
                thread_data[t].start = t*interval;
                thread_data[t].end = t*interval + interval;
                thread_data[t].param = param;
                thread_data[t].ab = 0;
        }
        thread_data[param->num_threads-1].end = numseq;

        if(ab && (mode == MODE_ARCH_COMP) ){
                for(t = 0;t < param->num_threads ;t++) {
                        MMALLOC(thread_data[t].ab,sizeof(struct arch_bag));
                        thread_data[t].ab->num_arch = ab->num_arch;
                        thread_data[t].ab->arch_posterior = 0;
                        thread_data[t].ab->archs = 0;
                        thread_data[t].ab->command_line = 0;
                        MMALLOC(thread_data[t].ab->archs,sizeof(struct model_bag*) * ab->num_arch );
                        MMALLOC(thread_data[t].ab->arch_posterior,sizeof(float) * ab->num_arch );
                        for(i = 0; i < ab->num_arch;i++){
                                thread_data[t].ab->archs[i] =copy_model_bag(ab->archs[i]);
                                thread_data[t].ab->arch_posterior[i] = prob2scaledprob(1.0);
                        }
                }
        }

        rc = pthread_attr_init(&attr);
        if(rc){
                sprintf(param->buffer,"ERROR; return code from pthread_attr_init() is %d\n", rc);
                //param->messages = append_message(param->messages, param->buffer);

                free_param(param);
                exit(EXIT_FAILURE);
        }
        pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

        for(t = 0;t < param->num_threads;t++) {
                switch (mode) {
                case MODE_GET_PROB:
                        rc = pthread_create(&threads[t], &attr, do_probability_estimation, (void *) &thread_data[t]);
                        break;
                case MODE_GET_LABEL:
                        rc = pthread_create(&threads[t], &attr, do_label_thread, (void *) &thread_data[t]);
                        break;
                case MODE_TRAIN:
                        rc = pthread_create(&threads[t], &attr, do_baum_welch_thread, (void *) &thread_data[t]);
                        break;

                case MODE_ARCH_COMP:
                        rc = pthread_create(&threads[t], &attr, do_arch_comparison, (void *) &thread_data[t]);
                        break;
                }

                if (rc) {
                        sprintf(param->buffer,"ERROR; return code from pthread_create() is %d\n", rc);
                        //param->messages = append_message(param->messages, param->buffer);
                        free_param(param);
                        exit(EXIT_FAILURE );
                }
        }

        pthread_attr_destroy(&attr);

        for (t = 0;t < param->num_threads;t++){
                rc = pthread_join(threads[t], NULL);
                if (rc){
                        sprintf(param->buffer,"ERROR; return code from pthread_join()is %d\n", rc);
                        //param->messages = append_message(param->messages, param->buffer);
                        free_param(param);
                        exit(EXIT_FAILURE );
                }
        }

        for (t = 0;t < param->num_threads;t++){
                for(i = 0; i < mb->num_models;i++){

                        mb->model[i] = copy_estimated_parameter(mb->model[i], thread_data[t].mb->model[i]);
                }
        }

        if(ab && (mode == MODE_ARCH_COMP) ){
                for(t = 0;t < param->num_threads ;t++) {
                        for(i = 0; i < ab->num_arch;i++){
                                free_model_bag(thread_data[t].ab->archs[i]);// =copy_model_bag(ab->archs[i]);
                                //here I sum the posteriors from the different runs!!!
                                //ab->arch_posterior[i] = logsum(ab->arch_posterior[i] , thread_data[t].ab->arch_posterior[i]);
                                ab->arch_posterior[i]  += thread_data[t].ab->arch_posterior[i];
                        }

                        MFREE(thread_data[t].ab->archs);// = malloc(sizeof(struct model_bag*) * ab->num_arch );
                        MFREE(thread_data[t].ab->arch_posterior);//  = malloc(sizeof(float) * ab->num_arch );
                        MFREE(thread_data[t].ab);// = malloc(struct arch_bag);

                }

                float sum = 0.0;
                sum = ab->arch_posterior[0];
                for(i = 1; i < ab->num_arch;i++){
                        sum = logsum(sum, ab->arch_posterior[i]);
                }
                for(i = 0; i < ab->num_arch;i++){
                        ab->arch_posterior[i] = ab->arch_posterior[i]  - sum;
                }
        }

        */

        for(t = 0;t < param->num_threads;t++) {
                free_model_bag(thread_data[t].mb);
        }

        MFREE(thread_data);
        return OK;
ERROR:
        return FAIL;
}



/** \fnstruct read_info** run_rna_dust(struct read_info** ri,struct parameters* param,struct fasta* reference_fasta ,int numseq)
    \brief Starts threads to run artifact matching  on subsets of sequences.


    \param ri @ref read_info.
    \param param @ref parameters.
    \param numseq Number of sequences.
    \param reference_fasta file to match against.

*/
int run_rna_dust(struct read_info** ri,struct parameters* param,struct fasta* reference_fasta ,int numseq)
{
        struct thread_data* thread_data = 0;
        int status;


        pthread_t threads[param->num_threads];
        pthread_attr_t attr;
        int t;

        int interval = 0;
        int rc;

        /*MMALLOC(thread_data,sizeof(struct thread_data)* param->num_threads);


        interval =  (int)((double)numseq /(double)param->num_threads);

        for(t = 0;t < param->num_threads ;t++) {
                thread_data[t].fasta = reference_fasta;
                thread_data[t].ri = ri;
                thread_data[t].start = t*interval;
                thread_data[t].end = t*interval + interval;
                thread_data[t].param = param;
        }
        thread_data[param->num_threads-1].end = numseq;
        //unsigned int seed = (unsigned int) (time(NULL) * ( 42));

        rc = pthread_attr_init(&attr);
        if(rc){
                sprintf(param->buffer,"ERROR; return code from pthread_attr_init() is %d\n", rc);
                //param->messages = append_message(param->messages, param->buffer);
                free_param(param);
                exit(EXIT_FAILURE);
        }
        pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

        for(t = 0;t < param->num_threads;t++) {

                rc = pthread_create(&threads[t], &attr, do_rna_dust, (void *) &thread_data[t]);
                if (rc) {
                        sprintf(param->buffer,"ERROR; return code from pthread_create() is %d\n", rc);
                        //param->messages = append_message(param->messages, param->buffer);
                        free_param(param);
                        exit(EXIT_FAILURE);
                }
        }

        pthread_attr_destroy(&attr);

        for (t = 0;t < param->num_threads;t++){
                rc = pthread_join(threads[t], NULL);
                if (rc){
                        sprintf(param->buffer,"ERROR; return code from pthread_join() is %d\n", rc);
                        //param->messages = append_message(param->messages, param->buffer);
                        free_param(param);
                        exit(EXIT_FAILURE);
                }
        }

        MFREE(thread_data);
        */
        return OK;
ERROR:
        return FAIL;
}



void* do_arch_comparison(void *threadarg)
{
        struct thread_data *data;
        data = (struct thread_data *) threadarg;

        struct read_info** ri  = data->ri;


        int start = data->start;
        int end = data->end;
        int i,j;

        //float raw[100];
        //float sum;

        for(i = start; i < end;i++){

                //fprintf(stderr,"%d\n", i);

                for(j = 0; j < data->ab->num_arch;j++){
                        data->ab->archs[j] = backward(data->ab->archs[j] , ri[i]->seq ,ri[i]->len);
                        //raw[j] = data->ab->archs[j]->b_score;
                        data->ab->arch_posterior[j]  += data->ab->archs[j]->b_score;
                        //	fprintf(stderr,"%d:  %f\n", j,   data->ab->archs[j]->b_score );
                }
                /*
                  sum = raw[0];
                  for(j = 1; j < data->ab->num_arch;j++){
                  sum = logsum(sum, raw[j]);
                  }
                  for(j = 0; j < data->ab->num_arch;j++){
                  data->ab->arch_posterior[j] = logsum(data->ab->arch_posterior[j], raw[j] - sum );
                  }*/
        }


        //pthread_exit((void *) 0);
}






/** \fn void* do_probability_estimation(void *threadarg)
    \brief Estimates the probability of a sequence but without labelling it.

    This function calculates the equivalent of a mapping Q value to reflext the confidence we have in a match to the user specified read architecture.


    \f[

    P_{wrong} = 1 - \frac{P(x|M_i) }{P(x|M) + P(x | R)}
    \f]

    where \f$P(x|M_i)\f$ is the probability of selecing one particular barcode, \f$P(x|M)\f$ is the total probability of the model and \f$P(x|R)\f$ is the probability of the random model. The resulting probability is converted into a phred scaled quality value:

    \f[
    Q = -10 * log_{10}( P_{wrong})
    \f]

    \param threadarg  A @ref thread_data object used to pass data / parameters to function.
*/
void* do_probability_estimation(void *threadarg)
{
        struct thread_data *data;
        data = (struct thread_data *) threadarg;

        struct read_info** ri  = data->ri;
        struct model_bag* mb = data->mb;

        int matchstart = data->param->matchstart;
        int matchend = data->param->matchend;

        int start = data->start;
        int end = data->end;
        int i;
        int tmp = 0;

        float pbest = 0;
        float Q = 0;
        for(i = start; i < end;i++){
                ri[i]->mapq = prob2scaledprob(0.0);
        }

        if(matchstart != -1 || matchend != -1){
                for(i = start; i < end;i++){
                        tmp = matchend - matchstart ;
                        mb = backward(mb, ri[i]->seq + matchstart , tmp);
                        mb = forward_max_posterior_decoding(mb, ri[i] , ri[i]->seq+matchstart ,tmp );
                        pbest = 1.0 - scaledprob2prob(  (ri[i]->bar_prob + mb->f_score) - logsum(mb->f_score,  mb->r_score));
                        if(!pbest){
                                Q = 40.0;
                        }else if(pbest == 1.0){
                                Q = 0.0;
                        }else{
                                Q = -10.0 * log10(pbest) ;
                        }
                        ri[i]->mapq = Q;
                }
        }else{
                for(i = start; i < end;i++){
                        //fprintf(stderr,"%d\n",i);
                        mb = backward(mb, ri[i]->seq ,ri[i]->len);
                        mb = forward_max_posterior_decoding(mb, ri[i], ri[i]->seq ,ri[i]->len);

                        pbest =ri[i]->mapq;

                        pbest = logsum(pbest, mb->f_score);
                        pbest = logsum(pbest, mb->r_score);

                        pbest = 1.0 - scaledprob2prob(  (ri[i]->bar_prob + mb->f_score ) - pbest);

                        if(!pbest){
                                Q = 40.0;
                        }else if(pbest == 1.0){
                                Q = 0.0;
                        }else{
                                Q = -10.0 * log10(pbest) ;
                        }

                        ri[i]->mapq = Q;
                }
        }
        //pthread_exit((void *) 0);
}












/** \fn void* do_label_thread(void *threadarg)
    \brief Estimates the probability of a sequence and labels it.

    This function rund the forward and backward algorithm on each sequence, obtails posterior label probabilities for each residue and calculates the optimal sequence labeling based on the posteriors. In addition the quality value to reflext the confidence we have in a matchis calculated:

    \f[

    P_{wrong} = 1 - \frac{P(x|M_i) }{P(x|M) + P(x | R)}
    \f]

    where \f$P(x|M_i)\f$ is the probability of selecing one particular barcode, \f$P(x|M)\f$ is the total probability of the model and \f$P(x|R)\f$ is the probability of the random model. The resulting probability is converted into a phred scaled quality value:

    \f[
    Q = -10 * log_{10}( P_{wrong})
    \f]

    After labelling each sequence is compared against user defined artifact sequences and low comlexity sequences are filtered out.

    \param threadarg  A @ref thread_data object used to pass data / parameters to function.
*/
void* do_label_thread(void *threadarg)
{
        struct thread_data *data;
        data = (struct thread_data *) threadarg;

        struct read_info** ri  = data->ri;
        struct model_bag* mb = data->mb;
        int status;
        int matchstart = data->param->matchstart;
        int matchend = data->param->matchend;

        int start = data->start;
        int end = data->end;
        int i;
        int tmp = 0;
        float pbest,Q;

        for(i = start; i < end;i++){
                ri[i]->mapq = prob2scaledprob(0.0);
        }

        if(matchstart != -1 || matchend != -1){
                for(i = start; i < end;i++){
                        tmp = matchend - matchstart ;
                        mb = backward(mb, ri[i]->seq + matchstart , tmp);
                        mb = forward_max_posterior_decoding(mb, ri[i] , ri[i]->seq+matchstart ,tmp );



                        pbest =ri[i]->mapq;

                        pbest = logsum(pbest, mb->f_score);
                        pbest = logsum(pbest, mb->r_score);

                        pbest = 1.0 - scaledprob2prob(  (ri[i]->bar_prob + mb->f_score) - pbest);
                        if(!pbest){
                                Q = 40.0;
                        }else if(pbest == 1.0){
                                Q = 0.0;

                        }else{
                                Q = -10.0 * log10(pbest) ;
                        }
                        ri[i]->mapq = Q;
                }
        }else{
                for(i = start; i < end;i++){
                        mb = backward(mb, ri[i]->seq ,ri[i]->len);
                        mb = forward_max_posterior_decoding(mb, ri[i], ri[i]->seq ,ri[i]->len);


                        pbest =ri[i]->mapq;

                        pbest = logsum(pbest, mb->f_score);
                        pbest = logsum(pbest, mb->r_score);



                        pbest = 1.0 - scaledprob2prob(  (ri[i]->bar_prob + mb->f_score ) - pbest);
                        if(!pbest){
                                Q = 40.0;
                        }else if(pbest == 1.0){
                                Q = 0.0;

                        }else{
                                Q = -10.0 * log10(pbest) ;
                        }
                        //fprintf(stderr,"%d	%f	%f %f	%f\n",i,mb->f_score,ri[i]->bar_prob,ri[i]->mapq,mb->r_score);
                        //fprintf(stderr,"%d	f:%f	%f	init:%f	r:%f:Q: %f\n",i,mb->f_score,scaledprob2prob( ri[i]->bar_prob),ri[i]->mapq,mb->r_score,Q);
                        ri[i]->mapq = Q;
                }
        }

        for(i = start; i < end;i++){
                ri[i]->bar_prob = 100;
                //print_labelled_reads(mb,data->param ,ri[i]);
                RUN(extract_reads(mb,data->param,ri[i]));
                //if((status =extract_reads(mb,data->param,ri[i])) != kslOK) KSLIB_XEXCEPTION(kslFAIL,"Extract reads failed.\n");
                //ri[i] = extract_reads(mb,data->param,ri[i]);
        }

        if(data->param->reference_fasta){
                ri = match_to_reference(data);
        }

        if(data->param->dust){
                ri = dust_sequences(data);
        }
        return NULL;
ERROR:
        return NULL;
        pthread_exit((void *) 0);
}




/** \fn void* do_rna_dust(void *threadarg)
    \brief Match all reads against a reference and remove low complexity sequences.

    \param threadarg  A @ref thread_data object used to pass data / parameters to function.
*/
void* do_rna_dust(void *threadarg)
{
        struct thread_data *data;
        data = (struct thread_data *) threadarg;

        struct read_info** ri  = data->ri;

        int start = data->start;
        int end = data->end;
        int i;

        for(i = start; i < end;i++){
                ri[i]->read_type =   EXTRACT_SUCCESS;
                //ri[i] = extract_reads(mb,data->param,ri[i]);
        }

        if(data->param->reference_fasta){
                ri = match_to_reference(data);
        }

        if(data->param->dust){
                ri = dust_sequences(data);
        }

        pthread_exit((void *) 0);
}




/** \fn struct read_info** dust_sequences(struct thread_data *data)
    \brief Removel low complexity sequences.

    Runs the DUST algorithm on the first 64 nucleotides of a read.
    \param threadarg  A @ref thread_data object used to pass data / parameters to function.
*/

struct read_info** dust_sequences(struct thread_data *data)
{
        struct read_info** ri = data->ri;
        const int start = data->start;
        const int end = data->end;

        int dust_cut = data->param->dust ;

        int i,j,c;
        int key = 0;
        double triplet[64];
        double s = 0.0;
        int len;
        for(j = 0;j < 64;j++){
                triplet[j] = 0.0;
        }

        for(i = start; i < end;i++){
                c = 0;
                while(ri[i]->seq[c] == 65){
                        c++;
                }


                key = ((ri[i]->seq[c] & 0x3 ) << 2 )|  (ri[i]->seq[c+1] & 0x3 );

                len = ri[i]->len;
                if(len > 64){
                        len = 64;
                }
                c+= 2;

                for(j = c;j < len;j++){
                        if(ri[i]->seq[j] == 65){
                                break;
                        }
                        key = key << 2 | (ri[i]->seq[j] & 0x3 );
                        triplet[key & 0x3F]++;
                        c++;
                }
                s = 0.0;
                for(j = 0;j < 64;j++){

                        s+= triplet[j] * (triplet[j] -1.0) / 2.0;
                        triplet[j] = 0.0;
                }
                s = s / (double)(c-3) *10.0; //should be number of triplets -2 : i.e. len -2 triples -1 = len -3;

                if(s > dust_cut){
#if DEBUG
                        char alphabet[] = "ACGTNN";
                        for(j = 0;j < len;j++){
                                fprintf(stderr,"%c",alphabet[(int)ri[i]->seq[j]]);
                        }
                        fprintf(stderr,"\tLOW:%f\n",s);
#endif
                        ri[i]->read_type = EXTRACT_FAIL_LOW_COMPLEXITY;
                }
        }
        return ri;
}


/** \fn struct read_info** match_to_reference(struct thread_data *data)
    \brief Matches reads to artifactual sequences.

    Exhaustively compares reads to a fasta file of known artifact sequences. Uses a SSE version of the Myers bit-parallel dynamic programming algorithm.

    \param threadarg  A @ref thread_data object used to pass data / parameters to function.
*/

struct read_info** match_to_reference(struct thread_data *data)
{
        struct read_info** ri = data->ri;
        struct fasta* reference = data->fasta ;
        const int start = data->start;
        const int end = data->end;

        int error_cut = data->param->filter_error ;

        int i,j,c;
        int test = 1;
        //int reverse = 0;
        unsigned char* seq[4];

        int _MM_ALIGN16 lengths[4];
        int _MM_ALIGN16 errors[4];
        int _MM_ALIGN16 sequence_id[4];

        for(i = start; i <= end-4;i+=4){
                test = 1;
                //reverse = 0;
                for(c = 0;c < 4;c++){
                        errors[c] = 100000;
                        sequence_id[c] = 0;

                }
                for(j =0; j < reference->numseq;j++){
                        seq[0] = (unsigned char* ) ri[i]->seq;
                        seq[1] = (unsigned char* ) ri[i+1]->seq;
                        seq[2] = (unsigned char* ) ri[i+2]->seq;
                        seq[3] = (unsigned char* ) ri[i+3]->seq;
                        lengths[0] =  ri[i]->len;
                        lengths[1] =  ri[i+1]->len;
                        lengths[2] =  ri[i+2]->len;
                        lengths[3] =  ri[i+3]->len;
                        validate_bpm_sse(seq,lengths,reference->string +  reference->s_index[j],reference->s_index[j+1] - reference->s_index[j],4);
                        for(c = 0;c < 4;c++){
                                if(lengths[c] < errors[c]){
                                        errors[c] = lengths[c];
                                        sequence_id[c] = j+1;
                                }
                        }


                        seq[0] = reverse_complement((unsigned char* ) ri[i]->seq,ri[i]->len);
                        seq[1] = reverse_complement((unsigned char* ) ri[i+1]->seq,ri[i+1]->len);
                        seq[2] = reverse_complement((unsigned char* ) ri[i+2]->seq,ri[i+2]->len);
                        seq[3] = reverse_complement((unsigned char* ) ri[i+3]->seq,ri[i+3]->len);
                        lengths[0] =  ri[i]->len;
                        lengths[1] =  ri[i+1]->len;
                        lengths[2] =  ri[i+2]->len;
                        lengths[3] =  ri[i+3]->len;
                        validate_bpm_sse(seq,lengths,reference->string +  reference->s_index[j],reference->s_index[j+1] - reference->s_index[j],4);
                        for(c = 0;c < 4;c++){
                                if(lengths[c] < errors[c]){
                                        errors[c] = lengths[c];
                                        sequence_id[c] = j+1;
                                }
                        }

                        seq[0] = reverse_complement((unsigned char* ) ri[i]->seq,ri[i]->len);
                        seq[1] = reverse_complement((unsigned char* ) ri[i+1]->seq,ri[i+1]->len);
                        seq[2] = reverse_complement((unsigned char* ) ri[i+2]->seq,ri[i+2]->len);
                        seq[3] = reverse_complement((unsigned char* ) ri[i+3]->seq,ri[i+3]->len);
                }
                for(c = 0;c < 4;c++){
                        if(errors[c] <= error_cut){
                                if(ri[i+c]->read_type == EXTRACT_SUCCESS){
                                        ri[i+c]->read_type  =  (sequence_id[c]  << 8) | EXTRACT_FAIL_MATCHES_ARTIFACTS;
                                }
                        }
                }
        }

        while(i < end){
                //fprintf(stderr,"Looking at %d	%d	%d\n",i,start,end);
                test = 1;
                //reverse = 0;
                for(j =0; j < reference->numseq;j++){
                        c = bpm_check_error(reference->string +  reference->s_index[j], (unsigned char* )ri[i]->seq,reference->s_index[j+1] - reference->s_index[j] , ri[i]->len);
                        if(c <= error_cut){
                                test = 0;
                                sequence_id[0] = j+1;
                                break;
                        }
                        ri[i]->seq = (char* )reverse_complement((unsigned char* ) ri[i]->seq,   ri[i]->len);
                        c = bpm_check_error(reference->string +  reference->s_index[j], (unsigned char* )ri[i]->seq,reference->s_index[j+1] - reference->s_index[j] , ri[i]->len);
                        if(c <= error_cut){
                                ri[i]->seq =(char* ) reverse_complement( (unsigned char* )ri[i]->seq,   ri[i]->len);
                                test = 0;
                                sequence_id[0] = j+1;
                                break;
                        }
                        ri[i]->seq = (char* )reverse_complement((unsigned char* ) ri[i]->seq,   ri[i]->len);
                }
                if(!test){
                        if(ri[i]->read_type == EXTRACT_SUCCESS){
                                ri[i]->read_type  = (sequence_id[0] << 8) |  EXTRACT_FAIL_MATCHES_ARTIFACTS;
                        }
                }
                i++;
        }


        return ri;
}


/** \fn struct read_info* emit_random_sequence(struct model_bag* mb, struct read_info* ri,int average_length,unsigned int* seed )
    \brief Emits sequences from random HMM model.

    \param mb  The model @ref model_bag .
    \param ri  @ref read_info - emitted sequences are written here.
    \param average_length Average length of the sequences.
    \param seed Seed used for randomization.


    \deprecated not used.
*/


int  emit_random_sequence(struct model_bag* mb, struct read_info* ri,int average_length, struct rng_state* rng)
{
        int status;
        int current_length = 0;
        int allocated_length = 100;
        double r;// = (float)rand()/(float)my_rand_max;
        double sum = prob2scaledprob(0.0f);
        //char alpha[] = "ACGTN";
        int nuc,i;
        MFREE(ri->seq);
        MFREE(ri->name);
        MFREE(ri->qual);
        MFREE(ri->labels);
        ri->seq = 0;
        ri->name = 0;
        ri->qual = 0;
        ri->labels = 0;
        ri->len = 0;
        ri->read_type = 0;
        MMALLOC(ri->seq,sizeof(char) * allocated_length);

        while(current_length < average_length){

                while(1){
                        //emission
                        r = tl_random_double(rng);
                        sum = prob2scaledprob(0.0f);
                        for(nuc = 0;nuc < 5;nuc++){
                                sum = logsum(sum, mb->model[0]->background_nuc_frequency[nuc] );
                                if(r <  scaledprob2prob(sum)){
                                        ri->seq[current_length] = nuc;
                                        //fprintf(stderr,"%c",alpha[(int)nuc]);
                                        current_length++;
                                        break;
                                }
                        }
                        if(current_length == allocated_length){
                                allocated_length = allocated_length*2;
                                MREALLOC(ri->seq,sizeof(char) * allocated_length );
                        }
                        //transition
                        r = tl_random_double(rng);
                        //r = (float)rand()/(float)my_rand_max;
                        // prob2scaledprob(1.0 - (1.0 / (float)len));
                        if(r > 1.0 - (1.0 / (float)average_length)){
                                break;
                        }
                }
                //fprintf(stderr,"	%d\n",current_length);
                if(current_length < average_length){
                        current_length = 0;
                }

                //if(current_length+2 >= MAX_HMM_SEQ_LEN){
                //	current_length = 0;
                //}
        }

        MREALLOC(ri->seq ,sizeof(char) * (current_length+1));
        ri->seq[current_length] = 0;
        ri->len = current_length;

        MMALLOC(ri->qual,sizeof(char) * (current_length+1));
        for(i = 0; i < current_length;i++){
                ri->qual[i] = 'B';
        }

        ri->qual[current_length] = 0;

        MMALLOC(ri->name,sizeof(char) *2);
        ri->name[0] = 'N';
        ri->name[1] = 0;

        MMALLOC(ri->labels,sizeof(char) * (current_length+1));
        return OK;
ERROR:
        return FAIL;
}



/** \fn struct model_bag* estimate_model_from_labels(struct model_bag* mb, struct parameters* param,  struct read_info** ri,int numseq)
    \brief Estimate model based on labelled sequences

    \bug not complete - is very buggy.
    \deprecated not used.
*/

struct model_bag* estimate_model_from_labels(struct model_bag* mb, struct parameters* param,  struct read_info** ri,int numseq)
{
        int i,j,c1,c2,c3,g;//,bar,mem;
        char alpha[6] = "ACGTNN";
        //set all counts to 1;

        mb = set_model_e_to_laplace(mb);


        char seq[100];

        int current_position = 0;
        int current_hmm = 0;
        int current_segment = -1;
        for(i =0; i < numseq;i++){
                if(ri[i]->read_type == EXTRACT_SUCCESS){
                        current_position = 0;
                        current_hmm = -1;
                        current_segment = -1;
                        // Do stuff with labels & sequences...
                        for(j = 0; j < ri[i]->len;j++){

                                c1 = mb->label[(int)ri[i]->labels[j]];
                                c2 = c1 & 0xFFFF; //which segment
                                c3 = (c1 >> 16) & 0x7FFF; // which HMM in a particular segment....

                                if(c2 != current_segment){ // new segment!
                                        if(current_segment != -1){
                                                fprintf(stderr,"%c segment:\n", param->read_structure->type[current_segment] );
                                                for(g = 0;g < current_position;g++){
                                                        fprintf(stderr,"%d ", seq[g]);

                                                }
                                                fprintf(stderr,"\n");

                                                fprintf(stderr,"%s\n",param->read_structure->sequence_matrix[current_segment][current_hmm]);
                                        }

                                        current_position = 0;
                                        current_segment = c2;
                                        current_hmm = c3;

                                        //update M->silent
                                        //update I->silent...
                                        //fprintf(stderr,"%d	%d\n", current_segment,current_hmm);



                                }else{ // move into next HMM column.....
                                        seq[current_position] = ri[i]->seq[j];

                                        current_position++;
                                }





                        }
                        int key = 0;
                        fprintf(stderr,"%s\n", ri[i]->name);
                        for(j = 0; j < ri[i]->len;j++){
                                c1 = mb->label[(int)ri[i]->labels[j+1]];
                                c2 = c1 & 0xFFFF;
                                c3 = (c1 >> 16) & 0x7FFF;
                                fprintf(stderr,"%c",  alpha[(int)ri[i]->seq[j]] );
                        }
                        fprintf(stderr,"\n");
                        for(j = 0; j < ri[i]->len;j++){
                                c1 = mb->label[(int)ri[i]->labels[j+1]];
                                c2 = c1 & 0xFFFF;
                                c3 = (c1 >> 16) & 0x7FFF;
                                fprintf(stderr,"%c",   param->read_structure->type[c2] );
                                if(param->read_structure->type[c2] == 'F'){
                                        key = (key << 2 )|  (ri[i]->seq[j] & 0x3);
                                }

                        }
                        fprintf(stderr,"	key = %d\n",key);

                        for(j = 0; j < ri[i]->len;j++){
                                c1 = mb->label[(int)ri[i]->labels[j+1]];
                                c2 = c1 & 0xFFFF;
                                c3 = (c1 >> 16) & 0x7FFF;
                                fprintf(stderr,"%d", c3   );
                                if(param->read_structure->type[c2] == 'B'){
                                        //				bar = c3;
                                        //				mem = c2;
                                }
                        }
                        //		fprintf(stderr,"	bar= %d (%s)\n",bar,   param->read_structure->sequence_matrix[mem][bar] );
                }



        }



        return mb;
}



/** \fn  struct read_info*  extract_reads(struct model_bag* mb, struct parameters* param,  struct read_info* ri)
    \brief Extracts reads from labelled raw sequences.

    This function extracts the mappable reads from the raw sequences. Barcodes and Fingerprint sequences are decoded and appended to re read names using the BC and FP tag.

    \param mb The HMM model @ref model_bag.
    \param param @ref parameters .
    \param ri The reads.

*/


int extract_reads(struct model_bag* mb, struct parameters* param,  struct read_info* ri)
{
        int j,c1,c2,c3,key,bar,mem,fingerlen,required_finger_len,status;//,ret;
        char* buffer = 0;

        MMALLOC(buffer,sizeof(char)* mb->current_dyn_length );
        //assert(buffer != NULL);
        int s_pos = 0;
        key = 0;
        bar = -1;
        mem = -1;
        //ret = 0;
        int offset = 0;
        int len;
        int hmm_has_barcode = 0;
        int too_short = 0;
        int in_read = 0;

        len = ri->len;
        if(param->matchstart != -1 || param->matchend != -1){
                offset = param->matchstart;
                len = param->matchend - param->matchstart;
        }
        required_finger_len = 0;
        for(j = 0; j < param->read_structure->num_segments;j++){
                if(param->read_structure->type[j] == 'F'){
                        required_finger_len += (int) strlen(param->read_structure->sequence_matrix[j][0]);
                }
        }
        //KSL_DPRINTF3(("Requiured_len: %d\n",required_finger_len ));

        if(param->confidence_threshold <=  ri->mapq){
                fingerlen = 0;
                for(j = 0; j < len;j++){
                        c1 = mb->label[(int)ri->labels[j+1]];
                        c2 = c1 & 0xFFFF; //which segment
                        c3 = (c1 >> 16) & 0x7FFF; // which HMM in a particular segment....
                        //KSL_DPRINTF3(("%c",   param->read_structure->type[c2] ));
                        if(param->read_structure->type[c2] == 'F'){
                                //	required_finger_len += (int) strlen(param->read_structure->sequence_matrix[c2][0]);
                                fingerlen++;

                                key = (key << 2 )|  (ri->seq[j+offset] & 0x3);
                        }
                        if(param->read_structure->type[c2] == 'B'){
                                hmm_has_barcode = 1;
                                bar = c3;

                                if(bar == param->read_structure->numseq_in_segment[c2]-1){
                                        //fprintf(stderr,"EXTRACTING N!!!!\n");
                                        hmm_has_barcode = -1;
                                }
                                mem = c2;
                        }

                        if(param->read_structure->type[c2] == 'R'){
                                s_pos++;
                                if(!in_read){
                                        in_read= 1;
                                }
                        }else{
                                if(in_read){
                                        if(s_pos <param->minlen){
                                                too_short = 1;
                                                break;
                                        }
                                }
                                in_read = 0;
                                s_pos = 0;
                        }
                }
                if(in_read){
                        if(s_pos <param->minlen){
                                too_short = 1;

                        }
                }
                //KSL_DPRINTF3(("\n"));
                //KSL_DPRINTF3(("len: %d\n",fingerlen ));
                if(!too_short){
                        if(hmm_has_barcode == -1){
                                ri->read_type  = EXTRACT_FAIL_BAR_FINGER_NOT_FOUND;
                        }else if(hmm_has_barcode && required_finger_len){
                                if(fingerlen == required_finger_len && bar != -1){
                                        ri = make_extracted_read(mb, param,ri);
                                        ri->barcode =  (mem << 16) |   bar;

                                        //ri->barcode_string = param->read_structure->sequence_matrix[mem][bar];
                                        if(required_finger_len <= 255){
                                                ri->fingerprint = (key <<  8) | required_finger_len;
                                        }else{
                                                ri->fingerprint = (key <<  8) | 255;
                                        }
                                        ri->read_type = EXTRACT_SUCCESS;
                                }else{
                                        ri->read_type  = EXTRACT_FAIL_BAR_FINGER_NOT_FOUND; // something wrong with the architecture
                                }
                        }else if(hmm_has_barcode){
                                if(bar != -1){
                                        ri = make_extracted_read(mb, param,ri);
                                        ri->barcode =  (mem << 16) |   bar;
                                        //ri->barcode_string = param->read_structure->sequence_matrix[mem][bar];
                                        ri->read_type = EXTRACT_SUCCESS;
                                }else{
                                        ri->read_type  = EXTRACT_FAIL_BAR_FINGER_NOT_FOUND;
                                }

                        }else if(required_finger_len){
                                if(fingerlen == required_finger_len){
                                        ri = make_extracted_read(mb, param,ri);
                                        if(required_finger_len <= 255){
                                                ri->fingerprint = (key <<  8) | required_finger_len;
                                        }else{
                                                ri->fingerprint = (key <<  8) | 255;
                                        }


                                        ri->read_type = EXTRACT_SUCCESS;
                                }else{
                                        ri->read_type  = EXTRACT_FAIL_BAR_FINGER_NOT_FOUND;
                                }
                        }else{
                                ri = make_extracted_read(mb, param,ri);

                                ri->read_type = EXTRACT_SUCCESS;
                        }
                }else{
                        ri->read_type = EXTRACT_FAIL_READ_TOO_SHORT;
                }
        }else{
                //fprintf(stderr,"UN\n");
                //print_labelled_reads(mb, param,ri);
                ri->read_type = EXTRACT_FAIL_ARCHITECTURE_MISMATCH;

        }

        ri->qual[ri->len] = 0;
        free (buffer);
        return OK;
ERROR:
        return FAIL;
}


/** \fn  struct read_info* make_extracted_read(struct model_bag* mb, struct parameters* param,  struct read_info* ri)
    \brief Reformats the read sequence.

    \param mb The HMM model @ref model_bag.
    \param param @ref parameters .
    \param ri The reads.

*/

struct read_info* make_extracted_read(struct model_bag* mb, struct parameters* param,  struct read_info* ri)
{

        //print_labelled_reads(mb, param,ri);

        int s_pos,j,c1,c2;
        //int multireadread = 0;
        s_pos = 0;
        for(j = 0; j < ri->len;j++){
                c1 = mb->label[(int)ri->labels[j+1]];
                c2 = c1 & 0xFFFF; //which segment

                //fprintf(stderr,"%d	%c	\n", ri->seq[j],param->read_structure->type[c2]);

                if(param->read_structure->type[c2] == 'R'){
                        //if(multireadread == 0){
                        //	multireadread = 1;
                        //}
                        ri->seq[s_pos] = ri->seq[j];
                        ri->qual[s_pos] = ri->qual[j];
                        s_pos++;
                        //}else if (param->multiread){
                }else{
                        ri->seq[s_pos] = 65; // 65 is the spacer! nucleotides are 0 -5....
                        ri->qual[s_pos] = 65;
                        s_pos++;
                }
        }
        ri->len = s_pos;
        //exit(0);
        return ri;
}

void print_labelled_reads(struct model_bag* mb, struct parameters* param,  struct read_info* ri)
{
        int j,c1,c2;
        char alphabet[] = "ACGTNN";
        //int multireadread = 0;
        fprintf(stderr,"%f RQ\n", ri->mapq);

        for(j = 0; j < ri->len;j++){
                fprintf(stderr,"%c",alphabet[(int)ri->seq[j]]);
        }
        fprintf(stderr,"\n");

        for(j = 0; j < ri->len;j++){
                c1 = mb->label[(int)ri->labels[j+1]];
                c2 = c1 & 0xFFFF; //which segment

                fprintf(stderr,"%c",param->read_structure->type[c2]);

        }
        fprintf(stderr,"\n\n");

}





/** \fn void* do_baum_welch_thread(void *threadarg)
    \brief Runs Baum-Welch procedure.

    We implemented the HMM training procedure to verify that the forward and backward recursion work as expected. In other words this is only used for testing.

    \param threadarg  A @ref thread_data object used to pass data / parameters to function.
*/
void* do_baum_welch_thread(void *threadarg)
{
        struct thread_data *data;
        data = (struct thread_data *) threadarg;

        struct read_info** ri  = data->ri;
        struct model_bag* mb = data->mb;

        int start = data->start;
        int end = data->end;
        int i;
        float pbest;
        float Q;
        for(i = start; i < end;i++){
                mb = backward(mb, ri[i]->seq ,ri[i]->len);
                mb = forward_max_posterior_decoding(mb, ri[i], ri[i]->seq ,ri[i]->len);

                pbest = 1.0 - scaledprob2prob(  (ri[i]->bar_prob + mb->f_score) - logsum(mb->f_score,  mb->r_score));
                if(!pbest){
                        Q = 40.0;
                }else if(pbest == 1.0){
                        Q = 0.0;

                }else{
                        Q = -10.0 * log10(pbest) ;
                }
                if(Q >= 10){
                        mb = forward_extract_posteriors(mb, ri[i]->seq,ri[i]->labels,ri[i]->len);
                }
        }
        pthread_exit((void *) 0);
}












/** \fn struct model* malloc_model(int main_length, int sub_length, int number_sub_models)
    \brief Allocates a HMM segment.

    \param main_length Length of the first HMM.
    \param sub_length Length of second... N HMM.
    \param number_sub_models Number of HMMs.
    \deprecated Not used anymore.
*/


struct model* malloc_model(int main_length, int sub_length, int number_sub_models)
{
        struct model* model = NULL;
        int status;
        int i = 0;
        int j = 0;

        ASSERT(number_sub_models  <=MAX_NUM_SUB_MODELS,"Oh dear");

        MMALLOC(model,sizeof(struct model));

        model->num_hmms =  (1+ number_sub_models);
        MMALLOC(model->hmms,sizeof(struct hmm*) * (1+ number_sub_models));
        for(i = 0; i < model->num_hmms;i++){
                MMALLOC(model->hmms[i],sizeof(struct hmm) );
        }


        model->hmms[0]->num_columns = main_length;
        MMALLOC(model->hmms[0]->hmm_column,sizeof(struct hmm_column*) * main_length);

        for(j = 0; j < main_length;j++){
                MMALLOC(model->hmms[0]->hmm_column[j],sizeof(struct hmm_column));
        }

        for(i = 1; i < model->num_hmms;i++){
                model->hmms[i]->num_columns = sub_length;
                MMALLOC(model->hmms[i]->hmm_column,sizeof(struct hmm_column*) * sub_length);
                for(j = 0; j < sub_length;j++){
                        MMALLOC(model->hmms[i]->hmm_column[j],sizeof(struct hmm_column));
                }
        }



        return model;
ERROR:
        //KSLIB_MESSAGE(status,"Something wrong in malloc_model.\n");
        return NULL;
}



/** \fn struct model* malloc_model_according_to_read_structure(int num_hmm, int length)
    \brief Allocates a HMM segment.

    \param length Length of all HMMs.
    \param num_hmm Number of HMMs.
*/



/** \fn struct model* init_model_according_to_read_structure(struct model* model,struct parameters* param , int key, double* background,int assumed_length)

    \brief Initialized whole HMM.

    \param model The model to be initialized.
    \param param @ref parameters.
    \param key The key to the HMM type.
    \param background Background nucleotide probabilities.
    \param assumed_length Length of the HMM.
*/




/** \fn void print_model(struct model* model)

    \brief Prints out all parameters of a HMM.

    \param model The model to be printed out.

*/


void print_model(struct model* model)
{
        int i,j,c;
        int len;
        float sum = 0;
        struct hmm_column* col =0;
        //fprintf(stderr,"Skip:%f Self:%f Next:%f\n",scaledprob2prob(model->skip), scaledprob2prob(model->random_self) , scaledprob2prob(model->random_next));
        for(i = 0; i < model->num_hmms;i++){
                fprintf(stderr,"HMM%d:		(get there:%f)\n",i, scaledprob2prob(model->silent_to_M[i][0]));



                len = model->hmms[i]->num_columns;
                //tmp = rs->sequence_matrix[key][i];
                for(j = 0; j < len;j++){
                        col = model->hmms[i]->hmm_column[j];
                        fprintf(stderr,"\t%d\t",j);
                        sum = 0.0;
                        for(c = 0; c < 5;c++){
                                sum += scaledprob2prob(col->m_emit[c]);
                                fprintf(stderr,"%0.4f ", scaledprob2prob(col->m_emit[c]));
                        }
                        fprintf(stderr,"%0.1fc ",sum);

                        sum = 0.0;
                        for(c = 0; c < 5;c++){
                                sum += scaledprob2prob(col->i_emit[c]);
                                fprintf(stderr,"%0.4f ", scaledprob2prob(col->i_emit[c]));

                        }
                        fprintf(stderr,"%0.1fc ",sum);

                        for(c = 0; c < 9;c++){
                                fprintf(stderr,"%0.4f ", scaledprob2prob(col->transition[c]));

                        }
                        fprintf(stderr,"%0.4f %0.4f %0.4f\n",scaledprob2prob(col->transition[MSKIP])+scaledprob2prob(col->transition[MM])+scaledprob2prob(col->transition[MI])+scaledprob2prob(col->transition[MD])  ,scaledprob2prob(col->transition[ISKIP])+scaledprob2prob(col->transition[II])+scaledprob2prob(col->transition[IM]),  scaledprob2prob(col->transition[DD]) + scaledprob2prob(col->transition[DM]) );
                }
        }
        /*fprintf(stderr," ESTIMATED::::: \n");
          for(i = 0; i < model->num_hmms;i++){
          fprintf(stderr,"HMM%d:		(get there:%f)\n",i, scaledprob2prob(model->silent_to_M[i][0]));



          len = model->hmms[i]->num_columns;
          //tmp = rs->sequence_matrix[key][i];
          for(j = 0; j < len;j++){
          col = model->hmms[i]->hmm_column[j];
          fprintf(stderr,"\t%d\t",j);
          sum = 0.0;
          for(c = 0; c < 5;c++){
          sum += scaledprob2prob(col->m_emit_e[c]);
          fprintf(stderr,"%0.4f ", scaledprob2prob(col->m_emit_e[c]));

          }
          fprintf(stderr,"%0.1fc ",sum);

          sum = 0.0;
          for(c = 0; c < 5;c++){
          sum += scaledprob2prob(col->i_emit_e[c]);
          fprintf(stderr,"%0.4f ", scaledprob2prob(col->i_emit_e[c]));

          }
          fprintf(stderr,"%0.1fc ",sum);

          for(c = 0; c < 9;c++){
          fprintf(stderr,"%0.4f ", scaledprob2prob(col->transition_e[c]));

          }
          fprintf(stderr,"%0.3f %0.3f %0.3f\n",scaledprob2prob(col->transition_e[MM])+scaledprob2prob(col->transition_e[MI])+scaledprob2prob(col->transition_e[MD])  ,scaledprob2prob(col->transition_e[II])+scaledprob2prob(col->transition_e[IM]),  scaledprob2prob(col->transition_e[DD]) + scaledprob2prob(col->transition_e[DM]) );
          }
          }*/
        fprintf(stderr,"Links:silent to\n");
        double sumM = prob2scaledprob( 0.0);
        double sumI =  prob2scaledprob( 0.0);

        for(i = 0; i < model->num_hmms;i++){
                len = model->hmms[i]->num_columns;
                //tmp = rs->sequence_matrix[key][i];
                for(j = 0; j < len;j++){
                        fprintf(stderr,"%d	%f	%f\n",i, scaledprob2prob(  model->silent_to_M[i][j]), scaledprob2prob(  model->silent_to_I[i][j]));
                        sumM = logsum(sumM, model->silent_to_M[i][j]);
                        sumI = logsum(sumI, model->silent_to_I[i][j]);
                }
                fprintf(stderr,"SANITY: %f\t%f\n",scaledprob2prob(sumM),scaledprob2prob(sumI));
                //fprintf(stderr,"%d	%f	%f	%f	%f\n",i, scaledprob2prob(  model->silent_to_M[i][0]), scaledprob2prob(  model->silent_to_I[i][0]),scaledprob2prob(   model->silent_to_M_e[i][0]),scaledprob2prob(  model->silent_to_I_e[i][0]));
        }
        fprintf(stderr,"Links:to silent \n");


        fprintf(stderr,"SKIP:\n");
        fprintf(stderr,"%f	%f\n", scaledprob2prob(model->skip) , scaledprob2prob(model->skip_e));



}


/** \fn void free_model(struct model* model)

    \brief Frees model.

    \param model The model to be freed .

*/







/** \fn struct model_bag* copy_model_bag(struct model_bag* org)

    \brief Copies HMM into new HMM.
    Used to copy HMM for each thread.

    \param org  The @ref model_bag to be copied.

*/





/** \fn struct model_bag* set_model_e_to_laplace(struct model_bag* mb)

    \brief Sets all estiamted probabilities to 1.

    \param mb  The model.

*/

struct model_bag* set_model_e_to_laplace(struct model_bag* mb)
{
        int g,i,j,c;

        struct model* m = 0;
        struct hmm_column* hmm_column = 0;


        for(g = 0;g < mb->num_models;g++){
                m = mb->model[g];


                for(i = 0; i < m->num_hmms;i++){


                        for(j = 0; j < m->hmms[i]->num_columns;j++){
                                hmm_column = m->hmms[i]->hmm_column[j];
                                //copy->silent_to_I[i][j] = org->silent_to_I[i][j];
                                if(m->silent_to_I[i][j] != prob2scaledprob(0.0)){
                                        m->silent_to_I_e[i][j] = prob2scaledprob( 1.0f); //org->silent_to_I_e[i][j];
                                }
                                //copy->silent_to_M[i][j] = org->silent_to_M[i][j];
                                if(m->silent_to_M[i][j] != prob2scaledprob(0.0)){
                                        m->silent_to_M_e[i][j] = prob2scaledprob( 1.0f); //org->silent_to_M_e[i][j];
                                }

                                for(c = 0; c< 5;c++){
                                        //copy_col->i_emit[c] = org_col->i_emit[c];
                                        hmm_column->i_emit_e[c] = prob2scaledprob(1.0f);

                                        //copy_col->m_emit[c] = org_col->m_emit[c];
                                        hmm_column->m_emit_e[c] = prob2scaledprob(1.0f);
                                }

                                for(c = 0; c < 9;c++){
                                        //c/opy_col->transition[c] = org_col->transition[c];
                                        if(hmm_column->transition[c] != prob2scaledprob(0.0)){
                                                hmm_column->transition_e[c] = prob2scaledprob(1.0f);
                                        }
                                }
                        }
                }
                //copy->skip = org->skip;
                if(m->skip != prob2scaledprob(0.0)){
                        m->skip_e = prob2scaledprob( 1.0f);//?//?org->skip_e;
                }
                //copy->num_hmms = org->num_hmms;
        }
        return mb;
}




/** \fn struct model* reestimate(struct model* m, int mode)

    \brief Sets new model parameters based on estiamted probabilities from data.

    \param mb  The model.
    \param mode Sets which parameters should be optimized.

*/






/** \fn struct model* copy_estimated_parameter(struct model* target, struct model* source )

    \brief Sums eatimated parameters from different models.
    Used to merge results from multiple threads.

    \param target  Model to hold the sums.
    \param source Model from which to copy parameters.

*/

struct model* copy_estimated_parameter(struct model* target, struct model* source )
{
        int i,j,c;

        struct hmm_column* target_col = 0;
        struct hmm_column* source_col = 0;




        for(i = 0; i < target->num_hmms;i++){

                //copy->I_to_silent[i] = org->I_to_silent[i];
                //target->I_to_silent_e[i] = logsum(target->I_to_silent_e[i] , source->I_to_silent_e[i]);//   org->I_to_silent_e[i];
                //copy->M_to_silent[i] = org->M_to_silent[i];
                //target->M_to_silent_e[i] = logsum(target->M_to_silent_e[i],source->M_to_silent_e[i]); // org->M_to_silent_e[i];

                for(j = 0; j < target->hmms[0]->num_columns;j++){

                        //copy->silent_to_I[i] = org->silent_to_I[i];
                        target->silent_to_I_e[i][j] = logsum(target->silent_to_I_e[i][j], source->silent_to_I_e[i][j]);
                        //copy->silent_to_M[i] = org->silent_to_M[i];
                        target->silent_to_M_e[i][j] = logsum(target->silent_to_M_e[i][j] ,source->silent_to_M_e[i][j]);//  org->silent_to_M_e[i];


                        target_col = target->hmms[i]->hmm_column[j];
                        source_col = source->hmms[i]->hmm_column[j];
                        for(c = 0; c< 5;c++){
                                //copy_col->i_emit[c] = org_col->i_emit[c];
                                target_col->i_emit_e[c] = logsum(target_col->i_emit_e[c] , source_col->i_emit_e[c] );  //org_col->i_emit_e[c];

                                //copy_col->m_emit[c] = org_col->m_emit[c];
                                target_col->m_emit_e[c] = logsum(target_col->m_emit_e[c], source_col->m_emit_e[c] );// org_col->m_emit_e[c];


                        }

                        for(c = 0; c < 9;c++){
                                //copy_col->transition[c] = org_col->transition[c];
                                target_col->transition_e[c] = logsum (target_col->transition_e[c], source_col->transition_e[c] ); // org_col->transition_e[c];

                        }
                }
        }
        //copy->skip = org->skip;
        target->skip_e = logsum(target->skip_e, source->skip_e );// copy->skip_e;



        return target;
}





/** \fn struct model_bag* init_model_bag(struct parameters* param,struct sequence_stats_info* ssi )

    \brief Initializes whole HMM model based on user input.


    \param param @ref parameters   Model to hold the sums.
    \param ssi sequence information including backgound nucleotide probabilities.

*/
