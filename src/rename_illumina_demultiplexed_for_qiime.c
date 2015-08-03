#if HAVE_CONFIG_H
#include "config.h"
#endif

#include "kslib.h"

#include "interface.h"
#include "nuc_code.h"
#include "io.h"
#include "misc.h"
#include "tagdust2.h"


struct sample_name_barcode_mapping{
	char* barcode;
	char* sample_name;
	int num;
};

int barcode_compare (const void *a, const void *b);

int binsearch(struct sample_name_barcode_mapping** sample_mapping, int num_entries, char* query);

int main (int argc,char * argv[]) {
	struct parameters* param = NULL;
	struct read_info** ri =NULL ;
	int status;
	int i,j,f,g,h;
	int num_samples = 0;
	init_nuc_code();
	char alphabet[] = "ACGTNN";
	char line[MAX_LINE];
	char c;
	param = interface(param,argc,argv);
	FILE* map_file = NULL;
	FILE* fastQ = NULL;
	int numseq;
	int  error;
	int min_error;
	
	struct sample_name_barcode_mapping** sample_mapping = NULL;
	int (*fp)(struct read_info** ,struct parameters*,FILE* , int* buffer_count) = 0;
	char* query = NULL;
	//char* tmp = NULL;
	char*new_name = NULL;
	fp = &read_fasta_fastq;
	
	if(param){
		if(param->infiles < 2){
			usage();
		}
		if((map_file = fopen(param->infile[0], "r")) == NULL) KSLIB_XEXCEPTION_SYS(kslEWRT,"Failed to open file:%s",param->buffer);
		num_samples = 0;
		while(fgets(line, MAX_LINE, map_file)){
			if(line[0] != '#'){
				j = 0;
				for(i = 0; i < MAX_LINE;i++){
					c = line[i];
					if(!c){
						break;
					}
					if( isspace((int)c)){
						j++;
					}else if(j == 1){
						if(nuc_code[(int)c] > 3){
							KSLIB_XEXCEPTION_SYS(kslEWRT,"Non Nucleotide letters in barcode:\n%s\n",line );
						}
					}
				}
				num_samples++;
				
			}
		}
		rewind(map_file);
		KSL_DPRINTF3(("%d samples found \n", num_samples));
		
		MMALLOC(sample_mapping, sizeof(struct sample_name_barcode_mapping*) * num_samples);
		for(i = 0; i < num_samples;i++){
			sample_mapping[i] = NULL;
			MMALLOC(sample_mapping[i], sizeof(struct sample_name_barcode_mapping ));
			sample_mapping[i]->barcode = NULL;
			sample_mapping[i]->sample_name = NULL;
			sample_mapping[i]->num = 0;
		}
		num_samples = 0;
		while(fgets(line, MAX_LINE, map_file)){
			if(line[0] != '#'){
				f = 0;
				j = 0;
				for(i = 0; i < MAX_LINE;i++){
					c = line[i];
					if(!c){
						KSL_DPRINTF3(("Quitting: %d (%d)\n",i,num_samples ));
						break;
					}
					KSL_DPRINTF3(("%d %c\n",i,c ));
					if(isspace((int)c)){
						j++;
						KSL_DPRINTF3(("Space found: %d\n",i ));
						if(j ==1){
							MMALLOC(sample_mapping[num_samples]->sample_name, sizeof(char)* ((i - f)+1) );
							for(g = f;g < i;g++) {
								sample_mapping[num_samples]->sample_name[g] = line[g];
							}
							sample_mapping[num_samples]->sample_name[g]  = 0;
							KSL_DPRINTF3(("Rrading:%d-%d\n",f,i ));
						}
						if(j ==2){
							MMALLOC(sample_mapping[num_samples]->barcode, sizeof(char)* ((i - f)+1));
							for(g = f;g < i;g++) {
								sample_mapping[num_samples]->barcode[g-f] = line[g];
							}
							sample_mapping[num_samples]->barcode[g-f]  = 0;
							KSL_DPRINTF3(("Rrading:%d-%d\n",f,i ));
						}
						f = i +1;
					}
				}
				num_samples++;
			}
		}
		fclose(map_file);
		
		qsort(sample_mapping,num_samples, sizeof(struct sample_name_barcode_mapping*), barcode_compare);

		
		for(i = 0;i < num_samples;i++){
			fprintf(stderr,"%d\t%s\t%s\n",i , sample_mapping[i]->sample_name,sample_mapping[i]->barcode);
		}
		MMALLOC(query, sizeof(char) * 1000);
		MMALLOC(new_name, sizeof(char) * 10000);
		
		
		fastQ = io_handler(fastQ , 1,param);
		//param->num_query = 10;
		ri = malloc_read_info(ri, param->num_query);
		while(1){
			if((status = fp(ri, param,fastQ,&numseq)) != kslOK)  exit(status);
			if(!numseq){
				break;
			}
			//while ((numseq = fp(ri, param,file)) != 0){
			for(j = 0;j < numseq;j++){
				g = (int)strlen(ri[j]->name);
				f = 0;
				for(i = 0; i < g;i++){
					c = ri[j]->name[i];
					if(c == ';' && f != 0){
						query[f] = 0;
						break;
					}else if(nuc_code[(int)c] <=3){
						query[f] = c;
						f++;
						if(f == 999){
							query[f] = 0;
							break;
						}
					}else{
						f = 0;
					}
				}
				//fprintf(stderr,"%s\n%s\n",ri[j]->name,query);
				f = binsearch(sample_mapping, num_samples, query);
				min_error = 0;
				if(f == -1){
					if(strlen(query) == strlen(sample_mapping[0]->barcode)){
						f = 0;
						min_error = 1000;
						for(i = 0;i < num_samples;i++){
							error = 0;
							for(g = 0; g < strlen(query);g++){
								if(query[g] != sample_mapping[i]->barcode[g]){
									error++;
								}
							}
							if(error < min_error){
								min_error = error;
								f = i;
							}
						}
					}

				}
				if(f != -1){
			
					h = 0;
					for(i = 0; i < g;i++){
						c = ri[j]->name[i];
						if(isalnum((int) c)){
							query[h] = c;
							h++;
						}
						if(isspace((int)c)){
							query[h] = 0;
							break;
						}
						
					}
					sample_mapping[f]->num++;
					snprintf(new_name, 10000,">%s_%d %s orig_bc=%s new_bc=%s bc_diffs=%d",sample_mapping[f]->sample_name,sample_mapping[f]->num, query,sample_mapping[f]->barcode,sample_mapping[f]->barcode ,min_error);
					for(i = 0; i < ri[j]->len ;i++){
						ri[j]->seq[i] = alphabet[(int) ri[j]->seq[i]];
					}
					
					fprintf(stdout,"%s\n%s\n", new_name, ri[j]->seq );
					
				}else{
					fprintf(stderr,"Warning: no barcode match for found for:\n%s\n",ri[j]->name);
				}
				//PC.634_1 FLP3FBN01ELBSX orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
			}
		}
		pclose(fastQ);
		
		for(i = 0;i < num_samples;i++){
			fprintf(stderr,"%d\t%d\t%s\t%s\n",i , sample_mapping[i]->num,sample_mapping[i]->sample_name,sample_mapping[i]->barcode);
		}
		
		MFREE(new_name);
		MFREE(query);
		for(i = 0; i < num_samples;i++){
			MFREE(sample_mapping[i]->sample_name);
			MFREE(sample_mapping[i]->barcode);
			MFREE(sample_mapping[i]);//, sizeof(struct sample_name_barcode_mapping ));
		}
		MFREE(sample_mapping);//, sizeof(struct sample_name_barcode_mapping*) * num_samples);
		
		
		
		free_read_info(ri, param->num_query);
		free_param(param);
		//merge(param);
	}
	
	return EXIT_SUCCESS;
ERROR:
	return EXIT_FAILURE;
	
}

int binsearch(struct sample_name_barcode_mapping** sample_mapping, int num_entries, char* query)
{
	int start, end,mid;
	int c;
	start = 0;
	end = num_entries -1;
	mid = (start+end)/2;
 
	while (start <= end) {
		
		c = strcmp(sample_mapping[mid]->barcode,query );
		KSL_DPRINTF3(("%s <-> %s %d\n",query,  sample_mapping[mid]->barcode,c ));
		if (c < 0 ){
			start = mid + 1;
		}else if (c == 0) {
			KSL_DPRINTF3(("FOUND!"));
			return mid;
		}else{
			end = mid - 1;
		}
		mid = (start+end)/2;
	}
	return -1;
}


int barcode_compare (const void *a, const void *b)
{
	
	//struct mys **a = (struct mys **)i1;
	//struct mys **b = (struct mys **)i2;
	//return (*b)->id - (*a)->id;
	
	const struct sample_name_barcode_mapping **elem1 = (const struct sample_name_barcode_mapping**) a;
	
	const struct sample_name_barcode_mapping **elem2 = (const struct sample_name_barcode_mapping**) b;
	
	return strcmp((*elem1)->barcode, (*elem2)-> barcode);
}




