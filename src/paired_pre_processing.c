#include "interface.h"
#include "io.h"
#include "misc.h"
#include  "nuc_code.h"

#include <ctype.h>

#include<sys/stat.h>

#define RED 1
#define BLACK 0


struct rb_node{
	struct rb_node* right;
	struct rb_node* left;
	//int value;
	void* p ;
	long int key;
	//unsigned int start;
	//unsigned int end;
	int color;
	//unsigned int strand:1;
	unsigned int num:31;
};

struct rb_node* insert_start(struct rb_node* n,long int key,void* p) ;
struct rb_node* insert (struct rb_node* n,long int key, void* p);

int rank(struct rb_node* n);

struct rb_node* set_num(struct rb_node* n);
int size(struct rb_node* n);

int isred(struct rb_node* n);
struct rb_node* rotateLeft(struct rb_node* n);
struct rb_node* rotateRight(struct rb_node* n);
struct rb_node* flipColors(struct rb_node* n) ;
void free_rbtree(struct rb_node* n);
void print_split_sequences(struct rb_node* n,char* out,struct read_info** ri, struct parameters* param,int (*fp)(struct read_info** ,struct parameters*,FILE* ) );

void concatenate_reads(struct parameters* param,int (*fp)(struct read_info** ,struct parameters*,FILE* ))
{
	int i = 0;
	int j = 0;
	param->sam = 0;
	if(!strcmp(".sam", param->infile[i] + (strlen(param->infile[i] ) - 4))){
		param->sam = 1;
	}else if (!strcmp(".bam", param->infile[i] + (strlen(param->infile[i] ) - 4))){
		param->sam = 2;
	}else if (!strcmp(".fa", param->infile[i] + (strlen(param->infile[i] ) - 3))){
		param->sam = 0;
		param->fasta = 1;
	}else if (!strcmp(".fq", param->infile[i] + (strlen(param->infile[i] ) - 3))){
		param->sam = 0;
	}else if (!strcmp(".fastq", param->infile[i] + (strlen(param->infile[i] ) - 6))){
		param->sam = 0;
	}else if (!strcmp(".fastaq", param->infile[i] + (strlen(param->infile[i] ) - 7))){
		param->sam = 0;
	}else if (!strcmp(".fasta", param->infile[i] + (strlen(param->infile[i] ) - 6))){
		param->sam = 0;
		param->fasta = 1;
	}else if(!strcmp(".sam.gz", param->infile[i] + (strlen(param->infile[i] ) - 7))){
		param->sam = 1;
		param->gzipped  = 1;
	}else if (!strcmp(".bam.gz", param->infile[i] + (strlen(param->infile[i] ) - 7))){
		param->sam = 2;
		param->gzipped  = 1;
	}else if (!strcmp(".fa.gz", param->infile[i] + (strlen(param->infile[i] ) - 6))){
		param->sam = 0;
		param->fasta = 1;
		param->gzipped  = 1;
	}else if (!strcmp(".fq.gz", param->infile[i] + (strlen(param->infile[i] ) - 6))){
		param->sam = 0;
		param->gzipped  = 1;
	}else if (!strcmp(".fastq.gz", param->infile[i] + (strlen(param->infile[i] ) - 9))){
		param->sam = 0;
		param->gzipped  = 1;
	}else if (!strcmp(".fastaq.gz", param->infile[i] + (strlen(param->infile[i] ) - 10))){
		param->sam = 0;
		param->gzipped  = 1;
	}else if (!strcmp(".fasta.gz", param->infile[i] + (strlen(param->infile[i] ) - 9))){
		param->sam = 0;
		param->gzipped  = 1;
	}else if (!strcmp(".fastq.bz2", param->infile[i] + (strlen(param->infile[i] ) - 10))){
		param->sam = 0;
		param->bzipped  = 1;
	}else if (!strcmp(".fq.bz2", param->infile[i] + (strlen(param->infile[i] ) - 7))){
		param->sam = 0;
		param->bzipped  = 1;
	}else{
		param->sam = -1;
	}
	
	struct read_info** ri1 = 0;
	struct read_info** ri2 = 0;
	
	FILE* file1 = 0;
	FILE* file2 = 0;

	
	ri1 = malloc(sizeof(struct read_info*) * param->num_query);
	ri2 = malloc(sizeof(struct read_info*) * param->num_query);
	
	//assert(ri1 !=0);
	
	for(i = 0; i < param->num_query;i++){
		ri1[i] = malloc(sizeof(struct read_info));
		ri1[i]->seq = 0;
		ri1[i]->name = 0;
		ri1[i]->qual = 0;
		ri1[i]->labels = 0;
		ri1[i]->len = 0;
		ri1[i]->cigar = 0;
		ri1[i]->bar_prob = 0;
		ri1[i]->md = 0;
		ri1[i]->mapq = -1.0;
		ri1[i]->strand = malloc(sizeof(unsigned int)* (LIST_STORE_SIZE+1));
		ri1[i]->hits = malloc(sizeof(unsigned int)* (LIST_STORE_SIZE+1));
		
		ri2[i] = malloc(sizeof(struct read_info));
		ri2[i]->seq = 0;
		ri2[i]->name = 0;
		ri2[i]->qual = 0;
		ri2[i]->labels = 0;
		ri2[i]->len = 0;
		ri2[i]->cigar = 0;
		ri2[i]->bar_prob = 0;
		ri2[i]->md = 0;
		ri2[i]->mapq = -1.0;
		ri2[i]->strand = malloc(sizeof(unsigned int)* (LIST_STORE_SIZE+1));
		ri2[i]->hits = malloc(sizeof(unsigned int)* (LIST_STORE_SIZE+1));
		
		
		
	}
	file1 =  io_handler(file1, 0,param);
	
	file2 =  io_handler(file2, 1,param);
	
	int numseq1,numseq2;
	
	char* spacer = 0;
	char* barcode=0;
	for(i = 0; i < param->read_structure->num_segments;i++){
		if(param->read_structure->type[i] == 'B'){
			barcode = param->read_structure->sequence_matrix[i][0];
		}
		
		if(param->read_structure->type[i] == 'S'){
			spacer = param->read_structure->sequence_matrix[i][0];
		}

	}
	char alpha[6] = "ACGTNN";
	FILE* outfile = 0;
	
	if(param->outfile){
		if ((outfile = fopen( param->outfile, "w")) == NULL){
			sprintf(param->buffer,"can't open output\n");
			fprintf(stderr,"%s",param->buffer);
			param->messages = append_message(param->messages, param->buffer);
			free_param(param);
			exit(EXIT_FAILURE);
			
			//fprintf(stderr,"can't open output\n");
			//exit(-1);
		}
	}else{
		outfile= stdout;
	}
	
	
	while ((numseq1 = fp(ri1, param,file1)) != 0){
		numseq2 = fp(ri2, param,file2);
		if(numseq1 != numseq2){
			sprintf(param->buffer,"Two files seem to be of different length.\n");
			fprintf(stderr,"%s",param->buffer);
			param->messages = append_message(param->messages, param->buffer);
			free_param(param);
			exit(EXIT_FAILURE);
			
		}
		for(i = 0; i < numseq1;i++){
			fprintf(outfile,"@");
			for(j = 0; j < strlen(ri1[i]->name);j++){
				if(isspace(ri1[i]->name[j])){
					ri1[i]->name[j] = 0;
					ri2[i]->name[j] = 0;
				}
				
			}
			if(strcmp(ri1[i]->name,ri2[i]->name)){
				sprintf(param->buffer,"Files seem to contain reads in different order:\n%s\n%s\n",ri1[i]->name,ri2[i]->name );
				fprintf(stderr,"%s",param->buffer);
				param->messages = append_message(param->messages, param->buffer);
				free_param(param);
				exit(EXIT_FAILURE);
			}
			
			fprintf(outfile,"%s",ri1[i]->name);
			
			fprintf(outfile,"\n");
			
			if(barcode){
				for(j = 0; j < strlen(barcode);j++){
					fprintf(outfile,"%c", alpha[ nuc_code[(int)  barcode[j]]]);
				}
			}
			
			for(j = 0; j < ri1[i]->len;j++){
				fprintf(outfile,"%c", alpha[(int)  ri1[i]->seq[j]]);
				
			}
			if(spacer){
				for(j = 0; j < strlen(spacer);j++){
					fprintf(outfile,"%c", alpha[nuc_code[(int) spacer[j]]]);
				}
			}
			for(j = 0; j < ri2[i]->len;j++){
				fprintf(outfile,"%c", alpha[(int) ri2[i]->seq[j]]);
				
			}
			fprintf(outfile,"\n");
			if(ri1[i]->qual){
				fprintf(outfile,"+\n");
				if(barcode){
					for(j = 0; j < strlen(barcode);j++){
						fprintf(outfile,"%c", alpha[nuc_code[(int) barcode[j]]]);
					}
				}
				for(j = 0; j < ri1[i]->len;j++){
					fprintf(outfile,"%c", ri1[i]->qual[j]);
					
				}
				
				if(spacer){
					for(j = 0; j < strlen(spacer);j++){
						fprintf(outfile,"%c", alpha[nuc_code[(int) spacer[j]]]);
					}
				}
				for(j = 0; j < ri2[i]->len;j++){
					fprintf(outfile,"%c", ri2[i]->qual[j]);
					
				}
				fprintf(outfile,"\n");
			}
		}
	}
	
	
	
	for(i = 0; i < param->num_query;i++){
		free(ri1[i]->strand);
		free(ri1[i]->hits);
		
		free(ri2[i]->strand);
		free(ri2[i]->hits);
		
		if(ri1[i]->cigar){
			free(ri1[i]->cigar);
		}
		if(ri2[i]->cigar){
			free(ri2[i]->cigar);
		}
		
		if(ri1[i]->md){
			free(ri1[i]->md);
		}
		if(ri2[i]->md){
			free(ri2[i]->md);
		}
		
		if(ri1[i]->name){
			free(ri1[i]->name);
		}
		if(ri2[i]->name){
			free(ri2[i]->name);
		}

		if(ri1[i]->seq){
			free(ri1[i]->seq);
		}
		if(ri2[i]->seq){
			free(ri2[i]->seq);
		}
		
		
		if(ri1[i]->qual){
			free(ri1[i]->qual );
		}
		
		if(ri2[i]->qual){
			free(ri2[i]->qual );
		}
		
		free(ri1[i]);
		free(ri2[i]);
	}
	free(ri1);
	free(ri2);
	
	if(param->sam == 2 || param->sam == 1 || param->gzipped ){
		pclose(file1);
		pclose(file2);

	}else{
		fclose(file1);
		fclose(file2);
	}
	if(param->outfile){
		fclose(outfile);
	}
}


void split(struct parameters* param,int (*fp)(struct read_info** ,struct parameters*,FILE* ))
{
	int i = 0;
	int j = 0;
	if(!param->outfile){
		sprintf(param->buffer,"You need to specify an output prefix.\n");
		fprintf(stderr,"%s",param->buffer);
		param->messages = append_message(param->messages, param->buffer);
		free_param(param);
		exit(EXIT_FAILURE);
		
	}
	
	
	param->sam = 0;
	if(!strcmp(".sam", param->infile[i] + (strlen(param->infile[i] ) - 4))){
		param->sam = 1;
	}else if (!strcmp(".bam", param->infile[i] + (strlen(param->infile[i] ) - 4))){
		param->sam = 2;
	}else if (!strcmp(".fa", param->infile[i] + (strlen(param->infile[i] ) - 3))){
		param->sam = 0;
		param->fasta = 1;
	}else if (!strcmp(".fq", param->infile[i] + (strlen(param->infile[i] ) - 3))){
		param->sam = 0;
	}else if (!strcmp(".fastq", param->infile[i] + (strlen(param->infile[i] ) - 6))){
		param->sam = 0;
	}else if (!strcmp(".fastaq", param->infile[i] + (strlen(param->infile[i] ) - 7))){
		param->sam = 0;
	}else if (!strcmp(".fasta", param->infile[i] + (strlen(param->infile[i] ) - 6))){
		param->sam = 0;
		param->fasta = 1;
	}else if(!strcmp(".sam.gz", param->infile[i] + (strlen(param->infile[i] ) - 7))){
		param->sam = 1;
		param->gzipped  = 1;
	}else if (!strcmp(".bam.gz", param->infile[i] + (strlen(param->infile[i] ) - 7))){
		param->sam = 2;
		param->gzipped  = 1;
	}else if (!strcmp(".fa.gz", param->infile[i] + (strlen(param->infile[i] ) - 6))){
		param->sam = 0;
		param->fasta = 1;
		param->gzipped  = 1;
	}else if (!strcmp(".fq.gz", param->infile[i] + (strlen(param->infile[i] ) - 6))){
		param->sam = 0;
		param->gzipped  = 1;
	}else if (!strcmp(".fastq.gz", param->infile[i] + (strlen(param->infile[i] ) - 9))){
		param->sam = 0;
		param->gzipped  = 1;
	}else if (!strcmp(".fastaq.gz", param->infile[i] + (strlen(param->infile[i] ) - 10))){
		param->sam = 0;
		param->gzipped  = 1;
	}else if (!strcmp(".fasta.gz", param->infile[i] + (strlen(param->infile[i] ) - 9))){
		param->sam = 0;
		param->gzipped  = 1;
	}else if (!strcmp(".fastq.bz2", param->infile[i] + (strlen(param->infile[i] ) - 10))){
		param->sam = 0;
		param->bzipped  = 1;
	}else if (!strcmp(".fq.bz2", param->infile[i] + (strlen(param->infile[i] ) - 7))){
		param->sam = 0;
		param->bzipped  = 1;
	}else{
		param->sam = -1;
	}
	
	struct read_info** ri = 0;
	
	FILE* file = 0;
	
	ri = malloc(sizeof(struct read_info*) * param->num_query);
	
	for(i = 0; i < param->num_query;i++){
		ri[i] = malloc(sizeof(struct read_info));
		ri[i]->seq = 0;
		ri[i]->name = 0;
		ri[i]->qual = 0;
		ri[i]->labels = 0;
		ri[i]->len = 0;
		ri[i]->cigar = 0;
		ri[i]->bar_prob = 0;
		ri[i]->md = 0;
		ri[i]->mapq = -1.0;
		ri[i]->strand = malloc(sizeof(unsigned int)* (LIST_STORE_SIZE+1));
		ri[i]->hits = malloc(sizeof(unsigned int)* (LIST_STORE_SIZE+1));
				
		
	}
	file =  io_handler(file, 0,param);
	struct rb_node* root = 0;
	
	
	
	//file2 =  io_handler(file2, 1,param);
	
	
	int numseq;
	
	long int key = 0;
	int c;
	int g;
	char* tmp = 0;
	j = 0;
	param->join = 0;
	while ((numseq = fp(ri, param,file)) != 0){
		for(i = 0;i < numseq;i++){
			key = 0;
			tmp = malloc(sizeof(char)*100);
			j =byg_end("BC:", ri[i]->name);
			g = 0;
			if(j){
				param->join |= 1;
				tmp[0]= 'B';
				tmp[1] = 'C';
				tmp[2] = ':';
				g = 3;
				for(c = j;c < strlen(ri[i]->name );c++){
					if(ri[i]->name[c] == ';'){
						tmp[g] = ri[i]->name[c];
						g++;
						break;
					}
					tmp[g] =ri[i]->name[c];
					g++;
					key = (key << 2l) | (long int) nuc_code[(int)ri[i]->name[c]];
				}
			}
			key = key << 8;
			j = byg_end("RS:", ri[i]->name);
			if(j){
				param->join |= 2;
				key = key | atoi(ri[i]->name+j);
				tmp[g]= 'R';
				g++;
				tmp[g] = 'S';
				g++;
				tmp[g] = ':';
				g++;
				for(c = j;c < strlen(ri[i]->name );c++){
					if(ri[i]->name[c] == ';'){
						tmp[g] = ri[i]->name[c];
						g++;
						break;
					}
					tmp[g] =ri[i]->name[c];
					g++;
					
				}
			}
			tmp[g] = 0;
			root = insert_start(root, key, tmp );
			
			
			
		}
	}
	
	
	
	if(param->sam == 2 || param->sam == 1 || param->gzipped ){
		pclose(file);
	}else{
		fclose(file);
	}
	
	
	
	print_split_sequences(root,param->outfile, ri,param,fp);
	free_rbtree(root);
	for(i = 0; i < param->num_query;i++){
		free(ri[i]->strand);
		free(ri[i]->hits);
		
		if(ri[i]->cigar){
			free(ri[i]->cigar);
		}
		if(ri[i]->md){
			free(ri[i]->md);
		}
		if(ri[i]->name){
			free(ri[i]->name);
		}
		if(ri[i]->seq){
			free(ri[i]->seq);
		}
		if(ri[i]->qual){
			free(ri[i]->qual );
		}
		if(ri[i]->labels){
			free(ri[i]->labels);
		}
		free(ri[i]);
	}
	free(ri);
		
}




int file_exists (char * name)
{
	struct stat buf;
	int ret,local_ret;
	ret = 0;
	local_ret= stat ( name, &buf );
	/* File found */
	if ( local_ret == 0 )
	{
		ret++;
		//return 1;
	}
	return ret;
}




void print_split_sequences(struct rb_node* n,char* out,struct read_info** ri, struct parameters* param,int (*fp)(struct read_info** ,struct parameters*,FILE* ) )
{
	
	
	if(n->left){
		print_split_sequences(n->left,out,ri,param,fp);
	}
	
	char* filename = 0;
	filename = malloc(sizeof(char) * 1000);
	filename[0] = 0;
	
	sprintf (filename, "%s_%s.fq" , out,(char*)n->p);
	int i;
	for(i = 0;i < strlen(filename );i++){
		if(filename[i] == ';'){
			filename[i] = '_';
		}
		if(filename[i] == ':'){
			filename[i] = '_';
		}
	}
	
	int go = 0;
	
	if(param->join == 3){ // both double reads and barcode
		if(byg_end("BC", filename)){
			if(byg_end("RS", filename)){
				go =1;
			}
		}
	}
	
	if(param->join == 1){ // just barcodes
		if(byg_end("BC", filename)){
			go = 1;
		}
	}
		
	
	if(param->join == 2){ // double reads...
		if(byg_end("RS", filename)){
			go =1;
		}
	}
	
	if(!go){
		sprintf (filename, "%s_%s.fq" , out,"not_extracted");
	}
	fprintf(stderr,"%s\n",filename );
	
	FILE* outfile = 0;
	FILE* infile = 0;
	if ((outfile = fopen( filename, "w")) == NULL){
		sprintf(param->buffer,"can't open output\n");
		fprintf(stderr,"%s",param->buffer);
		param->messages = append_message(param->messages, param->buffer);
		free_param(param);
		exit(EXIT_FAILURE);
	}
	int numseq;
	infile =  io_handler(infile,0,param);
	
	while ((numseq = fp(ri, param,infile)) != 0){
		for(i = 0;i < numseq;i++){
			if(byg_end((char*)n->p, ri[i]->name)){
				print_seq(ri[i],outfile);
			}
		}
	}
	
	fclose(outfile);
	
	if(param->sam == 2 || param->sam == 1 || param->gzipped ){
		pclose(infile);
	}else{
		fclose(infile);
	}
	
	
	
	free(filename);
	
	//fprintf(stderr,"%s\n", (char*)n->p);
	if(n->right){
		print_split_sequences(n->right,out,ri,param,fp);
	}
	//free(n);
}



struct rb_node* insert_start(struct rb_node* n,long int key,void* p)
{
	n = insert(n, key,p);
	n->color = BLACK;
	return n;
}

struct rb_node* insert (struct rb_node* n,long int key,void*p)
{
	if(n == 0){
		n = malloc(sizeof(struct rb_node));
		n->left = 0;
		n->right = 0;
		n->color = RED;
		n->p = p;
		n->key = key;
		n->num = 1;
		return n;
	}
	
	if(isred(n->left)&&isred(n->right)){
		n = flipColors(n);
	}
	
	if(key == n->key){
		free(p);
		//n->value = n->value  + value;
		//n->right = insert(n->right,key,p);
	}else if(key < n->key){
		n->left = insert(n->left,key,p);
	}else if(key >  n->key){
		n->right = insert(n->right,key,p);
	}
	
	if(isred(n->right) && !isred(n->left)){
		n = rotateLeft(n);
	}
	
	if(isred(n->left)){
		if(isred(n->left->left)){
			n = rotateRight(n);
		}
	}
	return set_num(n);
}



int rank(struct rb_node* n)
{
	if(n->left == 0){
		return 0;
	}
	return n->left->num;
}

struct rb_node* set_num(struct rb_node* n)
{
	n->num = size(n->left) + size(n->right) + 1;
	return n;
}

int size(struct rb_node* n)
{
	if (n == 0){
		return 0;
	}
	return n->num;
}


int isred(struct rb_node* n)
{
	if(!n){
		return 0;
	}
	return n->color == RED;
}

struct rb_node* rotateLeft(struct rb_node* n)
{
	struct rb_node* x;
	x = n->right;
	n->right = x->left;
	x->left = set_num(n);
	x->color = n->color;
	n->color = RED;
	return set_num(x);
}

struct rb_node* rotateRight(struct rb_node* n)
{
	struct rb_node* x;
	x = n->left;
	n->left = x->right;
	x->right = set_num(n);
	x->color = n->color;
	n->color = RED;
	return set_num(x);
}

struct rb_node* flipColors(struct rb_node* n)
{
	n->color = !(n->color);
	n->left->color = !(n->left->color);
	n->right->color = !(n->right->color);
	return n;
}

void free_rbtree(struct rb_node* n)
{
	if(n->left){
		free_rbtree(n->left);
	}
	free(n->p);
	//fprintf(stderr,"%s\n", (char*)n->p);
	if(n->right){
		free_rbtree(n->right);
	}
	free(n);
}




