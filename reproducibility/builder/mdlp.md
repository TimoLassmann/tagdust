% MDLP
% Timo Lassmann
% Today

# Introduction

__MDLP__  is a __C__ program to extract code from markdown files and stitch them together into compilable files. It is largely inspired by ocean[^1].


Literal programming

The entire program is written into one __C__ file [nstart.c](#file-nstart.c)

How to distribute?

 stand alonce C program
  then this has to be included in all future code...
 make an additional directory in default gnu style directory tree.
  A makefile will compile C code then check if it works by running it on itself.
 this resulting executable will then be used to extract code from the program markdown file..


[^1]: See http://ocean-lang.org


# Usage 

## MDLP

 Simply type: mdlp  markdown_file.md


## Pandoc 

	pandoc  -s -S   mdlp.md  --toc -o test.tex --template=template.tex  --listings --chapter


# The Program

There are three bits in MDLP that warrent some more detailed description:

1. data structures to keep track of code snippets
2. scanning of the input markdown file.
3. code to print snippets in the right order into files


## Data Structures

The input markdown file is read into one string stored in the `md_doc` structure. A red-black tree is used to index the various code snippets contained in the file. Each snippet has a name (the nearest heading above the code) and a linked list to store the actual code. Elements in the linked list contain a pointer to the copy of the file in memory and store the line numbers in the markdown file.

~~~~{.c}

struct line{
	char* content;
	int line_number;
	struct line* next;
	int len;
};

struct code_block{
	struct line** lines;
	char* name;
	int num_lines;
};

struct md_doc{
	char* file_string;
	struct rb_node* root;
};


~~~~


## Parsing

Markdown turns out to be pretty hard to parse in _C_. After three attepts I use a finite state machine consisting of two states: _code_ and _default markdown_. The program reads the contents of the input file character by character. At the beginning of new lines the code checks if there is some markdown syntax. Based on this variables are set defining the transition to the next state.

While parsing the file the latest header and line numbers are remembered. The name of the last header is used to name code snippets.


~~~~{.c}


	struct stat buf;
	int ret,i;
	int file_len = 0;
	
	ret= stat ( file_name, &buf );
	if (ret != 0){
		fprintf(stderr,"Error: file %s not found\n",file_name);
		exit(EXIT_FAILURE);
	}
	## Reading the file
~~~~

	
Here the fun starts:
	
~~~~{.c}
	
for(i = 1; i < file_len;i++){
	if(file_string[i] == '\n'){
		line_number++;
	}
	## Recognizing markdown
~~~~

	
### Recognizing markdown
Here I am at the start of a line. The code below checks if any special markdown syntax is present and sets the corresponding variables.
		
~~~~{.c}
	

if(file_string[i-1] == '\n'  ||  file_string[i-1] == 0 ){
	space = 0;
	wiggle = 0;
	newline = 0;
	switch (file_string[i]) {
		case ' ':
				space = 1;
				while(file_string[i+space] == ' '){
					space++;
				}
				break;
			case '\t':
				space = 4;
				break;
			case '~':
			case '`':
				wiggle = 1;
				while(file_string[i+wiggle] == '~' || file_string[i+wiggle] == '`' ){
					wiggle++;
				}
				break;
			case '\n':
				newline = 2;
				break;
			case '#':
				if(!fence){
					heading_pointer = i;
				}
				break;
				
			default:
				newline = 0;
				space = 0;
				wiggle = 0;
				break;
	}
}
~~~~
	
Depending on the current  `state` and the variables set above the code updates to the next state. There are only two possible states at the moment: code or normal markdown section. In future I should add a list state.
		
~~~~{.c}
switch (state) {
	case CODE_STATE:
		if(fence){
			if(wiggle >= 3){
				code_stop = i-1;
				while(file_string[i] != '\n'){
					i++;
				}
				line_number++;
				fence = 0;
				wiggle = 0;
				state = SECTION_STATE;
## enter code
				}
			}else{
				if(newline == 2 || space < 4){
					state = SECTION_STATE;
					code_stop = i;
## enter code
				}
			}
		break;
	case SECTION_STATE:
		if(space >= 4){
			state = CODE_STATE;
			code_start = i +1;
			code_start_line_number = line_number;
		}
		if(wiggle >= 3){
			while(file_string[i] != '\n'){
				i++;
			}
			line_number++;
			state = CODE_STATE;
			code_start = i +1;
			code_start_line_number = line_number;
			fence = 1;
			wiggle = 0;
		}
		break;
	default:
		break;
}
}

~~~~

### enter code

~~~~{.c}
md_doc->root = enter_code(md_doc->root,md_doc->file_string,code_start,code_stop,heading_pointer, code_start_line_number);

~~~~


## Reading the file
The call below loads the entire file into one string. To make the parsing easier I append a single newline character at the front.

~~~~{.c}

md_doc->file_string = get_input_into_string(md_doc->file_string, file_name, &file_len);
char* file_string = md_doc->file_string;
int state = SECTION_STATE;
int space = 0;
int wiggle = 0;

int heading_pointer = -1;
int newline = 0;
int fence = 0;
int code_start_line_number = 0;
int code_start = -1;
int code_stop = -1;
int line_number = 0;


~~~~



## Writing the code to file

Here is the rest of the code:

## The rest

~~~~{.c}

## Header

## Data Structures

## rbnode_def


struct rb_node* insert_start(struct rb_node* n, char* key,void* p) ;
struct rb_node* insert (struct rb_node* n, char*  key, void* p);

int rank(struct rb_node* n);

struct rb_node* set_num(struct rb_node* n);
int size(struct rb_node* n);

int isred(struct rb_node* n);
struct rb_node* rotateLeft(struct rb_node* n);
struct rb_node* rotateRight(struct rb_node* n);
struct rb_node* flipColors(struct rb_node* n) ;
void free_rbtree(struct rb_node* n);

struct rb_node* enter_code(struct rb_node* root,char* string,int start, int stop, int header, int line_number);
struct rb_node* find_entry(struct rb_node* n, char* key);
int traverse_tree(struct rb_node* root, struct rb_node* n,char* md_filename,int mode);

void print_code_tree(struct rb_node* root,char* target,char* md_filename,FILE* out, int type);
void usage();


struct md_doc*  parse_file(struct md_doc* md_doc,char* filename,int file_num);

char* shorten_pathname(char* p);
char* get_input_into_string(char* string,char* infile,int* len);


int main( int argc, char * argv[])
{
	//char* outfilename = 0;
	int c = 0;
	int help = 0;
	int infiles = 0;
	char** infile = 0;
	
	while (1){
		static struct option long_options[] ={
			{0, 0, 0, 0}
		};
		
		int option_index = 0;
		c = getopt_long_only (argc, argv,"h",long_options, &option_index);
		
		if (c == -1){
			break;
		}
		
		switch(c) {
			//case 'o':
			//	outfilename = optarg;
			//	break;
			case 'h':
				help = 1;
				break;
			case '?':
				exit(1);
				break;
		}
	}
	
	if(help){
		usage();
		exit(EXIT_SUCCESS);
	}
	
	infile = malloc(sizeof(char*)* (argc-optind));
	
	c = 0;
	while (optind < argc){
		infile[c] =  argv[optind++];
		c++;
	}
	infiles = c;
	
	if(!infiles){
		usage();
		exit(EXIT_SUCCESS);

	}
	
	struct md_doc* md_doc = malloc(sizeof(struct md_doc));
	md_doc->root = 0;
	md_doc->file_string = 0;
	md_doc = parse_file(md_doc,infile[0],0);
	
	traverse_tree(md_doc->root,md_doc->root, shorten_pathname(infile[0])  , 0);
	
	
	
	free_rbtree(md_doc->root);
	free(md_doc->file_string);
	free(md_doc);
	free(infile);
	exit(EXIT_SUCCESS);
}

int traverse_tree(struct rb_node* root, struct rb_node* n,char* md_filename,int mode)
{
	FILE* out = 0;
	//struct line* line_p = 0;
	struct stat buf;
	int add = 5;
	int type = 0;
	if(n->left){
		mode = traverse_tree(root,n->left,md_filename,mode);
	}
	//line_p = (struct line*) n->p;
	if(strncmp("FILE:",  n->key, 5)  == 0 || strncmp("File:",  n->key, 5)  == 0  || strncmp("file:",  n->key, 5)  == 0){
		add = 5;
		type = 0;
		while(!isalnum(n->key[add])){
			add++;
		}
		if((  n->key[(int)strlen(n->key)-1] == 'c' || n->key[(int)strlen(n->key)-1] == 'h' )  && n->key[(int)strlen(n->key)-2] == '.' ){
			type = 1;
		}
		if( n->key[(int)strlen(n->key)-1] == 'd' && n->key[(int)strlen(n->key)-2] == 'm'){
			type = -1;
		}
		
		fprintf(stderr,"FILE::::\n%s	%d\n",n->key+add,type );
		
		if(strcmp("mdlp_main.c",n->key+add) == 0){
			if(! stat ( n->key+add, &buf )){
				fprintf(stderr,"ERROR: I will not overwrite myself: %s\n",n->key+add);
				exit(EXIT_FAILURE);
			}
			
		}
		
		if(type >= 0){
			if (!(out = fopen(n->key+add, "w" ))){
				fprintf(stderr,"Cannot write to file '%s'\n",n->key+add);
				exit(EXIT_FAILURE);
			}
			print_code_tree(root, n->key,md_filename,out,type);
			fclose(out);
			chmod(n->key+add, S_IRWXU);
		}
		
		
		
		

		
	}
	if(n->right){
		mode = traverse_tree(root,n->right,md_filename,mode);
	}
	return mode;
}

void print_code_tree(struct rb_node* root,char* target, char* md_filename,FILE* out, int type)
{
	struct rb_node* n = 0;
	struct line* line_p = 0;
	
	int i,c;
	int line_number;
	int code_section_name = 0;
	n = root;
	n = find_entry(n,target );
	line_p = (struct line*) n->p;
	//fprintf(stderr,"TARGET:%s FOUND: %s\n", target,n->key);
	while (line_p) {
		line_number =line_p->line_number+1;
		if(type == 1){
			fprintf(out,"\n#line	%d	\"%s\"\n",line_number , md_filename );
		}
		for(i = 0; i < line_p->len;i++){
			if(line_p->content[i] == '\n'){
				line_number++;
			}
			if(line_p->content[i] == '#' && line_p->content[i+1] == '#'){
				line_number++;
				code_section_name = i;
				while(1){
					if(line_p->content[i] == '#' || !isalnum((int) line_p->content[i])){
						code_section_name++;
					}else{
						break;
					}
					i++;
				}
				
				i = code_section_name;
				while(1){
					
					if(line_p->content[i] =='\n' || line_p->content[i] == 0){
						break;
					}
					i++;
				}
				c = i;
				
				while(!isalnum((int)line_p->content[i])){
					line_p->content[i] = 0;
					i--;
				}
				//fprintf(stderr,"SUBSECTION:%s\n",line_p->content+code_section_name);
				i = c;
				print_code_tree( root,line_p->content+code_section_name,md_filename,out,type);
				if(type == 1){
					fprintf(out,"\n#line	%d	\"%s\"\n",line_number , md_filename );
				}
			}else{
				fprintf(out,"%c",line_p->content[i]);
			}
		}
		line_p = line_p->next;
	}
}

void usage()
{
	fprintf(stdout, "\n%s %s, Copyright (C) 2013 Timo Lassmann <%s>\n",PACKAGE_NAME, PACKAGE_VERSION,PACKAGE_BUGREPORT);
	fprintf(stdout, "\n");
	fprintf(stdout, "Usage:   lpcode <directory> <file1> <file2> .... \n\n");
	fprintf(stdout, "Options:\n");
	fprintf(stdout, "         -o         STR     output file name.\n");
	fprintf(stdout, "\n");
}

char* shorten_pathname(char* p)
{
	int i;
	char* tmp = p;
	for(i = 0; i< (int) strlen(p);i++){
		if(p[i] == '/'){
			tmp = p+i +1;
		}
	}
	return tmp;
}


struct md_doc* parse_file(struct md_doc* md_doc,char* file_name,int file_num)
{
	## Parsing
	return md_doc;
}



struct rb_node* enter_code(struct rb_node* root , char* string,int start, int stop, int header, int line_number)
{
	int i;
	i = header;
	while(1){
		if(string[i] == '#' || !isalnum((int) string[i])){
			header++;
		}else{
			break;
		}
		i++;
	}
	
	i = header;
	while(1){
		
		if(string[i] =='\n' || string[i] == 0){
			break;
		}
		i++;
	}
	
	
	
	while(!isalnum((int)string[i])){
		string[i] = 0;
		i--;
	}
	string[stop] = 0;
	
	//fprintf(stderr,"Insert:%s\n",string + header );
	
	struct line* line = malloc(sizeof(struct line));
	line->content = string+start;
	line->next = 0;
	line->line_number = line_number;
	line->len = stop - start;
	root =  insert_start(root,string + header, line );
	return root;
}

char* get_input_into_string(char* string,char* infile,int* len)
{
	long int i = 0;
	FILE *file = 0;
	size_t bytes_read;
	if (!(file = fopen( infile, "r" ))){
		return 0;
		fprintf(stderr,"Cannot open file '%s'\n", infile);
		exit(-1);
	}
	if (fseek(file,0,SEEK_END) != 0){
		(void)fprintf(stderr, "ERROR: fseek failed\n");
		(void)exit(EXIT_FAILURE);
	}
	i= ftell (file);
	if (fseek(file,0,SEEK_START) != 0){
		(void)fprintf(stderr, "ERROR: fseek failed\n");
		(void)exit(EXIT_FAILURE);
	}
	if(!string){
		string = malloc ((i+1+18)* sizeof(char));
	}
	string[0] = '\n';
	
	bytes_read = fread(string+1,sizeof(char), i, file);
	if(!bytes_read){
		fprintf(stderr,"Reading from file:%s failed \n", infile );
		exit(EXIT_FAILURE);
	}
	
	string[i+1] = 0;
	fclose(file);
	*len = (int) i+1;
	return string;
}


struct rb_node* find_entry(struct rb_node* n, char* key)
{
	int c;
	c =strcmp(key, n->key);
	
	/** I think that keys should never be the same in the lpcode program.*/
	if(!c){
		return n;
	}else if(c < 0){
		if(n->left){
			return find_entry(n->left,key);
		}else{
			fprintf(stderr,"%s not found in rbtree...\n",key);
			exit(EXIT_FAILURE);
		}
		
	}else{
		if(n->right){
			return find_entry(n->right,key);
		}else{
			fprintf(stderr,"%s not found in rbtree...\n",key);
			exit(EXIT_FAILURE);
		}
	}
	return n;
}

~~~~

## Red-Black Tree Implementation


Each node of the RB-tree contains a pointer to a `struct line`.

## rbnode_def

~~~~{.c}
struct rb_node{
	struct rb_node* right;
	struct rb_node* left;
	void* p ;
	char* key;
	int color;
	unsigned int num:31;
};
~~~~

As additional lines belonging to same code block are added a linked list is formed.


## Insert_code

~~~~{.c}


struct rb_node* insert_start(struct rb_node* n, char*  key,void* p)
{
	n = insert(n, key,p);
	n->color = BLACK;
	return n;
}

struct rb_node* insert (struct rb_node* n, char* key,void*p)
{
	int c;
	struct line* line_p = 0;
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
	c =strcmp(key, n->key);
	
	/** I think that keys should never be the same in the lpcode program.*/
	if(!c){
		line_p = (struct line*) n->p;
		while(line_p->next){
			line_p = line_p->next;
		}
		line_p->next = p;
	}else if(c < 0){
		n->left = insert(n->left,key,p);
	}else {
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
	struct line* line_p= 0;
	struct line* old_p= 0;
	
	if(n->left){
		free_rbtree(n->left);
	}
	
	
	if(n->right){
		free_rbtree(n->right);
	}
	line_p = (struct line*) n->p;
	old_p = (struct line*) n->p;
	
	while(line_p->next){
		line_p = line_p->next;
		free(old_p);
		old_p =line_p;
	}
	free(line_p);
	free(n);
}


~~~~





## File: nstart.c

	## The rest

	## Insert_code




## Header

Loads libraries. Not very interesting...

~~~~{.c}

#ifndef mdlp_main

#define mdlp_main

#define SEEK_START 0
#define SEEK_END 2

#ifdef DEBUG
#define DEBUG_TEST 1
#else
#define DEBUG_TEST 0
#endif

#define debug_print(fmt, ...) \
do { if (DEBUG_TEST) fprintf(stderr, "%s:%d:%s(): " fmt, __FILE__, \
__LINE__, __func__, __VA_ARGS__); } while (0)

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <getopt.h>
#include <time.h>
#include<sys/stat.h>
#include "config.h"


#define MAX_LINE 50000


#define RED 1
#define BLACK 0


#define CODE_STATE 1
#define SECTION_STATE 2
#define LIST_STATE 3

#endif


~~~~

