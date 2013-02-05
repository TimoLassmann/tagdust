/*
 
 Copyright (C) 2010 Timo Lassmann <timolassmann@gmail.com>
 
 This file is part of SAMstat.
 
 Delve is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 Delve is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with Delve.  If not, see <http://www.gnu.org/licenses/>.
 
 */

#include <stdio.h>
#include <getopt.h>
#include <stdlib.h>
#include <string.h>

#define OPT_UNMAPPED 1

struct parameters {
	char** infile;
	char* outfile;
	int kmer_size;
	int infiles;
	int quiet_flag;
	int num_query;
	char* print_unmapped;
	int solexa;
	int print_qual;
	int print_posteriors;
	int k_errors_allowed;
	char*  summary;
	char* format;
	char* filter;
	char* alt_lib_name;
	int gzipped;
	int sam;
};

struct parameters* interface(struct parameters* param,int argc, char *argv[]);

void free_param(struct parameters* param);
void usage(void);



