/*
 
 Copyright (C) 2013 Timo Lassmann <timolassmann@gmail.com>
 
 This file is part of TagDust.
 
 TagDust is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 TagDust is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with Tagdust.  If not, see <http://www.gnu.org/licenses/>.
 
 */

/*! \file exact.h
 \brief Functions to extract sequences without using errors.

 \author Timo Lassmann
 \bug No known bugs.
 */

#ifndef tagdust2_exact_h
#define tagdust2_exact_h
void exact_controller(struct parameters* param,int (*fp)(struct read_info** ,struct parameters*,FILE* ),int file_num);

int byg_end_barcode(const char* pattern,const char*text, int m, int n);// m is pattern length , n is text length....;
void print_seq_from_position_x(struct read_info* ri,FILE* out,int x);
#endif
