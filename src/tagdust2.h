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
 along with ;.  If not, see <http://www.gnu.org/licenses/>.
 
 */


/*! \file tagdust2.h
 \brief Defines some global macros 
*/



#ifndef tagdust2_tagdust2_h

#ifdef DEBUG
#define DEBUG_TEST 1
#else
#define DEBUG_TEST 0
#endif

#define debug_print(fmt, ...) \
do { if (DEBUG_TEST) fprintf(stderr, "%s:%d:%s(): " fmt, __FILE__, \
__LINE__, __func__, __VA_ARGS__); } while (0)


#define tagdust2_tagdust2_h

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <ctype.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>
#include <assert.h>

/* Safe string operation macros */
#define SAFE_STRCAT(dest, src, size) do { \
    size_t dest_len = strnlen(dest, size - 1); \
    if (dest_len < size - 1) { \
        size_t remaining = size - dest_len - 1; \
        size_t src_len = strnlen(src, remaining); \
        memcpy(dest + dest_len, src, src_len); \
        dest[dest_len + src_len] = '\0'; \
    } \
} while(0)

#define SAFE_STRCAT_CHECK(dest, src, size, retval) do { \
    size_t dest_len = strnlen(dest, size - 1); \
    if (dest_len >= size - 1) return retval; \
    size_t remaining = size - dest_len - 1; \
    size_t src_len = strnlen(src, remaining); \
    if (dest_len + src_len >= size) return retval; \
    memcpy(dest + dest_len, src, src_len); \
    dest[dest_len + src_len] = '\0'; \
} while(0)

#define SAFE_STRCPY(dest, src, size) do { \
    strncpy(dest, src, size - 1); \
    dest[size - 1] = '\0'; \
} while(0)

#define VALIDATE_STRING_LENGTH(str, max_len, error_msg) do { \
    if (strnlen(str, max_len + 1) > max_len) { \
        fprintf(stderr, "Error: %s (max %zu chars)\n", error_msg, (size_t)max_len); \
        return -1; \
    } \
} while(0)


/** \def MAX_SEQ_LEN 
 \brief Maximum length of sequences;
 */
//#define MAX_SEQ_LEN 500
/** \def MAX_LINE 
 \brief Maximum length of input line. 
 */
#define MAX_LINE 10000

/* Security-focused buffer size constants */
#define MSG_BUFFER_SIZE 1024        /* For param->buffer and error messages */
#define FILENAME_BUFFER_SIZE 1000   /* For file paths and names */
#define SEQUENCE_BUFFER_SIZE 1000   /* For biological sequences */
#define COMMAND_BUFFER_SIZE 10000   /* For command lines and large text */
#define HEADER_BUFFER_SIZE 512      /* For FASTQ headers */
#define BARCODE_BUFFER_SIZE 100     /* For barcode sequences */






#endif
