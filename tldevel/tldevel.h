#ifndef TLDEVEL_H
#define TLDEVEL_H

#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef TLDEVEL_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif

#define OK              0
#define FAIL            1

#define MESSAGE_MARGIN 22

#define MACRO_MIN(a,b)          (((a)<(b))?(a):(b))
#define MACRO_MAX(a,b)          (((a)>(b))?(a):(b))

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)
#define AT __FILE__ " line " TOSTRING(__LINE__)


#define ERROR_MSG(...) do {                     \
                error(AT, __VA_ARGS__ );        \
                goto ERROR;                     \
        }while (0)

#define WARNING_MSG(...) do {                   \
                warning(AT, __VA_ARGS__ );      \
        }while (0)


#define LOG_MSG(...) do {                       \
                log_message( __VA_ARGS__ );     \
        }while (0)

#define ASSERT(TEST,...)  if(!(TEST)) {         \
                error(AT,#TEST );               \
                error(AT, ##__VA_ARGS__);       \
                goto ERROR;                     \
        }

#if (DEBUGLEVEL >= 1)
#define DASSERT(TEST,...) if(!(TEST)) {           \
                        error(AT,#TEST );         \
                        error(AT, ##__VA_ARGS__); \
                        goto ERROR;               \
        }
#else
#define DASSERT(TEST,...)

#endif


#define ADDFAILED(x)  "Function \"" TOSTRING(x) "\" failed."

#define RUN(EXP) do {                               \
                if((EXP) != OK){                    \
                        ERROR_MSG(ADDFAILED(EXP));	\
                }                                   \
        }while (0)

#define RUNP(EXP) do {                              \
                if((EXP) == NULL){                  \
                        ERROR_MSG(ADDFAILED(EXP));	\
                }                                   \
        }while (0)

/* Functions to declare and use a timer */

#define DECLARE_TIMER(n) struct timespec ts1_##n; struct timespec ts2_##n;
#define START_TIMER(n) clock_gettime(CLOCK_MONOTONIC_RAW, &ts1_##n);
#define STOP_TIMER(n) clock_gettime(CLOCK_MONOTONIC_RAW, &ts2_##n);
#define GET_TIMING(n) (double)(ts2_##n.tv_sec - ts1_##n.tv_sec) + ((double)  ts2_##n.tv_nsec - ts1_##n.tv_nsec) / 1000000000.0

/* Memory functions  */

#define MFREE(p) do {                                           \
                if(p){                                          \
                        free(p);                                \
                        p = NULL;                               \
                }else{                                          \
                        WARNING_MSG("free on a null pointer");  \
                }                                               \
        } while (0)

#define MMALLOC(p,size) do {                                          \
                if (p != NULL){                                       \
                        ERROR_MSG( "malloc on a nun-null pointer");   \
                        goto ERROR;                                   \
                }                                                     \
                if(size == 0){                                        \
                        ERROR_MSG("malloc of size %d failed", size);	\
                        goto ERROR;                                   \
                }                                                     \
                if (((p) = malloc(size)) == NULL) {                   \
                        ERROR_MSG("malloc of size %d failed", size);	\
                        goto ERROR;                                   \
                }                                                     \
        } while (0)

#define MREALLOC(p, size) do {                                          \
                void *tmpp;                                             \
                if(size == 0){                                          \
                        ERROR_MSG("malloc of size %d failed", size);    \
                        goto ERROR;                                     \
                }                                                       \
                if ((p) == NULL) {                                      \
                        tmpp = malloc(size);                            \
                }else {                                                 \
                        tmpp = realloc((p), (size));                    \
                }                                                       \
                if (tmpp != NULL){                                      \
                        p = tmpp;                                       \
                }else {                                                 \
                        ERROR_MSG("realloc for size %d failed", size);  \
                        goto ERROR;                                     \
                }} while (0)

/* g memory functions */

EXTERN int get_dim1(void* ptr, int* d);
EXTERN int get_dim2(void* ptr, int* d);

#define CONCAT(X, Y) CONCAT_(X, Y)
#define CONCAT_(X, Y) X ## Y

#define ARGN(...) ARGN_(__VA_ARGS__)
#define ARGN_(_0, _1, _2, _3 , N, ...) N

#define NARG(...) ARGN(__VA_ARGS__ COMMA(__VA_ARGS__) 4, 3, 2, 1, 0)
#define HAS_COMMA(...) ARGN(__VA_ARGS__, 1, 1, 0)

#define SET_COMMA(...) ,

#define COMMA(...) SELECT_COMMA                     \
        (                                           \
                HAS_COMMA(__VA_ARGS__),             \
                HAS_COMMA(__VA_ARGS__ ()),          \
                HAS_COMMA(SET_COMMA __VA_ARGS__),   \
                HAS_COMMA(SET_COMMA __VA_ARGS__),   \
                HAS_COMMA(SET_COMMA __VA_ARGS__ ()) \
                )

#define SELECT_COMMA(_0, _1, _2, _3, _4) SELECT_COMMA_(_0, _1, _2, _3, _4)
#define SELECT_COMMA_(_0, _1, _2, _3, _4) COMMA_ ## _0 ## _1 ## _2 ## _3 ## _4

#define COMMA_00000 ,
#define COMMA_00001
#define COMMA_00010 ,
#define COMMA_00011 ,
#define COMMA_00100 ,
#define COMMA_00101 ,
#define COMMA_00110 ,
#define COMMA_00111 ,
#define COMMA_01000 ,
#define COMMA_01001 ,
#define COMMA_01010 ,
#define COMMA_01011 ,
#define COMMA_01100 ,
#define COMMA_01101 ,
#define COMMA_01110 ,
#define COMMA_01111 ,
#define COMMA_10000 ,
#define COMMA_10001 ,
#define COMMA_10010 ,
#define COMMA_10011 ,
#define COMMA_10100 ,
#define COMMA_10101 ,
#define COMMA_10110 ,
#define COMMA_10111 ,
#define COMMA_11000 ,
#define COMMA_11001 ,
#define COMMA_11010 ,
#define COMMA_11011 ,
#define COMMA_11100 ,
#define COMMA_11101 ,
#define COMMA_11110 ,
#define COMMA_11111 ,

#define ALLOC_1D_ARRAY_DEF(type)                                  \
        EXTERN int alloc_1D_array_size_ ##type (type **array, int dim1)

ALLOC_1D_ARRAY_DEF(char);
ALLOC_1D_ARRAY_DEF(int);
ALLOC_1D_ARRAY_DEF(float);
ALLOC_1D_ARRAY_DEF(double);

#define ALLOC_2D_ARRAY_DEF(type)                                        \
        EXTERN int alloc_2D_array_size_ ##type (type ***array, int dim1,int dim2)

ALLOC_2D_ARRAY_DEF(char);
ALLOC_2D_ARRAY_DEF(int);
ALLOC_2D_ARRAY_DEF(float);
ALLOC_2D_ARRAY_DEF(double);


#define FREE_VOID_DEF(type)                     \
        void gfree_void_ ##type(type *a)

#define FREE_1D_ARRAY_DEF(type)                 \
        void free_1d_array_ ##type(type **array)

#define FREE_2D_ARRAY_DEF(type)                   \
        void free_2d_array_ ##type(type ***array)


FREE_VOID_DEF(char);
FREE_VOID_DEF(int);
FREE_VOID_DEF(float);
FREE_VOID_DEF(double);

FREE_1D_ARRAY_DEF(char);
FREE_1D_ARRAY_DEF(int);
FREE_1D_ARRAY_DEF(float);
FREE_1D_ARRAY_DEF(double);

FREE_2D_ARRAY_DEF(char);
FREE_2D_ARRAY_DEF(int);
FREE_2D_ARRAY_DEF(float);
FREE_2D_ARRAY_DEF(double);



#define galloc(...) SELECTGALLOC(__VA_ARGS__)(__VA_ARGS__)
#define SELECTGALLOC(...) CONCAT(SELECTGALLOC_, NARG(__VA_ARGS__))(__VA_ARGS__)

#define SELECTGALLOC_0()

#define SELECTGALLOC_1(_1) _Generic ((_1),             \
                                     default: galloc_void \
                )

#define SELECTGALLOC_2(_1, _2) _Generic((_1),                        \
                                        char**: alloc_1D_array_size_char, \
                                        int**: alloc_1D_array_size_int,  \
                                        float**:  alloc_1D_array_size_float, \
                                        double**:alloc_1D_array_size_double \
                )

#define SELECTGALLOC_3(_1, _2, _3) _Generic((_1),                    \
                                                char***: _Generic((_2),  \
                                                                 int: alloc_2D_array_size_char \
                                                        ),              \
                                                int***: _Generic((_2),   \
                                                                int: alloc_2D_array_size_int \
                                                        ),              \
                                                float***: _Generic((_2), \
                                                                  int: alloc_2D_array_size_float \
                                                        ),              \
                                                double***: _Generic((_2), \
                                                                   int: alloc_2D_array_size_double \
                                                        )              \
                )




#define gfree(X) _Generic((&X),                            \
                          char*: gfree_void_char,          \
                          int*: gfree_void_int,            \
                          float*: gfree_void_float,        \
                          double*: gfree_void_double,      \
                          char**: free_1d_array_char,      \
                          int**: free_1d_array_int,        \
                          float**: free_1d_array_float,    \
                          double**: free_1d_array_double,  \
                          char***: free_2d_array_char,     \
                          int***: free_2d_array_int,       \
                          float***: free_2d_array_float,   \
                          double***: free_2d_array_double  \
                )(&X)

/* functions  */

EXTERN void error(const char *location, const char *format, ...);
EXTERN void warning(const char *location, const char *format, ...);
EXTERN void log_message( const char *format, ...);

#undef TLDEVEL_IMPORT
#undef EXTERN
#endif
