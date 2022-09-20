//Compile with C99 or above
// Written by Muralidhar Nalabothula July 2022
/**
 * @brief A simple library for ND arrays and perform blas operations on them. Read and write I/O is 
 * performed using netCDF library
 * Must link some blas and netCDF
 * Contains functions for float(s), double(d), float complex (c), double complex (z)
 * array object nd_arr(##TYPE) Ex: for double complex nd_arrz
 * all indices, counters must be in ND_indices
 */
#pragma once
#include <stdio.h>      
#include <complex.h>
#include <stddef.h>
#include <stdarg.h>
#include <stdlib.h>
#include <stdbool.h> 
#include <math.h>
#include <time.h>
#include <ctype.h>
#include <string.h>
#include <netcdf.h>
#include <tblis/tblis.h>
#include <stdint.h>
#include "nd_ulit.h"
//
#ifdef ND_MKL_BLAS
    #define lapack_complex_float    float _Complex
    #define lapack_complex_double   double _Complex
    #include <mkl.h>
#else
    #define LAPACK_COMPLEX_C99
    #include <cblas.h>
    #include <lapacke.h>
#endif
//

typedef int BLAS_INT; /* Internal type used for blas indices. set this according to the blas library (32 or 64 bit) */
typedef size_t ND_indices; /* type used for all indices internally. */
#define nd_idx (ND_indices []) /* macro for giving indices to functions*/
#define ND_ALL ((void *)0)


#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(EXIT_FAILURE);}


/*********************************************************************************************************************************/
/** Note the below line should always be starting line to generate header file */
/* PYTHON_HEADER_FILE START */
#define FUNCTION(FUN_NAME, TYPE_SMALL)         FUNCTION_HIDDEN(FUN_NAME, TYPE_SMALL)
#define FUNCTION_HIDDEN(FUN_NAME, TYPE_SMALL)  nd_ ## FUN_NAME ## _ ## TYPE_SMALL

#define ARRAY_T(TYPE_SMALL)         ARRAY_T_HIDDEN(TYPE_SMALL)
#define ARRAY_T_HIDDEN(TYPE_SMALL)  nd_arr ## _ ## TYPE_SMALL

#define ND_FUNCTION_CALL(FUN_NAME, TYPE_SMALL)         ND_FUNCTION_CALL_HIDDEN(FUN_NAME, TYPE_SMALL)
#define ND_FUNCTION_CALL_HIDDEN(FUN_NAME, TYPE_SMALL)  nd_ ## FUN_NAME ## _ ## TYPE_SMALL

#define BLAS_CALL(FUN_NAME, TYPE_SMALL)         BLAS_CALL_HIDDEN(FUN_NAME, TYPE_SMALL)
#define BLAS_CALL_HIDDEN(FUN_NAME, TYPE_SMALL)  cblas_ ## TYPE_SMALL ## FUN_NAME

#define TBLIS_CALL(FUN_NAME, TYPE_SMALL)         TBLIS_CALL_HIDDEN(FUN_NAME, TYPE_SMALL)
#define TBLIS_CALL_HIDDEN(FUN_NAME, TYPE_SMALL)  tblis_ ## FUN_NAME ## _ ## TYPE_SMALL 
/*********************************************************************************************************************************/

#if defined(COMPILE_ND_FLOAT)
    #define TYPE_S s    /* */
    #define TYPE_L float /* type */
    #define NetCDF_IO_TYPE NC_FLOAT
    #define NetCDF_FUN_TYPE float

#elif defined(COMPILE_ND_DOUBLE)
    #define TYPE_S d    /* */
    #define TYPE_L double /* type */
    #define NetCDF_IO_TYPE NC_DOUBLE
    #define NetCDF_FUN_TYPE double

#elif defined(COMPILE_ND_SINGLE_COMPLEX)
    #define TYPE_S c    /* */
    #define TYPE_L float complex /* type */
    #define NetCDF_IO_TYPE NC_FLOAT
    #define NetCDF_FUN_TYPE float

#elif defined(COMPILE_ND_DOUBLE_COMPLEX)
    #define TYPE_S z    /* */
    #define TYPE_L double complex /* type */
    #define NetCDF_IO_TYPE NC_DOUBLE
    #define NetCDF_FUN_TYPE double

#elif defined(COMPILE_ND_INT)
    #define TYPE_S i    /* */
    #define TYPE_L int /* type */
    #define NetCDF_IO_TYPE NC_INT
    #define NetCDF_FUN_TYPE int

#endif


#if defined(COMPILE_ND_DOUBLE_COMPLEX) || defined(COMPILE_ND_SINGLE_COMPLEX)
    #define BLAS_POINTER(VARIABLE_NAME) &VARIABLE_NAME    /* for C and Z, pointer is passed for constant in blas routines */
#else
    #define BLAS_POINTER(VARIABLE_NAME) VARIABLE_NAME    /*  for S and D, value is passed for constant in blas routines */
#endif

/*********************************************************************************************************************************/


typedef struct ARRAY_T(TYPE_S) {
    TYPE_L * data;        // store data pointer
    ND_indices * rank;    // rank of the tensor
    ND_indices * dims;    // pointer to dims array
    ND_indices * strides; // pointer to strides of an array (in elements and not in bytes)
    bool owner; /* set this to true if data belongs to this array. if it is referened then false. 
                when array is free data is not freed if owner is false */
} ARRAY_T(TYPE_S);



/**************************************************** alloc.c functions **********************************************************/
/*********************************************************************************************************************************/

void FUNCTION(init, TYPE_S) (ARRAY_T(TYPE_S) * nd_arr_in, const ND_indices rank, const ND_indices * dimensions);

void FUNCTION(uninit, TYPE_S) (ARRAY_T(TYPE_S) * nd_arr_in);

void FUNCTION(malloc, TYPE_S) (ARRAY_T(TYPE_S) * nd_arr_in);

void FUNCTION(calloc, TYPE_S) (ARRAY_T(TYPE_S) * nd_arr_in);

void FUNCTION(free, TYPE_S) (ARRAY_T(TYPE_S) * nd_arr_in);

void FUNCTION(init_tranpose, TYPE_S) (const ARRAY_T(TYPE_S) * nd_arr_in, const ND_indices * order, ARRAY_T(TYPE_S) * nd_arr_out);

void FUNCTION(init_slice, TYPE_S) (const ND_indices * start_idx, const ND_indices * end_idx, const ND_indices * step_idx, \
                            const ARRAY_T(TYPE_S) * nd_arr_in, ARRAY_T(TYPE_S) * nd_arr_out );

void FUNCTION(init_strip_dims, TYPE_S) (const ARRAY_T(TYPE_S) * nd_arr_in, const ND_indices n_dims_strip, ARRAY_T(TYPE_S) * nd_arr_out);

/**************************************************** array.c functions **********************************************************/
/*********************************************************************************************************************************/

TYPE_L * FUNCTION(ele, TYPE_S) (const ARRAY_T(TYPE_S) * nd_arr_in, const ND_indices * dimensions);

ND_indices FUNCTION(size, TYPE_S) (const ARRAY_T(TYPE_S) * nd_arr_in);

void FUNCTION(set_all, TYPE_S) (ARRAY_T(TYPE_S) * nd_arr_in, const TYPE_L set_constant );

void FUNCTION(reshape, TYPE_S) (const ARRAY_T(TYPE_S) * nd_arr_in, ARRAY_T(TYPE_S) * nd_arr_out);

void FUNCTION(strip_dims, TYPE_S) (const ARRAY_T(TYPE_S) * nd_arr_in, const ND_indices n_dims_strip, \
                                    const ND_indices * stripped_idxs, ARRAY_T(TYPE_S) * nd_arr_out);

void FUNCTION(slice, TYPE_S) (const ND_indices * start_idx, const ND_indices * end_idx, const ND_indices * step_idx, \
                            const ARRAY_T(TYPE_S) * nd_arr_in, ARRAY_T(TYPE_S) * nd_arr_out );

void FUNCTION(tranpose, TYPE_S) (const ARRAY_T(TYPE_S) * nd_arr_in, const ND_indices * order, ARRAY_T(TYPE_S) * nd_arr_out);

void FUNCTION(copy, TYPE_S) (const ARRAY_T(TYPE_S) * nd_arr_in, ARRAY_T(TYPE_S) * nd_arr_out);


/**************************************************** netcdf_io.c functions ******************************************************/
/*********************************************************************************************************************************/


void FUNCTION(read, TYPE_S) (const char* file_name, const char* var_name, ARRAY_T(TYPE_S) * nd_arr_in);

void FUNCTION(read_sub, TYPE_S) (const char* file_name, const char* var_name, ARRAY_T(TYPE_S) * nd_arr_in, ...);

void FUNCTION(write, TYPE_S) (const char* file_name, const char* var_name, const ARRAY_T(TYPE_S) * nd_arr_in, char ** dim_names, size_t * chunksize);

/*************************************************** linalg.c functions **********************************************************/
/*********************************************************************************************************************************/

#if defined(COMPILE_ND_DOUBLE_COMPLEX) || defined(COMPILE_ND_SINGLE_COMPLEX) || defined(COMPILE_ND_FLOAT) || defined(COMPILE_ND_DOUBLE)

void FUNCTION(matmul, TYPE_S) (const char TransA, const char TransB, const ARRAY_T(TYPE_S) * nd_arr_A, const ARRAY_T(TYPE_S) * nd_arr_B, \
                                ARRAY_T(TYPE_S) * nd_arr_C, const TYPE_L alpha, const TYPE_L beta, const ND_indices * A_idx,  \
                                                                        const ND_indices * B_idx,  const ND_indices * C_idx);

void FUNCTION(sum, TYPE_S) (char * str_A, char * str_C, ARRAY_T(TYPE_S) * nd_arrA, ARRAY_T(TYPE_S) * nd_arrC, const TYPE_L alpha, const TYPE_L beta);

void FUNCTION(einsum, TYPE_S) (char * einsum_indices, ARRAY_T(TYPE_S) * nd_arrA, ARRAY_T(TYPE_S) * nd_arrB, ARRAY_T(TYPE_S) * nd_arrC, \
                                                                                                const TYPE_L alpha, const TYPE_L beta);
#endif






