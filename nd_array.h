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

// #include <netcdf.h>

// //
// #ifdef MKL_BLAS
//     #define lapack_complex_float    float _Complex
//     #define lapack_complex_double   double _Complex
//     #include <mkl.h>
// #else
//     #define LAPACK_COMPLEX_C99
//     #include <cblas.h>
//     #include <lapacke.h>
// #endif

typedef int BLAS_INT; /* Internal type used for blas indices. set this according to the blas library (32 or 64 bit) */
typedef size_t ND_indices; /* type used for all indices internally. */
#define nd_idx (ND_indices []) /* macro for giving indices to functions*/
#define ND_ALL ((void *)0)
#define ND_END 0

// #define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(EXIT_FAILURE);}
#include "nd_ulit.h"
#include "nd_INT.h"
#include "nd_FLOAT.h"
#include "nd_DOUBLE.h"
#include "nd_SINGLE_COMPLEX.h"
#include "nd_DOUBLE_COMPLEX.h"

