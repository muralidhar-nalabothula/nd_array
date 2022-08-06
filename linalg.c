#include "nd_array_src.h"


/* Linear algebra function (blas, tblis, lapack ...)*/
/****************************************************************************************************/

#if defined(COMPILE_ND_DOUBLE_COMPLEX) || defined(COMPILE_ND_SINGLE_COMPLEX) || defined(COMPILE_ND_SINGLE_FLOAT) || defined(COMPILE_ND_DOUBLE)

static CBLAS_TRANSPOSE get_gemmn_T(char Trans);
static void nd_free_tblis(tblis_tensor * A_arr);
static void FUNCTION(to_tblis, TYPE_S)(const ARRAY_T(TYPE_S) * nd_arr_A, tblis_tensor * A_arr);
static void FUNCTION(to_tblis_scaled, TYPE_S)(const ARRAY_T(TYPE_S) * nd_arr_A, tblis_tensor * A_arr, const TYPE_L alpha);
static void check_string_dims(char * str_A, char * str_B, ND_indices * dims_A,  ND_indices * dims_B);
static void check_all_indices(char * str_A, char * str_B);




void FUNCTION(matmul, TYPE_S) (const char TransA, const char TransB, const ARRAY_T(TYPE_S) * nd_arr_A, const ARRAY_T(TYPE_S) * nd_arr_B, ARRAY_T(TYPE_S) * nd_arr_C,
                const TYPE_L alpha, const TYPE_L beta, const ND_indices * A_idx,  const ND_indices * B_idx,  const ND_indices * C_idx)

{
    /* THis will perform the following matrix multiplication.
    ----- A({A_idx},m,k)@ B({B_idx},k,n) = C({C_idx},m,n)  ------ */

    /* Some checks on data*/
    
    if ( (((nd_arr_A->rank) == NULL )  || ((nd_arr_B->rank) == NULL )) || ((nd_arr_C->rank) == NULL ) ) \
                                                    error_msg("Cannot accept uninitilized array") ; 
    
    if ( (((nd_arr_A->data) == NULL )  || ((nd_arr_B->data) == NULL )) || ((nd_arr_C->data) == NULL ) ) \
                                                    error_msg("Cannot accept NULL array. allocate them beforing passing to matmul");

    /*Find the max BLAS int to check for overflow of indices */
    intmax_t max_len_x;
    ND_indices max_len;
    
    for (max_len_x=INTMAX_MAX; (BLAS_INT)max_len_x!=max_len_x; max_len_x/=2);
    
    max_len = (ND_indices)max_len_x; 

    if ( ((*(nd_arr_A->rank) < (ND_indices) 2)  || (*(nd_arr_B->rank) < (ND_indices) 2)) || (*(nd_arr_C->rank) < (ND_indices) 2) ) \
                                                    error_msg("Matmul only accepts arrays with atleast rank 2") ;

    if ( (A_idx == NULL) && (*(nd_arr_A->rank)  != 2)) error_msg("Null is passed for A with dim > 2") ;

    if ( (B_idx == NULL) && (*(nd_arr_B->rank)  != 2)) error_msg("Null is passed for B with dim > 2") ;

    if ( (C_idx == NULL) && (*(nd_arr_C->rank)  != 2)) error_msg("Null is passed for C with dim > 2.") ;

    ND_indices ldA = (nd_arr_A->dims)[*(nd_arr_A->rank)-1] ;
    ND_indices ldB = (nd_arr_B->dims)[*(nd_arr_B->rank)-1];
    ND_indices ldC = (nd_arr_C->dims)[*(nd_arr_C->rank)-1];

    ND_indices m1, k1, k2, n2, m3,n3 ;

    if (TransA == 'C' || TransA == 'T')
    {
        m1 = (nd_arr_A->dims)[*(nd_arr_A->rank)-1]; /*  m,k--> k,m */
        k1 = (nd_arr_A->dims)[*(nd_arr_A->rank)-2];
    }

    else
    {
        k1 = (nd_arr_A->dims)[*(nd_arr_A->rank)-1]; /*  m,k--> k,m */
        m1 = (nd_arr_A->dims)[*(nd_arr_A->rank)-2];
    }

    if (TransB == 'C' || TransB == 'T')
    {
        k2 = (nd_arr_B->dims)[*(nd_arr_B->rank)-1]; /*  m,k--> k,m */
        n2 = (nd_arr_B->dims)[*(nd_arr_B->rank)-2];
    }

    else
    {
        n2 = (nd_arr_B->dims)[*(nd_arr_B->rank)-1]; /*  m,k--> k,m */
        k2 = (nd_arr_B->dims)[*(nd_arr_B->rank)-2];
    }

    n3 = (nd_arr_C->dims)[*(nd_arr_C->rank)-1]; /*  m,k--> k,m */
    m3 = (nd_arr_C->dims)[*(nd_arr_C->rank)-2];
    
    if ( ((k1 != k2) || (m1 !=m3)) || (n2 !=n3)  )
    {
        fprintf(stdout, "# [ Error !!!] : Incompatible matrix dimensions. Matrix product not possible (%zu, %zu) , (%zu, %zu) -> (%zu, %zu)! \n", 
        m1, k1, k2, n2, m3,n3);
        error_msg("Incompatible matrix dimensions. Matrix product not possible");
    }

    ND_indices idx_A = 0, idx_B = 0, idx_C = 0; 

    for (ND_indices i = 0; i< *(nd_arr_A->rank) - (ND_indices) 2; ++i)
    {
        if ( A_idx[i] >= (nd_arr_A->dims)[i] ) error_msg("Provided indices of A matrix for nd_matmul function are out of bounds");

        idx_A = idx_A + ((nd_arr_A->strides)[i] * A_idx[i]) ; 
    }

    // error_msg("");

    for (ND_indices i = 0; i< *(nd_arr_B->rank) - (ND_indices) 2; ++i)
    {
        if ( B_idx[i] >= (nd_arr_B->dims)[i] ) error_msg("Provided indices of B matrix for nd_matmul function are out of bounds");

        idx_B = idx_B + ((nd_arr_B->strides)[i] * B_idx[i]) ; 
    }

    for (ND_indices i = 0; i< *(nd_arr_C->rank) - (ND_indices) 2; ++i){
        if ( C_idx[i] >= (nd_arr_C->dims)[i] ) error_msg("Provided indices of C matrix for nd_matmul function are out of bounds");

        idx_C = idx_C + ((nd_arr_C->strides)[i] * C_idx[i]) ; 
    }
    /* call the gemm*/
    // 
    if ( (((((m1 >= max_len ||  n2 >= max_len) ||  k1 >= max_len) ||  ldA >= max_len) ||  ldB >= max_len) ||  ldC >= max_len)  ) \
                                            error_msg("BLAS indices Overflow. Compile the blas library with higher bit indices!");
    
    // Call blas gemm

    BLAS_CALL(gemm , TYPE_S) (CblasRowMajor, get_gemmn_T(TransA), get_gemmn_T(TransB), (BLAS_INT) m1 , (BLAS_INT) n2, (BLAS_INT) k1 ,\
                BLAS_POINTER(alpha), (nd_arr_A->data) + idx_A, (BLAS_INT) ldA,  (nd_arr_B->data) + idx_B,  (BLAS_INT) ldB,  BLAS_POINTER(beta), (nd_arr_C->data) + idx_C,  (BLAS_INT) ldC);
}

/****************************************************************************************************/


void FUNCTION(einsum, TYPE_S) (char * einsum_indices, ARRAY_T(TYPE_S) * nd_arrA, ARRAY_T(TYPE_S) * nd_arrB, ARRAY_T(TYPE_S) * nd_arrC, const TYPE_L alpha, const TYPE_L beta)
{
    //
    
    tblis_tensor A_arr, B_arr, C_arr ;

    size_t einsum_indices_len = strlen(einsum_indices) ; 
    
    if (einsum_indices_len>400) error_msg("Max length of Einsum indices string is 400");


    char * str_A = malloc(sizeof(char) * 5*400);
    char * str_B = str_A+400 ;
    char * str_C = str_A+800 ; 
    char * str_AB_cat = str_A+1200 ; // stores str_A + str_B
    //
    //
    size_t ldaA = 0, ldaB = 0, ldaC = 0 ; 
    size_t token_counter = 0;

    for (size_t i =0 ; i<einsum_indices_len; ++i)
    {
        if (isspace(einsum_indices[i]) == 0)
        {
            if (token_counter == 0)
            {
                if (einsum_indices[i] == '-' || einsum_indices[i+1] == '>' ) error_msg("Incompatible einsum indices provided") ;

                else if (einsum_indices[i] == ',') token_counter++ ; 

                else
                {
                    str_A[ldaA] = einsum_indices[i];
                    ldaA++ ;
                }
            }
            else if (token_counter == 1)
            {
                if (einsum_indices[i] == ',') error_msg("Incompatible einsum indices provided");
                //
                else if (einsum_indices[i] == '-')
                {
                    if (i <einsum_indices_len-1)
                    {
                        if (einsum_indices[i+1] != '>') error_msg("Incompatible einsum indices provided");
                    }
                    else error_msg("Incompatible einsum indices provided");
                }
                //
                else if (einsum_indices[i] == '>')
                {
                    if (i>0)
                    {
                        if (einsum_indices[i-1] != '-') error_msg("Incompatible einsum indices provided");

                        else token_counter++ ;
                    }
                    else error_msg("Incompatible einsum indices provided");
                }
                //
                else{
                    str_B[ldaB] = einsum_indices[i];
                    ldaB++ ; 
                }
            }
            else if (token_counter == 2)
            {
                if ((einsum_indices[i] == '>' || einsum_indices[i] == '-') || einsum_indices[i] == ',') error_msg("Incompatible einsum indices provided");

                else
                {
                    str_C[ldaC] = einsum_indices[i];
                    ldaC++;
                }
                
            }
            else error_msg("Incompatible einsum indices provided");

        }
    }

    str_C[ldaC] = '\0'; str_B[ldaB] = '\0'; str_A[ldaA] = '\0';

    if (token_counter != 2) error_msg("Incompatible einsum indices provided");

    if (  ((strlen(str_A) != (nd_arrA->rank)[0]) || (strlen(str_B) != (nd_arrB->rank)[0]) ) || (strlen(str_C) != (nd_arrC->rank)[0])) \
                            error_msg("Rank of input arrays doesn't match with given indices");

    if (  strlen(str_A) == (size_t)0 || strlen(str_B) == (size_t)0 ) error_msg("Rank of input arrays cannot be 0. Use scale/sum function instead");

    strcpy(str_AB_cat,str_A);
    strcat(str_AB_cat,str_B);

    check_all_indices(str_AB_cat,str_C);//check if the final indices are present in input arrays

    check_string_dims(str_A, str_A, nd_arrA->dims, nd_arrA->dims); // check is repeated indices are consistant in same array
    check_string_dims(str_B, str_B, nd_arrB->dims, nd_arrB->dims);
    check_string_dims(str_C, str_C, nd_arrC->dims, nd_arrC->dims);

    check_string_dims(str_A, str_B, nd_arrA->dims, nd_arrB->dims); // check is repeated indices are consistant in other arrays
    check_string_dims(str_A, str_C, nd_arrA->dims, nd_arrC->dims);
    check_string_dims(str_B, str_C, nd_arrB->dims, nd_arrC->dims);

    

    //to_tblis_scaled
    if (  ND_FUNCTION_CALL(size, TYPE_S) (nd_arrA)  < ND_FUNCTION_CALL(size, TYPE_S) (nd_arrB) )
    {
        ND_FUNCTION_CALL(to_tblis_scaled, TYPE_S)(nd_arrA, &A_arr, alpha);
        ND_FUNCTION_CALL(to_tblis, TYPE_S)(nd_arrB, &B_arr);
    }
    else
    {
        ND_FUNCTION_CALL(to_tblis, TYPE_S)(nd_arrA, &A_arr);
        ND_FUNCTION_CALL(to_tblis_scaled, TYPE_S)(nd_arrB, &B_arr, alpha);
    }
    //
    if (  strlen(str_C) == (size_t)0  )
    {
        tblis_scalar C_scalar ;
        TBLIS_CALL(init_scalar,TYPE_S)(&C_scalar, (TYPE_L)0.0f);
        tblis_tensor_dot(NULL, NULL,&A_arr, str_A, &B_arr, str_B, &C_scalar);
        (nd_arrC->data)[0] = beta*(nd_arrC->data)[0] + C_scalar.data.TYPE_S  ;
    }
    //
    else
    {
        ND_FUNCTION_CALL(to_tblis_scaled, TYPE_S)(nd_arrC, &C_arr, beta);
        tblis_tensor_mult(NULL, NULL, &A_arr, str_A, &B_arr, str_B, &C_arr, str_C);
    }

    //
    free(str_A);
    nd_free_tblis(&A_arr);
    nd_free_tblis(&B_arr);

    if ( strlen(str_C) != (size_t)0 ) nd_free_tblis(&C_arr);

}



/************************************* INTERNAL STATIC HELPER FUNCTIONS ******************************************************/

/* FUNCTION OF MATMUL to identify transpose/conjuate*/
static CBLAS_TRANSPOSE get_gemmn_T(char Trans)
{
    if (Trans == 'N')      return CblasNoTrans ;
    else if (Trans == 'T') return CblasTrans;
    else if (Trans == 'C') return CblasConjTrans;
    else 
    {
        error_msg("Can only take 'C' or 'T' or 'N' for Trans input");
        return CblasNoTrans ;
    }
}



static void check_string_dims(char * str_A, char * str_B, ND_indices * dims_A,  ND_indices * dims_B)
{   
    /* This Function checks if the string indices are consistant with the dims*/
    /* Before calling this function, you must check the ranks of A and B match with string lenghts*/
    if (str_A == NULL || str_B == NULL) error_msg("Cannot pass NULL strings to check_string_dims");

    char * str_found_ptr;
    size_t index_str_A;
    //
    size_t str_len_B ;
    //
    str_len_B = strlen(str_B);
    //
    for (size_t i = 0; i < str_len_B; ++i )
    {
        str_found_ptr = strchr(str_A, str_B[i]);
        if (str_found_ptr != NULL)
        {   
            index_str_A = (size_t)(str_found_ptr - str_A);
            if (dims_A[index_str_A] != dims_B[i]) error_msg("Dimensions are inconsistant with indices provided");
        }
    }

}


static void check_all_indices(char * str_A, char * str_B)
{   
    /* This Function checks if the all string indices of B are present in A are consistant with the dims*/
    /* To pass, mutiple strings, concatnate and pass to str_A*/
    if (str_A == NULL || str_B == NULL) error_msg("Cannot pass NULL strings to check_all_indices");

    char * str_found_ptr;
    size_t str_len_B ;

    str_len_B = strlen(str_B);

    for (size_t i = 0; i < str_len_B; ++i )
    {
        str_found_ptr = strchr(str_A, str_B[i]);
        if (str_found_ptr == NULL) error_msg("Final indices are not present in the initial indices");
    }

}

/* Convert nd_array to tblis array*/
static void FUNCTION(to_tblis, TYPE_S)(const ARRAY_T(TYPE_S) * nd_arr_A, tblis_tensor * A_arr)
{   
    /*Must be freed after using*/

    if (nd_arr_A->rank == NULL) error_msg("Cannot convert a nd_array to blis array");
    if (nd_arr_A->data == NULL) error_msg("Cannot convert a nd_array to blis array");

    /* Get max values to check for overflow*/
    intmax_t max_len, max_stride;

    for (max_len=INTMAX_MAX; (len_type)max_len!=max_len; max_len/=2);
    for (max_stride=INTMAX_MAX; (stride_type)max_stride!=max_stride; max_stride/=2);
    
    
    ND_indices rank = (nd_arr_A->rank)[0];
    len_type * blis_dims         =  malloc(sizeof(len_type)*rank) ;
    stride_type * blis_strides =  malloc(sizeof(stride_type)*rank) ;
    
    for (ND_indices i =0; i<rank; ++i )
    {
        if ( ((nd_arr_A->dims)[i] >= (ND_indices)max_len  ) || ((nd_arr_A->strides)[i] >= (ND_indices)max_stride ) ) \
                error_msg("Indices overflow. Complile tblis with larger bit indices for len_type and stride_type");

        else
        {
            blis_dims[i] =    (len_type) (nd_arr_A->dims)[i] ;
            blis_strides[i] = (stride_type) (nd_arr_A->strides)[i] ;
        }
    }

    TBLIS_CALL(init_tensor,TYPE_S) (A_arr, (unsigned) rank , blis_dims, nd_arr_A->data, blis_strides);

    // free(blis_dims);
    // free(blis_strides);
}



/* Convert nd_array to tblis array with scaling*/
static void FUNCTION(to_tblis_scaled, TYPE_S)(const ARRAY_T(TYPE_S) * nd_arr_A, tblis_tensor * A_arr, const TYPE_L alpha)
{
    /*Must be freed after using*/

    if (nd_arr_A->rank == NULL) error_msg("Cannot convert a nd_array to blis array");
    if (nd_arr_A->data == NULL) error_msg("Cannot convert a nd_array to blis array");

    /* Get max values to check for overflow*/
    intmax_t max_len, max_stride;

    for (max_len=INTMAX_MAX; (len_type)max_len!=max_len; max_len/=2);
    for (max_stride=INTMAX_MAX; (stride_type)max_stride!=max_stride; max_stride/=2);
    
    
    ND_indices rank = (nd_arr_A->rank)[0];
    len_type * blis_dims         =  malloc(sizeof(len_type)*rank) ;
    stride_type * blis_strides =  malloc(sizeof(stride_type)*rank) ;
    
    for (ND_indices i =0; i<rank; ++i )
    {
        if ( ((nd_arr_A->dims)[i] >= (ND_indices)max_len  ) || ((nd_arr_A->strides)[i] >= (ND_indices)max_stride ) ) \
                error_msg("Indices overflow. Complile tblis with larger bit indices for len_type and stride_type");
        
        else
        {
            blis_dims[i] =    (len_type) (nd_arr_A->dims)[i] ;
            blis_strides[i] = (stride_type) (nd_arr_A->strides)[i] ;
        }
    }

    TBLIS_CALL(init_tensor_scaled,TYPE_S) (A_arr, alpha, (unsigned) rank , blis_dims, nd_arr_A->data, blis_strides);

    // free(blis_dims);
    // free(blis_strides);
}



/* Free tblis arrays*/
static void nd_free_tblis(tblis_tensor * A_arr)
{   /*
    Note: This function only free, dimensions and strinds, It DOES NOT free the main data.
    */
    if ( (A_arr->len == NULL )|| (A_arr->stride == NULL)) error_msg("tblis_tensor already freed");
    //
    else
    {
        free(A_arr->len);
        free(A_arr->stride);
        A_arr->len = NULL ;
        A_arr->stride = NULL ;
    }
}


#endif


/* Level -1 Blas*/

/* Level -2 Blas*/

/* Level -3 Blas*/
/* Gemm. i-> j-->k  C(.....,i,k) = A(....,i,j)*B(.....,j,k) 
    C(...,m,n) = A(....,m,k)@B(....,k,n),
    "(2,3,:,:) (3,5,:,:) -> (4,5,:,:)", m, n, k
*/
