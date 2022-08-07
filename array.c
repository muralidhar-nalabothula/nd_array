#include "nd_array_src.h"

static void FUNCTION(slice_internal, TYPE_S) ( const ND_indices * restrict start_idx, const ND_indices * restrict end_idx, const ND_indices * restrict step_idx, \
        const ND_indices * restrict stride_F, const ND_indices * restrict stride_S, ND_indices idx_F,  ND_indices idx_S, const ND_indices ndim, \
        const ND_indices idim, const TYPE_L * restrict arrF, TYPE_L * restrict arrS);




/* Array operation function */
/****************************************************************************************************/

/* Function to get element pointer of an array */
TYPE_L * FUNCTION(ele, TYPE_S) (const ARRAY_T(TYPE_S) * nd_arr_in, const ND_indices * dimensions)
{
    /*returns the pointer to the particular element . so we can set and get the elements*/
    //
    if (nd_arr_in->rank == NULL) error_msg("uninitialized input array for nd_ele function");
    if (nd_arr_in->data == NULL) error_msg("NULL data array received for nd_ele function");

    ND_indices index_arr = 0, arg_val; 

    for(ND_indices i = 0; i < *(nd_arr_in->rank); ++i)
    {
        arg_val = dimensions[i];

        if (arg_val < ((nd_arr_in->dims)[i]) ) index_arr = index_arr + arg_val * ((nd_arr_in->strides)[i]);

        else error_msg("Array out of bound");
    }

    return (nd_arr_in->data) + index_arr;
}





/* Function to get Size of an array .*/
ND_indices FUNCTION(size, TYPE_S) (const ARRAY_T(TYPE_S) * nd_arr_in)
{   
    if (nd_arr_in->rank == NULL) error_msg("uninitialized input array for nd_size function");

    ND_indices size ; 

    if (*(nd_arr_in->rank) == (ND_indices) 0) size =  (ND_indices) 1 ;
    
    else size = *(nd_arr_in->dims) * *(nd_arr_in->strides);
    

    return size; 
}


/* Function to set all  elements of an array to constant*/
void FUNCTION(set_all, TYPE_S) (ARRAY_T(TYPE_S) * nd_arr_in, const TYPE_L set_constant )
{   
    if (nd_arr_in->rank == NULL) error_msg("uninitialized input array for nd_ele function");
    if (nd_arr_in->data == NULL) error_msg("NULL data array received for nd_ele function");

    const ND_indices arr_in_size = ND_FUNCTION_CALL(size, TYPE_S)(nd_arr_in);
    for (ND_indices i = 0 ; i < arr_in_size ; ++i ) 
    {
        (nd_arr_in->data)[i] = set_constant ;
    }

}

/* Function to get reshape an array .*/

void FUNCTION(reshape, TYPE_S) (const ARRAY_T(TYPE_S) * nd_arr_in, ARRAY_T(TYPE_S) * nd_arr_out)
{

    if (nd_arr_in->rank == NULL) error_msg("uninitialized input array for reshape function");
    if (nd_arr_out->rank == NULL) error_msg("uninitialized output array for reshape function");
    if (nd_arr_in->data == NULL) error_msg("NULL data array received for reshape function");
    if (nd_arr_out->data != NULL) error_msg("Cannot pass array with data to reshape function. This leads to memory leak");

    nd_arr_out->owner = false ;

    /*check size*/
    if ( ND_FUNCTION_CALL(size, TYPE_S) (nd_arr_in) !=  ND_FUNCTION_CALL(size, TYPE_S) (nd_arr_out) ) \
                                    error_msg("Size mismatch of reshaped array and original array");

    nd_arr_out->data = nd_arr_in->data;

}


/* Function to strip off the first n dimensions */
void FUNCTION(strip_dims, TYPE_S) (const ARRAY_T(TYPE_S) * nd_arr_in, const ND_indices n_dims_strip, const ND_indices * stripped_idxs, ARRAY_T(TYPE_S) * nd_arr_out)
{
    /* 
        Strips first n dims: 
        Ex: 
        for 5 dim array A, when we strip first two dims we get A[i,j,:,:,:]. n_dims_strip=2, stripped_idxs={i,j}
     */
    
    if (nd_arr_in->rank == NULL) error_msg("uninitialized input array for strip_dims function");
    if (nd_arr_in->data == NULL) error_msg("NULL data array received for strip_dims function");

    if (nd_arr_out->rank == NULL) error_msg("uninitialized output array for strip_dims function");
    if (nd_arr_out->data != NULL) error_msg("Cannot pass array with data to strip_dims function. This leads to memory leak");

    nd_arr_out->owner = false ;

    if (n_dims_strip >= (nd_arr_in->rank)[0]) error_msg("Number of dims to stip must be less than the rank of original array");

    ND_indices sub_tensor_rank = (nd_arr_in->rank)[0]-n_dims_strip ;
    
    if ( *(nd_arr_out->rank) != sub_tensor_rank) error_msg("Rank mismatch in nd_strip_dims function") ;

    for (ND_indices i = 1; i <= sub_tensor_rank; ++i)
    {
        if ((nd_arr_out->dims)[sub_tensor_rank - i] != (nd_arr_in->dims)[(nd_arr_in->rank)[0] - i]) error_msg("Dimension mismatch in nd_strip_dims function") ; 
    }

    ND_indices index_arr = 0, arg_val; 
    
    for(ND_indices i = 0; i < n_dims_strip; ++i)
    {
        arg_val = stripped_idxs[i];

        if (arg_val < ((nd_arr_in->dims)[i]) )
        {
            index_arr = index_arr + arg_val * ((nd_arr_in->strides)[i]);
        }
        else error_msg("Array out of bound");
    }

    nd_arr_out->data = nd_arr_in->data + index_arr;
}


/* Availble for user :) */
/* FUnction to slice a subtensor and copy to new tensor*/
void FUNCTION(slice, TYPE_S) (const ND_indices * start_idx, const ND_indices * end_idx, const ND_indices * step_idx, \
                            const ARRAY_T(TYPE_S) * nd_arr_in, ARRAY_T(TYPE_S) * nd_arr_out )
{
    if (nd_arr_in->rank == NULL) error_msg("uninitialized input array for slice function");
    if (nd_arr_in->data == NULL) error_msg("NULL data array received for slice function");

    if (nd_arr_out->rank == NULL) error_msg("uninitialized output array for slice function");
    if (nd_arr_out->data == NULL) error_msg("NULL data array received for slice function");

    ND_indices slice_rank = (ND_indices) 0;
    ND_indices arr_in_rank = *(nd_arr_in->rank);

    ND_indices * sliced_dims_temp = malloc(3*arr_in_rank * sizeof(ND_indices)); // create both temp and sliced and strides
    ND_indices * sliced_dims = sliced_dims_temp + arr_in_rank ; 
    ND_indices * out_strides = sliced_dims_temp + 2*arr_in_rank ;
    
    for (ND_indices i =0 ; i < arr_in_rank; i++ )
    {

        if ((end_idx[i] <= start_idx[i] || end_idx[i] > (nd_arr_in->dims)[i]) || start_idx[i] >= (nd_arr_in->dims)[i] ) \
                error_msg("Slicing index out of bound or ending index is larger than the starting index when slicing array.");

        else
        {
            sliced_dims_temp[i] = 1 + (end_idx[i]-start_idx[i] - 1)/step_idx[i];

            if (sliced_dims_temp[i] > 1)
            {
                sliced_dims[slice_rank] = sliced_dims_temp[i];
                slice_rank++;
            }
        }
    }

    ND_indices slice_size = 1;

    for (ND_indices i = 1; i < arr_in_rank+1; ++i)
    {
        out_strides[arr_in_rank - i] = slice_size;
        //printf("%zu \n", slice_size );
        slice_size = slice_size*sliced_dims_temp[arr_in_rank-i];
    }

    if (slice_rank == 0 || *(nd_arr_in->rank) == 0 ) error_msg("Ranks of sliced array or input array cannot be zero. Use nd_ele function");

    if (slice_size !=  ND_FUNCTION_CALL(size, TYPE_S) (nd_arr_out))
    {
        fprintf(stdout, "# [ Error !!!] : Sliced array size is %zu but the given sliced array size is %zu  \n",\
            slice_size, ND_FUNCTION_CALL(size, TYPE_S) (nd_arr_out) );
        error_msg("Incompatible sized in nd_slice function");
    }

    if (slice_rank !=  *(nd_arr_out->rank) ) error_msg("Ranks of inputed sliced array are wrong");

    for (ND_indices i = 0; i < slice_rank; ++i)
    {
        if (sliced_dims[i] != (nd_arr_out->dims)[i] ) error_msg("Dimensions of inputed sliced array are wrong");
    }
    
    ND_FUNCTION_CALL(slice_internal, TYPE_S) (start_idx, end_idx,step_idx,nd_arr_in->strides,out_strides, (ND_indices) 0,\
      (ND_indices) 0, *(nd_arr_in->rank), (ND_indices) 0, nd_arr_in->data, nd_arr_out->data);

    free(sliced_dims_temp);
}



// void FUNCTION(swap_owner, TYPE_S) ( ARRAY_T(TYPE_S) * nd_arr_new, ARRAY_T(TYPE_S) * nd_arr_old )
//{
//     nd_arr_new->owner = true;
//     nd_arr_old->owner = false;
// }

/* Function to transpose copy from one array to other. Same as like numpy.transpose()*/
void FUNCTION(tranpose, TYPE_S) (const ARRAY_T(TYPE_S) * nd_arr_in, const ND_indices * order, ARRAY_T(TYPE_S) * nd_arr_out)
{
    /*Note that you must input and out put array. No array is created. it just sets elements.
    This is not very cache friendly way !
    */
    if (nd_arr_in->rank == NULL) error_msg("uninitialized input array for transpose function");
    if (nd_arr_in->data == NULL) error_msg("NULL data array received for transpose function");

    if (nd_arr_out->rank == NULL) error_msg("uninitialized output array for transpose function");
    if (nd_arr_out->data == NULL) error_msg("NULL data array received for transpose function");

    if (ND_FUNCTION_CALL(size, TYPE_S) (nd_arr_in) != ND_FUNCTION_CALL(size, TYPE_S) (nd_arr_out)) \
                                error_msg("Incompatible size between input and transposed array"); //check dim

    if (*(nd_arr_in->rank) != *(nd_arr_out->rank)) error_msg("Incompatible size between input and transposed array"); // check rank

    // check if all transpose indices are in range
    for ( ND_indices i = 0 ; i < *(nd_arr_in->rank); ++i)
    {
        if ( order[i] >= *(nd_arr_in->rank) ) error_msg("Transpose indices must be less than rank.");
    }

    /* check dims*/
    for ( ND_indices i = 0 ; i < *(nd_arr_in->rank); ++i)
    {
        if (  (nd_arr_in->dims)[order[i]] != (nd_arr_out->dims)[i] ) error_msg("Incompatible transpose order.");
    }

    ND_indices rank = *(nd_arr_in->rank );
    ND_indices * idx_array = malloc(( 1 + rank ) * sizeof(ND_indices)); /*last value of index is remainder*/

    /* set the data*/
    ND_indices out_idx ;

    //
    for (ND_indices i = 0; i < ND_FUNCTION_CALL(size, TYPE_S) (nd_arr_in); ++i)
    {
        idx_array[rank] = i;
        out_idx = 0;

        for (ND_indices j = 0; j < rank; ++j)
        {
            idx_array[j] = idx_array[rank]/(nd_arr_in->strides)[j] ;
            idx_array[rank] = idx_array[rank] % (nd_arr_in->strides)[j] ;
        }
        
        for (ND_indices j = 0; j < rank; ++j)
        {
            out_idx = out_idx + idx_array[order[j]] * (nd_arr_out->strides)[j];
        }

        (nd_arr_out->data)[out_idx] = (nd_arr_in->data)[i] ; 
    }

    free(idx_array);

}




/* Function to deepcopy from one array to other*/
void FUNCTION(copy, TYPE_S) (const ARRAY_T(TYPE_S) * nd_arr_in, ARRAY_T(TYPE_S) * nd_arr_out)
{

    if (nd_arr_in->rank == NULL) error_msg("uninitialized input array for copy function");
    if (nd_arr_in->data == NULL) error_msg("NULL data array received for copy function");

    if (nd_arr_out->rank == NULL) error_msg("uninitialized output array for copy function");
    if (nd_arr_out->data == NULL) error_msg("NULL data array received for copy function");

    if (ND_FUNCTION_CALL(size, TYPE_S) (nd_arr_in) != ND_FUNCTION_CALL(size, TYPE_S) (nd_arr_out)) error_msg("Incompatible size between input and copy array");

    if (*(nd_arr_in->rank) != *(nd_arr_out->rank)) error_msg("Incompatible rank between input and copy array");


    /* check dims*/
    for ( ND_indices i = 0 ; i < *(nd_arr_in->rank); ++i)
    {
        if (  (nd_arr_in->dims)[i] != (nd_arr_out->dims)[i] ) error_msg("Incompatible dims for input and copy arrays.");

    }
    memcpy(nd_arr_out->data, nd_arr_in->data, sizeof(TYPE_L)* ND_FUNCTION_CALL(size, TYPE_S) (nd_arr_in));
    memcpy(nd_arr_out->rank, nd_arr_in->rank, sizeof(ND_indices) * ( (1 + 2* (nd_arr_in->rank)[0] ) ));
}


/************************************************************************ STATIC FUNCTIONS **********************************************************************/

/* Not availble for user */

static void FUNCTION(slice_internal, TYPE_S) ( const ND_indices * restrict start_idx, const ND_indices * restrict end_idx, const ND_indices * restrict step_idx, \
        const ND_indices * restrict stride_F, const ND_indices * restrict stride_S, ND_indices idx_F,  ND_indices idx_S, const ND_indices ndim, \
        const ND_indices idim, const TYPE_L * restrict arrF, TYPE_L * restrict arrS)
{
    
    //if (idim<ndim) printf("%zu %zu \n", idim, ndim);
    if (idim == ndim) arrS[idx_S] = arrF[idx_F] ;
    //
    else if ( (idim == ndim-1) && (step_idx[idim] == (ND_indices) 1) )
    {
        memcpy(arrS+idx_S, arrF + idx_F + (start_idx[idim] * stride_F[idim]), (end_idx[idim]-start_idx[idim])*sizeof(TYPE_L));
    }
    //
    else 
    {
        ND_indices idx_F1, idx_S1;

        for (ND_indices i = start_idx[idim] ; i < end_idx[idim] ; i = i + step_idx[idim] ){
            // if (idim<ndim) printf("%zu %zu %zu \n", i,start_idx[idim], end_idx[idim]);
            if (idim == (ND_indices) 0)
            {
                idx_F = (ND_indices) 0  ; 
                idx_S = (ND_indices) 0  ; 
            }
            idx_F1 = idx_F + (i * stride_F[idim]) ;
            idx_S1 = idx_S + (((i - start_idx[idim] )/step_idx[idim]) * stride_S[idim]) ;
            //printf("%zu %zu %zu \n",*idx_S, *idx_F , stride_S[idim] );
            ND_FUNCTION_CALL(slice_internal, TYPE_S)(start_idx, end_idx, step_idx, stride_F, stride_S, idx_F1, idx_S1, ndim, idim + (ND_indices) 1, arrF, arrS) ; 
        }
    }
}