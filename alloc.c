#include "nd_array_src.h"

/********************************************************************************************************************************************/
/********************************************************* INITIALIZATION  FUNCTIONS ********************************************************/


void FUNCTION(init, TYPE_S) (ARRAY_T(TYPE_S) * nd_arr_in, const ND_indices rank, const ND_indices * dimensions)
{
    /* first initiate the ND_indices array. Only stirdes, rank and dim arrays are created
    */
    ND_indices size_arr = 1;

    nd_arr_in->owner = false ; /* When initialized, set the data to false.*/ 

    ND_indices * rank_array = malloc((1 + 2*rank) * sizeof(ND_indices)); /* (2*rank + 1) elements*/

    if (rank_array == NULL) error_msg("Failed to allocate dimensions for nd_malloc array");
    
    nd_arr_in->rank = rank_array;
    //
    if (rank == (ND_indices)0)
    {
        nd_arr_in->dims = NULL;
        nd_arr_in->strides = NULL;
    }
    else
    {
        nd_arr_in->dims = rank_array+1;
        nd_arr_in->strides = rank_array+1+rank;
    }
    rank_array[0] = rank;
    /* Get all the args of type ND_indices*/
    
    if (rank != (ND_indices)0) memcpy(rank_array+1,dimensions,sizeof(ND_indices)*rank);
    /* set strides */ 
    for (ND_indices i = 0; i < rank; ++i)
    {
        rank_array[2*rank - i] = size_arr;
        size_arr = size_arr*rank_array[rank - i];
    }
    //
    nd_arr_in->data = NULL ; 
    //
}


void FUNCTION(uninit, TYPE_S) (ARRAY_T(TYPE_S) * nd_arr_in)
{
    /* Note that rank, dims, strides are all created with single malloc with rank being first pointer
        So Only free data and rank pointers
    */
    if (nd_arr_in->rank == NULL ) error_msg("Array already uninitialized");

    else if (nd_arr_in->owner) error_msg("Free the array before uninitialization");

    else 
    {
        free(nd_arr_in->rank);
        nd_arr_in->rank = NULL ; 
    }
}



/********************************************************************************************************************************************/
/********************************************************* MALLOC/CALLOC/FREE  FUNCTIONS ****************************************************/

/* malloc function*/
void FUNCTION(malloc, TYPE_S) (ARRAY_T(TYPE_S) * nd_arr_in)
{
    /* Allocated memory for the data. Must call free function.
    */
    if (nd_arr_in->owner == true) error_msg("Overwriting already allocated array (i.e, Owner array). Possible memory leak !") ;

    ND_indices size_arr = 1;
    //
    if (*(nd_arr_in->rank) == (ND_indices) 0)
    {
        size_arr =  (ND_indices) 1 ;
    }

    else
    {
        size_arr = *(nd_arr_in->dims) * *(nd_arr_in->strides);
    }
    //

    nd_arr_in->owner = true ; /* Set this to true as it owns the data created*/
    /* Create data pointer*/
    nd_arr_in->data = malloc(size_arr * sizeof(TYPE_L));

    if (nd_arr_in->data == NULL) error_msg("Failed to allocate data for nd_malloc array") ;
}



/* calloc function*/
void FUNCTION(calloc, TYPE_S) (ARRAY_T(TYPE_S) * nd_arr_in)
{
    /* Allocated memory for the data. Must call free function.
    */
    if (nd_arr_in->owner == true) error_msg("Overwriting already allocated array (i.e, Owner array). Possible memory leak !") ;

    ND_indices size_arr = 1;
    //
    if (*(nd_arr_in->rank) == (ND_indices) 0)
    {
        size_arr =  (ND_indices) 1 ;
    }

    else
    {
        size_arr = *(nd_arr_in->dims) * *(nd_arr_in->strides);
    }
    //

    nd_arr_in->owner = true ; /* Set this to true as it owns the data created*/
    /* Create data pointer*/
    nd_arr_in->data = calloc(size_arr, sizeof(TYPE_L));

    if (nd_arr_in->data == NULL) error_msg("Failed to allocate data for nd_malloc array") ;
}




/* Free function */
void FUNCTION(free, TYPE_S) (ARRAY_T(TYPE_S) * nd_arr_in)
{
    /* 
    Must be called when nd_malloc_? or nd_calloc_? functions are called, else memory leak !
    */
    if (nd_arr_in->owner)
    {   
        if (nd_arr_in->data == NULL) error_msg("Trying to free array with NULL data");
        free(nd_arr_in->data);
        nd_arr_in->data = NULL ; 
        nd_arr_in->owner = false;
    }
}

/********************************************************************************************************************************************/
/********************************************************************************************************************************************/


/********************************************************************************************************************************************/
/********************************************************* INITIALIZATION  FUNCTIONS ********************************************************/

/* These functions are to make life bit earier you can directly initialize these for transpose, init_slice, init_strip_dims*/

/* Function to initiate transpose array directly*/
void FUNCTION(init_tranpose, TYPE_S) (const ARRAY_T(TYPE_S) * nd_arr_in, const ND_indices * order, ARRAY_T(TYPE_S) * nd_arr_out)
{
    if (nd_arr_in->rank == NULL) error_msg("uninitialized input array for init_tranpose function");
    
    // check if all transpose indices are in range

    for ( ND_indices i = 0 ; i < *(nd_arr_in->rank); ++i)
    {
        if ( order[i] >= *(nd_arr_in->rank) ) error_msg("Transpose indices must be less than rank.");
    }

    ND_indices rank_out;

    rank_out = (nd_arr_in->rank)[0] ;
    //
    ND_indices * dimensions_out = malloc(sizeof(ND_indices) * rank_out);

    for ( ND_indices i = 0 ; i < rank_out; ++i)
    {
        dimensions_out[i] = (nd_arr_in->dims)[order[i]]  ;
    }
    //
    ND_FUNCTION_CALL(init, TYPE_S) (nd_arr_out, rank_out, dimensions_out);
    // free temp array
    free(dimensions_out);

}



/* Function to initiate slice array directly*/
void FUNCTION(init_slice, TYPE_S) (const ND_indices * start_idx, const ND_indices * end_idx, const ND_indices * step_idx, \
                            const ARRAY_T(TYPE_S) * nd_arr_in, ARRAY_T(TYPE_S) * nd_arr_out )
{
    if (nd_arr_in->rank == NULL) error_msg("uninitialized input array for init_slice function");

    if (*(nd_arr_in->rank) == 0 ) error_msg("Rank of input array cannot be zero. Use nd_ele function");

    ND_indices slice_rank = (ND_indices) 0;
    ND_indices arr_in_rank = *(nd_arr_in->rank);

    ND_indices * sliced_dims_temp = malloc(3*arr_in_rank * sizeof(ND_indices)); // create both temp and sliced and strides
    ND_indices * sliced_dims = sliced_dims_temp + arr_in_rank ; 
    
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
    //
    ND_FUNCTION_CALL(init, TYPE_S) (nd_arr_out, slice_rank, sliced_dims);
    // free temp array
    free(sliced_dims_temp);

}


/* Function to initiate strip_dims array directly*/
void FUNCTION(init_strip_dims, TYPE_S) (const ARRAY_T(TYPE_S) * nd_arr_in, const ND_indices n_dims_strip, ARRAY_T(TYPE_S) * nd_arr_out)
{
    if (nd_arr_in->rank == NULL) error_msg("uninitialized input array for init_strip_dims function");



    ND_indices rank_out;
    

    if (n_dims_strip >= (nd_arr_in->rank)[0]) error_msg("Number of dims to stip must be less than the rank of original array");

    rank_out = (nd_arr_in->rank)[0]-n_dims_strip ;

    ND_indices * dimensions_out = malloc(sizeof(ND_indices) * rank_out);
    
    for (ND_indices i = 1; i <= rank_out; ++i)
    {
        dimensions_out[rank_out - i] = (nd_arr_in->dims)[(nd_arr_in->rank)[0] - i] ; 
    }

    //
    ND_FUNCTION_CALL(init, TYPE_S) (nd_arr_out, rank_out, dimensions_out);
    // free temp array
    free(dimensions_out);

}