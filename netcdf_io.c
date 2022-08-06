#include "nd_array_src.h"


/* netCDF4 I/O function */
/****************************************************************************************************/


/* Function to read from netCDF */
void FUNCTION(read, TYPE_S) (const char* file_name, const char* var_name, ARRAY_T(TYPE_S) * nd_arr_in)
{   
    /* Note this function create the ARRAY_T(TYPE_S), so FUNCTION(free, TYPE_S) () must be called to free the memory
        DO NOT pass a pointer which already has data. this leads to memory leak
        Ex: FUNCTION(read, TYPE_S) ("ndb.BS_elph", "exc_elph", &temp_array);
    */
    //
    if (nd_arr_in->data != NULL) error_msg("Input array is uninitialized or has data."); // check if there is data or uninitialized

    if (nd_arr_in->rank != NULL) free(nd_arr_in->rank); // if not null, free the memeory

    //
    int ncid, var_id, retval, nc_rank; // variables 

    if ((retval = nc_open(file_name, NC_NOWRITE, &ncid))) ERR(retval); // open file

    if ((retval = nc_inq_varid(ncid, var_name, &var_id))) ERR(retval); // get the id of the req variable

    if ((retval = nc_inq_var(ncid, var_id, NULL, NULL, &nc_rank, NULL, NULL ))) ERR(retval); // get rank

    int * dim_ids             = malloc(((size_t) nc_rank)*sizeof(int));
    ND_indices * nd_dims      = malloc(((size_t) nc_rank)*sizeof(ND_indices));
    size_t * nd_dims_temp     = malloc(((size_t) nc_rank)*sizeof(size_t));

    if ((retval = nc_inq_var(ncid, var_id, NULL, NULL, NULL, dim_ids, NULL ))) ERR(retval); // get dims
    //
    for (ND_indices i = 0; i < ((ND_indices) nc_rank); ++i)
    {
        if ((retval = nc_inq_dimlen(ncid, dim_ids[i], nd_dims_temp + i))) ERR(retval);
        nd_dims[i] = (ND_indices) nd_dims_temp[i] ;
    }

    #if defined(COMPILE_ND_DOUBLE_COMPLEX) || defined(COMPILE_ND_SINGLE_COMPLEX)

        if ((nc_rank < 1) || nd_dims[(ND_indices) (nc_rank-1)] != 2 ) error_msg("Cannot convert a real to complex array.") ; 

        ND_FUNCTION_CALL(init, TYPE_S) (nd_arr_in, ((ND_indices) (nc_rank-1)), nd_dims); 
        
        ND_FUNCTION_CALL(malloc, TYPE_S) (nd_arr_in);// this must be free outside else memory leak

        NetCDF_FUN_TYPE * temp_array = malloc(ND_FUNCTION_CALL(size, TYPE_S)(nd_arr_in)*2 * sizeof(NetCDF_FUN_TYPE));

        if ((retval = NetCDF_FUNCTION_CALL(nc_get_var, NetCDF_FUN_TYPE) (ncid, var_id, temp_array))) ERR(retval); //get data in floats

        memcpy(nd_arr_in->data,temp_array, 2*sizeof(NetCDF_FUN_TYPE)* (ND_FUNCTION_CALL(size, TYPE_S)(nd_arr_in)) );

        free(temp_array);

    #else

        ND_FUNCTION_CALL(init, TYPE_S) (nd_arr_in, ((ND_indices) nc_rank), nd_dims); 
        
        ND_FUNCTION_CALL(malloc, TYPE_S) (nd_arr_in);  // this must be free outside else memory leak

        if ((retval = NetCDF_FUNCTION_CALL(nc_get_var, NetCDF_FUN_TYPE) (ncid, var_id, nd_arr_in->data))) ERR(retval); //get data in floats

    #endif

    if ((retval = nc_close(ncid))) ERR(retval); // close the file

    // free all temp arrays
    free(dim_ids);
    free(nd_dims);
    free(nd_dims_temp);
}



/* Function to write to netCDF */
void FUNCTION(write, TYPE_S) (const char* file_name, const char* var_name, const ARRAY_T(TYPE_S) * nd_arr_in, char ** dim_names)
{
    /* Dumps the ARRAY_T(TYPE_S) to file with filemane (file_name) and dimenstion name (dim_names)
        varible names
        Ex: FUNCTION(write, TYPE_S) ("nc.temp", "elph_ex", &temp_array, (char * [5]) {"nq", "modes", "Sf", "Si", "re_im"});
    */
    if (nd_arr_in->rank == NULL) error_msg("Cannot write a uninitialized array");
    if (nd_arr_in->data == NULL) error_msg("Cannot write a empty array");

    int retval; // error iD 
    int ncid, varid; // var iDs

    if ((retval = nc_create(file_name, NC_CLOBBER , &ncid))) ERR(retval); // create file

    #if defined(COMPILE_ND_DOUBLE_COMPLEX) || defined(COMPILE_ND_SINGLE_COMPLEX)

    int * dimids = malloc((*(nd_arr_in->rank) + 1) * sizeof(int));
    //
    for (ND_indices i =0 ; i< (*(nd_arr_in->rank)); ++i )
    {
        if ((retval = nc_def_dim(ncid, dim_names[i], (size_t) (nd_arr_in->dims)[i], dimids+i))) ERR(retval);   // get ids for the omega dimensions
    }

    if ((retval = nc_def_dim(ncid, "re_im", (size_t) 2, dimids+ (*(nd_arr_in->rank))  ))) ERR(retval);   // get ids for the omega dimensions
    //

    if ((retval = nc_def_var(ncid, var_name, NetCDF_IO_TYPE, (int) (1 + *(nd_arr_in->rank)), dimids, &varid))) ERR(retval); // Writing the data with name Raman_tensor
    
    if ((retval = nc_enddef(ncid))) ERR(retval);

    NetCDF_FUN_TYPE * temp_array = malloc( ND_FUNCTION_CALL(size, TYPE_S)(nd_arr_in) * 2 * sizeof(NetCDF_FUN_TYPE));
    
    memcpy(temp_array, nd_arr_in->data, 2*sizeof(NetCDF_FUN_TYPE)* (ND_FUNCTION_CALL(size, TYPE_S)(nd_arr_in)) );
    
    if ((retval =  NetCDF_FUNCTION_CALL(nc_put_var, NetCDF_FUN_TYPE) (ncid, varid, temp_array))) ERR(retval);
    
    free(temp_array);
    
    #else
    
    int * dimids = malloc((*(nd_arr_in->rank)) * sizeof(int));
    //
    for (ND_indices i =0 ; i< (*(nd_arr_in->rank)); ++i )
    {
        if ((retval = nc_def_dim(ncid, dim_names[i], (size_t) (nd_arr_in->dims)[i], dimids+i))) ERR(retval);   // get ids for the omega dimensions
    }
    //
    if ((retval = nc_def_var(ncid, var_name, NetCDF_IO_TYPE, (int) (*(nd_arr_in->rank)), dimids, &varid))) ERR(retval); // Writing the data with name Raman_tensor
    
    if ((retval = nc_enddef(ncid))) ERR(retval);
    
    if ((retval =  NetCDF_FUNCTION_CALL(nc_put_var, NetCDF_FUN_TYPE) (ncid, varid, nd_arr_in->data))) ERR(retval);
    
    #endif

    if ((retval = nc_close(ncid))) ERR(retval);
    free(dimids);
}
