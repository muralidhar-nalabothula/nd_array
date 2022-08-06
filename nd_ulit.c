/* This file contains functions common for all types*/

#include <stdio.h>
#include <stdlib.h>
#include "nd_ulit.h"

void error_msg(const char * error_str)
{
    fprintf(stdout, "# [ Error !!!] : %s \n",error_str);
    exit(EXIT_FAILURE);
}

