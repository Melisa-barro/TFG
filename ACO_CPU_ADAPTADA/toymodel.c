#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <limits.h>
#include <time.h>
#include "aco.h"

int *bs_optimum;  /* problem optimal solution (for toy model) */

float obj_function ( int *a_solution, int k )
{
    int     i, sc = 0;
    

    for ( i = 0 ; i < n ; i++ ) {
        if (bs_optimum[i] != a_solution[k*n + i]) sc++;
    }

   return ((float) sc);
}



