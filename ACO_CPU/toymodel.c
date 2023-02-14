/**
 * Ant Colony Optimization for the CellNopt
 *
 * @file toymodel.c
 * @author patricia.gonzalez@udc.es
 * @brief File contains a toy model objective function
 */
 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <limits.h>
#include <time.h>
#include "aco.h"

int *bs_optimum;  /* problem optimal solution (for toy model) */

double obj_function ( ant_struct *a )
/*    
      FUNCTION:       cost function that computes the distance to the known optimum
      INPUT:          optimum and ant-solution
      OUTPUT:         score

*/
{
    int     i, sc = 0;
  
    for ( i = 0 ; i < n ; i++ ) {
        if (bs_optimum[i] != a->solution[i]) sc++;
    }

   return ((double) sc);
}



