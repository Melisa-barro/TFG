/**
 * Ant Colony Optimization for the CellNopt
 *
 * @file ants.c
 * @author patricia.gonzalez@udc.es
 * @brief File contains procedures related with ants
 */
 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <limits.h>
#include <time.h>
#include "aco.h"

ant_struct *ant;
ant_struct *best_so_far_ant;

double   **pheromone;

int n_ants;      /* number of ants */

double rho;           /* parameter for evaporation */
double q_0;           /* probability of best choice in tour construction */

int as_flag;     /* ant system */
int eas_flag;    /* elitist ant system */
int mmas_flag;   /* MAX-MIN ant system */

double   trail_max;       /* maximum pheromone trail in MMAS */
double   trail_min;       /* minimum pheromone trail in MMAS */
double   trail_0;         /* initial pheromone level */
int     u_gb;

int n;		/* problem size */


void allocate_ants ( void )
/*    
      FUNCTION:       allocate the memory for the ant colony, the best-so-far and 
                      the iteration best ant
      INPUT:          none
      OUTPUT:         none
      (SIDE)EFFECTS:  allocation of memory for the ant colony and two ants that 
                      store intermediate tours

*/
{
    int     i;
  
    /* ANTS */
    if((ant = (ant_struct*) malloc(sizeof( ant_struct ) * n_ants +
		     sizeof(ant_struct *) * n_ants	 )) == NULL){
        printf("Out of memory, exit.");
        exit(1);
    }
    for ( i = 0 ; i < n_ants ; i++ ) {
        ant[i].solution        = (int*) calloc(n, sizeof(int));
    }

    /* BEST ANT */
    if((best_so_far_ant = (ant_struct*) malloc(sizeof( ant_struct ) )) == NULL){
        printf("Out of memory, exit.");
        exit(1);
    }
    best_so_far_ant->solution        = (int*) calloc(n, sizeof(int));
    

}



int find_best( void )
/*    
      FUNCTION:       find the best ant of the current iteration
      INPUT:          none
      OUTPUT:         index of struct containing the iteration best ant
      (SIDE)EFFECTS:  none
*/
{
    double   min;
    int   k, k_min;

    min = ant[0].score;
    k_min = 0;
    for( k = 1 ; k < n_ants ; k++ ) {
        if( ant[k].score < min ) {
            min = ant[k].score;
            k_min = k;
        }
    }
    return k_min;
}



int find_worst( void )
/*    
      FUNCTION:       find the worst ant of the current iteration
      INPUT:          none
      OUTPUT:         pointer to struct containing iteration best ant
      (SIDE)EFFECTS:  none
*/
{
    double   max;
    int   k, k_max;

    max = ant[0].score;
    k_max = 0;
    for( k = 1 ; k < n_ants ; k++ ) {
        if( ant[k].score > max ) {
            max = ant[k].score;
            k_max = k;
        }
    }
    return k_max;
}



/************************************************************
 ************************************************************
Procedures for pheromone manipulation 
 ************************************************************
 ************************************************************/


void check_pheromone_trail_limits( void )
/*
 FUNCTION:      MMAS keeps pheromone trails inside trail limits
 INPUT:         none
 OUTPUT:        none
 (SIDE)EFFECTS: pheromones are forced to interval [trail_min,trail_max]
 */
{
    int    i, j;
    
    for ( i = 0 ; i < n ; i++ ) {
        for ( j = 0 ; j < 2 ; j++ ) {
            if ( pheromone[i][j] < trail_min ) {
                pheromone[i][j] = trail_min;
            } else if ( pheromone[i][j] > trail_max ) {
                pheromone[i][j] = trail_max;            }
        }
    }
}


void init_pheromone_trails( double initial_trail )
/*
 FUNCTION:      initialize pheromone trails
 INPUT:         initial value of pheromone trails "initial_trail"
 OUTPUT:        none
 (SIDE)EFFECTS: pheromone matrix is reinitialized
 */
{
    int i, j;
    
    /* Initialize pheromone trails */
    for ( i = 0 ; i < n ; i++ ) {
        for ( j =0 ; j < 2 ; j++ ) {
            pheromone[i][j] = initial_trail;
        }
    }
}


void evaporation( void )
/*    
      FUNCTION:      implements the pheromone trail evaporation
      INPUT:         none
      OUTPUT:        none
      (SIDE)EFFECTS: pheromones are reduced by factor rho
*/
{ 
    int    i, j;

    for ( i = 0 ; i < n ; i++ ) {
        for ( j = 0 ; j < 2 ; j++ ) {
            pheromone[i][j] = (1 - rho) * pheromone[i][j];
        }
    }
}


void global_update_pheromone( ant_struct *a )
/*    
      FUNCTION:      reinforces edges used in ant k's solution
      INPUT:         pointer to ant that updates the pheromone trail 
      OUTPUT:        none
      (SIDE)EFFECTS: pheromones of arcs in ant k's tour are increased
*/
{  
    int i, j;
    double   d_tau;

    d_tau = 1.0 / a->score;
    for ( i = 0 ; i < n ; i++ ) {
        j = a->solution[i];
        pheromone[i][j] += d_tau;
    }
}



void global_update_pheromone_weighted( ant_struct *a, int weight )
/*    
      FUNCTION:      reinforces edges of the ant's tour with weight "weight"
      INPUT:         pointer to ant that updates pheromones and its weight  
      OUTPUT:        none
      (SIDE)EFFECTS: pheromones of arcs in the ant's tour are increased
*/
{  
    int      i, j;
    double        d_tau;

    d_tau = (double) weight / a->score;
    for ( i = 0 ; i < n ; i++ ) {
        j = a->solution[i];
        pheromone[i][j] += d_tau;
    }       
}


/****************************************************************
 ****************************************************************
Procedures implementing solution construction and related things
 ****************************************************************
 ****************************************************************/


void select_gate( ant_struct *a, int gate )
/*    
      FUNCTION:      chooses for an ant the next gate as the one with
                     maximal value of heuristic information times pheromone 
      INPUT:         pointer to ant and the construction step
      OUTPUT:        none 
      (SIDE)EFFECT:  ant moves to the next gate
*/
{ 
    double prob, prob0, prob1;
    prob = pheromone[gate][0]+pheromone[gate][1];
    prob0 = pheromone[gate][0]/prob;
    prob1 = pheromone[gate][1]/prob;

    
    if ( (q_0 > 0.0) && (ran01( &seed ) < q_0)  ) {
        /* with a probability q_0 make the best possible choice
         according to pheromone trails and heuristic information */
        /* we first check whether q_0 > 0.0, to avoid the very common case
         of q_0 = 0.0 to have to compute a random number, which is
         expensive computationally */
        if (prob1 < prob0 ) {
            a->solution[gate] = 0;
        }
        else {
            a->solution[gate] = 1;
        }
        return;
    }
    else {
        if (ran01(&seed) < prob0 ) {
            a->solution[gate] = 0;
        }
        else {
            a->solution[gate] = 1;
        }
    }
    
 }


/**************************************************************************
 **************************************************************************
Procedures specific to the ant's tour manipulation other than construction
***************************************************************************
 **************************************************************************/



void copy_from_to(ant_struct *a1, ant_struct *a2) 
{
/*    
      FUNCTION:       copy solution from ant a1 into ant a2
      INPUT:          pointers to the two ants a1 and a2 
      OUTPUT:         none
      (SIDE)EFFECTS:  a2 is copy of a1
*/
    int   i;
  
    a2->score = a1->score;
    for ( i = 0 ; i < n ; i++ ) {
        a2->solution[i] = a1->solution[i];
    }
}




