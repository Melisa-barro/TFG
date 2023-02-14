#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <limits.h>
#include <time.h>
#include "aco.h"

/* ------------------------------------------ CAMBIOS ------------------------------------------ */

int *ant_solution;
float *ant_score;
int *best_so_far_ant_solution;
float best_so_far_ant_score;
float *pheromone0;
float *pheromone1;

/* ---------------------------------------- FIN CAMBIOS ---------------------------------------- */

int n_ants;      /* number of ants */

float rho;           /* parameter for evaporation */
float q_0;           /* probability of best choice in tour construction */

int as_flag;     /* ant system */
int eas_flag;    /* elitist ant system */
int mmas_flag;   /* MAX-MIN ant system */

float   trail_max;       /* maximum pheromone trail in MMAS */
float   trail_min;       /* minimum pheromone trail in MMAS */
float   trail_0;         /* initial pheromone level */
int     u_gb;

int n;		/* problem size */


void allocate_ants ( void )
{

    if ((ant_solution = (int *) calloc((n_ants*n), sizeof(int))) == NULL) {
        printf("Out of memory, exit.");
        exit(1);
    }
    if ((ant_score = (float *) calloc(n_ants, sizeof(float))) == NULL) {
        printf("Out of memory, exit.");
        exit(1);
    }
    if ((best_so_far_ant_solution = (int *) calloc(n, sizeof(int))) == NULL) {
        printf("Out of memory, exit.");
        exit(1);
    }
    
}



int find_best( void )
{
    double min;
    int k, k_min;

    min = ant_score[0];
    k_min = 0;
    for( k = 1 ; k < n_ants ; k++ ) {
        if( ant_score[k] < min ) {
            min = ant_score[k];
            k_min = k;
        }
    }
    return k_min;
}



int find_worst( void )
    {
    double max;
    int k, k_max;

    max = ant_score[0];
    k_max = 0;
    for( k = 1 ; k < n_ants ; k++ ) {
        if( ant_score[k] > max ) {
            max = ant_score[k];
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
{
    int i, j;

    for ( i = 0 ; i < n ; i++ ) {
        if ( pheromone0[i] < trail_min ) {
            pheromone0[i] = trail_min;
        } else if ( pheromone0[i] > trail_max ) {
            pheromone0[i] = trail_max;
        }
    }
    for ( j = 0 ; j < n ; j++ ) {
        if ( pheromone1[j] < trail_min ) {
            pheromone1[j] = trail_min;
        } else if ( pheromone1[j] > trail_max ) {
            pheromone1[j] = trail_max;
        }
    }
}


void init_pheromone_trails( float initial_trail )
{
    int i;
    
    /* Initialize pheromone trails */
    for ( i = 0 ; i < n ; i++ ) {
        pheromone0[i] = initial_trail;
        pheromone1[i] = initial_trail;
    }
}


void evaporation( void )
{ 
    int i;

    for ( i = 0 ; i < n ; i++ ) {
        pheromone0[i]= (1 - rho) * pheromone0[i];
        pheromone1[i]= (1 - rho) * pheromone1[i];
    }
}


void global_update_pheromone( int *a_solution, float *a_score, int k )
{  
    int i;
    float d_tau;

    d_tau = 1.0 / a_score[k];
    for ( i = 0 ; i < n ; i++ ) {
        if(a_solution[k*n+i] == 0 )
            pheromone0[i] += d_tau;
        else
            pheromone1[i] += d_tau;
    }
}


void global_update_pheromone_best( int *a_solution, float a_score )
{  
    int i;
    float d_tau;

    d_tau = 1.0 / a_score;
    for ( i = 0 ; i < n ; i++ ) {
        if(a_solution[i] == 0 )
            pheromone0[i] += d_tau;
        else
            pheromone1[i] += d_tau;
    }
}


void global_update_pheromone_weighted( int *a_solution, float a_score, int weight )
{  
    int i;
    float d_tau;

    d_tau = (float) weight / a_score;
    for ( i = 0 ; i < n ; i++ ) {
        if(a_solution[i] == 0 )
            pheromone0[i] += d_tau;
        else
            pheromone1[i] += d_tau;
    }
}


/****************************************************************
 ****************************************************************
Procedures implementing solution construction and related things
 ****************************************************************
 ****************************************************************/

void select_gate( int *a_solution, int gate, int pos )
{ 
    float prob, prob0, prob1;
    prob = pheromone0[gate]+pheromone1[gate];
    prob0 = pheromone0[gate]/prob;
    prob1 = pheromone1[gate]/prob;


    if ( (q_0 > 0.0) && (ran01( &seedH ) < q_0)  ) {
    /* with a probability q_0 make the best possible choice
    according to pheromone trails and heuristic information */
    /* we first check whether q_0 > 0.0, to avoid the very common case
    of q_0 = 0.0 to have to compute a random number, which is
    expensive computationally */
        if (prob1 < prob0 ) {
            a_solution[pos] = 0;
        }
        else {
            a_solution[pos] = 1;
        }
        return;
    }
    else {
        if (ran01(&seedH) < prob0 ) {
            a_solution[pos] = 0;
        }
        else {
            a_solution[pos] = 1;
        }
    }
    
 }


/**************************************************************************
 **************************************************************************
Procedures specific to the ant's tour manipulation other than construction
***************************************************************************
 **************************************************************************/

void copy_from_to(int *a1_solution, float *a1_score, int *a2_solution, float *a2_score, int k) 
{
    int i;
    a2_score[0] = a1_score[k];
    for ( i = 0 ; i < n ; i++ )
        a2_solution[i] = a1_solution[k*n+i];
}
