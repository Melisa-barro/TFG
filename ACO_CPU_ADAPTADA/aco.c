#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>

#include "aco.h"

int ntry;


/*************************    ACO procedures  *****************************/

int termination_condition( void )
{
  return ( ((iteration >= max_iters) || (elapsed_time( REAL ) >= max_time)) ||
	  (best_so_far_ant_score <= optimal));
}



void construct_solutions()
{
    int k, j, pos;        /* counter variable */

    for ( k = 0 ; k < n_ants ; k++ ) {
       for ( j = 0 ; j < n ; j++ ) {
       		pos = k*n + j;
            select_gate( ant_solution, j, pos);
       }
        
        /* compute scores */
        ant_score[k] = obj_function(ant_solution, k);
    }
}


void init_ants()
{
    int k, j, rnd;   
    
    /* solution for ants initialized randomly */
    for ( k = 0 ; k < n_ants ; k++ ) {
        for ( j = 0 ; j < n ; j++ ) {
            rnd = (int) round( ran01( &seed ) ); /* random number 0 or 1 */
            ant_solution[k*n+j] = rnd;
        }
    	/* compute scores */
        ant_score[k] = obj_function(ant_solution, k);
    }
}


void init_aco( void )
{

    /* Allocate ants */
    allocate_ants();

    /* Initialize variables concerning statistics etc. */
    iteration    = 1;
    best_iteration = 1;
    restart_best = 1;
    n_restarts = 0;
    best_so_far_ant_score = INFTY;
    
    start_timers();
    best_time = 0.0;
    time_used = elapsed_time( REAL );
    time_passed = time_used;
   
    /* allocate pheromone matrix */
    if ((pheromone0 = (float*) calloc(n, sizeof(float))) == NULL) {
    	printf("Out of memory, exit.");
	exit(1);
    }
    if ((pheromone1 = (float*) calloc(n, sizeof(float))) == NULL) {
    	printf("Out of memory, exit.");
	exit(1);
    }

    /* Initialize pheromone trails */
    if (mmas_flag) {
        trail_max = 1. / ( (rho) * 0.5 );
        trail_min = trail_max / ( 2. * n );
        trail_0 = trail_max;
        init_pheromone_trails( trail_0 );
    }
    else {
        trail_0 = 0;
        init_pheromone_trails( trail_0 );
    }

    if (report) fprintf(report,"******** Try: %d **********\n",ntry);
    if (report_iter) fprintf(report_iter,"******** Try: %d **********\n",ntry);
   
 
}

void exit_aco( void )
{
    
    write_report ( );

    free( pheromone0 );
    free( pheromone1 );
    free( best_so_far_ant_solution );
    free( ant_solution );
    free( ant_score );
}
    
void update_statistics( void )
{

    int iteration_best_ant;

    iteration_best_ant = find_best(); /* iteration_best_ant is a global variable */
    if ( ant_score[iteration_best_ant] < best_so_far_ant_score ) {
        
        time_used = elapsed_time( REAL ); /* best sol found after time_used */
        copy_from_to( ant_solution, ant_score, best_so_far_ant_solution, &best_so_far_ant_score, iteration_best_ant );
        
        if ( report ) fprintf(report,"%f \t %f\n",best_so_far_ant_score,elapsed_time(REAL));
        if ( report_iter ) fprintf(report_iter,"%f \t %d\n",best_so_far_ant_score,iteration);
        
	    best_iteration = iteration;
        restart_best = iteration;
        best_time = time_used;

         if ( mmas_flag ) {
            trail_max = 1. / ( (rho) * best_so_far_ant_score );
            trail_min = trail_max / ( 2. * n );
            trail_0 = trail_max;
         }

    }
   
    if ( mmas_flag && (iteration - restart_best > restart_iters) ) {
        /* MAX-MIN Ant System was the first ACO algorithm to use
         pheromone trail re-initialisation as implemented
         here. Other ACO algorithms may also profit from this mechanism.
         */
	    n_restarts++;
        
	    init_pheromone_trails( trail_0 );
        restart_best = iteration;
        restart_time = elapsed_time( REAL );
    }
    
}


void pheromone_trail_update( void )  
{
  
    /* Simulate the pheromone evaporation of all pheromones */
    evaporation();

    /* Apply the pheromone deposit for different ACO variants */
    if (mmas_flag) mmas_update();
    else if (eas_flag) eas_update();
    else as_update();

    /* Check pheromone trail limits for MMAS */
    if ( mmas_flag )
        check_pheromone_trail_limits();

}

void as_update( void )
{
    int   k;
    
    for ( k = 0 ; k < n_ants ; k++ )
        global_update_pheromone( ant_solution, ant_score, k );
    
}



void eas_update( void )
{
    int   k;
    
    for ( k = 0 ; k < n_ants ; k++ )
        global_update_pheromone( ant_solution, ant_score, k );

    global_update_pheromone_weighted( best_so_far_ant_solution, best_so_far_ant_score, 2 );
    
}


void mmas_update( void )
{
    
    int iteration_best_ant;
   
    if ( iteration % u_gb ) {
        iteration_best_ant = find_best();
        global_update_pheromone( ant_solution, ant_score, iteration_best_ant );
    }
    else {
        global_update_pheromone_best( best_so_far_ant_solution, best_so_far_ant_score );
    }
    
    
    if ( ( iteration - restart_best ) < (int)(restart_iters/10) )
        u_gb = 10;
    else if ( (iteration - restart_best) < (int)(restart_iters/2) )
        u_gb = 5;
    else if ( (iteration - restart_best) < (int)(restart_iters/1.3) )
        u_gb = 3;
    else if ( (iteration - restart_best) < restart_iters )
        u_gb = 2;
    else
        u_gb = 1;

    
}


/*************************    ACO main    *********************************/

double aco_algorithm(){
    
    double  score;
   
    init_aco();
    
    /* iterations */
    while ( !termination_condition() ) {

        if ( iteration == 1 ) init_ants();
        else construct_solutions();

        update_statistics();

        pheromone_trail_update();

        iteration++;
    }
    
    score = best_so_far_ant_score;
    
    exit_aco();
    return(score);
}


/*************************    MAIN    *********************************/

int main(int argc, char **argv) {
    
   
    init_report ( );
    
    set_default_parameters ( );
    read_parameters ( );
    print_parameters ( );

    read_benchmark (argv[1]);
    
    for ( ntry = 0 ; ntry < max_tries ; ntry++ ) {
	    printf("try %d\n",ntry);
        aco_algorithm();

    }

    return (1);

}




