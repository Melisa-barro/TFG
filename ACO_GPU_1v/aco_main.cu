#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <cuda_runtime.h>

#include "aco.h"

int ntry;
int *ant_solution_dev, *bs_optimum_dev;
float *ant_score_dev;
float *pheromone0_dev, *pheromone1_dev;
long int *seed_dev;

__device__ double ran01_dev( long *idum )
{
    long k;
    double ans;

    k =(*idum)/IQ;
    *idum = IA * (*idum - k * IQ) - IR * k;
    if (*idum < 0 ) *idum += IM;
    ans = AM * (*idum);
    return ans;
}


__device__ float obj_function_dev ( int *a_solution_dev, int *bs_optimum_dev, int k, int n )
{
    int i, sc = 0;

    for ( i = 0 ; i < n ; i++ ) {
        if (bs_optimum_dev[i] != a_solution_dev[k*n + i]) sc++;
    }

    return ((float) sc);
}


__device__ void select_gate_dev( int *a_solution_dev, int gate, int pos, int my_ant, long int *seed_dev, float q_0, float *pheromone0_dev, float *pheromone1_dev )
{ 
    float prob, prob0, prob1;
    prob = pheromone0_dev[gate]+pheromone1_dev[gate];
    prob0 = pheromone0_dev[gate]/prob;
    prob1 = pheromone1_dev[gate]/prob;


    if ( (q_0 > 0.0) && (ran01_dev( &seed_dev[my_ant] ) < q_0)  ) {
    /* with a probability q_0 make the best possible choice
    according to pheromone trails and heuristic information */
    /* we first check whether q_0 > 0.0, to avoid the very common case
    of q_0 = 0.0 to have to compute a random number, which is
    expensive computationally */
        if (prob1 < prob0 ) {
            a_solution_dev[pos] = 0;
        }
        else {
            a_solution_dev[pos] = 1;
        }
        return;
    }
    else {
        if (ran01_dev(&seed_dev[my_ant]) < prob0 ) {
            a_solution_dev[pos] = 0;
        }
        else {
            a_solution_dev[pos] = 1;
        }
    }
    
 }
 
__global__ void construct_solutions_dev(int* ant_solution_dev, float* ant_score_dev, int* bs_optimum_dev, int n, long int *seed_dev, float q_0, float *pheromone0_dev, float *pheromone1_dev)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j, pos;        /* counter variable */

    for ( j = 0 ; j < n ; j++ ) {
        pos = i*n + j;
        select_gate_dev( ant_solution_dev, j, pos, i, seed_dev, q_0, pheromone0_dev, pheromone1_dev);
    }

    /* compute scores */
    ant_score_dev[i] = obj_function_dev(ant_solution_dev, bs_optimum_dev, i, n);

}

void construct_solutions_device()
{
    /* Copiar informacion necesaria a GPU */
    cudaMemcpy(pheromone0_dev, pheromone0, n * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(pheromone1_dev, pheromone1, n * sizeof(float), cudaMemcpyHostToDevice);

    construct_solutions_dev<<<1,n_ants>>>(ant_solution_dev, ant_score_dev, bs_optimum_dev, n, seed_dev, q_0, pheromone0_dev, pheromone1_dev);

    /* Esperamos a que termine la ejecucion en kernel */
    cudaError_t cudaerr = cudaDeviceSynchronize();

    /* Copiamos de vuelta los datos a CPU */
    cudaMemcpy(ant_solution, ant_solution_dev, (n_ants*n) * sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(ant_score, ant_score_dev, n_ants * sizeof(float), cudaMemcpyDeviceToHost);
}


/*************************    ACO procedures  *****************************/

int termination_condition( void )
{
  return ( ((iteration >= max_iters) || (elapsed_time( REAL ) >= max_time)) ||
	  (best_so_far_ant_score <= optimal));
}


void init_aco( void )
{

    int i;

    set_default_parameters ( );

    /* Allocate ants */
    allocate_ants();

    /* Allocate seed for each ant */

    if ((seed = (long int *) calloc((n_ants), sizeof(long int))) == NULL) {
        printf("Out of memory, exit.");
        exit(1);
    }
    for (i = 0 ; i < n_ants ; i++ ) seed[i]=seedH*(i+2);

    /* Reservamos memoria para los arrays que usaremos en GPU */
    cudaMalloc((void**)&seed_dev, (n_ants) * sizeof(long int));
    cudaMalloc((void**)&ant_solution_dev, (n_ants*n) * sizeof(int));
    cudaMalloc((void**)&bs_optimum_dev, n * sizeof(int));
    cudaMalloc((void**)&ant_score_dev, n_ants * sizeof(float));
    cudaMalloc((void**)&pheromone0_dev, n * sizeof(float));
    cudaMalloc((void**)&pheromone1_dev, n * sizeof(float));

    cudaMemcpy(bs_optimum_dev, bs_optimum, n * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(seed_dev, seed, n_ants * sizeof(long int), cudaMemcpyHostToDevice);

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

    /* allocate pheromone arrays */
    if ((pheromone0 = (float*) calloc(n, sizeof(float))) == NULL) {
        printf("Out of memory, exit.");
        exit(1);
    }
    if ((pheromone1 = (float*) calloc(n, sizeof(float))) == NULL) {
        printf("Out of memory, exit.");
        exit(1);
    }

    /* Initialize pheromone trails */
    trail_max = 1. / ( (rho) * 0.5 );
    trail_min = trail_max / ( 2. * n );
    trail_0 = trail_max;
    init_pheromone_trails( trail_0 );

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

    /* Liberar memoria de GPU */
    cudaFree(seed_dev);
    cudaFree(ant_solution_dev);
    cudaFree(ant_score_dev);
    cudaFree(bs_optimum_dev);
    cudaFree(pheromone0_dev);
    cudaFree(pheromone1_dev);
   
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

        trail_max = 1. / ( (rho) * best_so_far_ant_score );
        trail_min = trail_max / ( 2. * n );
        trail_0 = trail_max;

    }

    if ( iteration - restart_best > restart_iters ) {
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
    mmas_update();

    /* Check pheromone trail limits for MMAS */
    check_pheromone_trail_limits();

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

        construct_solutions_device();

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
    print_parameters ( );

    read_benchmark (argv[1]);

    for ( ntry = 0 ; ntry < max_tries ; ntry++ ) {
        printf("try %d\n",ntry);
        aco_algorithm();
    }

    return (1);

}

