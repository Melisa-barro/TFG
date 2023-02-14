/*******  Includes  *******/

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <cuda_runtime.h>


/********* Defines *********/

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define MASK 123459876
#define LINE_BUF_LEN     255
#define INFTY            LONG_MAX

typedef enum type_timer {REAL, VIRTUAL} TIMER_TYPE;
double elapsed_time(TIMER_TYPE type);


/********* Defines GPU *********/

#define  rho_dev 0.5
#define restart_iters_dev 100


/****** Global CPU variables ******/

int iteration;
int ntry;
int n_ants;
int n;
int *bs_optimum;
int *bs_optimum_device;
long int *seed, seedH;
long int *seed_device;
int *iteration_device;
int *ant_solution ;
float best_so_far_ant_score;
float *best_so_far_ant_score_device;
int *best_so_far_ant_solution_device;
FILE *final_report;
double time_passed;


/****** Global GPU variables ******/

__device__ float best_so_far_ant_score_dev;


/*******************    REPORTS    *******************/

void init_report( void )
{
    char temp_buffer[LINE_BUF_LEN];

    sprintf(temp_buffer,"final_report");
    final_report = fopen(temp_buffer, "w");

}

void write_report( void )
{

    if (final_report){
        fprintf(final_report,
                "Try %d:\t iters %d\t time %f\t best_score %f\t \n ",
                ntry,iteration,elapsed_time(REAL),best_so_far_ant_score);
    }
    fflush(final_report);

}

/************    BENCHMARKS  ***********/

void read_benchmark( char *bench_file_name )
{
    int num;
    FILE *sol_opt;

    if ((sol_opt = fopen(bench_file_name, "r"))==NULL){
        printf("No instance benchmark file specified, abort\n");
        exit(1);
    }
    else {
        int i = 0;
        if (fscanf(sol_opt, "%d ", &num) == 1)
        {
            n = num;
            i++;
        }

        if((bs_optimum = (int*) malloc(sizeof( int ) )) == NULL){
            printf("Out of memory, benchmark, exit.");
            exit(1);
        }
        bs_optimum = (int*) calloc(n, sizeof(int));

        while (fscanf(sol_opt, "%d ", &num) == 1)
        {
            bs_optimum[i]= num;
            i++;
        }
    }

}

/**************************   TIMERS  ***************************************/
static struct rusage res;
static struct timeval tp;
static double virtual_time, real_time;

void start_timers(void)
{
    getrusage( RUSAGE_SELF, &res );
    virtual_time = (double) res.ru_utime.tv_sec +
    (double) res.ru_stime.tv_sec +
    (double) res.ru_utime.tv_usec / 1000000.0 +
    (double) res.ru_stime.tv_usec / 1000000.0;

    gettimeofday( &tp, NULL );
    real_time = (double) tp.tv_sec +
    (double) tp.tv_usec / 1000000.0;
}



double elapsed_time(TIMER_TYPE type)
{
    if (type == REAL) {
        gettimeofday( &tp, NULL );
        return( (double) tp.tv_sec +
                (double) tp.tv_usec / 1000000.0
                - real_time );
    }
    else {
        getrusage( RUSAGE_SELF, &res );
        return( (double) res.ru_utime.tv_sec +
                (double) res.ru_stime.tv_sec +
                (double) res.ru_utime.tv_usec / 1000000.0 +
                (double) res.ru_stime.tv_usec / 1000000.0
                - virtual_time );
    }
    
}


/************    DEVICE ACO FUNCTIONS  *************/

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


__device__ float obj_function_dev ( int k, int *bs_optimum_device, int n_dev, int *ant_solution_dev )
{
    int i, sc = 0;

    for ( i = 0 ; i < n_dev ; i++ ) {
        if (bs_optimum_device[i] != ant_solution_dev[k*n_dev + i]) sc++;
    }

    return ((float) sc);
}


__device__ void select_gate_dev( int gate, int pos, int my_ant, long int *seed_device, float * pheromone0_dev, float *pheromone1_dev, int *ant_solution_dev )
{ 
    float prob, prob0;
    prob = pheromone0_dev[gate]+pheromone1_dev[gate];
    prob0 = pheromone0_dev[gate]/prob;

    if (ran01_dev(&seed_device[my_ant]) < prob0 ) {
        ant_solution_dev[pos] = 0;
    }
    else {
        ant_solution_dev[pos] = 1;
    }

 }
 
__device__ void construct_solutions_device( int k, int *bs_optimum_device, long int *seed_device, int n_dev, float * pheromone0_dev, float *pheromone1_dev, float *ant_score_dev,
 int *ant_solution_dev )
{
    int j, pos;        /* counter variable */

    for ( j = 0 ; j < n_dev ; j++ ) {
        pos = k*n_dev + j;
        select_gate_dev(  j, pos, k, seed_device,  pheromone0_dev, pheromone1_dev, ant_solution_dev);
    }

    /* compute scores */
    ant_score_dev[k] = obj_function_dev( k, bs_optimum_device, n_dev, ant_solution_dev);

}

__device__ int find_best_device( int n_ants_dev, float *ant_score_dev)
{
    float min;
    int k_min, k;

    min = ant_score_dev[0];
    k_min = 0;
    for( k = 1 ; k < n_ants_dev ; k++ ) {
        if( ant_score_dev[k] < min ) {
            min = ant_score_dev[k];
            k_min = k;
        }
    }
    return k_min;
}

__device__ void copy_from_to_device(int *a1_solution_dev, float *a1_score_dev, int *a2_solution_dev, float *a2_score_dev, int k, int n_dev) 
{
    int i;
    a2_score_dev[0] = a1_score_dev[k];
    for ( i = 0 ; i < n_dev ; i++ )
        a2_solution_dev[i] = a1_solution_dev[k*n_dev+i];
}

__device__ void init_pheromone_trails_device( float initial_trail, int n_dev, float * pheromone0_dev, float * pheromone1_dev)
{
    int i;

    for ( i = 0 ; i < n_dev ; i++ ) {
        pheromone0_dev[i] = initial_trail;
        pheromone1_dev[i] = initial_trail;
    }
}

__device__ void update_statistics_device( int *best_so_far_ant_solution_device, int n_dev, int n_ants_dev, int iteration_dev, float * pheromone0_dev, float * pheromone1_dev, float *ant_score_dev, int *ant_solution_dev, float *trail_max_dev, float *trail_min_dev, float *trail_0_dev, int *best_iteration_dev, int *restart_best_dev, int *n_restarts_dev )
{

    int iteration_best_ant_dev;

    iteration_best_ant_dev = find_best_device(n_ants_dev, ant_score_dev);

    if ( ant_score_dev[iteration_best_ant_dev] < best_so_far_ant_score_dev ) {

        copy_from_to_device( ant_solution_dev, ant_score_dev, best_so_far_ant_solution_device, &best_so_far_ant_score_dev, iteration_best_ant_dev, n_dev);

        *best_iteration_dev = iteration_dev;
        *restart_best_dev = iteration_dev;

        *trail_max_dev = 1. / ( (rho_dev) * best_so_far_ant_score_dev );
        *trail_min_dev = *trail_max_dev / ( 2. * n_dev );
        *trail_0_dev = *trail_max_dev;

    }

    if (iteration_dev - (*restart_best_dev) > restart_iters_dev)  {
        *n_restarts_dev++;
        init_pheromone_trails_device(*trail_0_dev, n_dev, pheromone0_dev,  pheromone1_dev);
        *restart_best_dev = iteration_dev;
    }

}

__device__ void evaporation_device( int n_dev, float * pheromone0_dev, float *pheromone1_dev )
{ 
    int i;

    for ( i = 0 ; i < n_dev ; i++ ) {
        pheromone0_dev[i]= (1 - rho_dev) * pheromone0_dev[i];
        pheromone1_dev[i]= (1 - rho_dev) * pheromone1_dev[i];
    }
}


__device__ void check_pheromone_trail_limits_device( int n_dev, float * pheromone0_dev, float *pheromone1_dev, float trail_max_dev, float trail_min_dev )
{
    int i, j;

    for ( i = 0 ; i < n_dev ; i++ ) {
        if ( pheromone0_dev[i] < trail_min_dev ) {
            pheromone0_dev[i] = trail_min_dev;
        } else if ( pheromone0_dev[i] > trail_max_dev ) {
            pheromone0_dev[i] = trail_max_dev;
        }
    }
    for ( j = 0 ; j < n_dev ; j++ ) {
    if ( pheromone1_dev[j] < trail_min_dev ) {
            pheromone1_dev[j] = trail_min_dev;
        } else if ( pheromone1_dev[j] > trail_max_dev ) {
            pheromone1_dev[j] = trail_max_dev;
        }
    }
}


__device__ void global_update_pheromone_device( int k, int n_dev, float * pheromone0_dev, float * pheromone1_dev, float *ant_score_dev, int *ant_solution_dev)
{  
    int i;
    float d_tau;

    d_tau = 1.0 / ant_score_dev[k];
    for ( i = 0 ; i < n_dev ; i++ ) {
        if( ant_solution_dev[k*n_dev+i] == 0 )
            pheromone0_dev[i] += d_tau;
        else
            pheromone1_dev[i] += d_tau;
    }
}

__device__ void global_update_pheromone_best_device(int *best_so_far_ant_solution_device, float best_so_far_ant_score_dev, int n_dev, float * pheromone0_dev, float * pheromone1_dev)
{  
    int i;
    float   d_tau;

    d_tau = 1.0 / best_so_far_ant_score_dev;
    for ( i = 0 ; i < n_dev ; i++ ) {
        if (best_so_far_ant_solution_device[i] == 0 )
            pheromone0_dev[i] += d_tau;
        else
            pheromone1_dev[i] += d_tau;
    }
}


__device__ void mmas_update_device( int *best_so_far_ant_solution_device, int iteration_dev, int n_ants_dev, int n_dev, int u_gb_dev, float * pheromone0_dev, float *pheromone1_dev, float *ant_score_dev, int *ant_solution_dev, int restart_best_dev )
{
    
    int iteration_best_ant_dev;

    if ( iteration_dev % u_gb_dev ) {
        iteration_best_ant_dev = find_best_device( n_ants_dev, ant_score_dev );
        global_update_pheromone_device( iteration_best_ant_dev, n_dev, pheromone0_dev, pheromone1_dev, ant_score_dev, ant_solution_dev);
    }
    else {
        global_update_pheromone_best_device(best_so_far_ant_solution_device, best_so_far_ant_score_dev, n_dev, pheromone0_dev, pheromone1_dev);
    }


    if ( ( iteration_dev - restart_best_dev ) < (int)(restart_iters_dev/10) )
        u_gb_dev = 10;
    else if ( (iteration_dev - restart_best_dev) < (int)(restart_iters_dev/2) )
        u_gb_dev = 5;
    else if ( (iteration_dev - restart_best_dev) < (int)(restart_iters_dev/1.3) )
        u_gb_dev = 3;
    else if ( (iteration_dev - restart_best_dev) < restart_iters_dev )
        u_gb_dev = 2;
    else
        u_gb_dev = 1;

    
}

__device__ void pheromone_trail_update_device(int *best_so_far_ant_solution_device, int iteration_dev, int n_ants_dev, int n_dev,  int u_gb_dev, float * pheromone0_dev, float *pheromone1_dev, float *ant_score_dev, int *ant_solution_dev, float trail_max_dev, float trail_min_dev, int restart_best_dev)  
{
  
    /* Simulate the pheromone evaporation of all pheromones */
    evaporation_device( n_dev, pheromone0_dev, pheromone1_dev  );
    mmas_update_device( best_so_far_ant_solution_device, iteration_dev, n_ants_dev, n_dev, u_gb_dev,  pheromone0_dev, pheromone1_dev, ant_score_dev, ant_solution_dev, restart_best_dev );
    check_pheromone_trail_limits_device( n_dev, pheromone0_dev, pheromone1_dev, trail_max_dev, trail_min_dev );

}


__device__ void init_aco_device( int n_dev, int n_ants_dev, float * pheromone0_dev, float * pheromone1_dev, float *trail_max_dev, float *trail_min_dev, float *trail_0_dev, int *best_iteration_dev,
				 int *restart_best_dev, int *n_restarts_dev ) {
	
    best_so_far_ant_score_dev = INFTY;
    *best_iteration_dev = 1;
    *restart_best_dev = 1;
    *n_restarts_dev = 0;

    *trail_max_dev 		= 1. / ( (rho_dev) * 0.5 );
    *trail_min_dev 		= *trail_max_dev / ( 2. * n_dev );
    *trail_0_dev 		= *trail_max_dev;

    init_pheromone_trails_device(*trail_0_dev, n_dev, pheromone0_dev, pheromone1_dev);
}


/****************      ACO DEVICE      ********************/

__global__ void aco_device( int *iteration_device, int *best_so_far_ant_solution_device, float *best_so_far_ant_score_device, int* bs_optimum_device, int n, int n_ants, long int *seed_device, int *ant_solution_dev ) {

    float trail_max_dev = 0;
    float trail_min_dev = 0;
    float trail_0_dev = 0;
    int best_iteration_dev = 0;
    int restart_best_dev = 0;
    int n_restarts_dev = 0;
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int iteration_dev =1;

    int	u_gb_dev = 20;
    int max_iters_dev = 10000;
    float optimal_dev = 0.0;

    extern __shared__ float pheromone0_dev[ ];
    float *pheromone1_dev = &pheromone0_dev[n];
    float *ant_score_dev = &pheromone1_dev[n];


    if ( i == 0 )
        init_aco_device(n, n_ants, pheromone0_dev, pheromone1_dev, &trail_max_dev, &trail_min_dev, &trail_0_dev, &best_iteration_dev, &restart_best_dev, &n_restarts_dev);

    __syncthreads();

    /* iterations */
    while ( !( (iteration_dev >= max_iters_dev) || (best_so_far_ant_score_dev <= optimal_dev) ) )
    {
        construct_solutions_device( i, bs_optimum_device, seed_device, n, pheromone0_dev, pheromone1_dev, ant_score_dev, ant_solution_dev );
        __syncthreads();

        if (i == 0) {
            update_statistics_device( best_so_far_ant_solution_device, n, n_ants, iteration_dev, pheromone0_dev, pheromone1_dev, ant_score_dev, ant_solution_dev, &trail_max_dev, &trail_min_dev, &trail_0_dev, &best_iteration_dev, &restart_best_dev, &n_restarts_dev );
            pheromone_trail_update_device( best_so_far_ant_solution_device, iteration_dev, n_ants, n, u_gb_dev, pheromone0_dev, pheromone1_dev, ant_score_dev, ant_solution_dev, trail_max_dev, trail_min_dev, restart_best_dev );
        }

        iteration_dev++;
        __syncthreads();

    }

    if ( i == 0) {
        iteration_device[0]=iteration_dev;
        best_so_far_ant_score_device[0]=best_so_far_ant_score_dev;
    }
}


/*************************    CPU ACO FUNCTIONS  *****************************/

void init_aco( void )
{

    int i;

    n_ants = 1024;
    seedH = (long int) time(NULL)*(ntry+1);

    best_so_far_ant_score = INFTY;
    iteration = 1;

    /* Allocate seed for each ant */
    if ((seed = (long int *) calloc((n_ants), sizeof(long int))) == NULL) {
        printf("Out of memory in seed vector, exit.");
        exit(1);
    }
    for (i = 0 ; i < n_ants ; i++) seed[i]=seedH*(i+2);

    /* Allocate memory of GPU arrays */
    cudaMalloc((void**)&seed_device, (n_ants) * sizeof(long int));
    cudaMalloc((void**)&bs_optimum_device, n * sizeof(int));
    cudaMalloc((void**)&best_so_far_ant_solution_device, n * sizeof(int));
    cudaMalloc((void**)&iteration_device, sizeof(int));
    cudaMalloc((void**)&best_so_far_ant_score_device, sizeof(float));
    cudaMalloc((void**)&ant_solution, n*n_ants*sizeof(int));


    cudaMemcpy(bs_optimum_device, bs_optimum, n * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(seed_device, seed, n_ants * sizeof(long int), cudaMemcpyHostToDevice);

    start_timers();
        
 }

void exit_aco( void )
{
    
    write_report( );

    free(seed);

    /* Free GPU memory */
    cudaFree(seed_device);
    cudaFree(bs_optimum_device);
    cudaFree(iteration_device);
    cudaFree(ant_solution);
    cudaFree(best_so_far_ant_score_device);
    cudaFree(best_so_far_ant_solution_device);
   
}


/*************************    ACO HOST    *********************************/

double aco_algorithm(){
   
    int 	sizeSM;

    /* Initialize and start time measure */
    init_aco();

    sizeSM = sizeof(float)*n + sizeof(float)*n_ants + sizeof(float)*n;


    cudaFuncSetCacheConfig(aco_device, cudaFuncCachePreferShared);
    aco_device<<<1, n_ants, sizeSM>>>( iteration_device, best_so_far_ant_solution_device, best_so_far_ant_score_device, bs_optimum_device, n, n_ants, seed_device, ant_solution);

    /* Wait for kernel to end */
    cudaError_t cudaerr = cudaDeviceSynchronize();

    /* End time measure */
    time_passed = elapsed_time( REAL );

    /* Copy score and number of iters back to host */
    cudaMemcpy(&best_so_far_ant_score, best_so_far_ant_score_device, sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(&iteration, iteration_device, sizeof(int), cudaMemcpyDeviceToHost);

    exit_aco();

    return(0.0);

}


/*************************    MAIN    *********************************/

int main(int argc, char **argv) {
    
    int max_tries = 30;

    init_report ( );
    read_benchmark (argv[1]);


    for ( ntry = 0 ; ntry < max_tries ; ntry++ ) {
        printf("try %d\n",ntry);
        aco_algorithm();
    }

    return (1);

}

