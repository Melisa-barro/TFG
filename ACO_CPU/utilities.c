/**
 * Ant Colony Optimization for the CellNopt
 *
 * @file utilities.c
 * @author patricia.gonzalez@udc.es
 * @brief File contains in-out and other misc. procedures
 */

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>

#include "aco.h"

int max_tries;

int iteration;         /* iteration counter */
int best_iteration;
int n_restarts;
int restart_best;
int restart_iters;

int max_iters;         /* maximum number of iterations */
long int seed;

double   max_time;          /* maximal allowed run time of a try  */
double   time_used;         /* time used until some given event */
double   time_passed;       /* time passed until some moment*/
double 	 best_time;
double 	 restart_time;

double optimal;           /* optimal solution or bound to find */


/* ------------------------------------------------------------------------ */

FILE *report_iter, *report, *final_report, *results_report;

char name_buf[LINE_BUF_LEN];
int  opt;


/**************************   TIMER  ***************************************/
static struct rusage res;
static struct timeval tp;
static double virtual_time, real_time;

void start_timers(void)
/*
 FUNCTION:       virtual and real time of day are computed and stored to
 allow at later time the computation of the elapsed time
 (virtual or real)
 INPUT:          none
 OUTPUT:         none
 (SIDE)EFFECTS:  virtual and real time are computed
 */
{
    getrusage( RUSAGE_SELF, &res );
    virtual_time = (double) res.ru_utime.tv_sec +
    (double) res.ru_stime.tv_sec +
    (double) res.ru_utime.tv_usec / 1000000.0 +
    (double) res.ru_stime.tv_usec / 1000000.0;
    
    gettimeofday( &tp, NULL );
    real_time =    (double) tp.tv_sec +
    (double) tp.tv_usec / 1000000.0;
}



double elapsed_time(TIMER_TYPE type)
/*
 FUNCTION:       return the time used in seconds (virtual or real, depending on type)
 INPUT:          TIMER_TYPE (virtual or real time)
 OUTPUT:         seconds since last call to start_timers (virtual or real)
 (SIDE)EFFECTS:  none
 */
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

/**************************   STATISTICS  ***************************************/

double ran01( long *idum )
/*    
      FUNCTION:       generate a random number that is uniformly distributed in [0,1]
      INPUT:          pointer to variable with the current seed
      OUTPUT:         random number uniformly distributed in [0,1]
      (SIDE)EFFECTS:  random number seed is modified (important, this has to be done!)
      ORIGIN:         numerical recipes in C
*/
{
  long k;
  double ans;

  k =(*idum)/IQ;
  *idum = IA * (*idum - k * IQ) - IR * k;
  if (*idum < 0 ) *idum += IM;
  ans = AM * (*idum);
  return ans;
}



long int random_number( long *idum )
/*    
      FUNCTION:       generate an integer random number
      INPUT:          pointer to variable containing random number seed
      OUTPUT:         integer random number uniformly distributed in {0,2147483647}
      (SIDE)EFFECTS:  random number seed is modified (important, has to be done!)
      ORIGIN:         numerical recipes in C
*/
{
  long int k;

  k =(*idum)/IQ;
  *idum = IA * (*idum - k * IQ) - IR * k;
  if (*idum < 0 ) *idum += IM;
  return *idum;
}



int ** generate_int_matrix( int n, int m)
/*    
      FUNCTION:       malloc a matrix and return pointer to it
      INPUT:          size of matrix as n x m 
      OUTPUT:         pointer to matrix
      (SIDE)EFFECTS:  
*/
{
  int i;
  int **matrix;

  if((matrix = (int**) malloc(sizeof(int) * n * m +
		      sizeof(int *) * n	 )) == NULL){
    printf("Out of memory, exit.");
    exit(1);
  }
  for ( i = 0 ; i < n ; i++ ) {
    matrix[i] = (int *)(matrix + n) + i*m;
  }

  return matrix;
}



double ** generate_double_matrix( int n, int m)
/*    
      FUNCTION:       malloc a matrix and return pointer to it
      INPUT:          size of matrix as n x m 
      OUTPUT:         pointer to matrix
      (SIDE)EFFECTS:  
*/
{

  int i;
  double **matrix;

  if((matrix = (double**) malloc(sizeof(double) * n * m +
		      sizeof(double *) * n	 )) == NULL){
    printf("Out of memory, exit.");
    exit(1);
  }
  for ( i = 0 ; i < n ; i++ ) {
    matrix[i] = (double *)(matrix + n) + i*m;
  }
  return matrix;
}



/**************************   IN-OUT  ***************************************/

void read_parameters( void )
/*
 FUNCTION:       read input file,
 INPUT:          none
 OUTPUT:         none
 COMMENTS:
 */
{
    char texto[20];
    double numero;

    FILE *params;
    if ((params = fopen("parameters.txt", "r"))==NULL)
    	printf("Without parameters file => default parameters...\n");
    else {
    	while (fscanf(params, "%s %lf", texto, &numero) > 1)
     	{
		if ( !strcmp(texto,"max_tries") ) max_tries = (int)numero;
		else if( !strcmp(texto,"n_ants") ) n_ants = (int)numero;
		else if( !strcmp(texto,"rho") ) rho = numero;
		else if( !strcmp(texto,"q_0") ) q_0 = numero;
		else if( !strcmp(texto,"max_iters") ) max_iters = (int)numero;
		else if( !strcmp(texto,"restart_iters") ) restart_iters = (int)numero;
		else if( !strcmp(texto,"max_time") ) max_time = numero;
		else if( !strcmp(texto,"as_flag") ) as_flag = (int)numero;
		else if( !strcmp(texto,"eas_flag") ) eas_flag = (int)numero;
		else if( !strcmp(texto,"mmas_flag") ) mmas_flag = (int)numero;
		else if( !strcmp(texto,"u_gb") ) u_gb = (int)numero;
		else if( !strcmp(texto,"optimal") ) optimal = numero;
		else printf(">>>>>>>>> Unknown parameter: %s\n",texto);
     	}
    
    fclose(params);
    }

}

void init_report( void )
/*
 FUNCTION:       prepare report file,
 INPUT:          none
 OUTPUT:         none
 COMMENTS:
 */
{
    char temp_buffer[LINE_BUF_LEN];

    sprintf(temp_buffer,"conv_report");
    report = fopen(temp_buffer, "w");
    sprintf(temp_buffer,"conv_report_iter");
    report_iter = fopen(temp_buffer, "w");
    sprintf(temp_buffer,"results_report");
    results_report = fopen(temp_buffer, "w");
    sprintf(temp_buffer,"final_report");
    final_report = fopen(temp_buffer, "w");

}


void set_default_parameters(void)
/*
 FUNCTION: set default parameter settings
 INPUT:    none
 OUTPUT:   none
 COMMENTS: none
 */
{
 
    max_tries	   = 30;
    n_ants         = 32;    /* number of ants */
    rho            = 0.5;
    q_0            = 0.0;
    max_iters      = 10000;
    seed           = (long int) time(NULL)*(ntry+1);
    max_time       = 200.0;
    optimal        = 0.0;
    as_flag        = 0;
    eas_flag       = 0;
    mmas_flag      = 1;
    u_gb	   = 20;
    restart_iters  = 100;
}


void print_parameters()
/*
 FUNCTION: print parameter settings
 INPUT:    none
 OUTPUT:   none
 COMMENTS: none
 */
{
    printf("\n Parameter settings are:\n");
    printf("max_tries\t\t %d\n", max_tries);
    printf("max_iters\t\t %d\n", max_iters);
    printf("max_time\t\t %.2f\n", max_time);
    printf("seed\t\t\t %ld\n", seed);
    printf("optimum\t\t\t %f\n", optimal);
    printf("n_ants\t\t\t %d\n", n_ants);
    printf("rho\t\t\t %.2f\n", rho);
    printf("q_0\t\t\t %.2f\n", q_0);
    printf("as_flag\t\t\t %d\n", as_flag);
    printf("eas_flag\t\t %d\n", eas_flag);
    printf("mmas_flag\t\t %d\n", mmas_flag);
    printf("restart_iters\t\t %d\n", restart_iters);
    printf("u_gb\t\t\t %d\n", u_gb);
}


void printSolution( int *t )
/*
 FUNCTION:       print the solution *t
 INPUT:          pointer to a solution
 OUTPUT:         none
 */
{
    int   i;
    
    printf("[ ");
    for( i = 0 ; i < n ; i++ ) {
        printf("%d ", t[i]);
    }
    printf(" ]\n");
}


void fprintSolution( int *t )
/*
 FUNCTION:       print the solution *t
 INPUT:          pointer to a solution
 OUTPUT:         none
 */
{
    int   i;

    if(results_report) {    
    	fprintf(results_report,"Try: %d, sol=[ ",ntry);
    	for( i = 0 ; i < n ; i++ ) {
        	fprintf(results_report,"%d ", t[i]);
    	}
    	fprintf(results_report," ]\n");
    }
}


void write_report( void )
/*
 FUNCTION:       write report file
 INPUT:          pointer to a solution
 OUTPUT:         none
 */
{

    if (final_report){
        fprintf(final_report,
                "Try %d:\t iters %d\t best_iter %d\t time %f\t best_time %f \t best_score %f\t restarts %d \n ",
                ntry,iteration,best_iteration,elapsed_time(REAL),best_time,best_so_far_ant->score,n_restarts);
    }
    fprintSolution(best_so_far_ant->solution);
    fflush(final_report);
    fflush(report_iter);
    fflush(report);
    fflush(results_report);

}

void read_benchmark( char *bench_file_name )
/*
 FUNCTION:       read the solution for the toy model
 INPUT:          file name
 OUTPUT:         none
 */
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
    	bs_optimum   = (int*) calloc(n, sizeof(int));
    
    	while (fscanf(sol_opt, "%d ", &num) == 1)
     		{
		bs_optimum[i]= num;
		i++;
		}
    }

}

    
