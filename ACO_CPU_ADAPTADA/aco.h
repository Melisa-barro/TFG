/***************************** ANTS **************************************/
#define HEURISTIC(m,n)     (1.0 / ((double) 0.1))
/* add a small constant to avoid division by zero if a distance is 
zero */

#define EPSILON            0.00000000000000000000000000000001

#define MAX_ANTS       1024     /* max no. of ants */


extern double   **total;     /* combination of pheromone times heuristic information */

extern int n_ants;      /* number of ants */
extern int n; 		/* problem size */

/* ------------------------------------------ CAMBIOS ------------------------------------------ */

extern int *ant_solution; /* array de tamaño numero hormigas * numero de caminos [n_ants*n] para guardar 0 o 1 dependiendo de si coge ese camino o no */
extern float *ant_score; /* array de tamaño numero de hormigas (tamaño de colonia) [n_ants] para guardar la puntuacion de cada una */

extern int *best_so_far_ant_solution; /* array con la solucion, que tiene tantos elementos como numero de caminos */
extern float best_so_far_ant_score;

extern float *pheromone0; /* arrays de feromonas con tamaño numero de caminos [n] */
extern float *pheromone1;

/* ---------------------------------------- FIN CAMBIOS ---------------------------------------- */

extern float rho;           /* parameter for evaporation */
extern float q_0;           /* probability of best choice in tour construction */

extern int as_flag;     /* = 1, run ant system */
extern int eas_flag;    /* = 1, run elitist ant system */
extern int mmas_flag;   /* = 1, run MAX-MIN ant system */

extern float   trail_max;       /* maximum pheromone trail in MMAS */
extern float   trail_min;       /* minimum pheromone trail in MMAS */
extern float   trail_0;         /* initial pheromone trail level */
extern int u_gb;            /* every u_gb iterations update with best-so-far ant */

extern int *bs_optimum;  /* problem optimal solution (for toy model) */

/***************************** IN-OUT **************************************/

#define LINE_BUF_LEN     255

extern int ntry;
extern int max_tries;

extern int iteration;    /* iteration counter */
extern int best_iteration;
extern int restart_best; 
extern int n_restarts; 
extern int max_iters;    /* maximum number of iterations */
extern int restart_iters;

extern double   max_time;     /* maximal allowed run time  */
extern double   time_used;    /* time used until some given event */
extern double   time_passed;  /* time passed until some moment*/
extern double 	best_time;
extern double 	restart_time;

extern double optimal;      /* optimal solution value or bound to find */

extern FILE *report_iter, *report, *final_report, *results_report;

extern char name_buf[LINE_BUF_LEN];
extern int  opt;

/***************************** TIMER **************************************/

typedef enum type_timer {REAL, VIRTUAL} TIMER_TYPE;

/***************************** UTILITIES **************************************/

#define INFTY                 LONG_MAX

#define TRUE  1
#define FALSE 0

/* general macros */

#define MAX(x,y)        ((x)>=(y)?(x):(y))
#define MIN(x,y)        ((x)<=(y)?(x):(y))

#define DEBUG( x )

#define TRACE( x )

/* constants for a random number generator, for details see numerical recipes in C */

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define MASK 123459876

extern long int seed;

double ran01 ( long *idum );

long int random_number ( long *idum );

int ** generate_int_matrix( int n, int m);

double ** generate_double_matrix( int n, int m);


/***************************** ANTS **************************************/

int termination_condition ( void );

void construct_solutions ( void );

void update_statistics ( void );

void pheromone_trail_update ( void );

void as_update ( void );

void eas_update ( void );

void mmas_update ( void );

void check_pheromone_trail_limits( void );

void init_pheromone_trails ( float initial_trail );

void evaporation ( void );

void global_update_pheromone ( int *a_solution, float *a_score, int pos );

void global_update_pheromone_best( int *a_solution, float a_score );

void global_update_pheromone_weighted ( int *a_solution, float a_score, int weight );

void compute_total_information( void );

void select_gate( int *a_solution, int gate, int pos );

int find_best ( void );

int find_worst( void );

void copy_from_to(int *a1_solution, float *a1_score, int *a2_solution, float *a2_score, int pos);

void allocate_ants ( void );

double compute_score ( int *s );

void pheromone_trail_update ( void );


/***************************** IN-OUT **************************************/

void set_default_parameters();

void read_parameters();

void init_report();

void write_report();

void print_parameters ( void );

void printSolution ( int *t);

void fprintSolution ( int *t);

/***************************** TIMER **************************************/

void start_timers(void);

double elapsed_time(TIMER_TYPE type);


/***************************** TOYMODEL **************************************/

float obj_function ( int *a_solution, int k );

void read_benchmark ( char *c );
