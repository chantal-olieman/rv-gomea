/**
 * AMaLGaM-Full.c
 *
 * Copyright (c) 1998-2010 Peter A.N. Bosman
 *
 * The software in this file is the proprietary information of
 * Peter A.N. Bosman.
 *
 * IN NO EVENT WILL THE AUTHOR OF THIS SOFTWARE BE LIABLE TO YOU FOR ANY
 * DAMAGES, INCLUDING BUT NOT LIMITED TO LOST PROFITS, LOST SAVINGS, OR OTHER
 * INCIDENTIAL OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OR THE INABILITY
 * TO USE SUCH PROGRAM, EVEN IF THE AUTHOR HAS BEEN ADVISED OF THE POSSIBILITY
 * OF SUCH DAMAGES, OR FOR ANY CLAIM BY ANY OTHER PARTY. THE AUTHOR MAKES NO
 * REPRESENTATIONS OR WARRANTIES ABOUT THE SUITABILITY OF THE SOFTWARE, EITHER
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, OR NON-INFRINGEMENT. THE
 * AUTHOR SHALL NOT BE LIABLE FOR ANY DAMAGES SUFFERED BY ANYONE AS A RESULT OF
 * USING, MODIFYING OR DISTRIBUTING THIS SOFTWARE OR ITS DERIVATIVES.
 *
 * Adapted Maximum-Likelihood Gaussian Model
 * Iterated Density-Estimation Evolutionary Algorithm
 * with a full covariance matrix
 *
 * In this implementation, minimization is assumed.
 *
 * The software in this file is the result of (ongoing) scientific research.
 * The following people have been actively involved in this research over
 * the years:
 * - Peter A.N. Bosman
 * - Dirk Thierens
 * - Jorn Grahl
 *
 * This is the most up-to-date literature reference regarding this software:
 *
 * P.A.N. Bosman. On Empirical Memory Design, Faster Selection of Bayesian
 * Factorizations and Parameter-Free Gaussian EDAs. In G. Raidl, E. Alba,
 * J. Bacardit, C. Bates Congdon, H.-G. Beyer, M. Birattari, C. Blum,
 * P.A.N. Bosman, D. Corne, C. Cotta, M. Di Penta, B. Doerr, R. Drechsler,
 * M. Ebner, J. Grahl, T. Jansen, J. Knowles, T. Lenaerts, M. Middendorf,
 * J.F. Miller, M. O'Neill, R. Poli, G. Squillero, K. Stanley, T. St?tzle
 * and J. van Hemert, editors, Proceedings of the Genetic and Evolutionary
 * Computation Conference - GECCO-2009, pages 389-396, ACM Press, New York,
 * New York, 2009.
 */

/*-=-=-=-=-=-=-=-=-=-=-=-=-=-= Section Includes -=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <sys/time.h>
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/

typedef struct FOS {
    int length;
    int **sets;
    int *set_length;
} FOS;

/*-=-=-=-=-=-=-=-=-=-=-=-= Section Header Functions -=-=-=-=-=-=-=-=-=-=-=-=*/
void *Malloc( long size );
double **matrixNew( int n, int m );
double vectorDotProduct( double *vector0, double *vector1, int n0 );
double *matrixVectorMultiplication( double **matrix, double *vector, int n0, int n1 );
double **matrixMatrixMultiplication( double **matrix0, double **matrix1, int n0, int n1, int n2 );
int blasDSWAP( int n, double *dx, int incx, double *dy, int incy );
int blasDAXPY(int n, double da, double *dx, int incx, double *dy, int incy);
void blasDSCAL( int n, double sa, double x[], int incx );
int linpackDCHDC( double a[], int lda, int p, double work[], int ipvt[] );
double **choleskyDecomposition( double **matrix, int n );
int linpackDTRDI( double t[], int ldt, int n );
double **matrixLowerTriangularInverse( double **matrix, int n );
void matrixWriteToFile( FILE *file, double **matrix, int n0, int n1 );
int *mergeSort( double *array, int array_size );
void mergeSortWithinBounds( double *array, int *sorted, int *tosort, int p, int q );
void mergeSortMerge( double *array, int *sorted, int *tosort, int p, int r, int q );
int *mergeSortFitness( double *objectives, double *constraints, int number_of_solutions );
void mergeSortFitnessWithinBounds( double *objectives, double *constraints, int *sorted, int *tosort, int p, int q );
void mergeSortFitnessMerge( double *objectives, double *constraints, int *sorted, int *tosort, int p, int r, int q );
int *mergeSortInt( int *array, int array_size );
void mergeSortWithinBoundsInt( int *array, int *sorted, int *tosort, int p, int q );
void mergeSortMergeInt( int *array, int *sorted, int *tosort, int p, int r, int q );
void startTimer( void );
double getTimer( void );
void printTimer( void );
void interpretCommandLine( int argc, char **argv );
void parseCommandLine( int argc, char **argv );
void parseOptions( int argc, char **argv, int *index );
void parseFOSElementSize( int *index, int argc, char** argv );
void printAllInstalledProblems( void );
void optionError( char **argv, int index );
void parseParameters( int argc, char **argv, int *index );
void printUsage( void );
void checkOptions( void );
void printVerboseOverview( void );
double randomRealUniform01( void );
int randomInt( int maximum );
double random1DNormalUnit( void );
int *randomPermutation( int n );
char *installedProblemName( int index );
int numberOfInstalledProblems( void );
double installedProblemLowerRangeBound( int index, int dimension );
double installedProblemUpperRangeBound( int index, int dimension );
short isParameterInRangeBounds( double parameter, int dimension );
void installedProblemEvaluation( int index, double *parameters, double *objective_value, double *constraint_value );
void installedProblemEvaluationWithoutRotation( int index, double *parameters, double *objective_value, double *constraint_value );
void sphereFunctionProblemEvaluation( double *parameters, double *objective_value, double *constraint_value );
double sphereFunctionProblemLowerRangeBound( int dimension );
double sphereFunctionProblemUpperRangeBound( int dimension );
void ellipsoidFunctionProblemEvaluation( double *parameters, double *objective_value, double *constraint_value );
double ellipsoidFunctionLowerRangeBound( int dimension );
double ellipsoidFunctionUpperRangeBound( int dimension );
void cigarFunctionProblemEvaluation( double *parameters, double *objective_value, double *constraint_value );
double cigarFunctionLowerRangeBound( int dimension );
double cigarFunctionUpperRangeBound( int dimension );
void tabletFunctionProblemEvaluation( double *parameters, double *objective_value, double *constraint_value );
double tabletFunctionLowerRangeBound( int dimension );
double tabletFunctionUpperRangeBound( int dimension );
void cigarTabletFunctionProblemEvaluation( double *parameters, double *objective_value, double *constraint_value );
double cigarTabletFunctionLowerRangeBound( int dimension );
double cigarTabletFunctionUpperRangeBound( int dimension );
void twoAxesFunctionProblemEvaluation( double *parameters, double *objective_value, double *constraint_value );
double twoAxesFunctionLowerRangeBound( int dimension );
double twoAxesFunctionUpperRangeBound( int dimension );
void differentPowersFunctionProblemEvaluation( double *parameters, double *objective_value, double *constraint_value );
double differentPowersFunctionLowerRangeBound( int dimension );
double differentPowersFunctionUpperRangeBound( int dimension );
void rosenbrockFunctionProblemEvaluation( double *parameters, double *objective_value, double *constraint_value );
double rosenbrockFunctionLowerRangeBound( int dimension );
double rosenbrockFunctionUpperRangeBound( int dimension );
void parabolicRidgeFunctionProblemEvaluation( double *parameters, double *objective_value, double *constraint_value );
double parabolicRidgeFunctionLowerRangeBound( int dimension );
double parabolicRidgeFunctionUpperRangeBound( int dimension );
void sharpRidgeFunctionProblemEvaluation( double *parameters, double *objective_value, double *constraint_value );
double sharpRidgeFunctionLowerRangeBound( int dimension );
double sharpRidgeFunctionUpperRangeBound( int dimension );
void griewankFunctionProblemEvaluation( double *parameters, double *objective_value, double *constraint_value );
double griewankFunctionLowerRangeBound( int dimension );
double griewankFunctionUpperRangeBound( int dimension );
void michalewiczFunctionProblemEvaluation( double *parameters, double *objective_value, double *constraint_value );
double michalewiczFunctionLowerRangeBound( int dimension );
double michalewiczFunctionUpperRangeBound( int dimension );
void rastriginFunctionProblemEvaluation( double *parameters, double *objective_value, double *constraint_value );
double rastriginFunctionLowerRangeBound( int dimension );
double rastriginFunctionUpperRangeBound( int dimension );
void sumOfEllipsoidsFunctionProblemEvaluation( double *parameters, double *objective_value, double *constraint_value );
double sumOfEllipsoidsFunctionLowerRangeBound( int dimension );
double sumOfEllipsoidsFunctionUpperRangeBound( int dimension );
void initialize( void );
void initializeMemory( void );
void initializeNewPopulation( int population_index );
void initializeNewPopulationHarikLobo(int population_index);
void initializeNewPopulationMemory( int population_index );
void initializeRandomNumberGenerator( void );
void initializeParameterRangeBounds( void );
void initializeFOS( int population_index );
void learnFOS(int population_index);
void initializeCovarianceMatrices( int population_index );
void initializeDistributionMultipliers( int population_index );
void initializePopulationAndFitnessValues( int population_index );
void initializeObjectiveRotationMatrix( void );
int determineNearestNeighbour( int index, double **S_matrix, int *mpm_number_of_indices, int mpm_length );
void computeRanks( void );
void computeRanksForOnePopulation( int population_index );
double distanceInParameterSpace( double *solution_a, double *solution_b );
void writeGenerationalStatisticsForOnePopulation( int population_index );
void writeGenerationalSolutions( short final );
void writeGenerationalSolutionsBest( short final );
short checkTerminationCondition( void );
short checkImmediateTerminationConditions( void );
short checkNumberOfEvaluationsTerminationCondition( void );
short checkVTRTerminationCondition( void );
short checkTimeLimitTerminationCondition( void );
void checkAverageFitnessTerminationCondition( void );
void determineBestSolutionInCurrentPopulations( int *population_of_best, int *index_of_best );
void checkFitnessVarianceTermination( void );
short checkFitnessVarianceTerminationSinglePopulation( int population_index );
void checkDistributionMultiplierTerminationCondition( void );
void makeSelections( void );
void makeSelectionsForOnePopulation( int population_index );
void makeSelectionsForOnePopulationUsingDiversityOnRank0( int population_index );
void makePopulation( int population_index );
void estimateParametersAllPopulations( void );
void estimateParameters( int population_index );
void estimateParametersML( int population_index );
void estimateMeanVectorML( int population_index );
void estimateFullCovarianceMatrixML( int population_index );
void estimateCovarianceMatrixML( int population_index );
void copyBestSolutionsToPopulation( int population_index );
void applyDistributionMultipliers( int population_index );
void generateAndEvaluateNewSolutionsToFillPopulation( int population_index );
void computeParametersForSampling( int population_index );
double *generateNewSolution( int population_index );
void adaptDistributionMultipliersForOnePopulation( int population_index );
short generationalImprovementForOnePopulation( int population_index, double *st_dev_ratio );
short betterFitness( double objective_value_x, double constraint_value_x, double objective_value_y, double constraint_value_y );
double getStDevRatio( int population_index, double *parameters );
void ezilaitini( void );
void ezilaitiniMemory( void );
void ezilaitiniMemoryOnePopulation( int population_index );
void ezilaitiniObjectiveRotationMatrix( void );
void ezilaitiniCovarianceMatrices( int population_index );
void ezilaitiniFOS( int population_index );
void ezilaitiniParametersForSampling( int population_index );
void run( void );
int main( int argc, char **argv );
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/





/*-=-=-=-=-=-=-=-=-=-=-=- Section Global Variables -=-=-=-=-=-=-=-=-=-=-=-=-*/
short      write_generational_statistics,    /* Whether to compute and write statistics every generation (0 = no). */
           write_generational_solutions,     /* Whether to write the population every generation (0 = no). */
           print_verbose_overview,           /* Whether to print a overview of settings (0 = no). */
           use_vtr,                          /* Whether to terminate at the value-to-reach (VTR) (0 = no). */
           vtr_hit_status,
           use_guidelines,                   /* Whether to override parameters with guidelines (for those that exist). */
           haveNextNextGaussian = 0,         /* Internally used variable for sampling the normal distribution. */
          *populations_terminated,           /* Which populations have terminated (array). */
           use_univariate_FOS,
           evaluations_for_statistics_hit;
FOS **linkage_model;
int        number_of_parameters,             /* The number of parameters to be optimized. */
          *population_sizes,                 /* The size of each population. */
           base_population_size,
          *selection_sizes,                  /* The size of the selection for each population. */
           ellipsoid_size,
           number_of_ellipsoids,
           FOS_element_size = -1,
            pruned_tree=1,sparse_tree=1,wait_with_pruning=0,
           maximum_number_of_evaluations,    /* The maximum number of evaluations. */
           number_of_evaluations,            /* The current number of times a function evaluation was performed. */
          *number_of_generations,            /* The current generation count. */
           total_number_of_generations,
           number_of_populations,            /* The number of parallel populations that initially partition the search space. */
           maximum_number_of_populations,
           number_of_subgenerations_per_population_factor,
           problem_index,                    /* The index of the optimization problem. */
          *samples_drawn_from_normal,        /* The number of samples drawn from the i-th normal in the last generation. */
          *out_of_bounds_draws,              /* The number of draws that resulted in an out-of-bounds sample. */
          *no_improvement_stretch,           /* The number of subsequent generations without an improvement while the distribution multiplier is <= 1.0, for each population separately. */
           maximum_no_improvement_stretch;   /* The maximum number of subsequent generations without an improvement while the distribution multiplier is <= 1.0. */
double     tau,                              /* The selection truncation percentile (in [1/population_size,1]). */
           alpha_AMS,                        /* The percentile of offspring to apply AMS (anticipated mean shift) to. */
           delta_AMS,                        /* The adaptation length for AMS (anticipated mean shift). */
        ***populations,                      /* The populations containing the solutions. */
         **objective_values,                 /* Objective values for population members. */
         **constraint_values,                /* Sum of all constraint violations for population members. */
          *population_best_obj_val,
          *population_best_constraint_val,
         **ranks,                            /* Ranks of population members. */
        ***selections,                       /* Selected solutions, one for each population. */
         **objective_values_selections,      /* Objective values of selected solutions. */
         **constraint_values_selections,     /* Sum of all constraint violations of selected solutions. */
          *lower_range_bounds,               /* The respected lower bounds on parameters. */
          *upper_range_bounds,               /* The respected upper bounds on parameters. */
          *lower_init_ranges,                /* The initialization range lower bound. */
          *upper_init_ranges,                /* The initialization range upper bound */
           lower_user_range,                 /* The initial lower range-bound indicated by the user (same for all dimensions). */
           upper_user_range,                 /* The initial upper range-bound indicated by the user (same for all dimensions). */
           rotation_angle,                   /* The angle of rotation to be applied to the problem. */
          *distribution_multipliers,         /* Distribution multipliers (AVS mechanism), one for each population. */
           distribution_multiplier_increase, /* The multiplicative distribution multiplier increase. */
           distribution_multiplier_decrease, /* The multiplicative distribution multiplier decrease. */
           st_dev_ratio_threshold,           /* The maximum ratio of the distance of the average improvement to the mean compared to the distance of one standard deviation before triggering AVS (SDR mechanism). */
           vtr,                              /* The value-to-reach (function value of best solution that is feasible). */
           fitness_variance_tolerance,       /* The minimum fitness variance level that is allowed. */
         **mean_vectors,                     /* The mean vectors, one for each population. */
         **mean_vectors_previous,            /* The mean vectors of the previous generation, one for each population. */
       ****covariance_matrices,              /* The covariance matrices to be used for sampling, one for each population. */
        ***full_covariance_matrix,
       ****cholesky_factors_lower_triangle,  /* The unique lower triangular matrix of the Cholesky factorization. */
           nextNextGaussian,                 /* Internally used variable for sampling the normal distribution. */
         **rotation_matrix,                  /* The rotation matrix to be applied before evaluating. */
           maximum_number_of_seconds;
int64_t    random_seed,                      /* The seed used for the random-number generator. */
           random_seed_changing;             /* Internally used variable for randomly setting a random seed. */
clock_t     start, end;
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/

void evolveDifferentialDependencies( int population_index );
   
double *S_vector;
double **S_matrix;
double **dependency_matrix;
int **checked_matrix;
double **checked_pairs;
double *first_individual;
double *second_individual;
double *fitness_of_first_individual;
int **dependency_pairs;
int current_waiting_position = 0;
int sparse_learning;
int differential_grouping_evals = 0;
int number_of_waiting_cycles = 0;
int number_of_pairs = 0;
int total_dependencies_found = 0;
int minimal_dependencies_per_run = 2;
int number_of_checked_pairs = 0;
int iteration = 0;
double epsilon = 0.0;
int FOS_element_ub;
short learn_linkage_tree = 0;
short random_linkage_tree = 0;
int static_linkage_tree;
int continued_learning;
int dependency_learning;
int block_start, block_size;
int pruning_ub;
int evolve_learning;
int pairs_per_run;

/*-=-=-=-=-=-=-=-=-=-=-=-=-=-= Section Constants -=-=-=-=-=-=-=-=-=-=-=-=-=-*/
#define PI 3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706798
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/






/*-=-=-=-=-=-=-=-=-=-=-= Section Elementary Operations -=-=-=-=-=-=-=-=-=-=-*/
/**
 * Allocates memory and exits the program in case of a memory allocation failure.
 */
void *Malloc( long size )
{
  void *result;

  result = (void *) malloc( size );
  if( !result )
  {
    printf("\n");
    printf("Error while allocating memory in Malloc( %ld ), aborting program.", size);
    printf("\n");

    exit( 0 );
  }

  return( result );
}

double nround (double n, int c) {
    double marge = pow(10, c);
    double up = n * marge;
    double ret = round(up) / marge;

    return ret;
}

void getMinMaxofPopulation(int variable, int population_index, double *min, double *max){
    *min = populations[population_index][0][variable];
    *max = populations[population_index][0][variable];
    for(int i = 0; i < population_sizes[population_index]; i++){
        if(populations[population_index][i][variable] < *min)
            *min = populations[population_index][i][variable];
        else if(populations[population_index][i][variable] > *max)
            *max = populations[population_index][i][variable];
    }
}

/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/





/*-=-=-=-=-=-=-=-=-=-=-=-=-=-= Section Matrix -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
/**
 * Creates a new matrix with dimensions n x m.
 */
double **matrixNew( int n, int m )
{
  int      i;
  double **result;

  result = (double **) malloc( n*( sizeof( double * ) ) );
  for( i = 0; i < n; i++ )
    result[i] = (double *) malloc( m*( sizeof( double ) ) );

  return( result );
}

/**
 * Computes the dot product of two vectors of the same dimensionality n0.
 */
double vectorDotProduct( double *vector0, double *vector1, int n0 )
{
  int    i;
  double result;

  result = 0.0;
  for( i = 0; i < n0; i++ )
    result += vector0[i]*vector1[i];

  return( result );
}

/**
 * Computes the multiplication Av of a matrix A and a vector v
 * where matrix A has dimensions n0 x n1 and vector v has
 * dimensionality n1.
 */
double *matrixVectorMultiplication( double **matrix, double *vector, int n0, int n1 )
{
  int     i;
  double *result;

  result = (double *) malloc( n0*sizeof( double ) );
  for( i = 0; i < n0; i++ )
    result[i] = vectorDotProduct( matrix[i], vector, n1 );

  return( result );
}

/**
 * Computes the matrix multiplication of two matrices A and B
 * of dimensions A: n0 x n1 and B: n1 x n2.
 */
double **matrixMatrixMultiplication( double **matrix0, double **matrix1, int n0, int n1, int n2 )
{
  int     i, j, k;
  double **result;

  result = (double **) malloc( n0*sizeof( double * ) );
  for( i = 0; i < n0; i++ )
    result[i] = (double *) malloc( n2*sizeof( double ) );

  for( i = 0; i < n0; i++ )
  {
    for( j = 0; j < n2; j++ )
    {
      result[i][j] = 0;
      for( k = 0; k < n1; k++ )
        result[i][j] += matrix0[i][k]*matrix1[k][j];
    }
  }

  return( result );
}

/**
 * BLAS subroutine.
 */
int blasDSWAP( int n, double *dx, int incx, double *dy, int incy )
{
  double dtmp;

  if (n > 0)
  {
    incx *= sizeof( double );
    incy *= sizeof( double );

    dtmp  = (*dx);
    *dx   = (*dy);
    *dy   = dtmp;

    while( (--n) > 0 )
    {
      dx = (double *) ((char *) dx + incx);
      dy = (double *) ((char *) dy + incy);
      dtmp = (*dx); *dx = (*dy); *dy = dtmp;
    }
  }

  return( 0 );
}

/**
 * BLAS subroutine.
 */
int blasDAXPY(int n, double da, double *dx, int incx, double *dy, int incy)
{
  double dtmp0, dtmp, *dx0, *dy0;

  if( n > 0 && da != 0. )
  {
    incx *= sizeof(double);
    incy *= sizeof(double);
    *dy  += da * (*dx);

    if( (n & 1) == 0 )
    {
      dx   = (double *) ((char *) dx + incx);
      dy   = (double *) ((char *) dy + incy);
      *dy += da * (*dx);
      --n;
    }
    n = n >> 1;
    while( n > 0 )
    {
      dy0   = (double *) ((char *) dy + incy);
      dy    = (double *) ((char *) dy0 + incy);
      dtmp0 = (*dy0);
      dtmp  = (*dy);
      dx0   = (double *) ((char *) dx + incx);
      dx    = (double *) ((char *) dx0 + incx);
      *dy0  = dtmp0 + da * (*dx0);
      *dy   = dtmp + da * (*dx);
      --n;
    }
  }

  return( 0 );
}

/**
 * BLAS subroutine.
 */
void blasDSCAL( int n, double sa, double x[], int incx )
{
  int i, ix, m;

  if( n <= 0 )
  {
  }
  else if( incx == 1 )
  {
    m = n % 5;

    for( i = 0; i < m; i++ )
    {
      x[i] = sa * x[i];
    }

    for( i = m; i < n; i = i + 5 )
    {
      x[i]   = sa * x[i];
      x[i+1] = sa * x[i+1];
      x[i+2] = sa * x[i+2];
      x[i+3] = sa * x[i+3];
      x[i+4] = sa * x[i+4];
    }
  }
  else
  {
    if( 0 <= incx )
    {
      ix = 0;
    }
    else
    {
      ix = ( - n + 1 ) * incx;
    }

    for( i = 0; i < n; i++ )
    {
      x[ix] = sa * x[ix];
      ix = ix + incx;
    }
  }
}

/**
 * LINPACK subroutine.
 */
int linpackDCHDC( double a[], int lda, int p, double work[], int ipvt[] )
{
  int    info, j, jp, k, l, maxl, pl, pu;
  double maxdia, temp;

  pl   = 1;
  pu   = 0;
  info = p;
  for( k = 1; k <= p; k++ )
  {
    maxdia = a[k-1+(k-1)*lda];
    maxl   = k;
    if( pl <= k && k < pu )
    {
      for( l = k+1; l <= pu; l++ )
      {
        if( maxdia < a[l-1+(l-1)*lda] )
        {
          maxdia = a[l-1+(l-1)*lda];
          maxl   = l;
        }
      }
    }

    if( maxdia <= 0.0 )
    {
      info = k - 1;

      return( info );
    }

    if( k != maxl )
    {
      blasDSWAP( k-1, a+0+(k-1)*lda, 1, a+0+(maxl-1)*lda, 1 );

      a[maxl-1+(maxl-1)*lda] = a[k-1+(k-1)*lda];
      a[k-1+(k-1)*lda]       = maxdia;
      jp                     = ipvt[maxl-1];
      ipvt[maxl-1]           = ipvt[k-1];
      ipvt[k-1]              = jp;
    }
    work[k-1]        = sqrt( a[k-1+(k-1)*lda] );
    a[k-1+(k-1)*lda] = work[k-1];

    for( j = k+1; j <= p; j++ )
    {
      if( k != maxl )
      {
        if( j < maxl )
        {
          temp                = a[k-1+(j-1)*lda];
          a[k-1+(j-1)*lda]    = a[j-1+(maxl-1)*lda];
          a[j-1+(maxl-1)*lda] = temp;
        }
        else if ( maxl < j )
        {
          temp                = a[k-1+(j-1)*lda];
          a[k-1+(j-1)*lda]    = a[maxl-1+(j-1)*lda];
          a[maxl-1+(j-1)*lda] = temp;
        }
      }
      a[k-1+(j-1)*lda] = a[k-1+(j-1)*lda] / work[k-1];
      work[j-1]        = a[k-1+(j-1)*lda];
      temp             = -a[k-1+(j-1)*lda];

      blasDAXPY( j-k, temp, work+k, 1, a+k+(j-1)*lda, 1 );
    }
  }

  return( info );
}

/**
 * Computes the lower-triangle Cholesky Decomposition
 * of a square, symmetric and positive-definite matrix.
 * Subroutines from LINPACK and BLAS are used.
 */
double **choleskyDecomposition( double **matrix, int n )
{
  int     i, j, k, info, *ipvt;
  double *a, *work, **result;

  a    = (double *) Malloc( n*n*sizeof( double ) );
  work = (double *) Malloc( n*sizeof( double ) );
  ipvt = (int *) Malloc( n*sizeof( int ) );

  k = 0;
  for( i = 0; i < n; i++ )
  {
    for( j = 0; j < n; j++ )
    {
      a[k] = matrix[i][j];
      k++;
    }
    ipvt[i] = 0;
  }

  info = linpackDCHDC( a, n, n, work, ipvt );

  result = matrixNew( n, n );
  if( info != n ) /* Matrix is not positive definite */
  {
    k = 0;
    for( i = 0; i < n; i++ )
    {
      for( j = 0; j < n; j++ )
      {
        result[i][j] = i != j ? 0.0 : sqrt( matrix[i][j] );
        k++;
      }
    }
  }
  else
  {
    k = 0;
    for( i = 0; i < n; i++ )
    {
      for( j = 0; j < n; j++ )
      {
        result[i][j] = i < j ? 0.0 : a[k];
        k++;
      }
    }
  }

  free( ipvt );
  free( work );
  free( a );

  return( result );
}

/**
 * LINPACK subroutine.
 */
int linpackDTRDI( double t[], int ldt, int n )
{
  int    j, k, info;
  double temp;

  info = 0;
  for( k = n; 1 <= k; k-- )
  {
    if ( t[k-1+(k-1)*ldt] == 0.0 )
    {
      info = k;
      break;
    }

    t[k-1+(k-1)*ldt] = 1.0 / t[k-1+(k-1)*ldt];
    temp = -t[k-1+(k-1)*ldt];

    if ( k != n )
    {
      blasDSCAL( n-k, temp, t+k+(k-1)*ldt, 1 );
    }

    for( j = 1; j <= k-1; j++ )
    {
      temp = t[k-1+(j-1)*ldt];
      t[k-1+(j-1)*ldt] = 0.0;
      blasDAXPY( n-k+1, temp, t+k-1+(k-1)*ldt, 1, t+k-1+(j-1)*ldt, 1 );
    }
  }

  return( info );
}

/**
 * Computes the inverse of a matrix that is of
 * lower triangular form.
 */
double **matrixLowerTriangularInverse( double **matrix, int n )
{
  int     i, j, k, info;
  double *t, **result;

  t = (double *) Malloc( n*n*sizeof( double ) );

  k = 0;
  for( i = 0; i < n; i++ )
  {
    for( j = 0; j < n; j++ )
    {
      t[k] = matrix[j][i];
      k++;
    }
  }

  info = linpackDTRDI( t, n, n );

  result = matrixNew( n, n );
  k = 0;
  for( i = 0; i < n; i++ )
  {
    for( j = 0; j < n; j++ )
    {
      result[j][i] = i > j ? 0.0 : t[k];
      k++;
    }
  }

  free( t );

  return( result );
}

/**
 * Writes the contents of a matrix of dimensions n0 x n1 to a file.
 */
void matrixWriteToFile( FILE *file, double **matrix, int n0, int n1 )
{
  int  i, j;
  char line_for_output[10000];

  sprintf( line_for_output, "[" );
  fputs( line_for_output, file );
  for( i = 0; i < n0; i++ )
  {
    sprintf( line_for_output, "[" );
    fputs( line_for_output, file );
    for( j = 0; j < n1; j++ )
    {
      sprintf( line_for_output, "%lf", matrix[i][j] );
      fputs( line_for_output, file );
      if( j < n1-1 )
      {
        sprintf( line_for_output, ", " );
        fputs( line_for_output, file );
      }
    }
    if( i == n0-1 )
      sprintf( line_for_output, "]" );
    else
      sprintf( line_for_output, "];" );
    fputs( line_for_output, file );
  }
  sprintf( line_for_output, "]\n" );
  fputs( line_for_output, file );
}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/





/*-=-=-=-=-=-=-=-=-=-=-=-=-= Section Merge Sort -=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
/**
 * Sorts an array of doubles and returns the sort-order (small to large).
 */
int *mergeSort( double *array, int array_size )
{
  int i, *sorted, *tosort;

  sorted = (int *) Malloc( array_size * sizeof( int ) );
  tosort = (int *) Malloc( array_size * sizeof( int ) );
  for( i = 0; i < array_size; i++ )
    tosort[i] = i;

  if( array_size == 1 )
    sorted[0] = 0;
  else
    mergeSortWithinBounds( array, sorted, tosort, 0, array_size-1 );

  free( tosort );

  return( sorted );
}

/**
 * Subroutine of merge sort, sorts the part of the array between p and q.
 */
void mergeSortWithinBounds( double *array, int *sorted, int *tosort, int p, int q )
{
  int r;

  if( p < q )
  {
    r = (p + q) / 2;
    mergeSortWithinBounds( array, sorted, tosort, p, r );
    mergeSortWithinBounds( array, sorted, tosort, r+1, q );
    mergeSortMerge( array, sorted, tosort, p, r+1, q );
  }
}

/**
 * Subroutine of merge sort, merges the results of two sorted parts.
 */
void mergeSortMerge( double *array, int *sorted, int *tosort, int p, int r, int q )
{
  int i, j, k, first;

  i = p;
  j = r;
  for( k = p; k <= q; k++ )
  {
    first = 0;
    if( j <= q )
    {
      if( i < r )
      {
        if( array[tosort[i]] < array[tosort[j]] )
          first = 1;
      }
    }
    else
      first = 1;

    if( first )
    {
      sorted[k] = tosort[i];
      i++;
    }
    else
    {
      sorted[k] = tosort[j];
      j++;
    }
  }

  for( k = p; k <= q; k++ )
    tosort[k] = sorted[k];
}

/**
 * Sorts an array of objectives and constraints
 * using constraint domination and returns the
 * sort-order (small to large).
 */
int *mergeSortFitness( double *objectives, double *constraints, int number_of_solutions )
{
  int i, *sorted, *tosort;

  sorted = (int *) Malloc( number_of_solutions * sizeof( int ) );
  tosort = (int *) Malloc( number_of_solutions * sizeof( int ) );
  for( i = 0; i < number_of_solutions; i++ )
    tosort[i] = i;

  if( number_of_solutions == 1 )
    sorted[0] = 0;
  else
    mergeSortFitnessWithinBounds( objectives, constraints, sorted, tosort, 0, number_of_solutions-1 );

  free( tosort );

  return( sorted );
}

/**
 * Subroutine of merge sort, sorts the part of the objectives and
 * constraints arrays between p and q.
 */
void mergeSortFitnessWithinBounds( double *objectives, double *constraints, int *sorted, int *tosort, int p, int q )
{
  int r;

  if( p < q )
  {
    r = (p + q) / 2;
    mergeSortFitnessWithinBounds( objectives, constraints, sorted, tosort, p, r );
    mergeSortFitnessWithinBounds( objectives, constraints, sorted, tosort, r+1, q );
    mergeSortFitnessMerge( objectives, constraints, sorted, tosort, p, r+1, q );
  }
}

/**
 * Subroutine of merge sort, merges the results of two sorted parts.
 */
void mergeSortFitnessMerge( double *objectives, double *constraints, int *sorted, int *tosort, int p, int r, int q )
{
  int i, j, k, first;

  i = p;
  j = r;
  for( k = p; k <= q; k++ )
  {
    first = 0;
    if( j <= q )
    {
      if( i < r )
      {
        if( betterFitness( objectives[tosort[i]], constraints[tosort[i]],
                           objectives[tosort[j]], constraints[tosort[j]] ) )
          first = 1;
      }
    }
    else
      first = 1;

    if( first )
    {
      sorted[k] = tosort[i];
      i++;
    }
    else
    {
      sorted[k] = tosort[j];
      j++;
    }
  }

  for( k = p; k <= q; k++ )
    tosort[k] = sorted[k];
}

int *mergeSortInt( int *array, int array_size )
{
    int i, *sorted, *tosort;

    sorted = (int *) Malloc( array_size * sizeof( int ) );
    tosort = (int *) Malloc( array_size * sizeof( int ) );
    for( i = 0; i < array_size; i++ )
        tosort[i] = i;

    if( array_size == 1 )
        sorted[0] = 0;
    else
        mergeSortWithinBoundsInt( array, sorted, tosort, 0, array_size-1 );

    free( tosort );

    return( sorted );
}

/**
 * Subroutine of merge sort, sorts the part of the array between p and q.
 */
void mergeSortWithinBoundsInt( int *array, int *sorted, int *tosort, int p, int q )
{
    int r;

    if( p < q )
    {
        r = (p + q) / 2;
        mergeSortWithinBoundsInt( array, sorted, tosort, p, r );
        mergeSortWithinBoundsInt( array, sorted, tosort, r+1, q );
        mergeSortMergeInt( array, sorted, tosort, p, r+1, q );
    }
}

void mergeSortMergeInt( int *array, int *sorted, int *tosort, int p, int r, int q )
{
    int i, j, k, first;

    i = p;
    j = r;
    for( k = p; k <= q; k++ )
    {
        first = 0;
        if( j <= q )
        {
            if( i < r )
            {
                if( array[tosort[i]] < array[tosort[j]] )
                    first = 1;
            }
        }
        else
            first = 1;

        if( first )
        {
            sorted[k] = tosort[i];
            i++;
        }
        else
        {
            sorted[k] = tosort[j];
            j++;
        }
    }

    for( k = p; k <= q; k++ )
        tosort[k] = sorted[k];
}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/


/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=- Section Timer -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/

void startTimer( void )
{
    start = clock();
}

double getTimer( void )
{
    end = clock();
    return ( ((double) (end - start)) / CLOCKS_PER_SEC );
}

void printTimer( void )
{
    double cpu_time_used;

    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("%.3f\n",cpu_time_used);
}



/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/


/*-=-=-=-=-=-=-=-=-=-=- Section Interpret Command Line -=-=-=-=-=-=-=-=-=-=-*/
/**
 * Parses and checks the command line.
 */
void interpretCommandLine( int argc, char **argv )
{
  startTimer();

  parseCommandLine( argc, argv );

  if( use_guidelines )
  {
    tau                              = 0.35;
    if( maximum_number_of_populations == 1 )
        base_population_size         = (int) (9.7 + 13.6*pow((double) number_of_parameters,0.5));
    else
        base_population_size         = 10;
    distribution_multiplier_decrease = 0.9;
    st_dev_ratio_threshold           = 1.0;
    maximum_no_improvement_stretch   = 25 + number_of_parameters;
  }
  
  ellipsoid_size = number_of_parameters;
  if( problem_index == 13 ) ellipsoid_size = 5;
  number_of_ellipsoids = (number_of_parameters + ellipsoid_size - 1) / ellipsoid_size;
  if( FOS_element_size == -1 )
    FOS_element_size = number_of_parameters;
  else if( FOS_element_size == 1 )
    use_univariate_FOS = 1;

  checkOptions();
}

/**
 * Parses the command line.
 * For options, see printUsage.
 */
void parseCommandLine( int argc, char **argv )
{
  int index;

  index = 1;

  parseOptions( argc, argv, &index );

  parseParameters( argc, argv, &index );
}

/**
 * Parses only the options from the command line.
 */
void parseOptions( int argc, char **argv, int *index )
{
  double dummy;

  write_generational_statistics = 0;
  write_generational_solutions  = 0;
  print_verbose_overview        = 0;
  use_vtr                       = 0;
  use_guidelines                = 0;

  for( ; (*index) < argc; (*index)++ )
  {
    if( argv[*index][0] == '-' )
    {
      /* If it is a negative number, the option part is over */
      if( sscanf( argv[*index], "%lf", &dummy ) && argv[*index][1] != '\0' )
        break;

      if( argv[*index][1] == '\0' )
        optionError( argv, *index );
      else if( argv[*index][2] != '\0' )
        optionError( argv, *index );
      else
      {
        switch( argv[*index][1] )
        {
          case '?': printUsage(); break;
          case 'P': printAllInstalledProblems(); break;
          case 's': write_generational_statistics = 1; break;
          case 'w': write_generational_solutions  = 1; break;
          case 'v': print_verbose_overview        = 1; break;
          case 'r': use_vtr                       = 1; break;
          case 'g': use_guidelines                = 1; break;
          case 'f': parseFOSElementSize( index, argc, argv ); break;
          default : optionError( argv, *index );
        }
      }
    }
    else /* Argument is not an option, so option part is over */
     break;
  }
}

void parseFOSElementSize( int *index, int argc, char** argv )
{
    short noError = 1;

    (*index)++;
    noError = noError && sscanf( argv[*index], "%d", &FOS_element_size );

    if( !noError )
    {
        printf("Error parsing parameters.\n\n");

        printUsage();
    }
}

/**
 * Writes the names of all installed problems to the standard output.
 */
void printAllInstalledProblems( void )
{
  int i, n;

  n = numberOfInstalledProblems();
  printf("Installed optimization problems:\n");
  for( i = 0; i < n; i++ )
    printf("%3d: %s\n", i, installedProblemName( i ));

  exit( 0 );
}

/**
 * Informs the user of an illegal option and exits the program.
 */
void optionError( char **argv, int index )
{
  printf("Illegal option: %s\n\n", argv[index]);

  printUsage();
}

/**
 * Parses only the EA parameters from the command line.
 */
void parseParameters( int argc, char **argv, int *index )
{
  int noError;

  if( (argc - *index) != 15 )
  {
    printf("Number of parameters is incorrect, require 15 parameters (you provided %d).\n\n", (argc - *index));

    printUsage();
  }

  noError = 1;
  noError = noError && sscanf( argv[*index+0], "%d", &problem_index );
  noError = noError && sscanf( argv[*index+1], "%d", &number_of_parameters );
  noError = noError && sscanf( argv[*index+2], "%lf", &lower_user_range );
  noError = noError && sscanf( argv[*index+3], "%lf", &upper_user_range );
  noError = noError && sscanf( argv[*index+4], "%lf", &rotation_angle );
  noError = noError && sscanf( argv[*index+5], "%lf", &tau );
  noError = noError && sscanf( argv[*index+6], "%d", &base_population_size );
  noError = noError && sscanf( argv[*index+7], "%d", &maximum_number_of_populations );
  noError = noError && sscanf( argv[*index+8], "%lf", &distribution_multiplier_decrease );
  noError = noError && sscanf( argv[*index+9], "%lf", &st_dev_ratio_threshold );
  noError = noError && sscanf( argv[*index+10], "%d", &maximum_number_of_evaluations );
  noError = noError && sscanf( argv[*index+11], "%lf", &vtr );
  noError = noError && sscanf( argv[*index+12], "%d", &maximum_no_improvement_stretch );
  noError = noError && sscanf( argv[*index+13], "%lf", &fitness_variance_tolerance );
  noError = noError && sscanf( argv[*index+14], "%lf", &maximum_number_of_seconds );

  if( !noError )
  {
    printf("Error parsing parameters.\n\n");

    printUsage();
  }
}

/**
 * Prints usage information and exits the program.
 */
void printUsage( void )
{
  printf("Usage: AMaLGaM-Full [-?] [-P] [-s] [-w] [-v] [-r] [-g] pro dim low upp rot tau pop nop dmd srt eva vtr imp tol\n");
  printf(" -?: Prints out this usage information.\n");
  printf(" -P: Prints out a list of all installed optimization problems.\n");
  printf(" -s: Enables computing and writing of statistics every generation.\n");
  printf(" -w: Enables writing of solutions and their fitnesses every generation.\n");
  printf(" -v: Enables verbose mode. Prints the settings before starting the run.\n");
  printf(" -r: Enables use of vtr in termination condition (value-to-reach).\n");
  printf(" -g: Uses guidelines to override parameter settings for those parameters\n");
  printf("     for which a guideline is known in literature. These parameters are:\n");
  printf("     tau pop dmd srt imp\n");
  printf("\n");
  printf("  pro: Index of optimization problem to be solved (minimization).\n");
  printf("  dim: Number of parameters.\n");
  printf("  low: Overall initialization lower bound.\n");
  printf("  upp: Overall initialization upper bound.\n");
  printf("  rot: The angle by which to rotate the problem.\n");
  printf("  tau: Selection percentile (tau in [1/pop,1], truncation selection).\n");
  printf("  pop: Population size per normal.\n");
  printf("  nop: The number of populations (parallel runs that initially partition the search space).\n");
  printf("  dmd: The distribution multiplier decreaser (in (0,1), increaser is always 1/dmd).\n");
  printf("  srt: The standard-devation ratio threshold for triggering variance-scaling.\n");
  printf("  eva: Maximum number of evaluations allowed.\n");
  printf("  vtr: The value to reach. If the objective value of the best feasible solution reaches\n");
  printf("       this value, termination is enforced (if -r is specified).\n");
  printf("  imp: Maximum number of subsequent generations without an improvement while the\n");
  printf("       the distribution multiplier is <= 1.0.\n");
  printf("  tol: The tolerance level for fitness variance (i.e. minimum fitness variance)\n");
  printf("  sec: The time limit in seconds.\n");
  exit( 0 );
}

/**
 * Checks whether the selected options are feasible.
 */
void checkOptions( void )
{
  if( number_of_parameters < 1 )
  {
    printf("\n");
    printf("Error: number of parameters < 1 (read: %d). Require number of parameters >= 1.", number_of_parameters);
    printf("\n\n");

    exit( 0 );
  }

  if( ((int) (tau*base_population_size)) <= 0 || tau >= 1 )
  {
    printf("\n");
    printf("Error: tau not in range (read: %e). Require tau in [1/pop,1] (read: [%e,%e]).", tau, 1.0/((double) base_population_size), 1.0);
    printf("\n\n");

    exit( 0 );
  }

  if( base_population_size < 1 )
  {
    printf("\n");
    printf("Error: population size < 1 (read: %d). Require population size >= 1.", base_population_size);
    printf("\n\n");

    exit( 0 );
  }

  if( maximum_number_of_populations < 1 )
  {
    printf("\n");
    printf("Error: number of populations < 1 (read: %d). Require number of populations >= 1.", maximum_number_of_populations);
    printf("\n\n");

    exit( 0 );
  }

  if( installedProblemName( problem_index ) == NULL )
  {
    printf("\n");
    printf("Error: unknown index for problem (read index %d).", problem_index );
    printf("\n\n");

    exit( 0 );
  }
}

/**
 * Prints the settings as read from the command line.
 */
void printVerboseOverview( void )
{
  int i;

  printf("### Settings ######################################\n");
  printf("#\n");
  printf("# Statistics writing every generation: %s\n", write_generational_statistics ? "enabled" : "disabled");
  printf("# Population file writing            : %s\n", write_generational_solutions ? "enabled" : "disabled");
  printf("# Use of value-to-reach (vtr)        : %s\n", use_vtr ? "enabled" : "disabled");
  printf("#\n");
  printf("###################################################\n");
  printf("#\n");
  printf("# Problem                 = %s\n", installedProblemName( problem_index ));
  printf("# Number of parameters    = %d\n", number_of_parameters);
  printf("# Initialization ranges   = ");
  for( i = 0; i < number_of_parameters; i++ )
  {
    printf("x_%d: [%e;%e]", i, lower_init_ranges[i], upper_init_ranges[i]);
    if( i < number_of_parameters-1 )
      printf("\n#                           ");
  }
  printf("\n");
  printf("# Boundary ranges         = ");
  for( i = 0; i < number_of_parameters; i++ )
  {
    printf("x_%d: [%e;%e]", i, lower_range_bounds[i], upper_range_bounds[i]);
    if( i < number_of_parameters-1 )
      printf("\n#                           ");
  }
  printf("\n");
  printf("# Rotation angle          = %e\n", rotation_angle);
  printf("# Tau                     = %e\n", tau);
  printf("# Population size/normal  = %d\n", base_population_size);
  printf("# FOS element size        = %d\n", FOS_element_size);
  printf("# Number of populations   = %d\n", maximum_number_of_populations);
  printf("# Dis. mult. decreaser    = %e\n", distribution_multiplier_decrease);
  printf("# St. dev. rat. threshold = %e\n", st_dev_ratio_threshold);
  printf("# Maximum numb. of eval.  = %d\n", maximum_number_of_evaluations);
  printf("# Value to reach (vtr)    = %e\n", vtr);
  printf("# Max. no improv. stretch = %d\n", maximum_no_improvement_stretch);
  printf("# Fitness var. tolerance  = %e\n", fitness_variance_tolerance);
  printf("# Random seed             = %ld\n", random_seed);
  printf("#\n");
  printf("###################################################\n");
}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/





/*-=-=-=-=-=-=-=-=-=-=-=-=-= Section Random Numbers -=-=-=-=-=-=-=-=-=-=-=-=*/
/**
 * Returns a random double, distributed uniformly between 0 and 1.
 */
double randomRealUniform01( void )
{
  int64_t n26, n27;
  double  result;

  random_seed_changing = (random_seed_changing * 0x5DEECE66DLLU + 0xBLLU) & ((1LLU << 48) - 1);
  n26                  = (int64_t)(random_seed_changing >> (48 - 26));
  random_seed_changing = (random_seed_changing * 0x5DEECE66DLLU + 0xBLLU) & ((1LLU << 48) - 1);
  n27                  = (int64_t)(random_seed_changing >> (48 - 27));
  result               = (((int64_t)n26 << 27) + n27) / ((double) (1LLU << 53));

  return( result );
}

/**
 * Returns a random integer, distributed uniformly between 0 and maximum.
 */
int randomInt( int maximum )
{
  int result;

  result = (int) (((double) maximum)*randomRealUniform01());

  return( result );
}

/**
 * Returns a random double, distributed normally with mean 0 and variance 1.
 */
double random1DNormalUnit( void )
{
  double v1, v2, s, multiplier, value;

  if( haveNextNextGaussian )
  {
    haveNextNextGaussian = 0;

    return( nextNextGaussian );
  }
  else
  {
    do
    {
      v1 = 2 * (randomRealUniform01()) - 1;
      v2 = 2 * (randomRealUniform01()) - 1;
      s = v1 * v1 + v2 * v2;
    } while (s >= 1);

    value                = -2 * log(s)/s;
    multiplier           = value <= 0.0 ? 0.0 : sqrt( value );
    nextNextGaussian     = v2 * multiplier;
    haveNextNextGaussian = 1;

    return( v1 * multiplier );
  }
}

/**
 * Returns a random compact (using integers 0,1,...,n-1) permutation
 * of length n using the Fisher-Yates shuffle.
 */
int *randomPermutation( int n )
{
    int i, j, dummy, *result;

    result = (int *) Malloc( n*sizeof( int ) );
    for( i = 0; i < n; i++ )
        result[i] = i;

    for( i = n-1; i > 0; i-- )
    {
        j         = randomInt( i+1 );
        dummy     = result[j];
        result[j] = result[i];
        result[i] = dummy;
    }

    return( result );
}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/





/*-=-=-=-=-=-=-=-=-=-=-=-=-=-= Section Problems -=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
/**
 * Returns the name of an installed problem.
 */
char *installedProblemName( int index )
{
  switch( index )
  {
    case  0: return( (char *) "Sphere" );
    case  1: return( (char *) "Ellipsoid" );
    case  2: return( (char *) "Cigar" );
    case  3: return( (char *) "Tablet" );
    case  4: return( (char *) "Cigar-tablet" );
    case  5: return( (char *) "Two axes" );
    case  6: return( (char *) "Different powers" );
    case  7: return( (char *) "Rosenbrock" );
    case  8: return( (char *) "Parabolic ridge" );
    case  9: return( (char *) "Sharp ridge" );
    case 10: return( (char *) "Griewank" );
    case 11: return( (char *) "Michalewicz" );
    case 12: return( (char *) "Rastrigin" );
    case 13: return( (char *) "Sum of ellipsoids" );
  }

  return( NULL );
}

/**
 * Returns the number of problems installed.
 */
int numberOfInstalledProblems( void )
{
  static int result = -1;

  if( result == -1 )
  {
    result = 0;
    while( installedProblemName( result ) != NULL )
      result++;
  }

  return( result );
}

/**
 * Returns the lower-range bound of an installed problem.
 */
double installedProblemLowerRangeBound( int index, int dimension )
{
  switch( index )
  {
    case  0: return( sphereFunctionProblemLowerRangeBound( dimension ) );
    case  1: return( ellipsoidFunctionLowerRangeBound( dimension ) );
    case  2: return( cigarFunctionLowerRangeBound( dimension ) );
    case  3: return( tabletFunctionLowerRangeBound( dimension ) );
    case  4: return( cigarTabletFunctionLowerRangeBound( dimension ) );
    case  5: return( twoAxesFunctionLowerRangeBound( dimension ) );
    case  6: return( differentPowersFunctionLowerRangeBound( dimension ) );
    case  7: return( rosenbrockFunctionLowerRangeBound( dimension ) );
    case  8: return( parabolicRidgeFunctionLowerRangeBound( dimension ) );
    case  9: return( sharpRidgeFunctionLowerRangeBound( dimension ) );
    case 10: return( griewankFunctionLowerRangeBound( dimension ) );
    case 11: return( michalewiczFunctionLowerRangeBound( dimension ) );
    case 12: return( rastriginFunctionLowerRangeBound( dimension ) );
    case 13: return( sumOfEllipsoidsFunctionLowerRangeBound( dimension ) );
  }

  return( 0.0 );
}

/**
 * Returns the upper-range bound of an installed problem.
 */
double installedProblemUpperRangeBound( int index, int dimension )
{
  switch( index )
  {
    case  0: return( sphereFunctionProblemUpperRangeBound( dimension ) );
    case  1: return( ellipsoidFunctionUpperRangeBound( dimension ) );
    case  2: return( cigarFunctionUpperRangeBound( dimension ) );
    case  3: return( tabletFunctionUpperRangeBound( dimension ) );
    case  4: return( cigarTabletFunctionUpperRangeBound( dimension ) );
    case  5: return( twoAxesFunctionUpperRangeBound( dimension ) );
    case  6: return( differentPowersFunctionUpperRangeBound( dimension ) );
    case  7: return( rosenbrockFunctionUpperRangeBound( dimension ) );
    case  8: return( parabolicRidgeFunctionUpperRangeBound( dimension ) );
    case  9: return( sharpRidgeFunctionUpperRangeBound( dimension ) );
    case 10: return( griewankFunctionUpperRangeBound( dimension ) );
    case 11: return( michalewiczFunctionUpperRangeBound( dimension ) );
    case 12: return( rastriginFunctionUpperRangeBound( dimension ) );
    case 13: return( sumOfEllipsoidsFunctionUpperRangeBound( dimension ) );
  }

  return( 0.0 );
}

/**
 * Returns whether a parameter is inside the range bound of
 * every problem.
 */
short isParameterInRangeBounds( double parameter, int dimension )
{
  if( parameter < installedProblemLowerRangeBound( problem_index, dimension ) ||
      parameter > installedProblemUpperRangeBound( problem_index, dimension ) ||
      isnan( parameter ) )
  {
    return( 0 );
  }

  return( 1 );
}

/**
 * Returns the value of the single objective
 * and the sum of all constraint violations
 * function after rotating the parameter vector.
 * Both are returned using pointer variables.
 */
void installedProblemEvaluation( int index, double *parameters, double *objective_value, double *constraint_value )
{
  int i, j;
  double *parameters_ellipsoid, *rotated_parameters, *rotated_parameters_ellipsoid;

  number_of_evaluations++;

  if( rotation_angle == 0.0 )
    installedProblemEvaluationWithoutRotation( index, parameters, objective_value, constraint_value );
  else
  {
    rotated_parameters = (double*) Malloc( number_of_parameters * sizeof( double ) );
    for( i = 0; i < number_of_ellipsoids; i++ )
    {
        parameters_ellipsoid = (double*) Malloc( ellipsoid_size * sizeof( double ) );
        for( j = 0; j < ellipsoid_size; j++ )
            parameters_ellipsoid[j] = parameters[ellipsoid_size*i+j];

        rotated_parameters_ellipsoid = matrixVectorMultiplication( rotation_matrix, parameters_ellipsoid, ellipsoid_size, ellipsoid_size );
        for( j = 0; j < ellipsoid_size; j++ )
            rotated_parameters[ellipsoid_size*i+j] = rotated_parameters_ellipsoid[j];

        free( parameters_ellipsoid );
        free( rotated_parameters_ellipsoid );
    }

    installedProblemEvaluationWithoutRotation( index, rotated_parameters, objective_value, constraint_value );

    free( rotated_parameters );
  }

  if( use_vtr && *constraint_value == 0 && *objective_value <= vtr  )
    vtr_hit_status = 1;

  if( (number_of_evaluations+1) % 20000 == 0 )
      evaluations_for_statistics_hit = 1;
}

/**
 * Returns the value of the single objective
 * and the sum of all constraint violations
 * without rotating the parameter vector.
 * Both are returned using pointer variables.
 */
void installedProblemEvaluationWithoutRotation( int index, double *parameters, double *objective_value, double *constraint_value )
{
  *objective_value  = 0.0;
  *constraint_value = 0.0;

  switch( index )
  {
    case  0: sphereFunctionProblemEvaluation( parameters, objective_value, constraint_value ); break;
    case  1: ellipsoidFunctionProblemEvaluation( parameters, objective_value, constraint_value ); break;
    case  2: cigarFunctionProblemEvaluation( parameters, objective_value, constraint_value ); break;
    case  3: tabletFunctionProblemEvaluation( parameters, objective_value, constraint_value ); break;
    case  4: cigarTabletFunctionProblemEvaluation( parameters, objective_value, constraint_value ); break;
    case  5: twoAxesFunctionProblemEvaluation( parameters, objective_value, constraint_value ); break;
    case  6: differentPowersFunctionProblemEvaluation( parameters, objective_value, constraint_value ); break;
    case  7: rosenbrockFunctionProblemEvaluation( parameters, objective_value, constraint_value ); break;
    case  8: parabolicRidgeFunctionProblemEvaluation( parameters, objective_value, constraint_value ); break;
    case  9: sharpRidgeFunctionProblemEvaluation( parameters, objective_value, constraint_value ); break;
    case 10: griewankFunctionProblemEvaluation( parameters, objective_value, constraint_value ); break;
    case 11: michalewiczFunctionProblemEvaluation( parameters, objective_value, constraint_value ); break;
    case 12: rastriginFunctionProblemEvaluation( parameters, objective_value, constraint_value ); break;
    case 13: sumOfEllipsoidsFunctionProblemEvaluation( parameters, objective_value, constraint_value ); break;
  }
}

void sphereFunctionProblemEvaluation( double *parameters, double *objective_value, double *constraint_value )
{
  int    i;
  double result;

  result = 0.0;
  for( i = 0; i < number_of_parameters; i++ )
    result += parameters[i]*parameters[i];

  *objective_value  = result;
  *constraint_value = 0;
}

double sphereFunctionProblemLowerRangeBound( int dimension )
{
  return( -1e+308 );
}

double sphereFunctionProblemUpperRangeBound( int dimension )
{
  return( 1e+308 );
}

void ellipsoidFunctionProblemEvaluation( double *parameters, double *objective_value, double *constraint_value )
{
  int    i;
  double result;

  result = 0.0;
  for( i = 0; i < number_of_parameters; i++ )
    result += pow( 10.0, 6.0*(((double) (i))/((double) (number_of_parameters-1))) )*parameters[i]*parameters[i];

  *objective_value  = result;
  *constraint_value = 0;
}

double ellipsoidFunctionLowerRangeBound( int dimension )
{
  return( -1e+308 );
}

double ellipsoidFunctionUpperRangeBound( int dimension )
{
  return( 1e+308 );
}

void cigarFunctionProblemEvaluation( double *parameters, double *objective_value, double *constraint_value )
{
  int    i;
  double result;

  result = parameters[0]*parameters[0];
  for( i = 1; i < number_of_parameters; i++ )
  {
    result += pow( 10.0, 6.0 )*parameters[i]*parameters[i];
  }

  *objective_value  = result;
  *constraint_value = 0;
}

double cigarFunctionLowerRangeBound( int dimension )
{
  return( -1e+308 );
}

double cigarFunctionUpperRangeBound( int dimension )
{
  return( 1e+308 );
}

void tabletFunctionProblemEvaluation( double *parameters, double *objective_value, double *constraint_value )
{
  int    i;
  double result;

  result = pow( 10.0, 6.0 )*parameters[0]*parameters[0];
  for( i = 1; i < number_of_parameters; i++ )
    result += parameters[i]*parameters[i];

  *objective_value  = result;
  *constraint_value = 0;
}

double tabletFunctionLowerRangeBound( int dimension )
{
  return( -1e+308 );
}

double tabletFunctionUpperRangeBound( int dimension )
{
  return( 1e+308 );
}

void cigarTabletFunctionProblemEvaluation( double *parameters, double *objective_value, double *constraint_value )
{
  int    i;
  double result;

  result = parameters[0]*parameters[0];
  for( i = 1; i < number_of_parameters-1; i++ )
    result += pow( 10.0, 4.0 )*parameters[i]*parameters[i];
  result += pow( 10.0, 8.0 )*parameters[number_of_parameters-1]*parameters[number_of_parameters-1];

  *objective_value  = result;
  *constraint_value = 0;
}

double cigarTabletFunctionLowerRangeBound( int dimension )
{
  return( -1e+308 );
}

double cigarTabletFunctionUpperRangeBound( int dimension )
{
  return( 1e+308 );
}

void twoAxesFunctionProblemEvaluation( double *parameters, double *objective_value, double *constraint_value )
{
  int    i;
  double result;

  result = 0.0;
  for( i = 0; i <= (number_of_parameters/2)-1; i++ )
    result += pow( 10.0, 6.0 )*parameters[i]*parameters[i];
  for( i = (number_of_parameters/2); i < number_of_parameters; i++ )
    result += parameters[i]*parameters[i];

  *objective_value  = result;
  *constraint_value = 0;
}

double twoAxesFunctionLowerRangeBound( int dimension )
{
  return( -1e+308 );
}

double twoAxesFunctionUpperRangeBound( int dimension )
{
  return( 1e+308 );
}

void differentPowersFunctionProblemEvaluation( double *parameters, double *objective_value, double *constraint_value )
{
  int    i;
  double result;

  result = 0.0;
  for( i = 0; i < number_of_parameters; i++ )
    result += pow( fabs(parameters[i]), 2.0 + 10.0*(((double) (i))/((double) (number_of_parameters-1))) );

  *objective_value  = result;
  *constraint_value = 0;
}

double differentPowersFunctionLowerRangeBound( int dimension )
{
  return( -1e+308 );
}

double differentPowersFunctionUpperRangeBound( int dimension )
{
  return( 1e+308 );
}

void rosenbrockFunctionProblemEvaluation( double *parameters, double *objective_value, double *constraint_value )
{
  int    i;
  double result;

  result = 0.0;
  for( i = 0; i < number_of_parameters-1; i++ )
    result += 100*(parameters[i+1]-parameters[i]*parameters[i])*(parameters[i+1]-parameters[i]*parameters[i]) + (1.0-parameters[i])*(1.0-parameters[i]);

  *objective_value  = result;
  *constraint_value = 0;
}

double rosenbrockFunctionLowerRangeBound( int dimension )
{
  return( -1e+308 );
}

double rosenbrockFunctionUpperRangeBound( int dimension )
{
  return( 1e+308 );
}

void parabolicRidgeFunctionProblemEvaluation( double *parameters, double *objective_value, double *constraint_value )
{
  int    i;
  double sum, result;

  sum = 0;
  for( i = 1; i < number_of_parameters; i++ )
    sum += parameters[i]*parameters[i];

  result = -parameters[0] + 100.0*sum;

  *objective_value  = result;
  *constraint_value = 0;
}

double parabolicRidgeFunctionLowerRangeBound( int dimension )
{
  return( -1e+308 );
}

double parabolicRidgeFunctionUpperRangeBound( int dimension )
{
  return( 1e+308 );
}

void sharpRidgeFunctionProblemEvaluation( double *parameters, double *objective_value, double *constraint_value )
{
  int    i;
  double sum, result;

  sum = 0;
  for( i = 1; i < number_of_parameters; i++ )
    sum += parameters[i]*parameters[i];

  result = -parameters[0] + 100.0*sqrt( sum );

  *objective_value  = result;
  *constraint_value = 0;
}

double sharpRidgeFunctionLowerRangeBound( int dimension )
{
  return( -1e+308 );
}

double sharpRidgeFunctionUpperRangeBound( int dimension )
{
  return( 1e+308 );
}

void griewankFunctionProblemEvaluation( double *parameters, double *objective_value, double *constraint_value )
{
  int    i;
  double yi, sum, prod, result;

  sum  = 0;
  prod = 1.0;
  for( i = 0; i < number_of_parameters; i++ )
  {
    yi    = parameters[i] - 100.0;
    sum  += yi*yi;

    yi    = (parameters[i] - 100.0)/sqrt((double) (i+1));
    prod *= cos( yi );
  }

  result = sum/4000.0 - prod + 1.0;

  *objective_value  = result;
  *constraint_value = 0;
}

double griewankFunctionLowerRangeBound( int dimension )
{
  return( -1e+308 );
}

double griewankFunctionUpperRangeBound( int dimension )
{
  return( 1e+308 );
}

void michalewiczFunctionProblemEvaluation( double *parameters, double *objective_value, double *constraint_value )
{
  int    i;
  double result;

  result  = 0.0;
  for( i = 0; i < number_of_parameters; i++ )
    result += -sin(parameters[i])*pow(sin(((i+1)*parameters[i]*parameters[i])/PI),20.0);

  *objective_value  = result;
  *constraint_value = 0;
}

double michalewiczFunctionLowerRangeBound( int dimension )
{
  return( 0.0 );
}

double michalewiczFunctionUpperRangeBound( int dimension )
{
  return( PI );
}

void rastriginFunctionProblemEvaluation( double *parameters, double *objective_value, double *constraint_value )
{
  int    i;
  double result;

  result  = 10*number_of_parameters;
  for( i = 0; i < number_of_parameters; i++ )
    result += parameters[i]*parameters[i] - 10.0*cos(2.0*PI*parameters[i]);

  *objective_value  = result;
  *constraint_value = 0;
}

double rastriginFunctionLowerRangeBound( int dimension )
{
  return( -1e+308 );
}

double rastriginFunctionUpperRangeBound( int dimension )
{
  return( 1e+308 );
}

void sumOfEllipsoidsFunctionProblemEvaluation( double *parameters, double *objective_value, double *constraint_value )
{
    int    i, j;
    double result;

    result = 0.0;
    for( i = 0; i < number_of_parameters; i++ )
    {
        j = i % ellipsoid_size;
        result += pow( 10.0, 6.0*(((double) (j))/((double) (ellipsoid_size-1))) )*parameters[i]*parameters[i];
    }

    *objective_value  = result;
    *constraint_value = 0;
}

double sumOfEllipsoidsFunctionLowerRangeBound( int dimension )
{
    return( -1e+308 );
}

double sumOfEllipsoidsFunctionUpperRangeBound( int dimension )
{
    return( 1e+308 );
}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/





/*-=-=-=-=-=-=-=-=-=-=-=-=- Section Initialize -=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
/**
 * Performs initializations that are required before starting a run.
 */
void initialize( void )
{
    total_number_of_generations = 0;
    number_of_evaluations = 0;
    number_of_populations = 0;
    number_of_subgenerations_per_population_factor = 8;
    vtr_hit_status = 0;
    distribution_multiplier_increase = 1.0/distribution_multiplier_decrease;
	block_start = 0;
	block_size = number_of_parameters;
	FOS_element_ub = number_of_parameters;

    initializeRandomNumberGenerator();

    initializeParameterRangeBounds();

    initializeMemory();

    initializeObjectiveRotationMatrix();
}

/**
 * Initializes the memory.
 */
void initializeMemory( void )
{
    populations                      = (double ***) Malloc( maximum_number_of_populations*sizeof( double ** ) );
    population_sizes                 = (int *) Malloc( maximum_number_of_populations*sizeof( int ) );
    selection_sizes                  = (int *) Malloc( maximum_number_of_populations*sizeof( int ) );
    populations_terminated           = (short *) Malloc( maximum_number_of_populations*sizeof( short ) );
    no_improvement_stretch           = (int *) Malloc( maximum_number_of_populations*sizeof( int ) );
    objective_values                 = (double **) Malloc( maximum_number_of_populations*sizeof( double * ) );
    population_best_obj_val          = (double *) Malloc( maximum_number_of_populations*sizeof( double ) );
    population_best_constraint_val   = (double *) Malloc( maximum_number_of_populations*sizeof( double ) );
    constraint_values                = (double **) Malloc( maximum_number_of_populations*sizeof( double * ) );
    ranks                            = (double **) Malloc( maximum_number_of_populations*sizeof( double * ) );
    selections                       = (double ***) Malloc( maximum_number_of_populations*sizeof( double ** ) );
    objective_values_selections      = (double **) Malloc( maximum_number_of_populations*sizeof( double * ) );
    constraint_values_selections     = (double **) Malloc( maximum_number_of_populations*sizeof( double * ) );
    cholesky_factors_lower_triangle  = (double ****) Malloc(maximum_number_of_populations*sizeof(double *** ) );
    mean_vectors                     = (double **) Malloc( maximum_number_of_populations*sizeof( double * ) );
    mean_vectors_previous            = (double **) Malloc( maximum_number_of_populations*sizeof( double * ) );
    covariance_matrices              = (double ****) Malloc( maximum_number_of_populations * sizeof( double ***) );
    full_covariance_matrix           = (double ***) Malloc( maximum_number_of_populations * sizeof( double **) );
    distribution_multipliers         = (double *) Malloc( maximum_number_of_populations*sizeof( double ) );
    samples_drawn_from_normal        = (int *) Malloc( maximum_number_of_populations*sizeof( int ) );
    out_of_bounds_draws              = (int *) Malloc( maximum_number_of_populations*sizeof( int ) );
    number_of_generations            = (int *) Malloc( maximum_number_of_populations*sizeof( int ) );
    linkage_model                    = (FOS **) Malloc( maximum_number_of_populations*sizeof( FOS *) );

	checked_matrix 					 = (int **) Malloc( number_of_parameters*sizeof( int * ) );
	dependency_matrix                = (double **) Malloc( number_of_parameters*sizeof( double * ) );
	checked_pairs                = (double **) Malloc( number_of_parameters*sizeof( double * ) );
	first_individual                 = (double *) Malloc( number_of_parameters*sizeof( double ) );
	second_individual                = (double *) Malloc( number_of_parameters*sizeof( double ) );
	fitness_of_first_individual      = (double *) Malloc( (number_of_parameters + 1)*sizeof( double * ) );
	dependency_pairs = (int **) Malloc( ((number_of_parameters*number_of_parameters)/2)*sizeof(int * ) );
    
	// rv gomea flags
    static_linkage_tree = 1;
	dependency_learning = 1;
	evolve_learning = 1;
	pruning_ub = 100;
	continued_learning=1;

    // to create a marginal product FOS we need to have a sparse tree, this flag should be one
    for(int j = 0; j < number_of_parameters; j++ ){
        dependency_matrix[j] = (double *) Malloc( number_of_parameters*sizeof( double ) );
        checked_matrix[j] = (int *) Malloc( number_of_parameters*sizeof( int ) );
    }
    pairs_per_run = number_of_parameters;
    number_of_checked_pairs = 0;
    int counter = 0;
    for (int i = 0; i < number_of_parameters; i++) {
        for (int j = i + 1; j < number_of_parameters; j++) {
            // add pairs to evaluate to the list
            dependency_pairs[counter] = (int *) Malloc(2 * sizeof(int));
            dependency_pairs[counter][0] = i;
            dependency_pairs[counter][1] = j;
            counter++;
        }
    }
    number_of_pairs = counter;
    for (int i = counter - 1; i >= 0; --i) {
        //generate a random number [0, n-1]
        int j = randomInt(i+1);

        //swap the last element with element at random index
        int *temp = dependency_pairs[i];
        dependency_pairs[i] = dependency_pairs[j];
        dependency_pairs[j] = temp;
    }

    // fill matrix with 0s
    for (int i = 0; i < number_of_parameters; i++) {
        for (int j = i; j < number_of_parameters; j++) {
            dependency_matrix[i][j] = 0.0;
            dependency_matrix[j][i] = 0.0;
            checked_matrix[i][j] = 0;
            checked_matrix[j][i] = 0;
        }
    }
    current_waiting_position = 0;
}

int *matchFOSElements( FOS *new_FOS, FOS *prev_FOS )
{
    int      i, j, a, b, matches, *permutation, 
            **FOS_element_similarity_matrix;

    permutation = (int *) Malloc( new_FOS->length*sizeof(int));
    FOS_element_similarity_matrix = (int**) Malloc((prev_FOS->length-number_of_parameters)*sizeof(int*));
    for( i = 0; i < prev_FOS->length-number_of_parameters; i++ )
        FOS_element_similarity_matrix[i] = (int*) Malloc((new_FOS->length-number_of_parameters)*sizeof(int));

    for( i = 0; i < number_of_parameters; i++ )
    {
        for( j = 0; j < number_of_parameters; j++ )
        {
            if( prev_FOS->sets[i][0] == new_FOS->sets[j][0] )
            {
                permutation[i] = j;
                break;
            }
        }
    }
    for( i = number_of_parameters; i < prev_FOS->length; i++ )
    {
        for( j = number_of_parameters; j < new_FOS->length; j++ )
        {
            a = 0; b = 0;
            matches = 0;
            while( a < prev_FOS->set_length[i] && b < new_FOS->set_length[j] )
            {
                if( prev_FOS->sets[i][a] < new_FOS->sets[j][b] )
                {
                    a++;
                }
                else if( prev_FOS->sets[i][a] > new_FOS->sets[j][b] )
                {
                    b++;
                }
                else
                {
                    a++;
                    b++;
                    matches++;
                }
            }
            FOS_element_similarity_matrix[i-number_of_parameters][j-number_of_parameters] = (int) 10000*(2.0*matches/(prev_FOS->set_length[i]+new_FOS->set_length[j]));
        }
    }

    for( i = 0; i < new_FOS->length; i++ )
    {
        int max_index = 0;
        int max_similarity = -1;
        for( j = number_of_parameters; j < prev_FOS->length; j++ )
        {
            if(FOS_element_similarity_matrix[j-number_of_parameters][i-number_of_parameters]>max_similarity){
                max_index = j;
                max_similarity = FOS_element_similarity_matrix[j-number_of_parameters][i-number_of_parameters];
            }
        }
        permutation[i] = max_index;
    }

    for( i = 0; i < prev_FOS->length-number_of_parameters; i++ )
        free( FOS_element_similarity_matrix[i] );
    free( FOS_element_similarity_matrix );

    return( permutation );
}

void initializeNewPopulationMemory( int population_index )
{
    int j;

    populations[population_index] = (double **) Malloc( population_sizes[population_index]*sizeof( double * ) );
    for( j = 0; j < population_sizes[population_index]; j++ )
        populations[population_index][j] = (double *) Malloc( number_of_parameters*sizeof( double ) );

    objective_values[population_index] = (double *) Malloc( population_sizes[population_index]*sizeof( double ) );

    constraint_values[population_index] = (double *) Malloc( population_sizes[population_index]*sizeof( double ) );

    ranks[population_index] = (double *) Malloc( population_sizes[population_index]*sizeof( double ) );

    selections[population_index] = (double **) Malloc( selection_sizes[population_index]*sizeof( double * ) );
    for( j = 0; j < selection_sizes[population_index]; j++ )
        selections[population_index][j] = (double *) Malloc( number_of_parameters*sizeof( double ) );

    objective_values_selections[population_index] = (double *) Malloc( selection_sizes[population_index]*sizeof( double ) );

    constraint_values_selections[population_index] = (double *) Malloc( selection_sizes[population_index]*sizeof( double ) );

    mean_vectors[population_index] = (double *) Malloc( number_of_parameters*sizeof( double ) );

    mean_vectors_previous[population_index] = (double *) Malloc( number_of_parameters*sizeof( double ) );
    
	full_covariance_matrix[population_index] = (double **) Malloc( number_of_parameters*sizeof( double * ) );
	for( j = 0; j < number_of_parameters; j++ )
		full_covariance_matrix[population_index][j] = (double *) Malloc( number_of_parameters*sizeof( double ) );

    populations_terminated[population_index] = 0;

    no_improvement_stretch[population_index] = 0;

    number_of_generations[population_index]  = 0;
}

void initializeNewPopulationHarikLobo( int population_index )
{
    if( population_index == 0 )
        population_sizes[population_index] = base_population_size;
    else
        population_sizes[population_index] = 2*population_sizes[population_index-1];
    selection_sizes[population_index] = (double) (tau * population_sizes[population_index] );
    number_of_populations++;

    initializeNewPopulation( population_index );
}

void initializeNewPopulation( int population_index )
{
    initializeNewPopulationMemory( population_index );

    initializePopulationAndFitnessValues( population_index );
    
    initializeDistributionMultipliers( population_index );

    computeRanksForOnePopulation( population_index );
}

/**
 * Initializes the random number generator.
 */
void initializeRandomNumberGenerator( void )
{
    //struct tm *timep;
    struct timeval tv;
    //time_t t;

    while( random_seed_changing == 0 )
    {
        //t                    = time( NULL );
        //timep                = localtime( &t );
    gettimeofday(&tv, NULL);
        random_seed_changing = (int64_t) tv.tv_usec; //((60*(long) timep->tm_min)) + (60*60*(long) timep->tm_hour) + ((long) timep->tm_sec);
    //printf("seed_changing: %30ld\n", random_seed_changing);
        random_seed_changing = (random_seed_changing/((int) (9.99*randomRealUniform01())+1))*(((int) (randomRealUniform01()*1000000.0))%10);
    }

    random_seed = random_seed_changing;
}

/**
 * Initializes the parameter range bounds.
 */
void initializeParameterRangeBounds( void )
{
  int i;

  lower_range_bounds = (double *) Malloc( number_of_parameters*sizeof( double ) );
  upper_range_bounds = (double *) Malloc( number_of_parameters*sizeof( double ) );
  lower_init_ranges  = (double *) Malloc( number_of_parameters*sizeof( double ) );
  upper_init_ranges  = (double *) Malloc( number_of_parameters*sizeof( double ) );

  for( i = 0; i < number_of_parameters; i++ )
  {
    lower_range_bounds[i] = installedProblemLowerRangeBound( problem_index, i );
    upper_range_bounds[i] = installedProblemUpperRangeBound( problem_index, i );
  }

  for( i = 0; i < number_of_parameters; i++ )
  {
    lower_init_ranges[i] = lower_user_range;
    if( lower_user_range < lower_range_bounds[i] )
      lower_init_ranges[i] = lower_range_bounds[i];

    upper_init_ranges[i] = upper_user_range;
    if( upper_user_range > upper_range_bounds[i] )
      upper_init_ranges[i] = upper_range_bounds[i];
  }
}

void initializeFOS( int population_index )
{
    int i;
        
	FOS *new_FOS = (FOS*) Malloc(sizeof(FOS));
	new_FOS->length      = (number_of_parameters + FOS_element_size - 1) / FOS_element_size;
	new_FOS->sets        = (int **) Malloc( new_FOS->length*sizeof( int * ) );
	new_FOS->set_length = (int *) Malloc( new_FOS->length*sizeof( int ) );
	for( i = 0; i < new_FOS->length; i++ )
	{
		new_FOS->sets[i] = (int *) Malloc( FOS_element_size*sizeof( int ) );
		new_FOS->set_length[i] = 0;
	}

	for( i = 0; i < number_of_parameters; i++ )
	{
		new_FOS->sets[i/FOS_element_size][i%FOS_element_size] = i;
		new_FOS->set_length[i/FOS_element_size]++;
	}
	linkage_model[population_index] = new_FOS;
}

/**
 * Initializes the distribution multipliers.
 */
void initializeDistributionMultipliers( int population_index )
{
    distribution_multipliers[population_index] = 1.0;
}

void initializeCovarianceMatrices( int population_index )
{
    int j, k;

    covariance_matrices[population_index] = (double ***) Malloc( linkage_model[population_index]->length*sizeof( double ** ) );
    for( j = 0; j < linkage_model[population_index]->length; j++ )
    {
      covariance_matrices[population_index][j] = (double **) Malloc( linkage_model[population_index]->set_length[j]*sizeof( double * ) );
      for( k = 0; k < linkage_model[population_index]->set_length[j]; k++ )
          covariance_matrices[population_index][j][k] = (double *) Malloc( linkage_model[population_index]->set_length[j]*sizeof( double ) );
    }
}

/**
 * Initializes the populations and the fitness values.
 */
void initializePopulationAndFitnessValues( int population_index )
{
    int     j, k;

    for( j = 0; j < population_sizes[population_index]; j++ )
    {
        for( k = 0; k < number_of_parameters; k++ )
            populations[population_index][j][k] = lower_init_ranges[k] + (upper_init_ranges[k] - lower_init_ranges[k])*randomRealUniform01();

        installedProblemEvaluation( problem_index, populations[population_index][j], &(objective_values[population_index][j]), &(constraint_values[population_index][j]) );
        if( j == 0 )
        {
            population_best_obj_val[population_index] = objective_values[population_index][j];
            population_best_constraint_val[population_index] = constraint_values[population_index][j];
        }
        else if( betterFitness(objective_values[population_index][j], constraint_values[population_index][j], population_best_obj_val[population_index], population_best_constraint_val[population_index] ) ){
            population_best_obj_val[population_index] = objective_values[population_index][j];
            population_best_constraint_val[population_index] = constraint_values[population_index][j];
        }
    }
}

/**
 * Computes the rotation matrix to be applied to any solution
 * before evaluating it (i.e. turns the evaluation functions
 * into rotated evaluation functions).
 */
void initializeObjectiveRotationMatrix( void )
{
  int      i, j, index0, index1;
  double **matrix, **product, theta, cos_theta, sin_theta;

  if( rotation_angle == 0.0 )
    return;

  matrix = (double **) Malloc( ellipsoid_size*sizeof( double * ) );
  for( i = 0; i < ellipsoid_size; i++ )
    matrix[i] = (double *) Malloc( ellipsoid_size*sizeof( double ) );

  rotation_matrix = (double **) Malloc( ellipsoid_size*sizeof( double * ) );
  for( i = 0; i < ellipsoid_size; i++ )
    rotation_matrix[i] = (double *) Malloc( ellipsoid_size*sizeof( double ) );

  /* Initialize the rotation matrix to the identity matrix */
  for( i = 0; i < ellipsoid_size; i++ )
  {
    for( j = 0; j < ellipsoid_size; j++ )
      rotation_matrix[i][j] = 0.0;
    rotation_matrix[i][i] = 1.0;
  }

  /* Construct all rotation matrices (quadratic number) and multiply */
  theta     = (rotation_angle/180.0)*PI;
  cos_theta = cos( theta );
  sin_theta = sin( theta );
  for( index0 = 0; index0 < ellipsoid_size-1; index0++ )
  {
    for( index1 = index0+1; index1 < ellipsoid_size; index1++ )
    {
      for( i = 0; i < ellipsoid_size; i++ )
      {
        for( j = 0; j < ellipsoid_size; j++ )
          matrix[i][j] = 0.0;
        matrix[i][i] = 1.0;
      }
      matrix[index0][index0] = cos_theta;
      matrix[index0][index1] = -sin_theta;
      matrix[index1][index0] = sin_theta;
      matrix[index1][index1] = cos_theta;

      product = matrixMatrixMultiplication( matrix, rotation_matrix, ellipsoid_size, ellipsoid_size, ellipsoid_size );
      for( i = 0; i < ellipsoid_size; i++ )
        for( j = 0; j < ellipsoid_size; j++ )
          rotation_matrix[i][j] = product[i][j];

      for( i = 0; i < ellipsoid_size; i++ )
        free( product[i] );
      free( product );
    }
  }

  for( i = 0; i < ellipsoid_size; i++ )
    free( matrix[i] );
  free( matrix );
}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/





/*-=-=-=-=-=-=-=-=-=-=-=-=-=- Section Ranking -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
/**
 * Computes the ranks of the solutions in all populations.
 */
void computeRanks( void )
{
  int i;

  for( i = 0; i < number_of_populations; i++ )
      computeRanksForOnePopulation( i );
}

/**
 * Computes the ranks of the solutions in one population.
 */
void computeRanksForOnePopulation( int population_index )
{
  int i, *sorted, rank;

  if( !populations_terminated[population_index] )
  {
      sorted = mergeSortFitness( objective_values[population_index], constraint_values[population_index], population_sizes[population_index] );

      rank                               = 0;
      ranks[population_index][sorted[0]] = rank;
      for( i = 1; i < population_sizes[population_index]; i++ )
      {
        if( objective_values[population_index][sorted[i]] != objective_values[population_index][sorted[i-1]] )
          rank++;

        ranks[population_index][sorted[i]] = rank;
      }

      free( sorted );
  }

}

/**
 * Computes the distance between two solutions a and b as
 * the Euclidean distance in parameter space.
 */
double distanceInParameterSpace( double *solution_a, double *solution_b )
{
  int    i;
  double value, result;

  result = 0.0;
  for( i = 0; i < number_of_parameters; i++ )
  {
    value   = solution_b[i] - solution_a[i];
    result += value*value;
  }
  result = sqrt( result );

  return( result );
}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/





/*-=-=-=-=-=-=-=-=-=-=-=-=-=- Section Output =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
/**
 * Writes (appends) statistics about the current generation to a
 * file named "statistics.dat".
 */
void writeGenerationalStatisticsForOnePopulation( int population_index )
{
    int     j;
    char    string[1000];
    double  population_objective_avg, population_objective_var, population_objective_best, population_objective_worst,
            population_constraint_avg, population_constraint_var, population_constraint_best, population_constraint_worst,
            overall_best_obj, overall_best_cons;
    FILE   *file;

    /* Average, best and worst */
    population_objective_avg    = 0.0;
    population_constraint_avg   = 0.0;
    population_objective_best   = population_best_obj_val[population_index];
    population_constraint_best  = population_best_constraint_val[population_index];
    population_objective_worst  = objective_values[population_index][0];
    population_constraint_worst = constraint_values[population_index][0];
    for( j = 0; j < population_sizes[population_index]; j++ )
    {
        population_objective_avg  += objective_values[population_index][j];
        population_constraint_avg += constraint_values[population_index][j];
        if( betterFitness( population_objective_worst, population_constraint_worst, objective_values[population_index][j], constraint_values[population_index][j] ) )
        {
            population_objective_worst = objective_values[population_index][j];
            population_constraint_worst = constraint_values[population_index][j];
        }
    }
    population_objective_avg  = population_objective_avg / ((double) population_sizes[population_index]);
    population_constraint_avg = population_constraint_avg / ((double) population_sizes[population_index]);

    overall_best_obj = population_objective_best;
    overall_best_cons = population_constraint_best;
    for( j = 0; j < number_of_populations; j++ )
    {
        if( betterFitness( population_best_obj_val[j], population_best_constraint_val[j], overall_best_obj, overall_best_cons ) )
        {
            overall_best_obj = population_best_obj_val[j];
            overall_best_cons = population_best_constraint_val[j];
        }
    }

    /* Variance */
    population_objective_var    = 0.0;
    population_constraint_var   = 0.0;
    for( j = 0; j < population_sizes[population_index]; j++ )
    {
        population_objective_var  += (objective_values[population_index][j] - population_objective_avg)*(objective_values[population_index][j] - population_objective_avg);
        population_constraint_var += (constraint_values[population_index][j] - population_constraint_avg)*(constraint_values[population_index][j] - population_constraint_avg);
    }
    population_objective_var  = population_objective_var / ((double) population_sizes[population_index]);
    population_constraint_var = population_constraint_var / ((double) population_sizes[population_index]);

    if( population_objective_var <= 0.0 )
        population_objective_var = 0.0;
    if( population_constraint_var <= 0.0 )
        population_constraint_var = 0.0;

    /* Then write them */
    file = NULL;
    if( total_number_of_generations == 0 )
    {
        file = fopen( "statistics.dat", "w" );

        sprintf( string, "# Generation  Evaluations  Time(s)  Best-obj. Best-cons. [Pop.index  Subgen.  Pop.size  Dis.mult.[0]  Pop.best.obj. Pop.avg.obj.  Pop.var.obj. Pop.worst.obj.  Pop.best.con. Pop.avg.con.  Pop.var.con. Pop.worst.con.]\n" );
        fputs( string, file );
    }
    else
        file = fopen( "statistics.dat", "a" );

    sprintf( string, "%10d %10d %11.3lf %15.10e %13e  ", total_number_of_generations, number_of_evaluations, getTimer(), overall_best_obj, overall_best_cons );
    fputs( string, file );

    sprintf( string, "[ %4d %6d %10d %13e %13e %13e %13e %13e %13e %13e %13e %13e ]\n", population_index, number_of_generations[population_index], population_sizes[population_index], distribution_multipliers[population_index], population_objective_best, population_objective_avg, population_objective_var, population_objective_worst, population_constraint_best, population_constraint_avg, population_constraint_var, population_constraint_worst );
    fputs( string, file );

    fclose( file );
}

/**
 * Writes the solutions to various files. The filenames
 * contain the generation. If the flag final is set
 * (final != 0), the generation number in the filename
 * is replaced with the word "final".
 *
 * all_populations_generation_xxxxx.dat : all populations combined
 * population_xxxxx_generation_xxxxx.dat: the individual populations
 * selection_xxxxx_generation_xxxxx.dat : the individual selections
 */
void writeGenerationalSolutions( short final )
{
  int   i, j, k;
  char  string[1000];
  FILE *file_all, *file_population, *file_selection;

  if( final )
    sprintf( string, "all_populations_generation_final.dat" );
  else
    sprintf( string, "all_populations_generation_%05d.dat", total_number_of_generations );
  file_all = fopen( string, "w" );

  for( i = 0; i < number_of_populations; i++ )
  {
    if( final )
      sprintf( string, "population_%05d_generation_final.dat", i );
    else
      sprintf( string, "population_%05d_generation_%05d.dat", i, number_of_generations[i] );
    file_population = fopen( string, "w" );

    if( number_of_generations[i] > 0 && !final )
    {
      sprintf( string, "selection_%05d_generation_%05d.dat", i, number_of_generations[i]-1 );
      file_selection = fopen( string, "w" );
    }

    /* Populations */
    for( j = 0; j < population_sizes[i]; j++ )
    {
      for( k = 0; k < number_of_parameters; k++ )
      {
        sprintf( string, "%13e", populations[i][j][k] );
        fputs( string, file_all );
        fputs( string, file_population );
        if( k < number_of_parameters-1 )
        {
          sprintf( string, " " );
          fputs( string, file_all );
          fputs( string, file_population );
        }
      }
      sprintf( string, "     " );
      fputs( string, file_all );
      fputs( string, file_population );
      sprintf( string, "%13e %13e", objective_values[i][j], constraint_values[i][j] );
      fputs( string, file_all );
      fputs( string, file_population );
      sprintf( string, "\n" );
      fputs( string, file_all );
      fputs( string, file_population );
    }

    fclose( file_population );

    /* Selections */
    if( number_of_generations[i] > 0 && !final )
    {
      for( j = 0; j < selection_sizes[i]; j++ )
      {
        for( k = 0; k < number_of_parameters; k++ )
        {
         sprintf( string, "%13e", selections[i][j][k] );
         fputs( string, file_selection );
         if( k < number_of_parameters-1 )
         {
           sprintf( string, " " );
           fputs( string, file_selection );
         }
         sprintf( string, "     " );
         fputs( string, file_selection );
        }
        sprintf( string, "%13e %13e", objective_values_selections[i][j], constraint_values_selections[i][j] );
        fputs( string, file_selection );
        sprintf( string, "\n" );
        fputs( string, file_selection );
      }
      fclose( file_selection );
    }
  }

  fclose( file_all );

  writeGenerationalSolutionsBest( final );
}

/**
 * Writes the best solution (measured in the single
 * available objective) to a file named
 * best_generation_xxxxx.dat where xxxxx is the
 * generation number. If the flag final is set
 * (final != 0), the generation number in the filename
 * is replaced with the word "final".The output
 * file contains the solution values with the
 * dimensions separated by a single white space,
 * followed by five white spaces and then the
 * single objective value for that solution
 * and its sum of constraint violations.
 */
void writeGenerationalSolutionsBest( short final )
{
  int   i, population_index_best, individual_index_best;
  char  string[1000];
  FILE *file;

  /* First find the best of all */
  determineBestSolutionInCurrentPopulations( &population_index_best, &individual_index_best );

  /* Then output it */
  if( final )
    sprintf( string, "best_generation_final.dat" );
  else
    sprintf( string, "best_generation_%05d.dat", total_number_of_generations );
  file = fopen( string, "w" );

  for( i = 0; i < number_of_parameters; i++ )
  {
    sprintf( string, "%13e", populations[population_index_best][individual_index_best][i] );
    fputs( string, file );
    if( i < number_of_parameters-1 )
    {
      sprintf( string, " " );
      fputs( string, file );
    }
  }
  sprintf( string, "     " );
  fputs( string, file );
  sprintf( string, "%13e %13e", objective_values[population_index_best][individual_index_best], constraint_values[population_index_best][individual_index_best] );
  fputs( string, file );
  sprintf( string, "\n" );
  fputs( string, file );

  fclose( file );
}

/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/





/*-=-=-=-=-=-=-=-=-=-=-=-=- Section Termination -=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
/**
 * Returns 1 if termination should be enforced, 0 otherwise.
 */
short checkTerminationCondition( void )
{
  short allTrue;
  int   i;

  if( total_number_of_generations == 0 )
      return( 0 );

  if( checkImmediateTerminationConditions() )
      return( 1 );

  checkAverageFitnessTerminationCondition();

  checkFitnessVarianceTermination();

  checkDistributionMultiplierTerminationCondition();

  if( number_of_populations < maximum_number_of_populations )
    return( 0 );

  allTrue = 1;
  for( i = 0; i < number_of_populations; i++ )
  {
    if( !populations_terminated[i] )
    {
      allTrue = 0;
      break;
    }
  }

  return( allTrue );
}

short checkImmediateTerminationConditions( void )
{
    if( total_number_of_generations == 0 )
        return( 0 );

    if( checkNumberOfEvaluationsTerminationCondition() )
        return( 1 );

    if( checkVTRTerminationCondition() )
        return( 1 );

    if( checkTimeLimitTerminationCondition() )
        return( 1 );

    return( 0 );
}

/**
 * Returns 1 if the maximum number of evaluations
 * has been reached, 0 otherwise.
 */
short checkNumberOfEvaluationsTerminationCondition( void )
{
  if( number_of_evaluations >= maximum_number_of_evaluations && maximum_number_of_evaluations > 0 )
    return( 1 );

  return( 0 );
}

/**
 * Returns 1 if the value-to-reach has been reached (in any population).
 */
short checkVTRTerminationCondition( void )
{
  return( use_vtr && vtr_hit_status );
}

short checkTimeLimitTerminationCondition( void )
{
    return( maximum_number_of_seconds > 0 && getTimer() > maximum_number_of_seconds );
}

void checkAverageFitnessTerminationCondition( void )
{
    int i, j;
    double *average_objective_values, *average_constraint_values;

    average_objective_values = (double*) Malloc( number_of_populations * sizeof(double) );
    average_constraint_values = (double*) Malloc( number_of_populations * sizeof(double) );
    for( i = number_of_populations-1; i >= 0; i-- )
    {
        average_objective_values[i] = 0;
        average_constraint_values[i] = 0;
        for( j = 0; j < population_sizes[i]; j++ )
        {
            average_objective_values[i] += objective_values[i][j];
            average_constraint_values[i] += constraint_values[i][j];
        }
        average_objective_values[i] /= population_sizes[i];
        average_constraint_values[i] /= population_sizes[i];
        if( i < number_of_populations-1 && betterFitness(average_objective_values[i+1], average_constraint_values[i+1], average_objective_values[i], average_constraint_values[i]) )
        {
            for( j = i; j >= 0; j-- )
                populations_terminated[j] = 1;
            break;
        }
    }
    free( average_objective_values );
    free( average_constraint_values );
}

/**
 * Determines which solution is the best of all solutions
 * in all current populations.
 */
void determineBestSolutionInCurrentPopulations( int *population_of_best, int *index_of_best )
{
  int i, j;

  (*population_of_best) = 0;
  (*index_of_best)      = 0;
  for( i = 0; i < number_of_populations; i++ )
  {
    for( j = 0; j < population_sizes[i]; j++ )
    {
      if( betterFitness( objective_values[i][j], constraint_values[i][j],
                         objective_values[(*population_of_best)][(*index_of_best)], constraint_values[(*population_of_best)][(*index_of_best)] ) )
      {
        (*population_of_best) = i;
        (*index_of_best)      = j;
      }
    }
  }
}

/**
 * Checks whether the fitness variance in any population
 * has become too small (user-defined tolerance).
 */
void checkFitnessVarianceTermination( void )
{
  int i;

  for( i = 0; i < number_of_populations; i++ )
  {
    if( !populations_terminated[i] )
      if( checkFitnessVarianceTerminationSinglePopulation( i ) )
        populations_terminated[i] = 1;
  }
}

/**
 * Returns 1 if the fitness variance in a specific population
 * has become too small (user-defined tolerance).
 */
short checkFitnessVarianceTerminationSinglePopulation( int population_index )
{
  int    i;
  double objective_avg, objective_var;

  objective_avg = 0.0;
  for( i = 0; i < population_sizes[population_index]; i++ )
    objective_avg  += objective_values[population_index][i];
  objective_avg = objective_avg / ((double) population_sizes[population_index]);

  objective_var = 0.0;
  for( i = 0; i < population_sizes[population_index]; i++ )
    objective_var  += (objective_values[population_index][i]-objective_avg)*(objective_values[population_index][i]-objective_avg);
  objective_var = objective_var / ((double) population_sizes[population_index]);

  if( objective_var <= 0.0 )
    objective_var = 0.0;

  if( objective_var <= fitness_variance_tolerance )
    return( 1 );

  return( 0 );
}

/**
 * Checks whether the distribution multiplier in any population
 * has become too small (1e-10).
 */
void checkDistributionMultiplierTerminationCondition( void )
{
  int i;
  short terminated;

  for( i = 0; i < number_of_populations; i++ )
  {
    if( !populations_terminated[i] )
    {
      terminated = 1;
      if( distribution_multipliers[i] > 1e-10 )
      {
          terminated = 0;
          break;
      }
      populations_terminated[i] = terminated;
    }
  }
}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/





/*-=-=-=-=-=-=-=-=-=-=-=-=-= Section Selection =-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/

/**
 * Performs truncation selection on a single population.
 */
void makeSelectionsForOnePopulation( int population_index )
{
  int i, j, *sorted;

  sorted = mergeSort( ranks[population_index], population_sizes[population_index] );

  if( ranks[population_index][sorted[selection_sizes[population_index]-1]] == 0 )
    makeSelectionsForOnePopulationUsingDiversityOnRank0( population_index );
  else
  {
    for( i = 0; i < selection_sizes[population_index]; i++ )
    {
      for( j = 0; j < number_of_parameters; j++ )
        selections[population_index][i][j] = populations[population_index][sorted[i]][j];

      objective_values_selections[population_index][i]  = objective_values[population_index][sorted[i]];
      constraint_values_selections[population_index][i] = constraint_values[population_index][sorted[i]];
    }
  }

  free( sorted );
}

/**
 * Performs selection from all solutions that have rank 0
 * based on diversity.
 */
void makeSelectionsForOnePopulationUsingDiversityOnRank0( int population_index )
{
  int     i, j, number_of_rank0_solutions, *preselection_indices,
         *selection_indices, index_of_farthest, number_selected_so_far;
  double *nn_distances, distance_of_farthest, value;

  number_of_rank0_solutions = 0;
  for( i = 0; i < population_sizes[population_index]; i++ )
  {
    if( ranks[population_index][i] == 0 )
      number_of_rank0_solutions++;
  }

  preselection_indices = (int *) Malloc( number_of_rank0_solutions*sizeof( int ) );
  j                    = 0;
  for( i = 0; i < population_sizes[population_index]; i++ )
  {
    if( ranks[population_index][i] == 0 )
    {
      preselection_indices[j] = i;
      j++;
    }
  }

  index_of_farthest    = 0;
  distance_of_farthest = objective_values[population_index][preselection_indices[0]];
  for( i = 1; i < number_of_rank0_solutions; i++ )
  {
    if( objective_values[population_index][preselection_indices[i]] > distance_of_farthest )
    {
      index_of_farthest    = i;
      distance_of_farthest = objective_values[population_index][preselection_indices[i]];
    }
  }

  number_selected_so_far                    = 0;
  selection_indices                         = (int *) Malloc( selection_sizes[population_index]*sizeof( int ) );
  selection_indices[number_selected_so_far] = preselection_indices[index_of_farthest];
  preselection_indices[index_of_farthest]   = preselection_indices[number_of_rank0_solutions-1];
  number_of_rank0_solutions--;
  number_selected_so_far++;

  nn_distances = (double *) Malloc( number_of_rank0_solutions*sizeof( double ) );
  for( i = 0; i < number_of_rank0_solutions; i++ )
    nn_distances[i] = distanceInParameterSpace( populations[population_index][preselection_indices[i]], populations[population_index][selection_indices[number_selected_so_far-1]] );

  while( number_selected_so_far < selection_sizes[population_index] )
  {
    index_of_farthest    = 0;
    distance_of_farthest = nn_distances[0];
    for( i = 1; i < number_of_rank0_solutions; i++ )
    {
      if( nn_distances[i] > distance_of_farthest )
      {
        index_of_farthest    = i;
        distance_of_farthest = nn_distances[i];
      }
    }

    selection_indices[number_selected_so_far] = preselection_indices[index_of_farthest];
    preselection_indices[index_of_farthest]   = preselection_indices[number_of_rank0_solutions-1];
    nn_distances[index_of_farthest]           = nn_distances[number_of_rank0_solutions-1];
    number_of_rank0_solutions--;
    number_selected_so_far++;

    for( i = 0; i < number_of_rank0_solutions; i++ )
    {
      value = distanceInParameterSpace( populations[population_index][preselection_indices[i]], populations[population_index][selection_indices[number_selected_so_far-1]] );
      if( value < nn_distances[i] )
        nn_distances[i] = value;
    }
  }

  for( i = 0; i < selection_sizes[population_index]; i++ )
  {
    for( j = 0; j < number_of_parameters; j++ )
      selections[population_index][i][j] = populations[population_index][selection_indices[i]][j];

    objective_values_selections[population_index][i]  = objective_values[population_index][selection_indices[i]];
    constraint_values_selections[population_index][i] = constraint_values[population_index][selection_indices[i]];
  }

  free( nn_distances );
  free( selection_indices );
  free( preselection_indices );
}

/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
int *mpm_number_of_indices;

void printFOS( FOS *fos )
{
    int i,j;
    printf("{");
    for( i = 0; i < fos->length; i++ )
    {
        printf("[");
        for( j = 0; j < fos->set_length[i]; j++ )
        {
            printf("%d", fos->sets[i][j]);
            if( j != fos->set_length[i]-1)
                printf(",");
        }
        printf("]");
        printf("\n");
    }
    printf("}\n");
}

double getSimilarity( int a, int b )
{
    if( pruning_ub < number_of_parameters && mpm_number_of_indices[a] + mpm_number_of_indices[b] > pruning_ub ) return( 0 );
    return( S_matrix[a][b] );
}

double **computeMIMatrix( double **cov_matrix, int n )
{
    int i, j;
    double si, sj, r;

	double **mat = (double**) Malloc( n * sizeof(double*) );
    for( i = 0; i < n; i++ )
	{
		mat[i] = (double*) Malloc( n * sizeof(double) );
	}

	for( i = 0; i < n; i++ )
    {
		mat[i][i] = 1e20;
        for( j = 0; j < i; j++ )
        {
            si = sqrt(cov_matrix[i][i]);
            sj = sqrt(cov_matrix[j][j]);
            r = cov_matrix[i][j]/(si*sj);
            if( r > 1 || r < -1 ) mat[i][j] = 1e10;
            mat[i][j] = log(sqrt(1/(1-r*r)));
            mat[j][i] = mat[i][j];
        }
    }
	return( mat );
}


FOS *learnLinkageTree( double **covariance_matrix , double **dependency_matrix)
{
  char     done;
  int      i, j, r0, r1, rswap, *indices, *order, *sorted,
          FOS_index, **mpm, mpm_length,
          **mpm_new, *mpm_new_number_of_indices, mpm_new_length,
          *NN_chain, NN_chain_length;
  double   mul0, mul1, **MI_matrix;
  FOS *new_FOS;

  /* Compute Mutual Information matrix */
  MI_matrix = NULL;
  if( learn_linkage_tree && !dependency_learning)
    MI_matrix = computeMIMatrix( covariance_matrix, number_of_parameters );


//    if(!dependency_learning)
//        printMIMatrix(MI_matrix,number_of_parameters, number_of_parameters);
  /* Initialize MPM to the univariate factorization */
  order                 = randomPermutation( number_of_parameters );
  mpm                   = (int **) Malloc( number_of_parameters*sizeof( int * ) );
  mpm_number_of_indices = (int *) Malloc( number_of_parameters*sizeof( int ) );
  mpm_length            = number_of_parameters;
  mpm_new               = NULL;
  for( i = 0; i < number_of_parameters; i++ )
  {
    indices                  = (int *) Malloc( 1*sizeof( int ) );
    indices[0]               = order[i];
    mpm[i]                   = indices;
    mpm_number_of_indices[i] = 1;
  }
  free( order );

  /* Initialize LT to the initial MPM */
  new_FOS                     = (FOS*) Malloc(sizeof(FOS));
  new_FOS->length             = number_of_parameters+number_of_parameters-1;
  new_FOS->sets               = (int **) Malloc( new_FOS->length*sizeof( int * ) );
  new_FOS->set_length         = (int *) Malloc( new_FOS->length*sizeof( int ) );
  FOS_index                                   = 0;
  for( i = 0; i < mpm_length; i++ )
  {
    new_FOS->sets[FOS_index]       = mpm[i];
    new_FOS->set_length[FOS_index] = mpm_number_of_indices[i];
    FOS_index++;
  }

  /* Initialize similarity matrix */
  S_matrix = NULL;
  if( !random_linkage_tree ){
    S_matrix = (double **) Malloc( number_of_parameters*sizeof( double * ) );
    for( i = 0; i < number_of_parameters; i++ )
      S_matrix[i] = (double *) Malloc( number_of_parameters*sizeof( double ) );
  }

  if( learn_linkage_tree )
  {
    if ( dependency_learning ) {
      for (i = 0; i < mpm_length; i++)
        for (j = 0; j < mpm_length; j++)
          S_matrix[i][j] = dependency_matrix[mpm[i][0]][mpm[j][0]];
      for (i = 0; i < mpm_length; i++)
        S_matrix[i][i] = 0;
    }
    else{
      for( i = 0; i < mpm_length; i++ )
        for( j = 0; j < mpm_length; j++ )
          S_matrix[i][j] = MI_matrix[mpm[i][0]][mpm[j][0]];
      for( i = 0; i < mpm_length; i++ )
        S_matrix[i][i] = 0;

      for( i = 0; i < number_of_parameters; i++ )
        free( MI_matrix[i] );
      free( MI_matrix );
    }

  }
  else if( random_linkage_tree )
  {
    S_vector = (double *) Malloc( number_of_parameters*sizeof(double));
    for( i = 0; i < number_of_parameters; i++ )
      S_vector[i] = randomRealUniform01();
  }
  else if( static_linkage_tree )
  {
    if( problem_index == 105 || problem_index == 106 )
    {
      for( i = 0; i < number_of_parameters-1; i++ )
      {
        for( j = i+1; j < number_of_parameters; j++ )
        {
          S_matrix[i][j] = 1.0 / covariance_matrix[mpm[i][0]][mpm[j][0]];
          S_matrix[j][i] = S_matrix[i][j];
        }
        S_matrix[i][i] = 0.0;
      }
    }
    else if (dependency_learning){
      for (i = 0; i < mpm_length; i++)
        for (j = 0; j < mpm_length; j++)
          S_matrix[i][j] = dependency_matrix[mpm[i][0]][mpm[j][0]];
      for (i = 0; i < mpm_length; i++)
        S_matrix[i][i] = 0;
    }
    else
    {
      for( i = 0; i < mpm_length; i++ )
      {
        for( j = 0; j < i; j++ )
        {
          if( mpm[i][0] < block_start || mpm[j][0] < block_start ) S_matrix[i][j] = randomRealUniform01();
          else if( (mpm[i][0]-block_start)/block_size == (mpm[j][0]-block_start)/block_size ) S_matrix[i][j] = randomRealUniform01() + 1e8;
          else S_matrix[i][j] = randomRealUniform01() + 1e3;
          S_matrix[j][i] = S_matrix[i][j];
        }
        S_matrix[i][i] = 0;
      }
    }
  }

  int *keep_FOS_element = (int *) Malloc( ((number_of_parameters*2))*sizeof( int ) );
  for( i = 0; i < number_of_parameters*2; i++ ){
    keep_FOS_element[i] = 1;
  }

  NN_chain        = (int *) Malloc( (number_of_parameters+2)*sizeof( int ) );
  NN_chain_length = 0;
  done            = 0;
  while( !done )
  {
    if( NN_chain_length == 0 )
    {
      NN_chain[NN_chain_length] = randomInt( mpm_length );
      NN_chain_length++;
    }

    if( NN_chain[NN_chain_length-1] >= mpm_length ) NN_chain[NN_chain_length-1] = mpm_length-1;

    while( NN_chain_length < 3 )
    {
      NN_chain[NN_chain_length] = determineNearestNeighbour( NN_chain[NN_chain_length-1], S_matrix, mpm_number_of_indices, mpm_length );
      NN_chain_length++;
    }

    while( NN_chain[NN_chain_length-3] != NN_chain[NN_chain_length-1] )
    {
      NN_chain[NN_chain_length] = determineNearestNeighbour( NN_chain[NN_chain_length-1], S_matrix, mpm_number_of_indices, mpm_length );
      if( ((getSimilarity(NN_chain[NN_chain_length-1],NN_chain[NN_chain_length]) == getSimilarity(NN_chain[NN_chain_length-1],NN_chain[NN_chain_length-2])))
          && (NN_chain[NN_chain_length] != NN_chain[NN_chain_length-2]) )
        NN_chain[NN_chain_length] = NN_chain[NN_chain_length-2];
      NN_chain_length++;
      if( NN_chain_length > number_of_parameters ){
        break;
      }
    }
    r0 = NN_chain[NN_chain_length-2];
    r1 = NN_chain[NN_chain_length-1];

    if( r1 >= mpm_length || r0 >= mpm_length || mpm_number_of_indices[r0]+mpm_number_of_indices[r1] > FOS_element_ub )
    {
      NN_chain_length = 1;
      NN_chain[0] = 0;
      if( FOS_element_ub < number_of_parameters )
      {
        done = 1;
        for( i = 1; i < mpm_length; i++ )
        {
          if( mpm_number_of_indices[i] + mpm_number_of_indices[NN_chain[0]] <= FOS_element_ub ) done = 0;
          if( mpm_number_of_indices[i] < mpm_number_of_indices[NN_chain[0]] ) NN_chain[0] = i;
        }
        if( done ) break;
      }
      continue;
    }

    if( r0 > r1 )
    {
      rswap = r0;
      r0    = r1;
      r1    = rswap;
    }
    NN_chain_length -= 3;



    if( r1 < mpm_length && r1 != r0 ) /* This test is required for exceptional cases in which the nearest-neighbor ordering has changed within the chain while merging within that chain */
    {
      indices = (int *) Malloc( (mpm_number_of_indices[r0]+mpm_number_of_indices[r1])*sizeof( int ) );

      i = 0;
      for( j = 0; j < mpm_number_of_indices[r0]; j++ )
      {
        indices[i] = mpm[r0][j];
        i++;
      }
      for( j = 0; j < mpm_number_of_indices[r1]; j++ )
      {
        indices[i] = mpm[r1][j];
        i++;
      }

      new_FOS->sets[FOS_index] = (int *) Malloc( (mpm_number_of_indices[r0]+mpm_number_of_indices[r1])*sizeof( int ) );
      new_FOS->set_length[FOS_index] = mpm_number_of_indices[r0]+mpm_number_of_indices[r1];
      sorted = mergeSortInt(indices, mpm_number_of_indices[r0]+mpm_number_of_indices[r1]);
      for( j = 0; j < mpm_number_of_indices[r0]+mpm_number_of_indices[r1]; j++ ){
        new_FOS->sets[FOS_index][j] = indices[sorted[j]];
      }

//            printf("similar: %f, epsilon: %f \n",getSimilarity(r0, r1), epsilon);
      if(getSimilarity(r0, r1) <= epsilon){
        keep_FOS_element[FOS_index] = 0;
//                printf("NotKeeping: %d \n",FOS_index);
      }
      if(FOS_index == (number_of_parameters*2)-2){
//                printf("dependency between biggest sets: %f\n", getSimilarity(r0, r1));
      }
      if( mpm_number_of_indices[r0]+mpm_number_of_indices[r1] > pruning_ub && keep_FOS_element[FOS_index] ){
        keep_FOS_element[FOS_index] = 0;
      }
      if( pruned_tree && keep_FOS_element[FOS_index] ){
        // we know we will merge r0 and r1, now lets check if they are all completely dependent
        int completely_dependent = 1;
        int all_checked = 1;
        for (i = 0; i < mpm_number_of_indices[r0]; i++){
          for (j = 0; j< mpm_number_of_indices[r1]; j++){
            if (dependency_matrix[mpm[r0][i]][mpm[r1][j]] <= 0.0){
              if(checked_matrix[mpm[r0][i]][mpm[r1][j]]==0){
                all_checked = 0;
              }
              completely_dependent = 0;
              break;
            }
          }
        }
        if ( completely_dependent ) { // remove subsets that build this set
          //remove r1
          int first_set_element = mpm[r0][0];
          int set_length = mpm_number_of_indices[r0];
          for (i = 0; i < FOS_index; i++) {
            if (new_FOS->set_length[i] == set_length && new_FOS->sets[i][0] == first_set_element) {
              keep_FOS_element[i] = 0;
            }
          }
          //remove r0
          first_set_element = mpm[r1][0];
          set_length = mpm_number_of_indices[r1];
          for (i = 0; i < FOS_index; i++) {
            if (new_FOS->set_length[i] == set_length && new_FOS->sets[i][0] == first_set_element) {
              keep_FOS_element[i] = 0;
            }
          }
        }
        else if(sparse_tree){ // dont add current set since it is not completely dependent
          if(!wait_with_pruning || all_checked){
            keep_FOS_element[FOS_index] = 0;
          }
        }
      }


      free( sorted );
      free( indices );

      mul0 = ((double) mpm_number_of_indices[r0])/((double) mpm_number_of_indices[r0]+mpm_number_of_indices[r1]);
      mul1 = ((double) mpm_number_of_indices[r1])/((double) mpm_number_of_indices[r0]+mpm_number_of_indices[r1]);
      if( random_linkage_tree )
      {
        S_vector[r0] = mul0*S_vector[r0]+mul1*S_vector[r1];
      }
      else
      {
        for( i = 0; i < mpm_length; i++ )
        {
          if( (i != r0) && (i != r1) )
          {
            S_matrix[i][r0] = mul0*S_matrix[i][r0] + mul1*S_matrix[i][r1];
            S_matrix[r0][i] = S_matrix[i][r0];
          }
        }
      }

      mpm_new                   = (int **) Malloc( (mpm_length-1)*sizeof( int * ) );
      mpm_new_number_of_indices = (int *) Malloc( (mpm_length-1)*sizeof( int ) );
      mpm_new_length            = mpm_length-1;
      for( i = 0; i < mpm_new_length; i++ )
      {
        mpm_new[i]                   = mpm[i];
        mpm_new_number_of_indices[i] = mpm_number_of_indices[i];
      }

      mpm_new[r0]                   = new_FOS->sets[FOS_index];
      mpm_new_number_of_indices[r0] = mpm_number_of_indices[r0]+mpm_number_of_indices[r1];
      if( r1 < mpm_length-1 )
      {
        mpm_new[r1]                   = mpm[mpm_length-1];
        mpm_new_number_of_indices[r1] = mpm_number_of_indices[mpm_length-1];

        if( random_linkage_tree )
        {
          S_vector[r1] = S_vector[mpm_length-1];
        }
        else
        {
          for( i = 0; i < r1; i++ )
          {
            S_matrix[i][r1] = S_matrix[i][mpm_length-1];
            S_matrix[r1][i] = S_matrix[i][r1];
          }

          for( j = r1+1; j < mpm_new_length; j++ )
          {
            S_matrix[r1][j] = S_matrix[j][mpm_length-1];
            S_matrix[j][r1] = S_matrix[r1][j];
          }
        }
      }

      for( i = 0; i < NN_chain_length; i++ )
      {
        if( NN_chain[i] == mpm_length-1 )
        {
          NN_chain[i] = r1;
          break;
        }
      }

      free( mpm );
      free( mpm_number_of_indices );
      mpm                   = mpm_new;
      mpm_number_of_indices = mpm_new_number_of_indices;
      mpm_length            = mpm_new_length;

      if( mpm_length == 1 ){
        done = 1;
      }

      FOS_index++;
    }
  }
//    printf("original FOS\n");
//    for( i =0; i < FOS_index; i++){
//        if (keep_FOS_element[i]){
//            printf("keep: ");
//        }
//        else{
//            printf("remove: ");
//        }
//        int setlenght = new_FOS->set_length[i];
//        for(int j = 0; j < setlenght; j++ ){
//            printf("%d, ", new_FOS->sets[i][j]);
//        }
//        printf("\n");
//    }

  if(( pruned_tree || epsilon > 0.0 || pruning_ub < number_of_parameters )){  //number of parasm > 100
    i = 0;
    j = 0;
    int new_lenght = 0;
    while (i < FOS_index){
      if (keep_FOS_element[i]) {
        if (i > j){
          while(keep_FOS_element[j] && j < FOS_index )
            j ++;
          new_FOS->set_length[j] = new_FOS->set_length[i];
          free(new_FOS->sets[j]);
          new_FOS->sets[j] = (int *) Malloc( (new_FOS->set_length[j])*sizeof( int ) );
          for(int k = 0; k < new_FOS->set_length[j]; k++ ){
            new_FOS->sets[j][k] = new_FOS->sets[i][k];
          }
          keep_FOS_element[j] = 1;
          keep_FOS_element[i] = 0;
        }
        if(i == j){
          j++;
        }
        new_lenght += 1;
      }
      i ++;
    }
    for(i = new_lenght; i<FOS_index; i++){
      free(new_FOS->sets[i]);
    }
    FOS_index = new_lenght;
  }

  new_FOS->length = FOS_index;
//    int *present =  (int *) Malloc( ((number_of_parameters)*sizeof( int ) ));
//    for(int i =0;i<number_of_parameters;i++ )present[i] = 0;
//    for(int i =0; i < new_FOS->length; i++){
//        for(int j = 0; j< new_FOS->set_length[i];j++){
//            if(present[new_FOS->sets[i][j]]){
//                printf("NOT MARGINAL %d \n ",new_FOS->sets[i][j] );
//                for(int k = 0; k < number_of_parameters;k++){
//                    printf("keepfos %d: %d\n",k, keep_FOS_element[k]);
//                }
//                printFOS(new_FOS);
//            }
//            else{
//                present[new_FOS->sets[i][j]]=1;
//            }
//        }
//    }
//    printf("marignal! \n");
//    printf("making Tree\n");
//    printFOS(new_FOS);
//    printf("NEW FOS\n");
//    for( i =0; i < FOS_index; i++){
//        int setlenght = new_FOS->set_length[i];
//        for(int j = 0; j < setlenght; j++ ){
//            printf("%d, ", new_FOS->sets[i][j]);
//        }
//        printf("\n");
//    }

  free( NN_chain );

  free( mpm_new );
  free( mpm_number_of_indices );
  free( keep_FOS_element );

  if( random_linkage_tree )
    free( S_vector );
  else
  {
    for( i = 0; i < number_of_parameters; i++ )
      free( S_matrix[i] );
    free( S_matrix );
  }

  return( new_FOS );
}



FOS *learnLinkageTreeRVGOMEA( int population_index )
{
    FOS *new_FOS;
    if( evolve_learning ){
        evolveDifferentialDependencies( population_index );
    }
    new_FOS = learnLinkageTree( full_covariance_matrix[population_index], dependency_matrix );
    return( new_FOS );
}

void learnFOS(int population_index)
{
    if( current_waiting_position == 0 ){
		if(number_of_generations[population_index] != 0){
			ezilaitiniCovarianceMatrices(population_index);
			ezilaitiniFOS(population_index);
		}
		estimateFullCovarianceMatrixML( population_index );
		linkage_model[population_index] = learnLinkageTreeRVGOMEA( population_index );
		initializeCovarianceMatrices( population_index );

        if( number_of_generations[population_index] == 0 ) {
            initializeDistributionMultipliers( population_index );
        }
		//printFOS( linkage_model[population_index] );
    }
    else if (current_waiting_position > 0){
        current_waiting_position--;
    }
}


void evolveDifferentialDependencies( int population_index ) {
    int i, j, k;
    double *individual_to_compare = (double *) Malloc(number_of_parameters * sizeof(double));
    double constraint_value;

    // initialize if no pairs are checked yet
    if (number_of_checked_pairs == 0) {
        double rand = randomRealUniform01();
        rand = 0.7;

        for (k = 0; k < number_of_parameters; k++) {
            double min = lower_init_ranges[k], max = upper_init_ranges[k];
            getMinMaxofPopulation(k, population_index, &min, &max);
            if (nround(min, 2) == nround(max, 2)) {
                max = upper_init_ranges[k];
            }
            first_individual[k] = min + ((max - min) * rand * 0.5);
            double parameter_diff = (max - min) * 0.5 * rand;
            second_individual[k] = parameter_diff + first_individual[k];
            individual_to_compare[k] = first_individual[k];
        }

        double objective_value, old_constraint, old_objective;
        // fill evaluation storage
        installedProblemEvaluation(problem_index, first_individual, &(old_objective), &(old_constraint) ); 
        //installedProblemEvaluation(problem_index, first_individual, &(old_objective), &(old_constraint), number_of_parameters, NULL, NULL, 0, 0);
        differential_grouping_evals = 1+ number_of_parameters;
        fitness_of_first_individual[number_of_parameters] = old_objective;
        fitness_of_first_individual[0] = old_objective;
        for (k = 0; k < number_of_parameters; k++) {
            individual_to_compare[k] = second_individual[k];
            installedProblemEvaluation(problem_index, individual_to_compare, &(objective_value), &(constraint_value) ); 
            //installedProblemEvaluation(problem_index, individual_to_compare, &(objective_value), &(constraint_value), 1, &(k), &(first_individual[k]), old_objective, old_constraint);

            fitness_of_first_individual[k] = objective_value;
            individual_to_compare[k] = first_individual[k];
        }
        int counter = number_of_pairs;
        for (int i = counter - 1; i >= 0; --i) {
            int j = randomInt(i+1);

            //swap the last element with element at random index
            int *temp = dependency_pairs[i];
            dependency_pairs[i] = dependency_pairs[j];
            dependency_pairs[j] = temp;
        }

    } else {
        for (k = 0; k < number_of_parameters; k++) {
            individual_to_compare[k] = first_individual[k];
        }
    }

    iteration += 1;
    int max_index = number_of_checked_pairs + pairs_per_run;
    if (max_index >= number_of_pairs) {
        max_index = number_of_pairs;
    }

    double original_objective = fitness_of_first_individual[number_of_parameters];

    for (k = 0; k < number_of_parameters; k++) {
        individual_to_compare[k] = first_individual[k];
    }
    int found_dependencies = 0;
    double max_dependency = 0.0;
    for (k = number_of_checked_pairs; k < max_index; k++) {
        i = dependency_pairs[k][0];
        j = dependency_pairs[k][1];

        double change_i, change_j, change_i_j;
        change_i = fitness_of_first_individual[i];
        change_j = fitness_of_first_individual[j];

        individual_to_compare[i] = second_individual[i];
        individual_to_compare[j] = second_individual[j];
        installedProblemEvaluation(problem_index, individual_to_compare, &(change_i_j), &(constraint_value) );
        //installedProblemEvaluation(problem_index, individual_to_compare, &(change_i_j), &(constraint_value), 1, &(j), &(first_individual[j]), fitness_of_first_individual[i], 0);
        differential_grouping_evals+=1;
        individual_to_compare[i] = first_individual[i];
        individual_to_compare[j] = first_individual[j];

        double delta_i, delta_j;

        change_i = change_i/original_objective;
        change_j = change_j/original_objective;


        delta_i = fabs(1.0 - change_i);
        delta_j = fabs(change_j - change_i_j);


        delta_i = nround(delta_i, 12);
        delta_j = nround(delta_j, 12);

        double dependency = 0.0;
        double inverted_difference;

        if(delta_j == 0.0) {
            double temp = delta_i;
            delta_i = delta_j;
            delta_j = temp;
        }
        if(delta_j != 0.0){
            inverted_difference = fabs(delta_i/delta_j);
            if(inverted_difference > 1.0){
                inverted_difference = fabs((double)delta_j/delta_i);
            }
        } else{
            inverted_difference = 1.0;
        }
        dependency = 1-inverted_difference;
        if (inverted_difference < 1) {
            found_dependencies += 1;
        } else{
            dependency = 0.0;
        }
        dependency_matrix[i][j] = dependency;
        dependency_matrix[j][i] = dependency;

        max_dependency = fmax(dependency, max_dependency);
        checked_matrix[i][j] = 1;
        checked_matrix[j][i] = 1;
    }
    total_dependencies_found += found_dependencies;
    number_of_checked_pairs += pairs_per_run;
    if (found_dependencies == 0) {
        int found_dependencies_per_run = total_dependencies_found / iteration;
        if (found_dependencies_per_run < minimal_dependencies_per_run) {
            current_waiting_position = number_of_waiting_cycles;
            number_of_waiting_cycles *= 2;
            iteration = 0; total_dependencies_found = 0;
        }
    }
    if (number_of_checked_pairs >= number_of_pairs){
        number_of_checked_pairs = 0;
        current_waiting_position = number_of_waiting_cycles;
        number_of_waiting_cycles *= 2;
        iteration = 0; total_dependencies_found = 0;
    }

    free(individual_to_compare);
}

int determineNearestNeighbour( int index, double **S_matrix, int *mpm_number_of_indices, int mpm_length )
{
  int i, result;

  result = 0;
  if( result == index )
    result++;
  for( i = 1; i < mpm_length; i++ )
  {
    if( ((S_matrix[index][i] > S_matrix[index][result]) || ((S_matrix[index][i] == S_matrix[index][result]) && (mpm_number_of_indices[i] < mpm_number_of_indices[result]))) && (i != index) )
      result = i;
  }

  return( result );
}

/*-=-=-=-=-=-=-=-=-=-=-=-=-= Section Variation -==-=-=-=-=-=-=-=-=-=-=-=-=-=*/
/**
 * First estimates the parameters of a normal distribution in the
 * parameter space from the selected sets of solutions (a separate
 * normal distribution for each population). Then copies the single
 * best selected solutions to their respective populations. Finally
 * fills up each population, after the variances have been scaled,
 * by drawing new samples from the normal distributions and applying
 * AMS to several of these new solutions. Then, the fitness ranks
 * are recomputed. Finally, the distribution multipliers are adapted
 * according to the SDR-AVS mechanism.
 */
void makePopulation( int population_index )
{
    if( populations_terminated[population_index] )
        return;

    estimateParameters( population_index );

    copyBestSolutionsToPopulation( population_index );

    applyDistributionMultipliers( population_index );

    generateAndEvaluateNewSolutionsToFillPopulation( population_index );

    computeRanksForOnePopulation( population_index );

    adaptDistributionMultipliersForOnePopulation( population_index );
    
    ezilaitiniParametersForSampling( population_index );
}

/**
 * Estimates the parameters of the multivariate normal
 * distribution for each population separately.
 */
void estimateParametersAllPopulations( void )
{
  int i;

  for( i = 0; i < number_of_populations; i++ )
    if( !populations_terminated[i] )
      estimateParameters( i );
}

/**
 * Estimates the paramaters of the multivariate
 * normal distribution for a specified population.
 */
void estimateParameters( int population_index )
{
    estimateMeanVectorML( population_index );
    
	learnFOS( population_index );
        
    estimateParametersML( population_index );
}

/**
 * Estimates (with maximum likelihood) the
 * parameters of a multivariate normal distribution
 * for a specified population.
 */
void estimateParametersML( int population_index )
{
  estimateMeanVectorML( population_index );

  estimateCovarianceMatrixML( population_index );
}

/**
 * Computes the sample mean for a specified population.
 */
void estimateMeanVectorML( int population_index )
{
  int i, j;

  if( number_of_generations[population_index] > 0 )
  {
    for( i = 0; i < number_of_parameters; i++ )
      mean_vectors_previous[population_index][i] = mean_vectors[population_index][i];
  }

  for( i = 0; i < number_of_parameters; i++ )
  {
    mean_vectors[population_index][i] = 0.0;

    for( j = 0; j < selection_sizes[population_index]; j++ )
      mean_vectors[population_index][i] += selections[population_index][j][i];

    mean_vectors[population_index][i] /= (double) selection_sizes[population_index];
  }

  /* Change the focus of the search to the best solution */
  if( distribution_multipliers[population_index] < 1.0 )
  {
      for( j = 0; j < number_of_parameters; j++ )
      {
         mean_vectors[population_index][j] = selections[population_index][0][j];
      }
  }
}

/**
 * Computes the matrix of sample covariances for
 * a specified population.
 *
 * It is important that the pre-condition must be satisified:
 * estimateMeanVector was called first.
 */
 void estimateFullCovarianceMatrixML( int population_index )
{
    int i, j, m;
    double cov;

    /* First do the maximum-likelihood estimate from data */
    for( i = 0; i < number_of_parameters; i++ )
    {
        for( j = 0; j < number_of_parameters; j++ )
        {
            cov = 0.0;
            for( m = 0; m < selection_sizes[population_index]; m++ )
                cov += (selections[population_index][m][i]-mean_vectors[population_index][i])*(selections[population_index][m][j]-mean_vectors[population_index][j]);

            cov /= (double) selection_sizes[population_index];
            full_covariance_matrix[population_index][i][j] = cov;
            full_covariance_matrix[population_index][j][i] = cov;
        }
    }
}

void estimateCovarianceMatrixML( int population_index )
{
  int i, j, k, m, a, b;
  double cov;

  /* First do the maximum-likelihood estimate from data */
    for( i = 0; i < linkage_model[population_index]->length; i++ )
    {
        for( j = 0; j < linkage_model[population_index]->set_length[i]; j++ )
        {
            for( k = j; k < linkage_model[population_index]->set_length[i]; k++ )
			{
              covariance_matrices[population_index][i][j][k] = 0.0;

              a = linkage_model[population_index]->sets[i][j];
              b = linkage_model[population_index]->sets[i][k];

			  cov = 0.0;
			  for( m = 0; m < selection_sizes[population_index]; m++ )
				  cov += (selections[population_index][m][a]-mean_vectors[population_index][a])*(selections[population_index][m][b]-mean_vectors[population_index][b]);

			  cov /= (double) selection_sizes[population_index];

              covariance_matrices[population_index][i][j][k] = cov;
              covariance_matrices[population_index][i][k][j] = cov;
          }
      }
  }
}

/**
 * Copies the single very best of the selected solutions
 * to their respective populations.
 */
void copyBestSolutionsToPopulation( int population_index )
{
  int k;

  if( !populations_terminated[population_index] )
  {
    for( k = 0; k < number_of_parameters; k++ )
      populations[population_index][0][k] = selections[population_index][0][k];

    objective_values[population_index][0]  = objective_values_selections[population_index][0];
    constraint_values[population_index][0] = constraint_values_selections[population_index][0];
  }
}

/**
 * Applies the distribution multipliers.
 */
void applyDistributionMultipliers( int population_index )
{
  int j, k, m;


  if( !populations_terminated[population_index] )
  {
    for( j = 0; j < linkage_model[population_index]->length; j++ )
      for( k = 0; k < linkage_model[population_index]->set_length[j]; k++ )
          for( m = 0; m < linkage_model[population_index]->set_length[j]; m++ )
              covariance_matrices[population_index][j][k][m] *= distribution_multipliers[population_index];
  }
}

/**
 * Generates new solutions for each
 * of the populations in turn.
 */
void generateAndEvaluateNewSolutionsToFillPopulation( int population_index )
{
  short   out_of_range;
  int     j, k, q, number_of_AMS_solutions;
  double *solution, *solution_AMS, shrink_factor, alpha_AMS, delta_AMS;

  solution_AMS = (double *) Malloc( number_of_parameters*sizeof( double ) );
  alpha_AMS = 0.5*tau*(((double) population_sizes[population_index])/((double) (population_sizes[population_index]-1)));
  delta_AMS = 2;

  computeParametersForSampling( population_index );

  if( !populations_terminated[population_index] )
  {
    number_of_AMS_solutions                     = (int) (alpha_AMS*(population_sizes[population_index]-1));
    samples_drawn_from_normal[population_index] = 0;
    out_of_bounds_draws[population_index]       = 0;
    q                                           = 0;
    for( j = 1; j < population_sizes[population_index]; j++ )
    {
      solution = generateNewSolution( population_index );

      for( k = 0; k < number_of_parameters; k++ )
        populations[population_index][j][k] = solution[k];

      if( (number_of_generations[population_index] > 0) && (q < number_of_AMS_solutions) )
      {
        out_of_range  = 1;
        shrink_factor = 2;
        while( (out_of_range == 1) && (shrink_factor > 1e-10) )
        {
          shrink_factor *= 0.5;
          out_of_range   = 0;
          for( k = 0; k < number_of_parameters; k++ )
          {
            solution_AMS[k] = solution[k] + shrink_factor*delta_AMS*distribution_multipliers[population_index]*(mean_vectors[population_index][k]-mean_vectors_previous[population_index][k]);
            if( !isParameterInRangeBounds( solution_AMS[k], k ) )
            {
              out_of_range = 1;
              break;
            }
          }
        }
        if( !out_of_range )
        {
          for( k = 0; k < number_of_parameters; k++ )
            populations[population_index][j][k] = solution_AMS[k];
        }
      }

      installedProblemEvaluation( problem_index, populations[population_index][j], &(objective_values[population_index][j]), &(constraint_values[population_index][j]) );

        if( betterFitness(objective_values[population_index][j], constraint_values[population_index][j], population_best_obj_val[population_index], population_best_constraint_val[population_index]))
        {
            population_best_obj_val[population_index] = objective_values[population_index][j];
            population_best_constraint_val[population_index] = constraint_values[population_index][j];
        }

      q++;

      free( solution );
    }
  }

  free( solution_AMS );
}

/**
 * Computes the Cholesky decompositions required for sampling
 * the multivariate normal distribution.
 */
void computeParametersForSampling( int population_index )
{
  int i;
  
  if( !use_univariate_FOS )
  {
    cholesky_factors_lower_triangle[population_index] = (double ***) Malloc( linkage_model[population_index]->length*sizeof( double ** ) );
    for( i = 0; i < linkage_model[population_index]->length; i++ )
    {
      cholesky_factors_lower_triangle[population_index][i] = choleskyDecomposition( covariance_matrices[population_index][i], linkage_model[population_index]->set_length[i] );
    }
  }
}

/**
 * Generates and returns a single new solution by drawing
 * a single sample from a specified model.
 */
double *generateNewSolution( int population_index )
{
  short   ready;
  int     i, j, times_not_in_bounds;
  double *result, *FOS_result, *z;

  times_not_in_bounds = -1;
  out_of_bounds_draws[population_index]--;

  ready = 0;
  do
  {
    times_not_in_bounds++;
    samples_drawn_from_normal[population_index]++;
    out_of_bounds_draws[population_index]++;
    if( times_not_in_bounds >= 100 )
    {
      result = (double *) Malloc( number_of_parameters*sizeof( double ) );
      for( i = 0; i < number_of_parameters; i++ )
        result[i] = lower_init_ranges[i] + (upper_init_ranges[i] - lower_init_ranges[i])*randomRealUniform01();
    }
    else
    {
      result = (double *) Malloc( number_of_parameters*sizeof( double ) );
      for( i = 0; i < linkage_model[population_index]->length; i++ )
      {
          z = (double *) Malloc( linkage_model[population_index]->set_length[i]*sizeof( double ) );

          for( j = 0; j < linkage_model[population_index]->set_length[i]; j++ )
             z[j] = random1DNormalUnit();
          
          if( use_univariate_FOS )
          {
              result[i] = z[0]*sqrt(covariance_matrices[population_index][i][0][0]) + mean_vectors[population_index][i];
          }
          else
          {
              FOS_result = matrixVectorMultiplication( cholesky_factors_lower_triangle[population_index][i], z, linkage_model[population_index]->set_length[i], linkage_model[population_index]->set_length[i] );
              for( j = 0; j < linkage_model[population_index]->set_length[i]; j++ )
                 result[linkage_model[population_index]->sets[i][j]] = FOS_result[j] + mean_vectors[population_index][linkage_model[population_index]->sets[i][j]];
              free( FOS_result );
          }
          
          free( z );
      }
    }

    ready = 1;
    for( i = 0; i < number_of_parameters; i++ )
    {
      if( !isParameterInRangeBounds( result[i], i ) )
      {
        ready = 0;
        break;
      }
    }
    if( !ready )
      free( result );
  }
  while( !ready );

  return( result );
}

/**
 * Adapts the distribution multipliers according to
 * the SDR-AVS mechanism.
 */
void adaptDistributionMultipliersForOnePopulation( int population_index )
{
  short  improvement;
  int    i;
  double st_dev_ratio;

  i = population_index;
    if( !populations_terminated[i] )
    {
      if( (((double) out_of_bounds_draws[i])/((double) samples_drawn_from_normal[i])) > 0.9 )
        distribution_multipliers[i] *= 0.5;
  
      improvement = generationalImprovementForOnePopulation( i, &st_dev_ratio );
  
      if( improvement )
      {
        no_improvement_stretch[i] = 0;

        if( distribution_multipliers[i] < 1.0 )
          distribution_multipliers[i] = 1.0;
  
        if( st_dev_ratio > st_dev_ratio_threshold )
          distribution_multipliers[i] *= distribution_multiplier_increase;
      }
      else
      {
        if( distribution_multipliers[i] <= 1.0 )
          (no_improvement_stretch[i])++;
  
        if( (distribution_multipliers[i] > 1.0) || (no_improvement_stretch[i] >= maximum_no_improvement_stretch) )
          distribution_multipliers[i] *= distribution_multiplier_decrease;
  
        if( (no_improvement_stretch[i] < maximum_no_improvement_stretch) && (distribution_multipliers[i] < 1.0) )
          distribution_multipliers[i] = 1.0;
      }
    }
}

/**
 * Determines whether an improvement is found for a specified
 * population. Returns 1 in case of an improvement, 0 otherwise.
 * The standard-deviation ratio required by the SDR-AVS
 * mechanism is computed and returned in the pointer variable.
 */
short generationalImprovementForOnePopulation( int population_index, double *st_dev_ratio )
{
  int     i, j, index_best_selected, index_best_population,
          number_of_improvements;
  double *average_parameters_of_improvements;

  /* Determine best selected solutions */
  index_best_selected = 0;
  for( i = 0; i < selection_sizes[population_index]; i++ )
  {
    if( betterFitness( objective_values_selections[population_index][i], constraint_values_selections[population_index][i],
                       objective_values_selections[population_index][index_best_selected], constraint_values_selections[population_index][index_best_selected] ) )
      index_best_selected = i;
  }

  /* Determine best in the population and the average improvement parameters */
  average_parameters_of_improvements = (double *) Malloc( number_of_parameters*sizeof( double ) );
  for( i = 0; i < number_of_parameters; i++ )
    average_parameters_of_improvements[i] = 0.0;

  index_best_population   = 0;
  number_of_improvements  = 0;
  for( i = 0; i < population_sizes[population_index]; i++ )
  {
    if( betterFitness( objective_values[population_index][i], constraint_values[population_index][i],
                       objective_values[population_index][index_best_population], constraint_values[population_index][index_best_population] ) )
      index_best_population = i;

    if( betterFitness( objective_values[population_index][i], constraint_values[population_index][i],
                       objective_values_selections[population_index][index_best_selected], constraint_values_selections[population_index][index_best_selected] ) )
    {
      number_of_improvements++;
      for( j = 0; j < number_of_parameters; j++ )
        average_parameters_of_improvements[j] += populations[population_index][i][j];
    }
  }

  /* Determine st.dev. ratio */
  *st_dev_ratio = 0.0;
  if( number_of_improvements > 0 )
  {
    for( i = 0; i < number_of_parameters; i++ )
      average_parameters_of_improvements[i] /= (double) number_of_improvements;

    *st_dev_ratio = getStDevRatio( population_index, average_parameters_of_improvements );
  }

  free( average_parameters_of_improvements );

  if( fabs( objective_values_selections[population_index][index_best_selected] - objective_values[population_index][index_best_population] ) == 0.0 )
    return( 0 );

  return( 1 );
}

/**
 * Returns 1 if x is better than y, 0 otherwise.
 * x is not better than y unless:
 * - x and y are both infeasible and x has a smaller sum of constraint violations, or
 * - x is feasible and y is not, or
 * - x and y are both feasible and x has a smaller objective value than y
 */
short betterFitness( double objective_value_x, double constraint_value_x, double objective_value_y, double constraint_value_y )
{
  short result;

  result = 0;

  if( constraint_value_x > 0 ) /* x is infeasible */
  {
    if( constraint_value_y > 0 ) /* Both are infeasible */
    {
      if( constraint_value_x < constraint_value_y )
       result = 1;
    }
  }
  else /* x is feasible */
  {
    if( constraint_value_y > 0 ) /* x is feasible and y is not */
      result = 1;
    else /* Both are feasible */
    {
      if( objective_value_x < objective_value_y )
        result = 1;
    }
  }

  return( result );
}

/**
 * Computes and returns the standard-deviation-ratio
 * of a given point for a given model.
 */
double getStDevRatio( int population_index, double *parameters )
{
  int      i, j;
  double **inverse, result, *x_min_mu, *z, z_uni;

  result = 0.0;
  
  if( use_univariate_FOS )
  {
    for( i = 0; i < number_of_parameters; i++ )
    {
        z_uni = (parameters[i]-mean_vectors[population_index][i])/sqrt(covariance_matrices[population_index][i][0][0]);
        if( fabs( z_uni ) > result )
            result = fabs( z_uni );
    }     
  }
  else
  {
      for( i = 0; i < linkage_model[population_index]->length; i++ )
      {
          inverse = matrixLowerTriangularInverse( cholesky_factors_lower_triangle[population_index][i], linkage_model[population_index]->set_length[i] );

          x_min_mu = (double *) Malloc( linkage_model[population_index]->set_length[i]*sizeof( double ) );

          for( j = 0; j < linkage_model[population_index]->set_length[i]; j++ )
            x_min_mu[j] = parameters[linkage_model[population_index]->sets[i][j]]-mean_vectors[population_index][linkage_model[population_index]->sets[i][j]];

          z = matrixVectorMultiplication( inverse, x_min_mu, linkage_model[population_index]->set_length[i], linkage_model[population_index]->set_length[i] );

          for( j = 0; j < linkage_model[population_index]->set_length[i]; j++ )
              if( fabs( z[j] ) > result )
                  result = fabs( z[j] );

          free( z );
          for( j = 0; j < linkage_model[population_index]->set_length[i]; j++ )
              free( inverse[j] );
          free( inverse );
          free( x_min_mu );
      }
  }

  return( result );
}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/





/*-=-=-=-=-=-=-=-=-=-=-=-=- Section Ezilaitini -=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
/**
 * Undoes initialization procedure by freeing up memory.
 */
void ezilaitini( void )
{
    ezilaitiniObjectiveRotationMatrix();

    ezilaitiniMemory();
}

void ezilaitiniMemoryOnePopulation( int population_index )
{
    int j;

    for( j = 0; j < population_sizes[population_index]; j++ )
      free( populations[population_index][j] );
    free( populations[population_index] );

    free( objective_values[population_index] );

    free( constraint_values[population_index] );

    free( ranks[population_index] );

    for( j = 0; j < selection_sizes[population_index]; j++ )
      free( selections[population_index][j] );
    free( selections[population_index] );

    free( objective_values_selections[population_index] );

    free( constraint_values_selections[population_index] );

    free( mean_vectors[population_index] );

    free( mean_vectors_previous[population_index] );
    
	for( j = 0; j < number_of_parameters; j++ )
		free( full_covariance_matrix[population_index][j]);
	free( full_covariance_matrix[population_index] );
}

/**
 * Undoes initialization procedure by freeing up memory.
 */
void ezilaitiniMemory( void )
{
  int i;

  for( i = 0; i < number_of_populations; i++ )
      ezilaitiniMemoryOnePopulation( i );

  free( linkage_model );
  free( full_covariance_matrix );
  free( covariance_matrices );
  free( cholesky_factors_lower_triangle );
  free( lower_range_bounds );
  free( upper_range_bounds );
  free( lower_init_ranges );
  free( upper_init_ranges );
  free( populations_terminated );
  free( no_improvement_stretch );
  free( populations );
  free( population_best_obj_val );
  free( population_best_constraint_val );
  free( objective_values );
  free( constraint_values );
  free( ranks );
  free( selections );
  free( objective_values_selections );
  free( constraint_values_selections );
  free( mean_vectors );
  free( mean_vectors_previous );
  free( population_sizes );
  free( selection_sizes );
  free( number_of_generations );
  free( distribution_multipliers );
  free( samples_drawn_from_normal );
  free( out_of_bounds_draws );
}

/**
 * Undoes initialization procedure by freeing up memory.
 */
void ezilaitiniObjectiveRotationMatrix( void )
{
  int i;

  if( rotation_angle == 0.0 )
    return;

  for( i = 0; i < ellipsoid_size; i++ )
    free( rotation_matrix[i] );
  free( rotation_matrix );
}

void ezilaitiniCovarianceMatrices( int population_index )
{
    int i,j,k;
    
    i = population_index;
    for( j = 0; j < linkage_model[i]->length; j++ )
    {
        for( k = 0; k < linkage_model[i]->set_length[j]; k++ )
            free( covariance_matrices[i][j][k] );
        free( covariance_matrices[i][j] );
    }
    free( covariance_matrices[i] );
}

void ezilaitiniFOS( int population_index )
{
    int i;
    
    for( i = 0; i < linkage_model[population_index]->length; i++ )
    {
        free( linkage_model[population_index]->sets[i] );
    }
    free( linkage_model[population_index]->sets );
    free( linkage_model[population_index]->set_length );
}

void ezilaitiniParametersForSampling( int population_index )
{
    int i, j;

    if( !use_univariate_FOS )
    {
        for( i = 0; i < linkage_model[population_index]->length; i++ )
        {
            for( j = 0; j < linkage_model[population_index]->set_length[i]; j++ )
                free( cholesky_factors_lower_triangle[population_index][i][j] );
            free( cholesky_factors_lower_triangle[population_index][i] );
        }
        free( cholesky_factors_lower_triangle[population_index] );
    }
}

/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/





/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=- Section Run -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/

void generationalStepAllPopulationsRecursiveFold(int population_index_smallest, int population_index_biggest );
void generationalStepAllPopulationsHarikLobo( void )
{
  int population_index_smallest, population_index_biggest;

  population_index_biggest  = number_of_populations-1;
  population_index_smallest = 0;
  while( population_index_smallest <= population_index_biggest )
  {
    if( !populations_terminated[population_index_smallest] )
      break;

    population_index_smallest++;
  }

  generationalStepAllPopulationsRecursiveFold( population_index_smallest, population_index_biggest );
}

void generationalStepAllPopulationsRecursiveFold( int population_index_smallest, int population_index_biggest )
{
  int i, j, population_index;

  for( i = 0; i < number_of_subgenerations_per_population_factor-1; i++ )
  {
    for( population_index = population_index_smallest; population_index <= population_index_biggest; population_index++ )
    {
      //printf("population[%d] - population size: %d\n", population_index, population_sizes[population_index]);
      if( !populations_terminated[population_index] )
      {
          makeSelectionsForOnePopulation( population_index );

          makePopulation( population_index );

          /*if( write_generational_statistics && evaluations_for_statistics_hit )
          {
            writeGenerationalStatisticsForOnePopulation( population_index );
            evaluations_for_statistics_hit = 0;
          }*/

          number_of_generations[population_index]++;

          if( checkImmediateTerminationConditions() )
          {
              for( j = 0; j < number_of_populations; j++ )
                  populations_terminated[j] = 1;
              return;
          }
      }
    }

    for( population_index = population_index_smallest; population_index < population_index_biggest; population_index++ )
      generationalStepAllPopulationsRecursiveFold( population_index_smallest, population_index );
  }
}

void runAllPopulationsHarikLobo( void )
{
    while( !checkTerminationCondition() )
    {
        if( number_of_populations < maximum_number_of_populations )
            initializeNewPopulationHarikLobo( number_of_populations );

        if( write_generational_statistics )
          writeGenerationalStatisticsForOnePopulation( number_of_populations-1 );

        if( write_generational_solutions )
            writeGenerationalSolutions( 0 );

        generationalStepAllPopulationsHarikLobo();

        total_number_of_generations++;
    }
}

/**
 * Runs the IDEA.
 */
void run( void )
{
  int population_of_best, index_of_best;

  initialize();

  if( print_verbose_overview )
    printVerboseOverview();

  runAllPopulationsHarikLobo();

  if( write_generational_solutions )
      writeGenerationalSolutions( 1 );

  if( write_generational_statistics )
    writeGenerationalStatisticsForOnePopulation( number_of_populations-1 );

  printf("evals %d ", number_of_evaluations);

  determineBestSolutionInCurrentPopulations( &population_of_best, &index_of_best );
  printf("obj_val %6.2e ", objective_values[population_of_best][index_of_best]);

  ezilaitini();

  printf("time %f ", getTimer());
  printf("generations %d\n", total_number_of_generations);

}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/





/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=- Section Main -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
/**
 * The main function:
 * - interpret parameters on the command line
 * - run the algorithm with the interpreted parameters
 */
int main( int argc, char **argv )
{
  interpretCommandLine( argc, argv );

  run();

  return( 0 );
}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
