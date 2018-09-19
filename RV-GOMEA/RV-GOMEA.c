/**
 *
 * RV-GOMEA
 *
 * If you use this software for any purpose, please cite the most recent publication:
 * A. Bouter, C. Witteveen, T. Alderliesten, P.A.N. Bosman. 2017.
 * Exploiting Linkage Information in Real-Valued Optimization with the Real-Valued
 * Gene-pool Optimal Mixing Evolutionary Algorithm. In Proceedings of the Genetic 
 * and Evolutionary Computation Conference (GECCO 2017).
 * DOI: 10.1145/3071178.3071272
 *
 * Copyright (c) 1998-2017 Peter A.N. Bosman
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
 * The software in this file is the result of (ongoing) scientific research.
 * The following people have been actively involved in this research over
 * the years:
 * - Peter A.N. Bosman
 * - Dirk Thierens
 * - Jörn Grahl
 * - Anton Bouter
 * 
 */

/*-=-=-=-=-=-=-=-=-=-=-=-=-=-= Section Includes -=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <sys/time.h>
#include "../util/Tools.h"
#include "../util/SO_optimization.h"
#include "../util/Optimization.h"
#include "../util/FOS.h"
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/


/*-=-=-=-=-=-=-=-=-=-=-=-= Section Header Functions -=-=-=-=-=-=-=-=-=-=-=-=*/
void interpretCommandLine( int argc, char **argv );
void run( void );
int *mergeSortFitness( double *objectives, double *constraints, int number_of_solutions );
void mergeSortFitnessWithinBounds( double *objectives, double *constraints, int *sorted, int *tosort, int p, int q );
void mergeSortFitnessMerge( double *objectives, double *constraints, int *sorted, int *tosort, int p, int r, int q );
void parseCommandLine( int argc, char **argv );
void parseOptions( int argc, char **argv, int *index );
void parseFOSElementSize( int *index, int argc, char** argv );
void optionError( char **argv, int index );
void parseParameters( int argc, char **argv, int *index );
void printUsage( void );
void checkOptions( void );
void printVerboseOverview( void );
void initialize( void );
void initializeMemory( void );
void initializeNewPopulation( void );
void initializeNewPopulationMemory( int population_index );
void initializeFOS( int population_index );
void initializeProblem( void );
void initializeDistributionMultipliers(int population_index );
void initializePopulationAndFitnessValues( int population_index );
void inheritDistributionMultipliers( FOS *new_FOS, FOS *prev_FOS, double *multipliers );
FOS *learnLinkageTreeRVGOMEA( int population_index );
void computeRanksForAllPopulations( void );
void computeRanksForOnePopulation( int population_index );
void writeGenerationalStatisticsForOnePopulation( int population_index );
void writeGenerationalSolutions( short final );
void writeGenerationalSolutionsBest( short final );
short checkTerminationCondition( void );
short checkSubgenerationTerminationConditions( void );
short checkTimeLimitTerminationCondition( void );
short checkNumberOfEvaluationsTerminationCondition( void );
short checkVTRTerminationCondition( void );
void checkAverageFitnessTerminationCondition( void );
void determineBestSolutionInCurrentPopulations( int *population_of_best, int *index_of_best );
void checkFitnessVarianceTermination( void );
short checkFitnessVarianceTerminationSinglePopulation( int population_index );
void checkDistributionMultiplierTerminationCondition( void );
void makeSelections( void );
void makeSelectionsForOnePopulation( int population_index );
void makeSelectionsForOnePopulationUsingDiversityOnRank0( int population_index );
void estimateParameters( int population_index );
void estimateMeanVectorML( int population_index );
void estimateFullCovarianceMatrixML( int population_index );
void estimateParametersML( int population_index );
void estimateCovarianceMatricesML( int population_index );
void initializeCovarianceMatrices( int population_index );
void copyBestSolutionsToAllPopulations( void );
void copyBestSolutionsToPopulation( int population_index );
void getBestInPopulation( int population_index, int *individual_index );
void getOverallBest( int *population_index, int *individual_index );
void evaluateCompletePopulation( int population_index );
void applyDistributionMultipliersToAllPopulations( void );
void applyDistributionMultipliers(int population_index );
void generateAndEvaluateNewSolutionsToFillAllPopulations( void );
void generateAndEvaluateNewSolutionsToFillPopulation( int population_index );
void computeParametersForSampling(int population_index);
short generateNewSolutionFromFOSElement(int population_index, int FOS_index, int individual_index, short apply_AMS );
short applyAMS( int population_index, int individual_index );
void applyForcedImprovements(int population_index, int individual_index, int donor_index );
double *generateNewPartialSolutionFromFOSElement( int population_index, int FOS_index );
short adaptDistributionMultipliers(int population_index , int FOS_index);
short generationalImprovementForOnePopulationForFOSElement(int population_index, int FOS_index, double *st_dev_ratio );
double getStDevRatioForFOSElement( int population_index, double *parameters, int FOS_index );
void ezilaitini( void );
void ezilaitiniMemory( void );
void ezilaitiniDistributionMultipliers( int population_index );
void ezilaitiniParametersForSampling( int population_index );
void ezilaitiniParametersAllPopulations( void );
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/



/*-=-=-=-=-=-=-=-=-=-=-=- Section Global Variables -=-=-=-=-=-=-=-=-=-=-=-=-*/
short  print_verbose_overview,                              /* Whether to print a overview of settings (0 = no). */
       use_guidelines,                                      /* Whether to override parameters with guidelines (for those that exist). */
       fix_seed;                                            /* Whether a fixed seed is used. */
int    base_population_size,                                /* The size of the first population in the multi-start scheme. */
      *selection_sizes,                                     /* The size of the selection for each population. */
       total_number_of_writes,                              /* Total number of times a statistics file has been written. */
       maximum_number_of_populations,                       /* The maximum number of populations in the multi-start scheme. */
       number_of_subgenerations_per_population_factor,      /* The subgeneration factor in the multi-start scheme. */
     **samples_drawn_from_normal,                           /* The number of samples drawn from the i-th normal in the last generation. */
     **out_of_bounds_draws,                                 /* The number of draws that resulted in an out-of-bounds sample. */
      *no_improvement_stretch,                              /* The number of subsequent generations without an improvement while the distribution multiplier is <= 1.0, for each population separately. */
       maximum_no_improvement_stretch,                      /* The maximum number of subsequent generations without an improvement while the distribution multiplier is <= 1.0. */
     **individual_NIS;                                      /* The number of generations a solution has not improved. */
double maximum_number_of_evaluations,                       /* The maximum number of evaluations. */
       maximum_number_of_seconds,                           /* The maximum number of seconds. */
       tau,                                                 /* The selection truncation percentile (in [1/population_size,1]). */
    ***populations,                                         /* The populations containing the solutions. */
     **objective_values,                                    /* Objective values for population members. */
     **constraint_values,                                   /* Sum of all constraint violations for population members. */
     **ranks,                                               /* Ranks of population members. */
    ***selections,                                          /* Selected solutions, one for each population. */
     **objective_values_selections,                         /* Objective values of selected solutions. */
     **constraint_values_selections,                        /* Sum of all constraint violations of selected solutions. */
     **distribution_multipliers,                            /* Distribution multipliers of each FOS element of each population. */
       distribution_multiplier_increase,                    /* The multiplicative distribution multiplier increase. */
       distribution_multiplier_decrease,                    /* The multiplicative distribution multiplier decrease. */
       st_dev_ratio_threshold,                              /* The maximum ratio of the distance of the average improvement to the mean compared to the distance of one standard deviation before triggering AVS (SDR mechanism). */
       fitness_variance_tolerance,                          /* The minimum fitness variance level that is allowed. */
     **mean_vectors,                                        /* The mean vectors, one for each population. */
     **mean_shift_vector,                                   /* The mean vectors of the previous generation, one for each population. */
   ****decomposed_covariance_matrices,                      /* The covariance matrices to be used for the sampling. */
   ****decomposed_cholesky_factors_lower_triangle,          /* The unique lower triangular matrix of the Cholesky factorization for every linkage tree element. */
    ***full_covariance_matrix,
       eta_ams = 1.0,
       eta_cov = 1.0;
FOS  **linkage_model;
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/

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
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/


/*-=-=-=-=-=-=-=-=-=-=- Section Interpret Command Line -=-=-=-=-=-=-=-=-=-=-*/
/**
 * Parses and checks the command line.
 */
void interpretCommandLine( int argc, char **argv )
{
    startTimer();

    use_univariate_FOS = 0;
    learn_linkage_tree = 0;
    static_linkage_tree = 0;
    random_linkage_tree = 0;
    FOS_element_size = -1;
    block_start = 0;
    evaluations_for_statistics_hit = 0;
    haveNextNextGaussian = 0;

    parseCommandLine( argc, argv );

    if( use_guidelines )
    {
        tau                              = 0.35;
        if( maximum_number_of_populations == 1 )
            base_population_size         = (int) (36.1 + 7.58*log2((double) number_of_parameters));
        else
            base_population_size         = 10;
        //base_population_size           = (int) (10.0*pow((double) number_of_parameters,0.5));
        //base_population_size           = (int) (17.0 + 3.0*pow((double) number_of_parameters,1.5));
        //base_population_size           = (int) (4.0*pow((double) number_of_parameters,0.5));
        distribution_multiplier_decrease = 0.9;
        st_dev_ratio_threshold           = 1.0;
        maximum_no_improvement_stretch   = 25 + number_of_parameters;
    }
    FOS_element_ub = number_of_parameters;
    block_size = number_of_parameters;
    if( problem_index == 13 || problem_index == 15 ) block_size = 5;
    number_of_blocks = (number_of_parameters + block_size - 1) / block_size;
    if( FOS_element_size == -1 ) FOS_element_size = number_of_parameters;
    if( FOS_element_size == -2 ) learn_linkage_tree = 1;
    if( FOS_element_size == -3 ) static_linkage_tree = 1;
    if( FOS_element_size == -4 ) {static_linkage_tree = 1; FOS_element_ub = 100;}
    if( FOS_element_size == -5 ) {random_linkage_tree = 1; static_linkage_tree = 1; FOS_element_ub = 100;}
    if( FOS_element_size == 1 ) use_univariate_FOS = 1;

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
    black_box_evaluations         = 0;

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
                case 'h': printUsage(); break;
                case 'P': printAllInstalledProblems(); break;
                case 's': write_generational_statistics = 1; break;
                case 'w': write_generational_solutions  = 1; break;
                case 'v': print_verbose_overview        = 1; break;
                case 'r': use_vtr                       = 1; break;
                case 'g': use_guidelines                = 1; break;
                case 'b': black_box_evaluations         = 1; break;
                case 'f': parseFOSElementSize( index, argc, argv ); break;
                case 'S': fix_seed                      = 1; break;
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
    noError = noError && sscanf( argv[*index+10], "%lf", &maximum_number_of_evaluations );
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
    printf("Usage: RV-GOMEA [-?] pro dim low upp rot tau pop nop dmd srt eva vtr imp tol\n");
    printf(" -h: Prints out this usage information.\n");
    printf(" -P: Prints out a list of all installed optimization problems.\n");
    printf(" -s: Enables computing and writing of statistics every generation.\n");
    printf(" -w: Enables writing of solutions and their fitnesses every generation.\n");
    printf(" -v: Enables verbose mode. Prints the settings before starting the run.\n");
    printf(" -r: Enables use of vtr in termination condition (value-to-reach).\n");
    printf(" -b: Enables counting every partial evaluation as a full evaluation.\n");
    printf(" -f %%d: Sets linkage model that is used. Positive: Use a FOS with elements of %%d consecutive variables.\n");
    printf("     Use -1 for full linkage model, -2 for dynamic linkage tree learned from the population, -3 for fixed linkage tree learned from distance measure,\n");
    printf("     -4 for bounded fixed linkage tree learned from distance measure, -5 for fixed bounded linkage tree learned from random distance measure.\n");
    printf(" -g: Uses guidelines to override parameter settings for those parameters\n");
    printf("     for which a guideline is known in literature. These parameters are:\n");
    printf("     tau pop dmd srt imp\n");
    printf(" -S: A fixed random seed is used.\n");

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
        printf("Error: number of populations < 1 (read: %d). Require number of populations >= 1.", number_of_populations);
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
    
    if( rotation_angle > 0 && ( !learn_linkage_tree && FOS_element_size > 1 && FOS_element_size != block_size && FOS_element_size != number_of_parameters) )
    {
        printf("\n");
        printf("Error: invalid FOS element size (read %d). Must be %d, %d or %d.", FOS_element_size, 1, block_size, number_of_parameters );
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
    printf("# Max num of populations  = %d\n", maximum_number_of_populations);
    printf("# Dis. mult. decreaser    = %e\n", distribution_multiplier_decrease);
    printf("# St. dev. rat. threshold = %e\n", st_dev_ratio_threshold);
    printf("# Maximum numb. of eval.  = %lf\n", maximum_number_of_evaluations);
    printf("# Value to reach (vtr)    = %e\n", vtr);
    printf("# Max. no improv. stretch = %d\n", maximum_no_improvement_stretch);
    printf("# Fitness var. tolerance  = %e\n", fitness_variance_tolerance);
    printf("# Random seed             = %ld\n", random_seed);
    printf("#\n");
    printf("###################################################\n");
}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/


/*-=-=-=-=-=-=-=-=-=-=-=-=- Section Initialize -=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
/**
 * Performs initializations that are required before starting a run.
 */
void initialize( void )
{
    total_number_of_generations = 0;
    total_number_of_writes = 0;
    number_of_evaluations = 0;
    number_of_populations = 0;
    number_of_subgenerations_per_population_factor = 8;
    vtr_hit_status = 0;
    elitist_objective_value = 1e308;
    elitist_constraint_value = 1e308;

    initializeRandomNumberGenerator();
    if( fix_seed ) random_seed_changing = 611270;

    initializeProblem();

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
    constraint_values                = (double **) Malloc( maximum_number_of_populations*sizeof( double * ) );
    ranks                            = (double **) Malloc( maximum_number_of_populations*sizeof( double * ) );
    selections                       = (double ***) Malloc( maximum_number_of_populations*sizeof( double ** ) );
    objective_values_selections      = (double **) Malloc( maximum_number_of_populations*sizeof( double * ) );
    constraint_values_selections     = (double **) Malloc( maximum_number_of_populations*sizeof( double * ) );
    decomposed_cholesky_factors_lower_triangle = (double ****) Malloc(maximum_number_of_populations*sizeof(double *** ) );
    mean_vectors                     = (double **) Malloc( maximum_number_of_populations*sizeof( double * ) );
    mean_shift_vector                = (double **) Malloc( maximum_number_of_populations*sizeof( double * ) );
    decomposed_covariance_matrices   = (double ****) Malloc( maximum_number_of_populations * sizeof( double ***) );
    full_covariance_matrix           = (double ***) Malloc( maximum_number_of_populations * sizeof( double **) );
    distribution_multipliers         = (double **) Malloc( maximum_number_of_populations*sizeof( double * ) );
    samples_drawn_from_normal        = (int **) Malloc( maximum_number_of_populations*sizeof( int * ) );
    out_of_bounds_draws              = (int **) Malloc( maximum_number_of_populations*sizeof( int * ) );
    number_of_generations            = (int *) Malloc( maximum_number_of_populations*sizeof( int ) );
    linkage_model                    = (FOS **) Malloc( maximum_number_of_populations*sizeof( FOS *) );
    individual_NIS                   = ( int ** ) Malloc( maximum_number_of_populations*sizeof( int * ) );
}

void initializeNewPopulationMemory( int population_index )
{
    int i,j;

    if( population_index == 0 )
        population_sizes[population_index] = base_population_size;
    else
        population_sizes[population_index] = 2*population_sizes[population_index-1];
    selection_sizes[population_index] = (double) (tau * population_sizes[population_index] );

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

    mean_shift_vector[population_index] = (double *) Malloc( number_of_parameters*sizeof( double ) );
    
    individual_NIS[population_index] = (int*) Malloc( population_sizes[population_index]*sizeof(int));
    
    if( learn_linkage_tree )
    {
        distribution_multipliers[population_index]  = (double *) Malloc( 1*sizeof( double ) );
        samples_drawn_from_normal[population_index] = (int *) Malloc( 1*sizeof( int ) );
        out_of_bounds_draws[population_index]       = (int *) Malloc( 1*sizeof( int ) );
        linkage_model[population_index]             = (FOS *) Malloc( sizeof( FOS ) );
        linkage_model[population_index]->length     = 1;
        linkage_model[population_index]->sets       = (int **) Malloc( linkage_model[population_index]->length*sizeof( int * ) );
        linkage_model[population_index]->set_length = (int *) Malloc( linkage_model[population_index]->length*sizeof( int ) );
        for( i = 0; i < linkage_model[population_index]->length; i++ )
            linkage_model[population_index]->sets[i] = (int *) Malloc( 1*sizeof( int ) );
    }
    else initializeFOS( population_index );

    populations_terminated[population_index] = 0;

    no_improvement_stretch[population_index] = 0;

    number_of_generations[population_index] = 0;
}

void initializeNewPopulation()
{
    initializeNewPopulationMemory( number_of_populations );

    initializePopulationAndFitnessValues( number_of_populations );
    
    if( !learn_linkage_tree )
    {
        initializeCovarianceMatrices( number_of_populations );

        initializeDistributionMultipliers( number_of_populations );
    }
    
    computeRanksForOnePopulation( number_of_populations );

    number_of_populations++;
}

/**
 * Initializes the linkage tree
 */
void initializeFOS( int population_index )
{
    int      i;
    FILE    *file;
    FOS     *new_FOS;

    fflush( stdout );
    file = fopen( "FOS.in", "r" );
    if( file != NULL )
    {
        if( population_index == 0 )
            new_FOS = readFOSFromFile( file );
        else
            new_FOS = copyFOS( linkage_model[0] );
    }
    else if( static_linkage_tree )
    {
        if( population_index == 0 )
            new_FOS = learnLinkageTreeRVGOMEA( population_index );
        else
            new_FOS = copyFOS( linkage_model[0] );
    }
    else
    {
        new_FOS = (FOS*) Malloc(sizeof(FOS));
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
    }
    linkage_model[population_index] = new_FOS;
}

void initializeProblem( void )
{
    switch( problem_index )
    {
    default: break;
    }
}

/**
 * Initializes the distribution multipliers.
 */
void initializeDistributionMultipliers( int population_index )
{
    int j;

    if( learn_linkage_tree )
    {
        free( distribution_multipliers[population_index] );
        free( samples_drawn_from_normal[population_index] );
        free( out_of_bounds_draws[population_index] );
    }

    distribution_multipliers[population_index] = (double *) Malloc( linkage_model[population_index]->length*sizeof( double ) );
    for( j = 0; j < linkage_model[population_index]->length; j++ )
        distribution_multipliers[population_index][j] = 1.0;

    samples_drawn_from_normal[population_index] = (int *) Malloc( linkage_model[population_index]->length*sizeof( int ) );
    out_of_bounds_draws[population_index]       = (int *) Malloc( linkage_model[population_index]->length*sizeof( int ) );

    distribution_multiplier_increase = 1.0/distribution_multiplier_decrease;
}

/**
 * Initializes the populations and the fitness values.
 */
void initializePopulationAndFitnessValues( int population_index )
{
    int     j, k;

    for( j = 0; j < population_sizes[population_index]; j++ )
    {
        individual_NIS[population_index][j] = 0;
        for( k = 0; k < number_of_parameters; k++ )
            populations[population_index][j][k] = lower_init_ranges[k] + (upper_init_ranges[k] - lower_init_ranges[k])*randomRealUniform01();

        installedProblemEvaluation( problem_index, populations[population_index][j], &(objective_values[population_index][j]), &(constraint_values[population_index][j]), number_of_parameters, NULL, NULL, 0, 0 );
    }
}


/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/

FOS *learnLinkageTreeRVGOMEA( int population_index )
{
    int i;
    FOS *new_FOS;

    new_FOS = learnLinkageTree( full_covariance_matrix[population_index] );
    if( learn_linkage_tree && number_of_generations[population_index] > 0 )
        inheritDistributionMultipliers( new_FOS, linkage_model[population_index], distribution_multipliers[population_index] );

    if( learn_linkage_tree )
    {
        for( i = 0; i < linkage_model[population_index]->length; i++ )
            free( linkage_model[population_index]->sets[i] );
        free( linkage_model[population_index]->sets );
        free( linkage_model[population_index]->set_length );
        free( linkage_model[population_index]);
    }
    return( new_FOS );
}

void inheritDistributionMultipliers( FOS *new_FOS, FOS *prev_FOS, double *multipliers )
{
    int      i, *permutation;
    double   *multipliers_copy;

    multipliers_copy = (double*) Malloc(new_FOS->length*sizeof(double));
    for( i = 0; i < new_FOS->length; i++ )
        multipliers_copy[i] = multipliers[i];

    permutation = matchFOSElements( new_FOS, prev_FOS );

    for( i = 0; i < new_FOS->length; i++ )
        multipliers[permutation[i]] = multipliers_copy[i];

    free( multipliers_copy );
    free( permutation );
}

/*-=-=-=-=-=-=-=-=-=-=-=-=-=- Section Ranking -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
/**
 * Computes the ranks of the solutions in all populations.
 */
void computeRanksForAllPopulations( void )
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
            population_constraint_avg, population_constraint_var, population_constraint_best, population_constraint_worst;
    FILE   *file;

    /* Average, best and worst */
    population_objective_avg    = 0.0;
    population_constraint_avg   = 0.0;
    population_objective_best   = objective_values[population_index][0];
    population_constraint_best  = constraint_values[population_index][0];
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
        if( betterFitness( objective_values[population_index][j], constraint_values[population_index][j], population_objective_best, population_constraint_best ) )
        {
            population_objective_best = objective_values[population_index][j];
            population_constraint_best = constraint_values[population_index][j];
        }
    }
    population_objective_avg  = population_objective_avg / ((double) population_sizes[population_index]);
    population_constraint_avg = population_constraint_avg / ((double) population_sizes[population_index]);

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
    if( total_number_of_writes == 0 )
    {
        file = fopen( "statistics.dat", "w" );

        sprintf( string, "# Generation  Evaluations  Time(s)  Best-obj. Best-cons. [Pop.index  Subgen.  Pop.size  Dis.mult.[0]  Pop.best.obj. Pop.avg.obj.  Pop.var.obj. Pop.worst.obj.  Pop.best.con. Pop.avg.con.  Pop.var.con. Pop.worst.con.]\n" );
        fputs( string, file );
    }
    else
        file = fopen( "statistics.dat", "a" );

    sprintf( string, "%10d %11lf %11.3lf %15.10e %13e  ", total_number_of_generations, number_of_evaluations, getTimer(), elitist_objective_value, elitist_constraint_value );
    fputs( string, file );

    //sprintf( string, "[ %4d %6d %10d %13e %13e %13e %13e %13e %13e %13e %13e %13e ]", population_index, number_of_generations[population_index], population_sizes[population_index], distribution_multipliers[population_index][0], population_objective_best, population_objective_avg, population_objective_var, population_objective_worst, population_constraint_best, population_constraint_avg, population_constraint_var, population_constraint_worst );
    //fputs( string, file );
    
    sprintf( string, "\n");
    fputs( string, file );

    fclose( file );

    total_number_of_writes++;
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

    file_selection = NULL;
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
    static int *c = NULL;
    if( c == NULL ){c = (int *) Malloc( sizeof( int ) ); c[0]=0;}

    /* First find the best of all */
    determineBestSolutionInCurrentPopulations( &population_index_best, &individual_index_best );

    /* Then output it */
    if( final )
        sprintf( string, "best_generation_final.dat" );
    else
        sprintf( string, "best_generation_%05d.dat", c[0] );
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
    c[0]++;

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

    if( checkSubgenerationTerminationConditions() )
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

short checkSubgenerationTerminationConditions( void )
{
    if( checkNumberOfEvaluationsTerminationCondition() )
        return( 1 );

    if( checkVTRTerminationCondition() )
        return( 1 );

    if( checkTimeLimitTerminationCondition() )
        return( 1 );

    return( 0 );
}

short checkTimeLimitTerminationCondition( void )
{
    return( maximum_number_of_seconds > 0 && getTimer() > maximum_number_of_seconds );
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
    int i, j;
    short converged;

    for( i = 0; i < number_of_populations; i++ )
    {
        if( !populations_terminated[i] )
        {
            converged = 1;
            for( j = 0; j < linkage_model[i]->length; j++ )
            {
                if( distribution_multipliers[i][j] > 1e-10 )
                {
                    converged = 0;
                    break;
                }
            }

            if( converged )
                populations_terminated[i] = 1;
        }
    }
}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/





/*-=-=-=-=-=-=-=-=-=-=-=-=-= Section Selection =-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
/**
 * Makes a set of selected solutions for each population.
 */
void makeSelections( void )
{
    int i;

    for( i = 0; i < number_of_populations; i++ )
        if( !populations_terminated[i] )
            makeSelectionsForOnePopulation( i );
}

/**
 * Performs truncation selection on a single population.
 */
void makeSelectionsForOnePopulation( int population_index )
{
    int i, j, *sorted;

    sorted = mergeSort( ranks[population_index], population_sizes[population_index] );

    if( ranks[population_index][sorted[selection_sizes[population_index]-1]] == 0 ){
        makeSelectionsForOnePopulationUsingDiversityOnRank0( population_index );}
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
        nn_distances[i] = distanceEuclidean( populations[population_index][preselection_indices[i]], populations[population_index][selection_indices[number_selected_so_far-1]], number_of_parameters );

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
            value = distanceEuclidean( populations[population_index][preselection_indices[i]], populations[population_index][selection_indices[number_selected_so_far-1]], number_of_parameters );
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

    ezilaitiniParametersForSampling( population_index );
}

/**
 * Estimates the paramaters of the multivariate
 * normal distribution for a specified population.
 */
void estimateParameters( int population_index )
{
    if( !populations_terminated[population_index] )
    {
        estimateMeanVectorML( population_index );
        
        if( learn_linkage_tree )
        {
            estimateFullCovarianceMatrixML( population_index );
            
            linkage_model[population_index] = learnLinkageTreeRVGOMEA( population_index );

            initializeCovarianceMatrices( population_index );

            if( number_of_generations[population_index] == 0 )
                initializeDistributionMultipliers( population_index );
        }
        
        estimateParametersML( population_index );
    }
}

void estimateParametersML( int population_index )
{
    int i, j;
    
    /* Change the focus of the search to the best solution */
    for( i = 0; i < linkage_model[population_index]->length; i++ )
        if( distribution_multipliers[population_index][i] < 1.0 )
            for( j = 0; j < linkage_model[population_index]->set_length[i]; j++)
                mean_vectors[population_index][linkage_model[population_index]->sets[i][j]] = selections[population_index][0][linkage_model[population_index]->sets[i][j]];

    estimateCovarianceMatricesML( population_index );
}

/**
 * Computes the sample mean for a specified population.
 */
void estimateMeanVectorML( int population_index )
{
    int i, j;
    double new_mean;
    
    for( i = 0; i < number_of_parameters; i++ )
    {
        new_mean = 0.0;
        for( j = 0; j < selection_sizes[population_index]; j++ )
            new_mean += selections[population_index][j][i];
        new_mean /= (double) selection_sizes[population_index];

        if( number_of_generations[population_index] > 0 )
            mean_shift_vector[population_index][i] = new_mean - mean_vectors[population_index][i];

        mean_vectors[population_index][i] = new_mean;
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

    full_covariance_matrix[population_index] = (double **) Malloc( number_of_parameters*sizeof( double * ) );
    for( j = 0; j < number_of_parameters; j++ )
        full_covariance_matrix[population_index][j] = (double *) Malloc( number_of_parameters*sizeof( double ) );
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

void estimateCovarianceMatricesML( int population_index )
{
    int i, j, k, m, vara, varb;
    double cov;

    /* First do the maximum-likelihood estimate from data */
    for( i = 0; i < linkage_model[population_index]->length; i++ )
    {
        for( j = 0; j < linkage_model[population_index]->set_length[i]; j++ )
        {
            vara = linkage_model[population_index]->sets[i][j];
            for( k = j; k < linkage_model[population_index]->set_length[i]; k++ )
            {
                varb = linkage_model[population_index]->sets[i][k];
                
                if( learn_linkage_tree )
                {
                    cov = full_covariance_matrix[population_index][vara][varb];
                }
                else
                {
                    cov = 0.0;
                    for( m = 0; m < selection_sizes[population_index]; m++ )
                        cov += (selections[population_index][m][vara]-mean_vectors[population_index][vara])*(selections[population_index][m][varb]-mean_vectors[population_index][varb]);

                    cov /= (double) selection_sizes[population_index];
                }
                decomposed_covariance_matrices[population_index][i][j][k] = (1-eta_cov)*decomposed_covariance_matrices[population_index][i][j][k]+ eta_cov*cov;
                decomposed_covariance_matrices[population_index][i][k][j] = decomposed_covariance_matrices[population_index][i][j][k];
            }
        }
    }
}

void initializeCovarianceMatrices( int population_index )
{
    int j, k, m;

    decomposed_covariance_matrices[population_index] = (double ***) Malloc( linkage_model[population_index]->length * sizeof( double **) );
    for( j = 0; j < linkage_model[population_index]->length; j++ )
    {
        decomposed_covariance_matrices[population_index][j] = (double **) Malloc( linkage_model[population_index]->set_length[j]*sizeof( double * ) );
        for( k = 0; k < linkage_model[population_index]->set_length[j]; k++)
        {
            decomposed_covariance_matrices[population_index][j][k] = (double *) Malloc( linkage_model[population_index]->set_length[j]*sizeof( double ) );
            for( m = 0; m < linkage_model[population_index]->set_length[j]; m++)
            {
                decomposed_covariance_matrices[population_index][j][k][m] = 1.0;
            }
        }
    }
}

void copyBestSolutionsToAllPopulations( void )
{
    int i;

    for( i = 0; i < number_of_populations; i++ )
        copyBestSolutionsToPopulation( i );
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

void getBestInPopulation( int population_index, int *individual_index )
{
    int i;

    *individual_index = 0;
    for( i = 0; i < population_sizes[population_index]; i++ )
        if( betterFitness(objective_values[population_index][i], constraint_values[population_index][i], objective_values[population_index][*individual_index], constraint_values[population_index][*individual_index]))
            *individual_index = i;
}

void getOverallBest( int *population_index, int *individual_index )
{
    int i, best_individual_index;

    *population_index = 0;
    getBestInPopulation( 0, &best_individual_index );
    *individual_index = best_individual_index;
    for( i = 0; i < number_of_populations; i++ )
    {
        getBestInPopulation( i, &best_individual_index );
        if( betterFitness(objective_values[i][best_individual_index], constraint_values[i][best_individual_index], objective_values[*population_index][*individual_index], constraint_values[*population_index][*individual_index]))
        {
            *population_index = i;
            *individual_index = best_individual_index;
        }
    }
}

void evaluateCompletePopulation( int population_index )
{
    int j;

    for( j = 0; j < population_sizes[population_index]; j++ )
        installedProblemEvaluation( problem_index, populations[population_index][j], &(objective_values[population_index][j]), &(constraint_values[population_index][j]), number_of_parameters, NULL, NULL, 0, 0 );
}

/**
 * Applies the distribution multipliers.
 */
void applyDistributionMultipliersToAllPopulations( void )
{
    int i;

    for( i = 0; i < number_of_populations; i++ )
        applyDistributionMultipliers( i );
}

void applyDistributionMultipliers( int population_index )
{
    int j, k, m;

    if( !populations_terminated[population_index] )
    {
        for( j = 0; j < linkage_model[population_index]->length; j++ )
            for( k = 0; k < linkage_model[population_index]->set_length[j]; k++ )
                for( m = 0; m < linkage_model[population_index]->set_length[j]; m++ )
                    decomposed_covariance_matrices[population_index][j][k][m] *= distribution_multipliers[population_index][j];
    }
}

void generateAndEvaluateNewSolutionsToFillAllPopulations( void )
{
    int i;

    for( i = 0; i < number_of_populations; i++ )
        generateAndEvaluateNewSolutionsToFillPopulation( i );
}

/**
 * Generates new solutions for each
 * of the populations in turn.
 */
void generateAndEvaluateNewSolutionsToFillPopulation( int population_index )
{
    short   generationalImprovement, *FOS_element_caused_improvement, all_multipliers_leq_one, *individual_improved, apply_AMS;
    int     oj, i, j, k, *fos_order, number_of_AMS_solutions,
            best_individual_index;
    double  alpha_AMS;

    computeParametersForSampling( population_index );

    if( !populations_terminated[population_index] )
    {
        if( !black_box_evaluations && (number_of_generations[population_index]+1) % 50 == 0 )
            evaluateCompletePopulation( population_index );

        FOS_element_caused_improvement = (short *) Malloc( linkage_model[population_index]->length*sizeof( short ) );
        individual_improved = (short *) Malloc( population_sizes[population_index]*sizeof( short ) );
        for( k = 1; k < population_sizes[population_index]; k++ )
            individual_improved[k] = 0;

        alpha_AMS = 0.5*tau*(((double) population_sizes[population_index])/((double) (population_sizes[population_index]-1)));
        number_of_AMS_solutions = (int) (alpha_AMS*(population_sizes[population_index]-1));
        fos_order = randomPermutation(linkage_model[population_index]->length);
        for( oj = 0; oj < linkage_model[population_index]->length; oj++ )
        {
            j = fos_order[oj];

            samples_drawn_from_normal[population_index][j] = 0;
            out_of_bounds_draws[population_index][j]       = 0;
            FOS_element_caused_improvement[j] = 0;

            apply_AMS = 1;
            for( k = 1; k < population_sizes[population_index]; k++ )
            {
                if( k > number_of_AMS_solutions ) apply_AMS = 0;
                individual_improved[k] |= generateNewSolutionFromFOSElement( population_index, j, k, apply_AMS );
            }

            FOS_element_caused_improvement[j] = adaptDistributionMultipliers( population_index, j );
        }
        free( fos_order );

        if( number_of_generations[population_index] > 0 )
        {
            for( k = 1; k <= number_of_AMS_solutions; k++ )
                individual_improved[k] |= applyAMS(population_index, k);
        }
        
        for( i = 1; i < population_sizes[population_index]; i++ )
            if( !individual_improved[i] ) individual_NIS[population_index][i]++;
            else individual_NIS[population_index][i] = 0;
        
        getBestInPopulation( population_index, &best_individual_index );
        for( k = 1; k < population_sizes[population_index]; k++ )
            if( individual_NIS[population_index][k] > maximum_no_improvement_stretch )
                applyForcedImprovements( population_index, k, best_individual_index );


        generationalImprovement = 0;
        for( j = 0; j < linkage_model[population_index]->length; j++ )
            if( FOS_element_caused_improvement[j] )
                generationalImprovement = 1;

        if( generationalImprovement )
            no_improvement_stretch[population_index] = 0;
        else
        {
            all_multipliers_leq_one = 1;
            for( j = 0; j < linkage_model[population_index]->length; j++ )
                if( distribution_multipliers[population_index][j] > 1.0 )
                {all_multipliers_leq_one = 0; break;}

            if( all_multipliers_leq_one )
                (no_improvement_stretch[population_index])++;
        }

        free( individual_improved );
        free( FOS_element_caused_improvement );
    }
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
        decomposed_cholesky_factors_lower_triangle[population_index] = (double ***) Malloc(linkage_model[population_index]->length * sizeof(double**));
        for( i = 0; i < linkage_model[population_index]->length; i++ )
            decomposed_cholesky_factors_lower_triangle[population_index][i] = choleskyDecomposition( decomposed_covariance_matrices[population_index][i], linkage_model[population_index]->set_length[i] );
    }
}

/**
 * Generates and returns a single new solution by drawing
 * a sample for the variables in the selected FOS element
 * and inserting this into the population.
 */
double *generateNewPartialSolutionFromFOSElement( int population_index, int FOS_index )
{

    short   ready;
    int     i, times_not_in_bounds, num_indices, *indices;
    double *result, *z;

    num_indices = linkage_model[population_index]->set_length[FOS_index];
    indices = linkage_model[population_index]->sets[FOS_index];

    times_not_in_bounds = -1;
    out_of_bounds_draws[population_index][FOS_index]--;

    ready = 0;
    do
    {
        times_not_in_bounds++;
        samples_drawn_from_normal[population_index][FOS_index]++;
        out_of_bounds_draws[population_index][FOS_index]++;
        if( times_not_in_bounds >= 100 )
        {
            result = (double *) Malloc( num_indices*sizeof( double ) );
            for( i = 0; i < num_indices; i++ )
                result[i] = lower_init_ranges[indices[i]] + (upper_init_ranges[indices[i]] - lower_init_ranges[indices[i]])*randomRealUniform01();
        }
        else
        {
            z = (double *) Malloc( num_indices*sizeof( double ) );

            for( i = 0; i < num_indices; i++ )
                z[i] = random1DNormalUnit();

            if( use_univariate_FOS )
            {
                result = (double*) Malloc(1*sizeof(double));
                result[0] = z[0]*sqrt(decomposed_covariance_matrices[population_index][FOS_index][0][0]) + mean_vectors[population_index][indices[0]];
            }
            else
            {
                result = matrixVectorMultiplication( decomposed_cholesky_factors_lower_triangle[population_index][FOS_index], z, num_indices, num_indices );
                for( i = 0; i < num_indices; i++ )
                    result[i] += mean_vectors[population_index][indices[i]];
            }

            free( z );
        }

        ready = 1;
        for( i = 0; i < num_indices; i++ )
        {
            if( !isParameterInRangeBounds( result[i], indices[i] ) )
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
 * Generates and returns a single new solution by drawing
 * a single sample from a specified model.
 */
short generateNewSolutionFromFOSElement( int population_index, int FOS_index, int individual_index, short apply_AMS )
{
    int j, m, im, *indices, num_indices, *touched_indices, num_touched_indices;
    double *result, *individual_backup, obj_val, cons_val, delta_AMS, shrink_factor;
    short improvement, any_improvement, out_of_range;

    delta_AMS = 2.0;
    improvement = 0;
    any_improvement = 0;
    num_indices = linkage_model[population_index]->set_length[FOS_index];
    indices = linkage_model[population_index]->sets[FOS_index];
    num_touched_indices = num_indices;
    touched_indices = indices;
    individual_backup = (double*) Malloc( num_touched_indices * sizeof( double ) );

    for( j = 0; j < num_touched_indices; j++ )
        individual_backup[j] = populations[population_index][individual_index][touched_indices[j]];

    result = generateNewPartialSolutionFromFOSElement(population_index, FOS_index);

    for( j = 0; j < num_indices; j++ )
        populations[population_index][individual_index][indices[j]] = result[j];

    if( apply_AMS && (number_of_generations[population_index] > 0) )
    {
        out_of_range  = 1;
        shrink_factor = 2;
        while( (out_of_range == 1) && (shrink_factor > 1e-10) )
        {
            shrink_factor *= 0.5;
            out_of_range   = 0;
            for( m = 0; m < num_indices; m++ )
            {
                im = indices[m];
                result[m] = populations[population_index][individual_index][im]+shrink_factor*delta_AMS*distribution_multipliers[population_index][FOS_index]*mean_shift_vector[population_index][im];
                if( !isParameterInRangeBounds( result[m], im ) )
                {
                    out_of_range = 1;
                    break;
                }
            }
        }
        if( !out_of_range )
        {
            for( m = 0; m < num_indices; m++ )
            {
                populations[population_index][individual_index][indices[m]] = result[m];
            }
        }
    }

    installedProblemEvaluation( problem_index, populations[population_index][individual_index], &obj_val, &cons_val, num_touched_indices, touched_indices, individual_backup, objective_values[population_index][individual_index], constraint_values[population_index][individual_index] );
    improvement = betterFitness(obj_val, cons_val, objective_values[population_index][individual_index], constraint_values[population_index][individual_index]);
    if( improvement )
    {
        any_improvement = 1;
        objective_values[population_index][individual_index] = obj_val;
        constraint_values[population_index][individual_index] = cons_val;
        for( j = 0; j < num_touched_indices; j++ )
            individual_backup[j] = populations[population_index][individual_index][touched_indices[j]];
    }
    free( result );

    if( !any_improvement && randomRealUniform01() >= 0.05 )
    {
        for( j = 0; j < num_touched_indices; j++ )
            populations[population_index][individual_index][touched_indices[j]] = individual_backup[j];
    }
    else
    {
        objective_values[population_index][individual_index] = obj_val;
        constraint_values[population_index][individual_index] = cons_val;
    }

    free( individual_backup );
    return( any_improvement );
}

short applyAMS( int population_index, int individual_index )
{
    short out_of_range, improvement;
    double shrink_factor, delta_AMS, *solution_AMS, obj_val, cons_val;
    int m;

    delta_AMS     = 2;
    out_of_range  = 1;
    shrink_factor = 2;
    improvement   = 0;
    solution_AMS = (double *) Malloc( number_of_parameters*sizeof( double ) );
    while( (out_of_range == 1) && (shrink_factor > 1e-10) )
    {
        shrink_factor *= 0.5;
        out_of_range   = 0;
        for( m = 0; m < number_of_parameters; m++ )
        {
            solution_AMS[m] = populations[population_index][individual_index][m]
                    + shrink_factor*delta_AMS*mean_shift_vector[population_index][m]; //*distribution_multipliers[population_index][FOS_index]
            if( !isParameterInRangeBounds( solution_AMS[m], m ) )
            {
                out_of_range = 1;
                break;
            }
        }
    }
    if( !out_of_range )
    {
        installedProblemEvaluation( problem_index, solution_AMS, &obj_val, &cons_val, number_of_parameters, NULL, NULL, 0, 0 );
        if( randomRealUniform01() < 0.05 || betterFitness(obj_val, cons_val, objective_values[population_index][individual_index], constraint_values[population_index][individual_index]))
        {
            objective_values[population_index][individual_index] = obj_val;
            constraint_values[population_index][individual_index] = cons_val;
            for( m = 0; m < number_of_parameters; m++ )
                populations[population_index][individual_index][m] = solution_AMS[m];
            improvement = 1;
        }
    }
    free( solution_AMS );
    return( improvement );
}

void applyForcedImprovements( int population_index, int individual_index, int donor_index )
{
    int i, io, j, *order, *touched_indices, num_touched_indices, FOS_element_index;
    double *FI_backup, obj_val, cons_val, alpha;
    short improvement;

    improvement = 0;
    alpha = 1.0;
    
    while( alpha >= 0.01 )
    {
        alpha *= 0.5;
        order = randomPermutation( linkage_model[population_index]->length );
        for( io = 0; io < linkage_model[population_index]->length; io++ )
        {
            i = order[io];
            touched_indices = linkage_model[population_index]->sets[i];
            num_touched_indices = linkage_model[population_index]->set_length[i];
            if( rotation_angle > 0 && linkage_model[population_index]->set_length[i] < block_size )
            {
                num_touched_indices = block_size;
                FOS_element_index = linkage_model[population_index]->sets[i][0]/block_size;
                touched_indices = (int*) Malloc( block_size*sizeof( int ) );
                for( j = 0; j < block_size; j++ )
                    touched_indices[j] = FOS_element_index*block_size + j;
            }
            FI_backup = (double*) Malloc( num_touched_indices*sizeof(double));
            for( j = 0; j < num_touched_indices; j++ )
            {
                FI_backup[j] = populations[population_index][individual_index][touched_indices[j]];
                populations[population_index][individual_index][touched_indices[j]] = alpha*populations[population_index][individual_index][touched_indices[j]] + (1-alpha)*populations[population_index][donor_index][touched_indices[j]];
            }
            installedProblemEvaluation( problem_index, populations[population_index][individual_index], &obj_val, &cons_val, num_touched_indices, touched_indices, FI_backup, objective_values[population_index][individual_index], constraint_values[population_index][individual_index] );
            improvement = betterFitness( obj_val, cons_val, objective_values[population_index][individual_index], constraint_values[population_index][individual_index] );
            //printf("alpha=%.1e\tf=%.30e\n",alpha,obj_val);

            if( !improvement )
                for( j = 0; j < num_touched_indices; j++ )
                    populations[population_index][individual_index][touched_indices[j]] = FI_backup[j];
            else
            {
                objective_values[population_index][individual_index] = obj_val;
                constraint_values[population_index][individual_index] = cons_val;
            }

            free( FI_backup );

            if( rotation_angle > 0 && linkage_model[population_index]->set_length[i] < block_size )
                free( touched_indices );

            if( improvement )
                break;
        }
        free( order );
        if( improvement )
            break;
    }

    if( improvement )
    {
        objective_values[population_index][individual_index] = obj_val;
        constraint_values[population_index][individual_index] = cons_val;
    }
    else
    {
        for( i = 0; i < number_of_parameters; i++ )
            populations[population_index][individual_index][i] = populations[population_index][donor_index][i];
        objective_values[population_index][individual_index] = objective_values[population_index][donor_index];
        constraint_values[population_index][individual_index] = constraint_values[population_index][donor_index];
    }
}

/**
 * Adapts distribution multipliers according to SDR-AVS mechanism.
 * Returns whether the FOS element with index FOS_index has caused
 * an improvement in population_index.
 */
short adaptDistributionMultipliers( int population_index, int FOS_index )
{
    short  improvementForFOSElement;
    int    i, j;
    double st_dev_ratio, increase_for_FOS_element, decrease_for_FOS_element;

    i = population_index;
    j = FOS_index;
    improvementForFOSElement = 0;
    increase_for_FOS_element = distribution_multiplier_increase;
    decrease_for_FOS_element = 1.0/increase_for_FOS_element;
    if( !populations_terminated[i] )
    {
        if( (((double) out_of_bounds_draws[i][j])/((double) samples_drawn_from_normal[i][j])) > 0.9 )
            distribution_multipliers[i][j] *= 0.5;

        improvementForFOSElement = generationalImprovementForOnePopulationForFOSElement( i, j, &st_dev_ratio );

        if( improvementForFOSElement )
        {
            if( distribution_multipliers[i][j] < 1.0 )
                distribution_multipliers[i][j] = 1.0;

            if( st_dev_ratio > st_dev_ratio_threshold )
                distribution_multipliers[i][j] *= increase_for_FOS_element;
        }
        else
        {
            if( (distribution_multipliers[i][j] > 1.0) || (no_improvement_stretch[i] >= maximum_no_improvement_stretch) )
                distribution_multipliers[i][j] *= decrease_for_FOS_element;

            if( no_improvement_stretch[i] < maximum_no_improvement_stretch && distribution_multipliers[i][j] < 1.0)
                distribution_multipliers[i][j] = 1.0;
        }
    }
    return( improvementForFOSElement );
}


/**
 * Determines whether an improvement is found for a specified
 * population. Returns 1 in case of an improvement, 0 otherwise.
 * The standard-deviation ratio required by the SDR-AVS
 * mechanism is computed and returned in the pointer variable.
 */
short generationalImprovementForOnePopulationForFOSElement( int population_index, int FOS_index, double *st_dev_ratio )
{
    int     i, j, index_best_population,
            number_of_improvements, *indices, num_indices;
    double *average_parameters_of_improvements;
    short generationalImprovement;

    generationalImprovement = 0;
    indices = linkage_model[population_index]->sets[FOS_index];
    num_indices = linkage_model[population_index]->set_length[FOS_index];

    // Determine best in the population and the average improvement parameters
    average_parameters_of_improvements = (double *) Malloc( num_indices*sizeof( double ) );
    for( i = 0; i < num_indices; i++ )
        average_parameters_of_improvements[i] = 0.0;

    index_best_population   = 0;
    number_of_improvements  = 0;
    for( i = 0; i < population_sizes[population_index]; i++ )
    {
        if( betterFitness( objective_values[population_index][i], constraint_values[population_index][i],
                           objective_values[population_index][index_best_population], constraint_values[population_index][index_best_population] ) )
            index_best_population = i;

        if( betterFitness( objective_values[population_index][i], constraint_values[population_index][i],
                           objective_values_selections[population_index][0], constraint_values_selections[population_index][0] ) )
        {
            number_of_improvements++;
            for( j = 0; j < num_indices; j++ )
                average_parameters_of_improvements[j] += populations[population_index][i][indices[j]];
        }
    }

    // Determine st.dev. ratio
    *st_dev_ratio = 0.0;
    if( number_of_improvements > 0 )
    {
        for( i = 0; i < num_indices; i++ )
            average_parameters_of_improvements[i] /= (double) number_of_improvements;

        *st_dev_ratio = getStDevRatioForFOSElement( population_index, average_parameters_of_improvements, FOS_index );
        generationalImprovement = 1;
    }

    free( average_parameters_of_improvements );

    return( generationalImprovement );
}

/**
 * Computes and returns the standard-deviation-ratio
 * of a given point for a given model.
 */
double getStDevRatioForFOSElement( int population_index, double *parameters, int FOS_index )
{
    int      i, *indices, num_indices;
    double **inverse, result, *x_min_mu, *z;

    indices = linkage_model[population_index]->sets[FOS_index];
    num_indices = linkage_model[population_index]->set_length[FOS_index];
    x_min_mu = (double *) Malloc( num_indices*sizeof( double ) );

    for( i = 0; i < num_indices; i++ )
        x_min_mu[i] = parameters[i]-mean_vectors[population_index][indices[i]];
    result = 0.0;
    
    if( use_univariate_FOS )
    {
        result = fabs( x_min_mu[0]/sqrt(decomposed_covariance_matrices[population_index][FOS_index][0][0]) );
    }
    else
    {
        inverse = matrixLowerTriangularInverse( decomposed_cholesky_factors_lower_triangle[population_index][FOS_index], num_indices );
        z = matrixVectorMultiplication( inverse, x_min_mu, num_indices, num_indices );

        for( i = 0; i < num_indices; i++ )
        {
            if( fabs( z[i] ) > result )
                result = fabs( z[i] );
        }
        
        free( z );
        for( i = 0; i < num_indices; i++ )
            free( inverse[i] );
        free( inverse );
    }

    free( x_min_mu );

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

/**
 * Undoes initialization procedure by freeing up memory.
 */
void ezilaitiniMemory( void )
{
    int i, j, k;

    for( i = 0; i < number_of_populations; i++ )
    {
        for( j = 0; j < population_sizes[i]; j++ )
            free( populations[i][j] );
        free( populations[i] );

        free( objective_values[i] );

        free( constraint_values[i] );

        free( ranks[i] );

        for( j = 0; j < selection_sizes[i]; j++ )
            free( selections[i][j] );
        free( selections[i] );

        free( objective_values_selections[i] );

        free( constraint_values_selections[i] );

        free( mean_vectors[i] );

        free( mean_shift_vector[i] );
        
        if( !learn_linkage_tree )
        {
            for( j = 0; j < linkage_model[i]->length; j++)
            {
                for( k = 0; k < linkage_model[i]->set_length[j]; k++ )
                    free( decomposed_covariance_matrices[i][j][k] );
                free( decomposed_covariance_matrices[i][j] );
            }
            free( decomposed_covariance_matrices[i] );
        }
        
        free( individual_NIS[i] );
        
        ezilaitiniDistributionMultipliers( i );

        ezilaitiniFOS( linkage_model[i] );
    }

    free( distribution_multipliers );
    free( samples_drawn_from_normal );
    free( out_of_bounds_draws );
    free( individual_NIS );
    free( full_covariance_matrix );
    free( decomposed_covariance_matrices );
    free( decomposed_cholesky_factors_lower_triangle );
    free( lower_range_bounds );
    free( upper_range_bounds );
    free( lower_init_ranges );
    free( upper_init_ranges );
    free( populations_terminated );
    free( no_improvement_stretch );
    free( populations );
    free( objective_values );
    free( constraint_values );
    free( ranks );
    free( selections );
    free( objective_values_selections );
    free( constraint_values_selections );
    free( mean_vectors );
    free( mean_shift_vector );
    free( population_sizes );
    free( selection_sizes );
    free( number_of_generations );
    free( linkage_model );
}

/**
 * Undoes initialization procedure by freeing up memory.
 */
void ezilaitiniDistributionMultipliers( int population_index )
{
    free( distribution_multipliers[population_index] );
    free( samples_drawn_from_normal[population_index] );
    free( out_of_bounds_draws[population_index] );
}

void ezilaitiniCovarianceMatrices( int population_index )
{
    int i,j,k;
    
    i = population_index;
    for( j = 0; j < linkage_model[i]->length; j++ )
    {
        for( k = 0; k < linkage_model[i]->set_length[j]; k++ )
            free( decomposed_covariance_matrices[i][j][k] );
        free( decomposed_covariance_matrices[i][j] );
    }
    free( decomposed_covariance_matrices[i] );
}

void ezilaitiniParametersAllPopulations( void )
{
    int i;

    for( i = 0; i < number_of_populations; i++)
    {
        ezilaitiniParametersForSampling(i);
    }
}

/**
 * Frees memory of the Cholesky decompositions required for sampling.
 */
void ezilaitiniParametersForSampling( int population_index )
{
    int i, j;

    if( !use_univariate_FOS )
    {
        for( i = 0; i < linkage_model[population_index]->length; i++ )
        {
            for( j = 0; j < linkage_model[population_index]->set_length[i]; j++ )
                free( decomposed_cholesky_factors_lower_triangle[population_index][i][j] );
            free( decomposed_cholesky_factors_lower_triangle[population_index][i] );
        }
        free( decomposed_cholesky_factors_lower_triangle[population_index] );
    }
    if( learn_linkage_tree )
    {
        ezilaitiniCovarianceMatrices( population_index );
        
        for( i = 0; i < number_of_parameters; i++ )
            free( full_covariance_matrix[population_index][i] );
        free( full_covariance_matrix[population_index] );
    }
}

/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/





/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=- Section Run -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/

void generationalStepAllPopulationsRecursiveFold(int population_index_smallest, int population_index_biggest );
void generationalStepAllPopulations()
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
            if( !populations_terminated[population_index] )
            {
                makeSelectionsForOnePopulation( population_index );

                makePopulation( population_index );

                number_of_generations[population_index]++;

                if( checkSubgenerationTerminationConditions() )
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

void runAllPopulations()
{
    while( !checkTerminationCondition() )
    {
        if( number_of_populations < maximum_number_of_populations )
        {
            initializeNewPopulation();
            if( total_number_of_generations == 0 && write_generational_statistics )
                writeGenerationalStatisticsForOnePopulation( number_of_populations-1 );

            if( total_number_of_generations == 0 && write_generational_solutions )
                writeGenerationalSolutions( 0 );
        }

        generationalStepAllPopulations();

        if( write_generational_statistics )
            writeGenerationalStatisticsForOnePopulation( number_of_populations-1 );

        if( write_generational_solutions )
            writeGenerationalSolutions( 0 );

        total_number_of_generations++;
    }
}

/**
 * Runs the IDEA.
 */
void run( void )
{
    initialize();

    if( print_verbose_overview )
        printVerboseOverview();

    runAllPopulations();

    printf("evals %f ", number_of_evaluations);

    printf("obj_val %6.2e ", elitist_objective_value);

    ezilaitini();

    printf("time %lf ", getTimer());
    printf("generations %d\n", total_number_of_generations);
}

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
