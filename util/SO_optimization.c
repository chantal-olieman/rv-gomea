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
#include "SO_optimization.h"
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/

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
    case 14: return( (char *) "Circles in a Square" );
    case 15: return( (char *) "TrapSphere" );
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
    case 14: return( ciasBRFunctionLowerRangeBound( dimension ) );
    case 15: return( trapSphereFunctionLowerRangeBound( dimension ) );
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
    case 14: return( ciasBRFunctionUpperRangeBound( dimension ) );
    case 15: return( trapSphereFunctionUpperRangeBound( dimension ) );
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
        if( lower_user_range > upper_range_bounds[i] )
            lower_init_ranges[i] = lower_range_bounds[i];

        upper_init_ranges[i] = upper_user_range;
        if( upper_user_range > upper_range_bounds[i] )
            upper_init_ranges[i] = upper_range_bounds[i];
        if( upper_user_range < lower_range_bounds[i] )
            upper_init_ranges[i] = upper_range_bounds[i];
    }
}

/**
 * Returns the value of the single objective
 * and the sum of all constraint violations
 * function after rotating the parameter vector.
 * Both are returned using pointer variables.
 * Number of evaluations is increased by the
 * ratio ([0..1]) of new parameters that have
 * been changed.
 */
void installedProblemEvaluation( int index, double *parameters, double *objective_value, double *constraint_value, int number_of_touched_parameters, int *touched_parameters_indices, double *parameters_before, double objective_value_before, double constraint_value_before )
{
    double *rotated_parameters, *rotated_parameters_before, *touched_parameters, *touched_parameters_before, last_obj, last_cons;
    int i, j, c, prev_block, cur_block, *block_indices;

    touched_parameters = NULL;
    rotated_parameters = NULL;
    if( !(touched_parameters_indices == NULL || black_box_evaluations) )
    {
        touched_parameters = (double*) Malloc( number_of_touched_parameters*sizeof( double ) );
        for( i = 0; i < number_of_touched_parameters; i++ )
            touched_parameters[i] = parameters[touched_parameters_indices[i]];
    }

    if( rotation_angle == 0.0 )
    {
        if( touched_parameters_indices == NULL || black_box_evaluations ) installedProblemEvaluationWithoutRotation( index, parameters, objective_value, constraint_value, number_of_parameters, NULL, NULL, NULL, 0, 0 );
        else installedProblemEvaluationWithoutRotation( index, parameters, objective_value, constraint_value, number_of_touched_parameters, touched_parameters_indices, touched_parameters, parameters_before, objective_value_before, constraint_value_before );
    }
    else
    {
        if( black_box_evaluations || touched_parameters_indices == NULL  )
        {
            rotated_parameters = rotateAllParameters( parameters );
            installedProblemEvaluationWithoutRotation( index, rotated_parameters, objective_value, constraint_value, number_of_parameters, NULL, NULL, NULL, 0, 0 );
        }
        else
        {
            if( problem_index == 13 && (learn_linkage_tree || static_linkage_tree || use_univariate_FOS) ) // indices must be sorted
            {
                if( touched_parameters_indices != NULL )
                    free( touched_parameters );
                prev_block = -1;
                touched_parameters = (double*) Malloc( block_size*sizeof( double ) );
                touched_parameters_before = (double*) Malloc( block_size*sizeof( double ) );
                block_indices = (int*) Malloc( block_size*sizeof( int ) );
                last_obj = objective_value_before;
                last_cons = constraint_value_before;
                for( i = 0; i < number_of_touched_parameters; i++ )
                {
                    cur_block = touched_parameters_indices[i]/block_size;
                    if( cur_block != prev_block )
                    {
                        if( rotated_parameters != NULL ) free( rotated_parameters );

                        for( j = 0; j < block_size; j++ )
                        {
                            block_indices[j] = cur_block*block_size+j;
                            touched_parameters[j] = parameters[cur_block*block_size+j];
                            touched_parameters_before[j] = parameters[cur_block*block_size+j];
                        }
                        c = 0;
                        while( i+c < number_of_touched_parameters && touched_parameters_indices[i+c]/block_size == cur_block )
                        {
                            touched_parameters_before[touched_parameters_indices[i+c]%block_size] = parameters_before[i+c];
                            c++;
                        }

                        rotated_parameters = matrixVectorMultiplication( rotation_matrix, touched_parameters, block_size, block_size );
                        rotated_parameters_before = matrixVectorMultiplication( rotation_matrix, touched_parameters_before, block_size, block_size );
                        installedProblemEvaluationWithoutRotation( index, parameters, objective_value, constraint_value, block_size, block_indices, rotated_parameters, rotated_parameters_before, last_obj, last_cons );
                        last_obj = *objective_value;
                        last_cons = *constraint_value;
                        free( rotated_parameters_before );
                    }
                    prev_block = cur_block;
                }
                free( block_indices );
                free( touched_parameters_before );
            }
            else
            {
                rotated_parameters = matrixVectorMultiplication( rotation_matrix, touched_parameters, number_of_touched_parameters, number_of_touched_parameters );
                rotated_parameters_before = matrixVectorMultiplication( rotation_matrix, parameters_before, number_of_touched_parameters, number_of_touched_parameters );
                installedProblemEvaluationWithoutRotation( index, parameters, objective_value, constraint_value, number_of_touched_parameters, touched_parameters_indices, rotated_parameters, rotated_parameters_before, objective_value_before, constraint_value_before );
                free( rotated_parameters_before );
            }
        }

        free( rotated_parameters );
    }

    if( use_vtr && !vtr_hit_status && *constraint_value == 0 && *objective_value <= vtr  )
    {
        if( touched_parameters_indices != NULL )
            installedProblemEvaluation( index, parameters, objective_value, constraint_value, number_of_parameters, NULL, NULL, 0, 0 );
        if( *constraint_value == 0 && *objective_value <= vtr  )
        {
            vtr_hit_status = 1;
            elitist_objective_value = *objective_value;
            elitist_constraint_value = *constraint_value;
        }
    }

    if( !vtr_hit_status && betterFitness(*objective_value, *constraint_value, elitist_objective_value, elitist_constraint_value) )
    {
        elitist_objective_value = *objective_value;
        elitist_constraint_value = *constraint_value;
    }

    if( touched_parameters_indices != NULL )
        free( touched_parameters );
}

/**
 * Returns the value of the single objective
 * and the sum of all constraint violations
 * without rotating the parameter vector.
 * Both are returned using pointer variables.
 */
void installedProblemEvaluationWithoutRotation( int index, double *parameters, double *objective_value, double *constraint_value, int number_of_touched_parameters, int *touched_parameters_indices, double *touched_parameters, double *parameters_before, double objective_value_before, double constraint_value_before )
{
    *objective_value  = 0.0;
    *constraint_value = 0.0;
    if( black_box_evaluations || touched_parameters_indices == NULL )
    {
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
        case 14: ciasBRFunctionProblemEvaluation( parameters, objective_value, constraint_value ); break;
        case 15: trapSphereFunctionProblemEvaluation( parameters, objective_value, constraint_value ); break;
        }
        number_of_evaluations++;
    }
    else
    {
        switch( index )
        {
        case  0: sphereFunctionPartialProblemEvaluation( parameters, objective_value, constraint_value, number_of_touched_parameters, touched_parameters_indices, touched_parameters, parameters_before, objective_value_before, constraint_value_before ); break;
        case  1: ellipsoidFunctionPartialProblemEvaluation( parameters, objective_value, constraint_value, number_of_touched_parameters, touched_parameters_indices, touched_parameters, parameters_before, objective_value_before, constraint_value_before ); break;
        case  2: cigarFunctionProblemEvaluation( parameters, objective_value, constraint_value ); break;
        case  3: tabletFunctionProblemEvaluation( parameters, objective_value, constraint_value ); break;
        case  4: cigarTabletFunctionProblemEvaluation( parameters, objective_value, constraint_value ); break;
        case  5: twoAxesFunctionProblemEvaluation( parameters, objective_value, constraint_value ); break;
        case  6: differentPowersFunctionProblemEvaluation( parameters, objective_value, constraint_value ); break;
        case  7: rosenbrockFunctionPartialProblemEvaluation( parameters, objective_value, constraint_value, number_of_touched_parameters, touched_parameters_indices, touched_parameters, parameters_before, objective_value_before, constraint_value_before ); break;
        case  8: parabolicRidgeFunctionProblemEvaluation( parameters, objective_value, constraint_value ); break;
        case  9: sharpRidgeFunctionProblemEvaluation( parameters, objective_value, constraint_value ); break;
        case 10: griewankFunctionProblemEvaluation( parameters, objective_value, constraint_value ); break;
        case 11: michalewiczFunctionPartialProblemEvaluation( parameters, objective_value, constraint_value, number_of_touched_parameters, touched_parameters_indices, touched_parameters, parameters_before, objective_value_before, constraint_value_before ); break;
        case 12: rastriginFunctionPartialProblemEvaluation( parameters, objective_value, constraint_value, number_of_touched_parameters, touched_parameters_indices, touched_parameters, parameters_before, objective_value_before, constraint_value_before ); break;
        case 13: sumOfEllipsoidsFunctionPartialProblemEvaluation( parameters, objective_value, constraint_value, number_of_touched_parameters, touched_parameters_indices, touched_parameters, parameters_before, objective_value_before, constraint_value_before ); break;
        case 14: ciasBRFunctionProblemEvaluation( parameters, objective_value, constraint_value ); break;
        case 15: trapSphereFunctionProblemEvaluation( parameters, objective_value, constraint_value ); break;
        }
        number_of_evaluations += number_of_touched_parameters/(double)number_of_parameters;
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

void sphereFunctionPartialProblemEvaluation( double *parameters, double *objective_value, double *constraint_value, int number_of_touched_parameters, int *touched_parameters_indices, double *touched_parameters, double *parameters_before, double objective_value_before, double constraint_value_before )
{
    int    i;
    double result;

    result = objective_value_before;
    for( i = 0; i < number_of_touched_parameters; i++ )
    {
        result += touched_parameters[i]*touched_parameters[i];
        result -= parameters_before[i]*parameters_before[i];
    }

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

void ellipsoidFunctionPartialProblemEvaluation( double *parameters, double *objective_value, double *constraint_value, int number_of_touched_parameters, int *touched_parameters_indices, double *touched_parameters, double *parameters_before, double objective_value_before, double constraint_value_before )
{
    int    i;
    double result;

    result = objective_value_before;
    for( i = 0; i < number_of_touched_parameters; i++ )
    {
        result += pow( 10.0, 6.0*(((double) (touched_parameters_indices[i]))/((double) (number_of_parameters-1))) )*touched_parameters[i]*touched_parameters[i];
        result -= pow( 10.0, 6.0*(((double) (touched_parameters_indices[i]))/((double) (number_of_parameters-1))) )*parameters_before[i]*parameters_before[i];
    }

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

void rosenbrockFunctionPartialProblemEvaluation( double *parameters, double *objective_value, double *constraint_value, int number_of_touched_parameters, int *touched_parameters_indices, double *touched_parameters, double *parameters_before, double objective_value_before, double constraint_value_before )
{
    int    i, j;
    double result;

    result = objective_value_before;
    if( number_of_touched_parameters == 1 )
    {
        i = touched_parameters_indices[0];
        if( i > 0 )
        {
            result += 100*(parameters[i]-parameters[i-1]*parameters[i-1])*(parameters[i]-parameters[i-1]*parameters[i-1]) + (1.0-parameters[i-1])*(1.0-parameters[i-1]);
            result -= 100*(parameters_before[0]-parameters[i-1]*parameters[i-1])*(parameters_before[0]-parameters[i-1]*parameters[i-1]) + (1.0-parameters[i-1])*(1.0-parameters[i-1]);
        }
        if( i < number_of_parameters-1 )
        {
            result += 100*(parameters[i+1]-parameters[i]*parameters[i])*(parameters[i+1]-parameters[i]*parameters[i]) + (1.0-parameters[i])*(1.0-parameters[i]);
            result -= 100*(parameters[i+1]-parameters_before[0]*parameters_before[0])*(parameters[i+1]-parameters_before[0]*parameters_before[0]) + (1.0-parameters_before[0])*(1.0-parameters_before[0]);
        }
    }
    else if( number_of_touched_parameters == 2 && touched_parameters_indices[1]-touched_parameters_indices[0] == 1 )
    {
        i = touched_parameters_indices[0];
        j = touched_parameters_indices[1];

        result += 100*(parameters[j]-parameters[i]*parameters[i])*(parameters[j]-parameters[i]*parameters[i]) + (1.0-parameters[i])*(1.0-parameters[i]);
        result -= 100*(parameters_before[1]-parameters_before[0]*parameters_before[0])*(parameters_before[1]-parameters_before[0]*parameters_before[0]) + (1.0-parameters_before[0])*(1.0-parameters_before[0]);
        if( i > 0 )
        {
            result += 100*(parameters[i]-parameters[i-1]*parameters[i-1])*(parameters[i]-parameters[i-1]*parameters[i-1]) + (1.0-parameters[i-1])*(1.0-parameters[i-1]);
            result -= 100*(parameters_before[0]-parameters[i-1]*parameters[i-1])*(parameters_before[0]-parameters[i-1]*parameters[i-1]) + (1.0-parameters[i-1])*(1.0-parameters[i-1]);
        }
        if( j < number_of_parameters-1 )
        {
            result += 100*(parameters[j+1]-parameters[j]*parameters[j])*(parameters[j+1]-parameters[j]*parameters[j]) + (1.0-parameters[j])*(1.0-parameters[j]);
            result -= 100*(parameters[j+1]-parameters_before[1]*parameters_before[1])*(parameters[j+1]-parameters_before[1]*parameters_before[1]) + (1.0-parameters_before[1])*(1.0-parameters_before[1]);
        }
    }
    else
        rosenbrockFunctionProblemEvaluation( parameters, &result, constraint_value );

    *objective_value = result;
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

void michalewiczFunctionPartialProblemEvaluation( double *parameters, double *objective_value, double *constraint_value, int number_of_touched_parameters, int *touched_parameters_indices, double *touched_parameters, double *parameters_before, double objective_value_before, double constraint_value_before )
{
    int    i;
    double result;

    result  = objective_value_before;
    for( i = 0; i < number_of_touched_parameters; i++ )
    {
        result += -sin(touched_parameters[i])*pow(sin(((touched_parameters_indices[i]+1)*touched_parameters[i]*touched_parameters[i])/PI),20.0);
        result -= -sin(parameters_before[i])*pow(sin(((touched_parameters_indices[i]+1)*parameters_before[i]*parameters_before[i])/PI),20.0);
    }

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

void rastriginFunctionPartialProblemEvaluation( double *parameters, double *objective_value, double *constraint_value, int number_of_touched_parameters, int *touched_parameters_indices, double *touched_parameters, double *parameters_before, double objective_value_before, double constraint_value_before )
{
    int    i;
    double result;

    result  = objective_value_before;
    for( i = 0; i < number_of_touched_parameters; i++ )
    {
        result -= parameters_before[i]*parameters_before[i] - 10.0*cos(2.0*PI*parameters_before[i]);
        result += touched_parameters[i]*touched_parameters[i] - 10.0*cos(2.0*PI*touched_parameters[i]);
    }

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
        j = i % block_size;
        result += pow( 10.0, 6.0*(((double) (j))/((double) (block_size-1))) )*parameters[i]*parameters[i];
    }

    *objective_value  = result;
    *constraint_value = 0;
}

void sumOfEllipsoidsFunctionPartialProblemEvaluation( double *parameters, double *objective_value, double *constraint_value, int number_of_touched_parameters, int *touched_parameters_indices, double *touched_parameters, double *parameters_before, double objective_value_before, double constraint_value_before )
{
    int    i, j;
    double result;

    result = objective_value_before;

    for( i = 0; i < number_of_touched_parameters; i++ )
    {
        j = touched_parameters_indices[i] % block_size;
        result -= pow( 10.0, 6.0*(((double) (j))/((double) (block_size-1))) )*parameters_before[i]*parameters_before[i];
        result += pow( 10.0, 6.0*(((double) (j))/((double) (block_size-1))) )*touched_parameters[i]*touched_parameters[i];
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

void ciasBRFunctionProblemEvaluation( double *parameters, double *objective_value, double *constraint_value )
{
    int    i, j, nc;
    double result, xi0, xi1, xj0, xj1, distance;

    nc = number_of_parameters/2;

    for( i = 0; i < number_of_parameters; i++ )
    {
        if( parameters[i] < 0 )
            parameters[i] = 0;

        if( parameters[i] > 1 )
            parameters[i] = 1;
    }

    result = -1.0;
    for( i = 0; i < nc; i++ )
        for( j = i+1; j < nc; j++ )
        {
            xi0      = parameters[2*i];
            xi1      = parameters[2*i+1];
            xj0      = parameters[2*j];
            xj1      = parameters[2*j+1];
            distance = (xi0-xj0)*(xi0-xj0) + (xi1-xj1)*(xi1-xj1);
            if( result < 0 || distance < result )
                result = distance;
        }
    result = sqrt( result );

    *objective_value  = -result;
    *constraint_value = 0;
}

double ciasBRFunctionLowerRangeBound( int dimension )
{
    return( -1e+308 );
}

double ciasBRFunctionUpperRangeBound( int dimension )
{
    return( 1e+308 );
}

void trapSphereFunctionProblemEvaluation( double *parameters, double *objective_value, double *constraint_value )
{
    int    i, j, m, k, u;
    double result;

    *objective_value = 0;
    *constraint_value = 0;

    k = 5;
    m = (number_of_parameters / 2) / k;
    for( i = 0; i < m; i++ )
    {
        u = 0;
        for( j = 0; j < k; j++ )
            u += parameters[i*k+j]<0.5?0:1;

        if( u == k )
            result = 1.0;
        else
            result = ((double) (k-1-u))/((double) k);
        *objective_value += (1.0-result);
    }
    for( i = number_of_parameters/2; i < number_of_parameters; i++ )
    {
        *objective_value += parameters[i]*parameters[i];
    }
}

double trapSphereFunctionLowerRangeBound( int dimension )
{
    return( -1e+308 );
}

double trapSphereFunctionUpperRangeBound( int dimension )
{
    return( 1e+308 );
}

/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
