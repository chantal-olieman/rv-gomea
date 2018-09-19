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

#ifndef SO_OPTIMIZATION_H
#define SO_OPTIMIZATION_H

/*-=-=-=-=-=-=-=-=-=-=-=-=-=-= Section Includes -=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
#include "Optimization.h"
#include "FOS.h"
#include "Tools.h"
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/



/*-=-=-=-=-=-=-=-=-=-=-=-= Section Header Functions -=-=-=-=-=-=-=-=-=-=-=-=*/
short betterFitness( double objective_value_x, double constraint_value_x, double objective_value_y, double constraint_value_y );
char *installedProblemName( int index );
int numberOfInstalledProblems( void );
void printAllInstalledProblems( void );
double installedProblemLowerRangeBound( int index, int dimension );
double installedProblemUpperRangeBound( int index, int dimension );
void initializeParameterRangeBounds( void );
short isParameterInRangeBounds( double parameter, int dimension );
void installedProblemEvaluation( int index, double *parameters, double *objective_value, double *constraint_value, int number_of_touched_parameters, int *touched_parameters_indices, double *parameters_before, double objective_value_before, double constraint_value_before );
void installedProblemEvaluationWithoutRotation( int index, double *parameters, double *objective_value, double *constraint_value, int number_of_touched_parameters, int *touched_parameters_indices, double *touched_parameters, double *parameters_before, double objective_value_before, double constraint_value_before );
void sphereFunctionProblemEvaluation( double *parameters, double *objective_value, double *constraint_value );
void sphereFunctionPartialProblemEvaluation( double *parameters, double *objective_value, double *constraint_value, int number_of_touched_parameters, int *touched_parameters_indices, double *touched_parameters, double *parameters_before, double objective_value_before, double constraint_value_before );
double sphereFunctionProblemLowerRangeBound( int dimension );
double sphereFunctionProblemUpperRangeBound( int dimension );
void ellipsoidFunctionProblemEvaluation( double *parameters, double *objective_value, double *constraint_value );
void ellipsoidFunctionPartialProblemEvaluation( double *parameters, double *objective_value, double *constraint_value, int number_of_touched_parameters, int *touched_parameters_indices, double *touched_parameters, double *parameters_before, double objective_value_before, double constraint_value_before );
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
void rosenbrockFunctionPartialProblemEvaluation( double *parameters, double *objective_value, double *constraint_value, int number_of_touched_parameters, int *touched_parameters_indices, double *touched_parameters, double *parameters_before, double objective_value_before, double constraint_value_before );
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
void michalewiczFunctionPartialProblemEvaluation( double *parameters, double *objective_value, double *constraint_value, int number_of_touched_parameters, int *touched_parameters_indices, double *touched_parameters, double *parameters_before, double objective_value_before, double constraint_value_before );
double michalewiczFunctionLowerRangeBound( int dimension );
double michalewiczFunctionUpperRangeBound( int dimension );
void rastriginFunctionProblemEvaluation( double *parameters, double *objective_value, double *constraint_value );
void rastriginFunctionPartialProblemEvaluation( double *parameters, double *objective_value, double *constraint_value, int number_of_touched_parameters, int *touched_parameters_indices, double *touched_parameters, double *parameters_before, double objective_value_before, double constraint_value_before );
double rastriginFunctionLowerRangeBound( int dimension );
double rastriginFunctionUpperRangeBound( int dimension );
void sumOfEllipsoidsFunctionProblemEvaluation( double *parameters, double *objective_value, double *constraint_value );
void sumOfEllipsoidsFunctionPartialProblemEvaluation( double *parameters, double *objective_value, double *constraint_value, int number_of_touched_parameters, int *touched_parameters_indices, double *touched_parameters, double *parameters_before, double objective_value_before, double constraint_value_before );
double sumOfEllipsoidsFunctionLowerRangeBound( int dimension );
double sumOfEllipsoidsFunctionUpperRangeBound( int dimension );
void ciasBRFunctionProblemEvaluation( double *parameters, double *objective_value, double *constraint_value );
double ciasBRFunctionLowerRangeBound( int dimension );
double ciasBRFunctionUpperRangeBound( int dimension );
void trapSphereFunctionProblemEvaluation( double *parameters, double *objective_value, double *constraint_value );
double trapSphereFunctionLowerRangeBound( int dimension );
double trapSphereFunctionUpperRangeBound( int dimension );
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/

/*-=-=-=-=-=-=-=-=-=-=-=- Section Global Variables -=-=-=-=-=-=-=-=-=-=-=-=-*/
double elitist_objective_value,
       elitist_constraint_value;
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/

#endif
