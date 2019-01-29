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
#include "link.h"




void EvaluateCPP(){
//    double *temp_parameters = (double*)Malloc(sizeof(double)* number_of_parameters);
//    for(int i = 0; i< number_of_parameters; i++) temp_parameters[i] = parameters[i];
//    double result = F14(temp_parameters, number_of_parameters);
//    *objective_value = schwefel(temp_parameters, number_of_parameters);
//    Benchmarks *fp;

//    *objective_value = 0.4;
//    if( *objective_value )
//    printf("restul: %f \n", result);
//    *constraint_value = 0.0;
}

double EvaluateBenchmark(  double *parameters, int benchmark_index ) {
    if( !initialized ) {
        fp = generateFuncObj(benchmark_index - 21);
        initialized = 1;
        printf("Function index: %d\n", fp->getID());
    }
    double objective = fp->compute(parameters);
//    double res = parameters[0];
//    printf("One evaluation: %f  (i = %f)\n", objective, parameters[0]);
//    parameters[0] = res+100;
//    objective = fp->compute(parameters);
//    printf("Delta i: %f \n", objective);
//    parameters[0] = res;
//    parameters[1] = parameters[1]+100;
//    objective = fp->compute(parameters);
//    printf("Delta j: %f \n", objective);
//    parameters[0] = res+100;
//    objective = fp->compute(parameters);
//    printf("Delta ij: %f \n", objective);


    return objective;
}