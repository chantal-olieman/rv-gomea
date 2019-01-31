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
#include "CEC_benchmark.h"
//#include "../temp-dir/aaa_c_connector.h"
//#include "../cpp/F14.h"
//#include "../cpp/Header.h"
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/

////
////int sign(double x)
////{
////    if (x > 0) return 1;
////    if (x < 0) return -1;
////    return 0;
////}
////
////double hat(double x)
////{
////    if (x==0) return 0; else return log(abs(x));
////}
////
////double c1(double x)
////{
////    if (x>0) return 10; else return 5.5;
////}
////
////double c2(double x)
////{
////    if (x>0) return 7.9; else return 3.1;
////}
//
//
//int sign(double x)
//{
//    if (x > 0) return 1;
//    if (x < 0) return -1;
//    return 0;
//}
//
//double hat(double x)
//{
//    if (x==0)
//    {
//        return 0;
//    }
//    else
//    {
//        return log(abs(x));
//    }
//}
//
//double c1(double x)
//{
//    if (x>0)
//    {
//        return 10;
//    }
//    else
//    {
//        return 5.5;
//    }
//}
//
//double c2(double x)
//{
//    if (x>0)
//    {
//        return 7.9;
//    }
//    else
//    {
//        return 3.1;
//    }
//}
//
//
//double** readOvectorVec( int *s)
//{
//    // read O vector from file in csv format, seperated by s_size groups
//    double** d = (double**) malloc(s_size*sizeof(double*));
//    FILE *fpt;
//    char FileName[300];
//    int vector_lenght = 1000;
//    sprintf(FileName, "../cdatafiles/F%d-xopt.txt", ID);
//    fpt = fopen(FileName,"r");
//
//    string value;
//    string line;
//    int c = 0;                      // index over 1 to dim
//    int i = -1;                      // index over 1 to s_size
//    int up = 0;                   // current upper bound for one group
//    double read;
//    if (fpt==NULL)
//    {
//        while ( fscanf(fpt,"%lf",read))
//        {
//            if (c==up)             // out (start) of one group
//            {
//                // printf("=\n");
//                i++;
//                d[i] =  (double*) malloc(s[i]*sizeof(double));
//                up += s[i];
//            }
//            // printf("c=%d\ts=%d\ti=%d\tup=%d\tindex=%d\n",c,s[i],i,up,c-(up-s[i]));
//            d[i][c-(up-s[i])] = read;
//            // printf("1\n");
//            c++;
//            // printf("2\n");
//        }
//        file.close();
//    }
//    else
//    {
//        cout<<"Cannot open datafiles"<<endl;
//    }
//    fclose(fpt);
//    return d;
//}
//
//
//void transform_osz(double* z, int dim)
//{
//    // apply osz transformation to z
//    for (int i = 0; i < dim; ++i)
//    {
//        z[i] = sign(z[i]) * exp( hat(z[i]) + 0.049 * ( sin( c1(z[i]) * hat(z[i]) ) + sin( c2(z[i])* hat(z[i]) )  ) ) ;
//    }
//}
//
//void transform_asy(double* z, double beta, int dim)
//{
//    for (int i = 0; i < dim; ++i)
//    {
//        if (z[i]>0)
//        {
//            z[i] = pow(z[i], 1 + beta * i/((double) (dim-1)) * sqrt(z[i]) );
//        }
//    }
//}
//
//
//// for single group non-separable function
//double schwefel( double *x , int dim ){
//    int    j;
//    double s1 = 0;
//    double s2 = 0;
//
//    // T_{osz}
//    transform_osz(x,dim);
//
//    // T_{asy}^{0.2}
//    transform_asy(x, 0.2, dim);
//
//    for (j = 0; j < dim; j++) {
//        s1 += x[j];
//        s2 += (s1 * s1);
//    }
//
//    return(s2);
//}
//
//int* readS(int num)
//{
//    s =  (int*)Malloc(sizeof(int)* num);
//
//    sprintf(FileName, "../cdatafiles/F%d-s.txt", ID);
//    fpt = fopen(FileName,"r");
//    int read;
//    int c = 0;                      // index over for all s
//    double read;
//    if (fpt==NULL)
//    {
//        while ( fscanf(fpt,"%d",read))
//        {
//            s[c++] = read;
//        }
//        file.close();
//    }
//    else
//    {
//        cout<<"Cannot open datafiles"<<endl;
//    }
//    fclose(fpt);
//
//    return s;
//}
//
//double F14(double*x){
//    OvectorVec = NULL;
//    minX = -100;
//    maxX = 100;
//    ID = 14;
//    s_size = 20;
//    dimension = 905;        // because of overlapping
//    overlap = 5;
//    int i;
//    double result=0.0;
//
//    if(OvectorVec == NULL)
//    {
//        s = readS(s_size);
//        OvectorVec = readOvectorVec(s);
//        // // inspect OvectorVec
//        // for (int i = 0; i < s_size; ++i)
//        //   {
//        //     for (int j=0; j< s[i]; j++)
//        //       {
//        //         printf("%.1f\t",OvectorVec[i][j]);
//        //       }
//        //     printf("\n");
//        //   }
//        Pvector = readPermVector();
//        r25 = readR(25);
//        r50 = readR(50);
//        r100 = readR(100);
//        w = readW(s_size);
//    }
//
//    // s_size non-separable part with rotation
//    int c = 0;
//    for (i = 0; i < s_size; i++)
//    {
//        // cout<<"c="<<c<<", i="<<i<<endl;
//        anotherz1 = rotateVectorConflict(i, c, x);
//        // cout<<"done rot"<<endl;
//        result += w[i] * schwefel(anotherz1, s[i]);
//        delete []anotherz1;
//        // cout<<result<<endl;
//    }
//
//    return(result);
//}
//
//
//
//
void CECProblemEvaluation( double *parameters, double *objective_value, double *constraint_value, int problem_index ){
    for(int i = 0; i < 1000; i++) if(parameters[i] < -100) parameters[i] = -100;
    for(int i = 0; i < 1000; i++) if(parameters[i] > 100) parameters[i] = 100;
//    double result = 0.0;
    double result = Evaluate(parameters, problem_index);

    *objective_value = result;
    *constraint_value = 0.0;
}


double CECFunctionLowerRangeBound( int dimension )
{
    return( -1e+308 );
}

double CECFunctionUpperRangeBound( int dimension )
{
    return( 1e+308 );
}
