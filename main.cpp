#include <iostream>
#include <cstdlib>
#include "SimplexAbstract.hpp"
#include <math.h>
#include <list>
#include <time.h>
#include "SimplexAlpha.hpp"
#include "Covering_CORE.hpp"
#include "OptimizedCovering.hpp"
#include "Lex.hpp"
#include <mpi.h>
#include <stdio.h>
#include <cstdio>
#include <iomanip>     
#include <functional>

using namespace std;

Complex komplex;
VectorInt v0, v1, v2, v3, v4;
MatrixVectorInt matV;
SimplexAlpha simplex, simplex0;
MatrixSimplexAlpha com;

Simplex s0, s1, s2, s3, s4;
SimplexProd subComplex;

SimplicialMap map0, copyMap, map1, map2;
ComplexProduct KxK;
SubComplexJ L, JsubL, JsubL1, JsubL2;
LocalSearch suchen, suchenCore1, suchen2; 
MatrixIntList ll;
MatrixInt Test;
SimplexAlpha testSimplex;

Covering c;
OptimizedCovering O;

double fRand(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

int main(int argc, char **argv) {

    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);
    ios_base::sync_with_stdio(false);


        srand(time(NULL));

        int M = 40000;
        double r = 0.1;
        int Ns = 100;

        int maxInt = 3;
        int numOfMaxSimplex = 4;
        komplex.initComplex(numOfMaxSimplex, maxInt + 1);

        komplex.K.A[0][0].initSimplex(3, 0, 1, 2);
        komplex.K.A[0][1].initSimplex(3, 0, 1, 3);
        komplex.K.A[0][2].initSimplex(3, 0, 2, 3);
        komplex.K.A[0][3].initSimplex(3, 1, 2, 3);

        KxK.initComplexProduct(komplex);
        KxK.escMaximalSimplices();
        map1.projection1(KxK);
        map0.projection2(KxK);

        L.initSubComplexJ(96);
        int counting = 0;

        for (int i = 0; i < KxK.listOfFacets.m; i++) {
                for (int j = 0; j < KxK.listOfFacets.rowLength.getA(0, i); j++) {
                        L.initA(counting, KxK.listOfFacets.A[i][j]);
                        counting += 1;
                }
        }

	 L.initZero_Skeleton();
	
         komplex.initAdjMat();
         komplex.graph.addWeight(0, 1, 1);
         komplex.graph.addWeight(0, 2, 1);
         komplex.graph.addWeight(0, 3, 1);
         komplex.graph.addWeight(1, 2, 1);
         komplex.graph.addWeight(1, 3, 1);
         komplex.graph.addWeight(2, 3, 1);


        int ccount = 0;

        O.initOptimizedCovering_CORE(L, komplex, map1, map0, M, r,  argc, argv);
        O.runOptimizedCovering_CORE(L, komplex, map1, map0, M, Ns, r, argc, argv);


        cout << "\n\n\n";
        O.closeWriterO();

    
//    clock_gettime(CLOCK_MONOTONIC, &end); 
//    double time_taken;
//    time_taken = (end.tv_sec - start.tv_sec) * 1e9;
//    time_taken = (time_taken + (end.tv_nsec - start.tv_nsec)) * 1e-9;
//  
//    cout << "El tiempo que el programa tomó es: " << fixed
//         << time_taken << setprecision(9);
//    cout << " sec" << endl;
    return 0;
}
