//Micro-benchmark: cost of one LocalSearch random-walk step
//(one contiguous() over |J| facets + two d() evaluations with Dijkstra),
//which is what the C++ code pays 40,000 times per AddFacet attempt.
#include <iostream>
#include <ctime>
#include "SimplexAbstract.hpp"
#include "SimplexAlpha.hpp"

using namespace std;

Complex komplex;
ComplexProduct KxK;
SubComplexJ L;
SimplicialMap map1, map0, probe;

int main() {
    komplex.initComplex(4, 4);
    komplex.K.A[0][0].initSimplex(3, 0, 1, 2);
    komplex.K.A[0][1].initSimplex(3, 0, 1, 3);
    komplex.K.A[0][2].initSimplex(3, 0, 2, 3);
    komplex.K.A[0][3].initSimplex(3, 1, 2, 3);

    KxK.initComplexProduct(komplex);
    map1.projection1(KxK);
    map0.projection2(KxK);
    komplex.initAdjMat();

    L.initSubComplexJ(96);
    int counting = 0;
    for (int i = 0; i < KxK.listOfFacets.m; i++)
        for (int j = 0; j < KxK.listOfFacets.rowLength.getA(0, i); j++) {
            L.initA(counting, KxK.listOfFacets.A[i][j]);
            counting += 1;
        }
    L.initZero_Skeleton();

    probe.initSimplicialMapCopy(map1);
    //mutate the probe map so the union test is non-trivial
    probe.updateSimplicialMapImageA(0, 1, 2);

    struct timespec t0, t1;
    int N = 2000;

    //1) contiguous() over all 96 facets of L
    clock_gettime(CLOCK_MONOTONIC, &t0);
    for (int i = 0; i < N; i++)
        probe.contiguous(map0, L, komplex);
    clock_gettime(CLOCK_MONOTONIC, &t1);
    double us = ((t1.tv_sec - t0.tv_sec) * 1e9 + (t1.tv_nsec - t0.tv_nsec)) * 1e-3 / N;
    printf("contiguous(pi2, L=all 96 facets): %.1f us/call\n", us);

    //2) d() with 16 Dijkstra runs
    clock_gettime(CLOCK_MONOTONIC, &t0);
    for (int i = 0; i < N; i++)
        probe.d(map0, L, komplex);
    clock_gettime(CLOCK_MONOTONIC, &t1);
    us = ((t1.tv_sec - t0.tv_sec) * 1e9 + (t1.tv_nsec - t0.tv_nsec)) * 1e-3 / N;
    printf("d(pi2, L) [16 Dijkstra runs]:     %.1f us/call\n", us);

    //3) estimated cost of one failed AddFacet attempt, M = 40000
    double step = us * 2.0 + 200.0; //d() is evaluated twice per accepted check; contiguous dominates
    printf("\nM=40000-step walk over |J|=96:    ~%.1f s per AddFacet attempt\n",
           40000 * (us + 2 * us) / 1e6);
    return 0;
}