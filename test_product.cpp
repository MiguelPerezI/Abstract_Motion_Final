//Verification harness: builds K x K through the same pipeline as main.cpp
//for K = dDelta^2 (the paper's circle, Section 5.1) and K = dDelta^3 (S^2).
#include <iostream>
#include "SimplexAbstract.hpp"
#include "SimplexAlpha.hpp"

using namespace std;

Complex komplex;
ComplexProduct KxK;
SimplicialMap map1, map0;

void testComplex(int numFacets, int maxInt) {
    komplex.initComplex(numFacets, maxInt);
    if (numFacets == 3) {
        komplex.K.A[0][0].initSimplex(2, 0, 1);
        komplex.K.A[0][1].initSimplex(2, 0, 2);
        komplex.K.A[0][2].initSimplex(2, 1, 2);
    } else {
        komplex.K.A[0][0].initSimplex(3, 0, 1, 2);
        komplex.K.A[0][1].initSimplex(3, 0, 1, 3);
        komplex.K.A[0][2].initSimplex(3, 0, 2, 3);
        komplex.K.A[0][3].initSimplex(3, 1, 2, 3);
    }

    komplex.setN_Vertices();
    printf("K: %d facets, %d vertices\n", komplex.numSimplex, komplex.n);

    KxK.initComplexProduct(komplex);

    int sum = 0;
    for (int i = 0; i < KxK.listOfFacets.m; i++)
        for (int j = 0; j < KxK.listOfFacets.rowLength.getA(0, i); j++)
            sum += 1;
    printf("KxK: %d maximal facets, zero skeleton %dx%d\n\n",
           sum, KxK.zero_skeleton.getM(), KxK.zero_skeleton.getN());

    map1.projection1(KxK);
    map0.projection2(KxK);
}

int main() {
    printf("=== Circle (paper Section 5.1, K = dDelta^2) ===\n");
    testComplex(3, 2);

    printf("=== S^2 (repo's current target, K = dDelta^3) ===\n");
    testComplex(4, 3);

    return 0;
}