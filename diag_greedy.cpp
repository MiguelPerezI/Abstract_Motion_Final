//diag_greedy.cpp — Empirical diagnosis on the REPO's own code path.
//Replays Covering's first RCC cycles (AddFacet over the 96 tetrahedra of
//dDelta^3 x dDelta^3 with the repo's ghostLocalSearch at the real M and r)
//and logs, for every simplex, how often adding it to a random maximal J
//succeeded vs failed. Answers: which (sigma,tau) pairs does the randomized
//greedy actually fail to place, and how big do its domains stay?
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <vector>
#include "SimplexAbstract.hpp"
#include "SimplexAlpha.hpp"
#include "LocalSearch.hpp"

using namespace std;

Complex komplex;
ComplexProduct KxK;
SubComplexJ L, J;
SimplicialMap map1, map0;
LocalSearch suchen;

int main(int argc, char **argv) {
    int M = (argc > 1) ? atoi(argv[1]) : 4000;     // walk bound per attempt
    int cycles = (argc > 2) ? atoi(argv[2]) : 30;  // number of RCC cycles

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
    suchen.initLocalSearch(L, komplex, map1, map0, M, 0.1);

    // product-pair id of each of the 96 tetrahedra: tetra index -> (i,j)
    vector<pair<int,int>> tetPair;
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            for (int k = 0; k < 6; k++)
                tetPair.push_back({i, j});

    vector<int> success(96, 0), fail(96, 0);

    srand(time(NULL));
    for (int c = 0; c < cycles; c++) {
        // random maximal starting subcomplex J (single random tetrahedron)
        J.initSubComplexJ(1);
        J.initA(0, L.listOfFacets.A[0][0]);   // initialize the slot's inner matrix
        J.initZero_Skeleton();
        J.resetSubComplexJ(L.listOfFacets.A[0][rand() % 96]);

        // try to add each candidate (AddFacet semantics, one at a time)
        for (int attempt = 0; attempt < 96; attempt++) {
            int t = rand() % 96;
            // skip if already in J
            bool inJ = false;
            for (int f = 0; f < J.listOfFacets.n; f++)
                if (J.listOfFacets.A[0][f].compareSimplexAlpha(L.listOfFacets.A[0][t]) == 1) { inJ = true; break; }
            if (inJ) continue;

            J.pushSimplexAlpha(L.listOfFacets.A[0][t]);
            suchen.ghostLocalSearch(J, komplex, map1, map0, M, 0.1);
            if (suchen.ghost == 420) success[t] += 1;
            else { fail[t] += 1; J.popSimplexAlpha(); }
        }
        if ((c + 1) % 5 == 0)
            printf("cycle %d/%d done\n", c + 1, cycles);
    }

    // aggregate per product pair (sigma_i x sigma_j)
    vector<int> sucPair(16, 0), faiPair(16, 0);
    for (int t = 0; t < 96; t++) {
        sucPair[tetPair[t].first * 4 + tetPair[t].second] += success[t];
        faiPair[tetPair[t].first * 4 + tetPair[t].second] += fail[t];
    }

    printf("\n=== success / failure per tetrahedron (by product pair) ===\n");
    printf("      ");
    for (int j = 0; j < 4; j++) printf("  s%d      ", j);
    printf("\n");
    for (int i = 0; i < 4; i++) {
        printf("s%d   ", i);
        for (int j = 0; j < 4; j++) {
            int p = i * 4 + j;
            printf("%3d/%-3d  ", sucPair[p], faiPair[p]);
        }
        printf("\n");
    }
    printf("\n=== individual tetrahedra with failures (tet : fail/success : pair) ===\n");
    for (int t = 0; t < 96; t++)
        if (fail[t] > 0)
            printf("  tet %2d: fail %3d  ok %3d   pair (s%d,s%d)\n",
                   t, fail[t], success[t], tetPair[t].first, tetPair[t].second);
    int totS = 0, totF = 0;
    for (int t = 0; t < 96; t++) { totS += success[t]; totF += fail[t]; }
    printf("\ntotal attempts: %d ok, %d failed (%.1f%% ok)\n",
           totS, totF, 100.0 * totS / (totS + totF));
    return 0;
}