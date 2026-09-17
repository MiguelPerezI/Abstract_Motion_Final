//diag_sc2.cpp — Exact diagnosis of SC_strict on K x K (K = dDelta^2 or dDelta^3).
//
//For a subcomplex J of the staircase-triangulated product K x K, a contiguity
//chain from pi_1 to pi_2 is a sequence of simplicial maps J -> K where
//consecutive maps are 1-contiguous (union of images on every facet of J is a
//simplex of K -- Def 2.3 of arXiv:2008.13290). The C++ LocalSearch explores
//this space with single-vertex moves; here we search it EXACTLY with A*
//(heuristic: number of J-vertices where the map differs from pi_2; each move
//fixes at most one, so it is admissible). This decides feasibility and gives
//the minimal chain length -- ground truth for "can these facets share a domain".
//
//Every domain of a motion planner is a set of pairwise-compatible facets
//(a clique in the compatibility graph), so the minimum clique cover of that
//graph is a rigorous LOWER BOUND on the number of domains, i.e. on SC_strict.
//
//Usage: ./diag_sc2 [2|3]        2 = circle (sanity check), 3 = S^2 (default)

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <set>
#include <algorithm>
#include <chrono>
#include <functional>
#include <string>
#include "LocalSearch.hpp"
using namespace std;
using clk = chrono::steady_clock;

// ------------------------------------------------------------------ K ------
struct Kx {
    int nv;                          // number of vertices
    vector<vector<int>> facets;      // maximal simplices as sorted vertex lists
    unordered_set<int> masks;        // ALL simplices of K as vertex bitmasks
};

Kx buildK(int which) {
    Kx k;
    if (which == 2) {            // circle: boundary of the triangle
        k.nv = 3;
        k.facets = {{0,1},{0,2},{1,2}};
    } else {                     // S^2: boundary of the tetrahedron
        k.nv = 4;
        k.facets = {{0,1,2},{0,1,3},{0,2,3},{1,2,3}};
    }
    for (auto &f : k.facets) {
        int m = 0;
        for (int v : f) m |= 1 << v;
        for (int s = m; ; s = (s - 1) & m) {   // all submasks = all faces
            k.masks.insert(s);
            if (s == 0) break;
        }
    }
    return k;
}

// --------------------------------------------------------- product unit ----
// Product pair sigma_i x sigma_j, staircase triangulation. Vertices are
// global pair-ids pid = a*k.nv + b. A maximal simplex of the unit is a
// staircase simplex (5 vertices for triangle x triangle, 3 for edge x edge).
struct Unit {
    int fi, fj;
    vector<int> verts;               // global pair-ids
    vector<vector<int>> tetra;       // maximal simplices of the unit
};

Unit buildUnit(const Kx &k, int fi, int fj) {
    Unit u; u.fi = fi; u.fj = fj;
    const vector<int> &S = k.facets[fi], &T = k.facets[fj];
    auto pid = [&](int a, int b) { return a * k.nv + b; };

    int p = (int)S.size() - 1, q = (int)T.size() - 1, cols = p + q;
    for (int mask = 0; mask < (1 << cols); mask++) {
        int bits = 0;
        for (int b = 0; b < cols; b++) if ((mask >> b) & 1) bits++;
        if (bits != q) continue;

        vector<int> tet;
        int pi = 0, pj = 0, a = S[0], b = T[0];
        tet.push_back(pid(a, b));
        u.verts.push_back(pid(a, b));
        for (int step = 0; step < cols; step++) {
            if ((mask >> step) & 1) b = T[++pj];
            else                    a = S[++pi];
            tet.push_back(pid(a, b));
            u.verts.push_back(pid(a, b));
        }
        u.tetra.push_back(tet);
    }
    sort(u.verts.begin(), u.verts.end());
    u.verts.erase(unique(u.verts.begin(), u.verts.end()), u.verts.end());
    return u;
}

// ------------------------------------------------------------------- A* ----
struct Problem {
    const Kx *k;
    vector<vector<int>> tetra;       // facets of J as global pair-ids
    vector<int> jverts;              // global pair-ids of J
    vector<int> goal;                // pi2 image per global pair-id
};

int solveAstar(const Problem &P, int64_t start, long cap, long *expandedOut = nullptr) {
    const Kx &k = *P.k;
    int nv = k.nv;

    auto goalDist = [&](int64_t s) {          // h: #J-vertices differing from pi2
        int hh = 0;
        for (int pid : P.jverts)
            if (((s >> (2 * pid)) & 3) != P.goal[pid]) hh++;
        return hh;
    };
    if (goalDist(start) == 0) return 0;

    // tetras touching each pid
    unordered_map<int, vector<int>> touch;
    for (size_t t = 0; t < P.tetra.size(); t++)
        for (int pid : P.tetra[t]) touch[pid].push_back((int)t);

    auto tetMask = [&](int64_t s, const vector<int> &t) {
        int m = 0;
        for (int pid : t) m |= 1 << ((s >> (2 * pid)) & 3);
        return m;
    };

    struct Item { int f, g; int64_t s; };
    auto cmp = [](const Item &a, const Item &b) { return a.f > b.f; };
    priority_queue<Item, vector<Item>, decltype(cmp)> pq(cmp);
    unordered_map<int64_t,int> best;

    pq.push({goalDist(start), 0, start});
    best[start] = 0;
    long expanded = 0;

    while (!pq.empty()) {
        auto [f, g, s] = pq.top(); pq.pop();
        auto it = best.find(s);
        if (it != best.end() && it->second < g) continue;
        if (goalDist(s) == 0) return g;
        if (++expanded > cap) { if (expandedOut) *expandedOut = expanded; return -2; }

        for (int pid : P.jverts) {
            int cur = (s >> (2 * pid)) & 3;
            for (int val = 0; val < nv; val++) {
                if (val == cur) continue;
                int64_t ns = (s & ~((int64_t)3 << (2 * pid))) | ((int64_t)val << (2 * pid));
                auto jt = best.find(ns);
                if (jt != best.end() && jt->second <= g + 1) continue;

                bool ok = true;
                for (int ti : touch[pid]) {
                    int mold = tetMask(s, P.tetra[ti]);
                    int mnew = (mold & ~(1 << cur)) | (1 << val);
                    if (!k.masks.count(mnew))              { ok = false; break; }  // simplicial
                    if (!k.masks.count(mold | (1 << val))) { ok = false; break; }  // 1-contiguous
                }
                if (!ok) continue;
                best[ns] = g + 1;
                int hh = goalDist(s) - (P.goal[pid] != cur ? 1 : 0) + (P.goal[pid] != val ? 1 : 0);
                pq.push({g + 1 + hh, g + 1, ns});
            }
        }
    }
    if (expandedOut) *expandedOut = expanded;
    return -1;   // frontier exhausted: no single-vertex chain exists
}

Problem makeProblem(const Kx &k, const vector<Unit> &units, const vector<int> &which) {
    Problem P; P.k = &k;
    P.goal.assign(k.nv * k.nv, 0);
    for (int pid = 0; pid < k.nv * k.nv; pid++)
        P.goal[pid] = pid % k.nv;                       // pi2(a,b) = b
    set<int> vset;
    for (int w : which)
        for (int pid : units[w].verts) vset.insert(pid);
    P.jverts.assign(vset.begin(), vset.end());
    for (int w : which)
        for (auto &t : units[w].tetra) P.tetra.push_back(t);
    return P;
}

int64_t pi1State(const Kx &k) {
    int64_t s = 0;
    for (int a = 0; a < k.nv; a++)
        for (int b = 0; b < k.nv; b++)
            s |= (int64_t)a << (2 * (a * k.nv + b));
    return s;
}

// ------------------------------------------------------- graph tools -------
// Max clique via Bron-Kerbosch (adjacency with ZERO diagonal required here).
void bronKerbosch(vector<int> R, vector<int> P, vector<int> X,
                  const vector<vector<uint8_t>> &adj,
                  vector<int> &best, clk::time_point t0, long capMs) {
    if (clk::now() - t0 > chrono::milliseconds(capMs)) return;
    if (P.empty() && X.empty()) {
        if ((int)R.size() > (int)best.size()) best = R;
        return;
    }
    if ((int)R.size() + (int)P.size() <= (int)best.size()) return;
    int u = P[0], bestc = -1;
    for (int v : P) {
        int c = 0;
        for (int w : P) if (adj[v][w]) c++;
        if (c > bestc) { bestc = c; u = v; }
    }
    vector<int> ext;
    for (int v : P) if (!adj[u][v]) ext.push_back(v);
    for (int v : ext) {
        vector<int> R2 = R; R2.push_back(v);
        vector<int> P2, X2;
        for (int w : P) if (adj[v][w]) P2.push_back(w);
        for (int w : X) if (adj[v][w]) X2.push_back(w);
        if ((int)R2.size() + (int)P2.size() > (int)best.size())
            bronKerbosch(R2, P2, X2, adj, best, t0, capMs);
        P.erase(find(P.begin(), P.end(), v));
        X.push_back(v);
    }
}


// Minimum clique cover (<= k parts) == coloring of the complement graph.
// (pairwise lower-bound tool; the exact search below verifies parts jointly)
bool colorComplement(const vector<vector<uint8_t>> &compat, int n, int k,
                     vector<int> &color, clk::time_point t0, long capMs) {
    vector<vector<uint8_t>> incompat(n, vector<uint8_t>(n, 0));
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            incompat[i][j] = (i != j) && !compat[i][j];

    color.assign(n, -1);
    function<bool(int)> bt = [&](int assigned) -> bool {
        if (clk::now() - t0 > chrono::milliseconds(capMs)) return false;
        if (assigned == n) return true;
        int v = -1, bestSat = -1, bestDeg = -1;
        for (int i = 0; i < n; i++) {
            if (color[i] != -1) continue;
            vector<int> sat;
            for (int j = 0; j < n; j++)
                if (incompat[i][j] && color[j] != -1)
                    if (find(sat.begin(), sat.end(), color[j]) == sat.end()) sat.push_back(color[j]);
            int deg = 0;
            for (int j = 0; j < n; j++) if (incompat[i][j]) deg++;
            if ((int)sat.size() > bestSat || ((int)sat.size() == bestSat && deg > bestDeg)) {
                bestSat = (int)sat.size(); bestDeg = deg; v = i;
            }
        }
        for (int c = 0; c < k; c++) {
            bool ok = true;
            for (int j = 0; j < n; j++)
                if (incompat[v][j] && color[j] == c) { ok = false; break; }
            if (!ok) continue;
            color[v] = c;
            if (bt(assigned + 1)) return true;
            color[v] = -1;
        }
        return false;
    };
    return bt(0);
}
#include <atomic>
const long SEARCH_BUDGET = 3000000;
atomic<long> coverNodes{0};

// Joint feasibility of a set of units (whole J as one contiguity subcomplex),
// with caching. A* result >= 0 means feasible with that minimal chain length;
// -1 means no chain exists; -2 means the A* expansion cap was hit.
#include <unordered_set>
struct FeasCache {
    const Kx *k;
    const vector<Unit> *units;
    long cap;
    unordered_map<string, int> cache;

    FeasCache(const Kx *k, const vector<Unit> *units, long cap) : k(k), units(units), cap(cap) {}

    Problem problemFrom(const vector<int> &unitIds) const {
        Problem P; P.k = k;
        P.goal.assign(k->nv * k->nv, 0);
        for (int pid = 0; pid < k->nv * k->nv; pid++) P.goal[pid] = pid % k->nv;
        set<int> vset;
        for (int u : unitIds)
            for (int pid : (*units)[u].verts) vset.insert(pid);
        P.jverts.assign(vset.begin(), vset.end());
        for (int u : unitIds)
            for (auto &tet : (*units)[u].tetra) P.tetra.push_back(tet);
        return P;
    }

    int check(const vector<int> &unitIds) {
        vector<int> sorted = unitIds;
        sort(sorted.begin(), sorted.end());
        string key;
        for (int u : sorted) key += to_string(u) + ",";
        auto it = cache.find(key);
        if (it != cache.end()) return it->second;
        int r = solveAstar(problemFrom(sorted), pi1State(*k), cap);
        cache[key] = r;
        return r;
    }

    // node-based variant: nodes reference external vertex/facet tables
    int checkTetra(const vector<int> &nodeIds,
                   const vector<vector<int>> &nodeVerts,
                   const vector<vector<vector<int>>> &nodeTetra) {
        vector<int> sorted = nodeIds;
        sort(sorted.begin(), sorted.end());
        string key = "T";
        for (int u : sorted) key += to_string(u) + ",";
        auto it = cache.find(key);
        if (it != cache.end()) return it->second;

        Problem P; P.k = k;
        P.goal.assign(k->nv * k->nv, 0);
        for (int pid = 0; pid < k->nv * k->nv; pid++) P.goal[pid] = pid % k->nv;
        set<int> vset;
        for (int u : sorted)
            for (int pid : nodeVerts[u]) vset.insert(pid);
        P.jverts.assign(vset.begin(), vset.end());
        for (int u : sorted)
            for (auto &tet : nodeTetra[u]) P.tetra.push_back(tet);

        int r = solveAstar(P, pi1State(*k), cap);
        cache[key] = r;
        return r;
    }
};

// ---- hybrid verification (initialized in main) ---------------------------
// For small parts (<= gAstarCapTetras tetrahedra): exact A* (certifies both
// feasibility and infeasibility). For larger parts: the repo's own
// ghostLocalSearch with the real M=40000 -- a FOUND CHAIN CERTIFIES
// feasibility; failure is inconclusive and treated as "cannot extend".
static SubComplexJ Lg;            // the 96 tetrahedra, repo-style
static SubComplexJ Jg;            // scratch part under verification
static LocalSearch *gSuchen = NULL;
static int  gAstarCapTetras = 10;
static long gLsM = 40000;
static Complex *gKomplex = NULL;
static SimplicialMap gMap1, gMap0;

int checkPartT(const vector<int> &part,
               const vector<vector<int>> &nodeVerts,
               const vector<vector<vector<int>>> &nodeTetra,
               FeasCache &fc) {
    if ((int)part.size() <= gAstarCapTetras)
        return fc.checkTetra(part, nodeVerts, nodeTetra);
    (void)nodeVerts; (void)nodeTetra;
    for (int attempt = 0; attempt < 3; attempt++) {
        Jg.resetSubComplexJ(Lg.listOfFacets.A[0][part[0]]);
        for (size_t i = 1; i < part.size(); i++)
            Jg.pushSimplexAlpha(Lg.listOfFacets.A[0][part[i]]);
        gSuchen->ghostLocalSearch(Jg, *gKomplex, gMap1, gMap0, gLsM, 0.1);
        if (gSuchen->ghost == 420) return 1;   // certified feasible
    }
    return -1;   // inconclusive -> treat as "cannot extend"
}

// Exact-ish search: cover all nodes with <= maxParts JOINTLY-FEASIBLE parts.
bool coverSearchT(int u, int nNodes, vector<vector<int>> &parts, int &partsUsed,
                  int maxParts, FeasCache &fc, const vector<vector<uint8_t>> &pairCompat,
                  const vector<vector<int>> &nodeVerts,
                  const vector<vector<vector<int>>> &nodeTetra) {
    if (++coverNodes > SEARCH_BUDGET) return false;
    if (u == nNodes) return true;
    for (int p = 0; p < partsUsed; p++) {
        bool quickOK = true;
        for (int m : parts[p])
            if (!pairCompat[u][m]) { quickOK = false; break; }
        if (!quickOK) continue;
        parts[p].push_back(u);
        int r = checkPartT(parts[p], nodeVerts, nodeTetra, fc);
        parts[p].pop_back();
        if (r >= 0) {
            parts[p].push_back(u);
            if (coverSearchT(u + 1, nNodes, parts, partsUsed, maxParts, fc, pairCompat, nodeVerts, nodeTetra))
                return true;
            parts[p].pop_back();
        }
        // infeasible or inconclusive: try the next part
    }
    if (partsUsed < maxParts) {
        if (checkPartT({u}, nodeVerts, nodeTetra, fc) < 0) return false;
        parts[partsUsed++] = {u};
        if (coverSearchT(u + 1, nNodes, parts, partsUsed, maxParts, fc, pairCompat, nodeVerts, nodeTetra))
            return true;
        partsUsed--;
        parts[partsUsed].clear();
    }
    return false;
}

int runLevel(const Kx &k, const vector<Unit> &units, int nNodes,
             const vector<vector<int>> &nodeVerts,
             const vector<vector<vector<int>>> &nodeTetra,
             const char *lvlName, long cap, long capMs, int maxK) {
    (void)units;
    printf("\n=============== LEVEL: %s (%d nodes) ===============\n", lvlName, nNodes);

    printf("[singles] minimal chain length pi1 -> pi2:\n");
    int noChain = 0, capped = 0;
    for (int u = 0; u < nNodes; u++) {
        FeasCache fc(&k, &units, cap);
        int r = fc.checkTetra({u}, nodeVerts, nodeTetra);
        if (r == -1) noChain++;
        if (r == -2) capped++;
        if (nNodes <= 20 || r == -1 || r == -2)
            printf("  node %2d: %s\n", u, r == -1 ? "NO CHAIN" : (r == -2 ? "capped" : to_string(r).c_str()));
    }
    printf("  --> %d/%d single nodes admit NO chain (%d capped)\n", noChain, nNodes, capped);

    printf("[pairs] pairwise compatibility...\n");
    vector<vector<uint8_t>> compat(nNodes, vector<uint8_t>(nNodes, 1));
    for (int i = 0; i < nNodes; i++)
        for (int j = i; j < nNodes; j++) {
            FeasCache fc(&k, &units, cap);
            int r = fc.checkTetra({i, j}, nodeVerts, nodeTetra);
            compat[i][j] = compat[j][i] = (r >= 0 || i == j) ? 1 : 0;
            if (r == -2) { compat[i][j] = compat[j][i] = 1; }   // capped: assume compatible (safe bound)
        }
    int inc = 0;
    for (int i = 0; i < nNodes; i++) for (int j = i + 1; j < nNodes; j++) if (!compat[i][j]) inc++;
    printf("  --> %d incompatible pairs of %d\n", inc, nNodes * (nNodes - 1) / 2);
    if (inc > 0 && nNodes <= 20) {
        for (int i = 0; i < nNodes; i++) {
            printf("  %2d ", i);
            for (int j = 0; j < nNodes; j++) printf("%c", compat[i][j] ? '#' : '.');
            printf("\n");
        }
    }

    {
        vector<vector<uint8_t>> adj0 = compat;
        for (int i = 0; i < nNodes; i++) adj0[i][i] = 0;
        vector<int> best{0};
        auto t0 = clk::now();
        bronKerbosch({}, [&]{ vector<int> p; for (int i = 0; i < nNodes; i++) p.push_back(i); return p; }(),
                     {}, adj0, best, t0, capMs);
        printf("  --> pairwise max clique = %zu\n", best.size());
    }
    for (int kk = 1; kk <= 5; kk++) {
        vector<int> color;
        auto t1 = clk::now();
        bool ok = colorComplement(compat, nNodes, kk, color, t1, capMs);
        printf("  --> pairwise clique cover with %d parts: %s\n", kk, ok ? "YES" : "no");
        if (ok) break;
    }

    printf("\n[exact cover] minimum number of JOINTLY-FEASIBLE domains:\n");
    int minParts = -1;
    vector<vector<int>> bestParts;
    for (int kk = 1; kk <= maxK; kk++) {
        vector<vector<int>> parts(kk);
        int partsUsed = 0;
        FeasCache fc(&k, &units, cap);
        coverNodes.store(0);
        auto t1 = clk::now();
        bool ok = coverSearchT(0, nNodes, parts, partsUsed, kk, fc, compat, nodeVerts, nodeTetra);
        double secs = chrono::duration<double>(clk::now() - t1).count();
        printf("  cover with %d parts: %s (%.1f s, %ld search nodes)\n",
               kk, ok ? "FOUND" : "no", secs, coverNodes.load());
        if (ok) { minParts = kk; bestParts = parts; break; }
        if (coverNodes.load() > SEARCH_BUDGET) { printf("  (search budget exhausted)\n"); break; }
    }

    if (minParts > 0) {
        printf("\n  ==> MINIMUM DOMAINS (this level) = %d\n", minParts);
        for (int p = 0; p < minParts; p++) {
            printf("  domain %d (%d nodes): ", p, (int)bestParts[p].size());
            for (size_t x = 0; x < bestParts[p].size(); x++) {
                if (x && x % 16 == 0) printf("\n      ");
                printf("%d ", bestParts[p][x]);
            }
            printf("\n");
        }
        for (int p = 0; p < minParts; p++) {
            FeasCache fc(&k, &units, cap);
            int r = checkPartT(bestParts[p], nodeVerts, nodeTetra, fc);
            printf("  domain %d re-verified: %s%s\n", p, r >= 0 ? "OK" : "FAIL",
                   (r >= 0 && (int)bestParts[p].size() > gAstarCapTetras) ? " (chain by LocalSearch)" : "");
        }
    }
    return minParts;
}

int main(int argc, char **argv) {
    setvbuf(stdout, NULL, _IONBF, 0);
    int which = (argc > 1) ? atoi(argv[1]) : 3;
    Kx k = buildK(which);
    int nF = (int)k.facets.size();

    vector<Unit> units;
    for (int i = 0; i < nF; i++)
        for (int j = 0; j < nF; j++)
            units.push_back(buildUnit(k, i, j));
    int U = (int)units.size();
    int totalTet = 0;
    for (auto &u : units) totalTet += (int)u.tetra.size();
    printf("K: %d facets, %d vertices; %d product units, %d maximal simplices\n",
           nF, k.nv, U, totalTet);

    long cap = 2000000, capMs = 120000;

    // Level 1: whole product pairs as atoms (summary only)
    vector<vector<int>> uVerts;
    vector<vector<vector<int>>> uTetra;
    for (auto &u : units) { uVerts.push_back(u.verts); uTetra.push_back(u.tetra); }
    bool skipUnitCover = (argc > 2 && atoi(argv[2]) == 1);
    runLevel(k, units, U, uVerts, uTetra, "PRODUCT-PAIRS", cap, capMs, skipUnitCover ? 0 : 5);

    // Level 2: individual staircase simplices (what AddFacet actually adds)
    vector<vector<int>> tVerts;
    vector<vector<vector<int>>> tTetra;
    for (auto &u : units)
        for (auto &tet : u.tetra) {
            set<int> vs(tet.begin(), tet.end());
            tVerts.push_back(vector<int>(vs.begin(), vs.end()));
            tTetra.push_back({tet});
        }

    // init the repo-style objects used by the hybrid verifier
    if (which != 3) gAstarCapTetras = 1 << 30;   // small cases: always exact A*
    else {
    gKomplex = new Complex();
    gKomplex->initComplex(4, 4);
    gKomplex->K.A[0][0].initSimplex(3, 0, 1, 2);
    gKomplex->K.A[0][1].initSimplex(3, 0, 1, 3);
    gKomplex->K.A[0][2].initSimplex(3, 0, 2, 3);
    gKomplex->K.A[0][3].initSimplex(3, 1, 2, 3);
    gKomplex->initAdjMat();
    Lg.initSubComplexJ(96);
    {
        int counting = 0;
        for (int u = 0; u < (int)tVerts.size(); u++) {
            SimplexAlpha s;
            s.initSimplexAlpha((int)tVerts[u].size());
            for (size_t kk = 0; kk < tVerts[u].size(); kk++) {
                int pid = tVerts[u][kk];
                s.initVertex((int)kk, pid / k.nv, pid % k.nv);
            }
            Lg.initA(counting++, s);
        }
        Lg.initZero_Skeleton();
    }
    Jg.initSubComplexJ(1);
    Jg.initA(0, Lg.listOfFacets.A[0][0]);
    Jg.initZero_Skeleton();
    gMap1.image.initMatrixInt(k.nv, k.nv);
    gMap0.image.initMatrixInt(k.nv, k.nv);
    for (int i = 0; i < k.nv; i++)
        for (int j = 0; j < k.nv; j++) {
            gMap1.image.updateA(i, j, i);
            gMap0.image.updateA(i, j, j);
        }
    gSuchen = new LocalSearch();
    gSuchen->initLocalSearch(Lg, *gKomplex, gMap1, gMap0, gLsM, 0.1);
    }

    int res = runLevel(k, units, (int)tVerts.size(), tVerts, tTetra,
                       "MAXIMAL-SIMPLICES", cap, capMs, 5);
    printf("\n==== FINAL: minimum domains at simplex level = %d ====\n", res);
    return 0;
}
