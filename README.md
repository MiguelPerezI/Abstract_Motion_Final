# Abstract_Motion

C++/MPI implementation of the **randomized greedy algorithms for simplicial
complexity** from the paper [*A randomized greedy algorithm for motion planning
in simplicial complexes*](https://arxiv.org/abs/2008.13290)
(C. Ortiz, J. Lara, A. González, G. Borat — arXiv:2008.13290).

The code computes, for a finite simplicial complex `K`, an increasingly good
approximation of **SC_strict(K)** — the strict simplicial complexity, the
minimum number of *domains* (subcomplexes of `K×K` on which a contiguity chain
between the two projections exists) needed to cover `K×K`. This is a
combinatorial upper bound on Farber's **topological complexity TC(K)**, which
is the classical robot motion-planning invariant: `TC(K) ≤ SC_strict(K)`.

Current target (set in `main.cpp`): `K = ∂Δ³` (the 2-sphere `S²`), whose
product `K×K` has **96 maximal facets**.

---

## 1. Background in one page

* A **motion planner** on a domain `J ⊆ K×K` is a simplicial map
  `π₁|_J → π₂|_J` (or a chain of such maps). Two simplicial maps
  `φ, φ′ : J → K` are **contiguous** (Def. 2.3 of the paper) when
  `φ(σ) ∪ φ′(σ)` is a simplex of `K` for every facet `σ` of `J`.
* A **contiguity chain** `π₁|_J → φ₁ → ... → π₂|_J` is a planner on `J`;
  the number of maps minus one is its length.
* `SC_strict(K)` = (minimum number of domains covering `K×K`) − 1, where each
  domain admits a chain. `TC(K)` is the same number where each domain admits a
  single contiguous pair.
* Products `σ×τ` are triangulated with the **staircase triangulation**
  (maximal simplices correspond to interleavings; triangle×triangle gives
  6 tetrahedra, edge×edge gives 2 triangles).
* The algorithms search for these domains **randomized and greedily**:
  * **Alg. 1 – LocalSearch** (`LocalSearch.hpp: ghostLocalSearch`): random walk
    in the graph of simplicial maps `J → K`, accepting a step when the map
    stays contiguous to `π₂` and it is closer to it (with an escape
    probability `r` and a step budget `M`).
  * **Alg. 2 – Reduce** (`LocalSearch.hpp: reduceChain`): shrinks a chain by
    dropping intermediate maps.
  * **Alg. 3 – AddFacet / RCC** (`Covering_CORE.hpp: simplexScrutiny`): tries
    to extend a domain by one facet; this is the workhorse.
  * **Alg. 4 – Covering** (`Covering_CORE.hpp: runCovering`): builds a whole
    domain cover by repeated RCC from random seeds.
  * **Alg. 5 – OptimizedCovering** (`OptimizedCovering.hpp`): Alg. 4 plus a
    "bubble" refinement that merges/re-balances large domains between cycles.

Because the search is randomized, a *failed* run proves nothing — it only
means the walk didn't find a chain within `M` steps. For **exact** yes/no
answers use the diagnostic tools in §5.

---

## 2. Repository layout

| File | Role |
|---|---|
| `main.cpp` | Experiment driver: builds `K`, `K×K`, calls OptimizedCovering |
| `SimplexAlpha.{hpp,cpp}` | Complex, product complex, subcomplex `J`, simplicial maps, contiguity & distance |
| `SimplexAbstract.{hpp,cpp}` | Basic containers (Simplex, Matrix, VectorInt, ...) |
| `LocalSearch.{hpp,cpp}` | Alg. 1 (LocalSearch) and Alg. 2 (Reduce) |
| `Covering_CORE.hpp` / `Covering.hpp` | Alg. 3 (AddFacet/RCC) and Alg. 4 (Covering), MPI consensus |
| `OptimizedCovering.hpp` | Alg. 5 (the one `main.cpp` runs) |
| `RCC.hpp`, `Lex.hpp` | Alternative/auxiliary code paths |
| `Matrix.{hpp,cpp}` | Adjacency matrix + Dijkstra (graph distance `d`) |
| `job-motion-*.pbs`, `job-motion-LNS.sh`, `job-script.sh` | Cluster submission scripts (SLURM/PBS, 96 ranks) |
| `compileComand.txt` | Original one-line build + run command |
| `test_product.cpp`, `bench_step.cpp` | Verification harness / micro-benchmark |
| `diag_sc2.cpp`, `diag_greedy.cpp` | Exact S² diagnosis tools (see §5) |

**Known limitation:** the code writes the *domains* it finds to `coreN.txt`
(one file per MPI rank) but currently discards the *planner chains*
(`pMap` population is commented out), so `core0.txt` is empty. Reactivating
that output is on the todo list.

---

## 3. Build and run

Requirements: a C++ compiler and an MPI implementation (OpenMPI or Intel MPI).

```bash
# build (matches compileComand.txt; -O2 recommended over the original -O0)
mpic++ main.cpp SimplexAbstract.cpp SimplexAlpha.cpp LocalSearch.cpp Matrix.cpp -o absMtion666 -O2 -pthread

# run on a laptop / desktop (4 ranks)
mpirun -np 4 ./absMtion666

# run on the cluster (96 ranks, 120 h wall clock)
sbatch job-motion-LNS.sh        # SLURM variant
# or qsub job-motion.CAMD.pbs   # PBS variant
```

Every rank writes its own output file `coreN.txt` (`N` = rank id). The
OptimizedCovering driver runs `Ns = 100` refinement cycles over the partition
produced by Covering.

Tunable parameters live at the top of `main.cpp`:

| Parameter | Default | Meaning |
|---|---|---|
| `M`  | `40000` | step budget of each random walk (LocalSearch bound) |
| `r`  | `0.1`   | escape probability of the greedy acceptance rule |
| `Ns` | `100`   | number of OptimizedCovering refinement cycles |

---

## 4. Running the experiments

### 4.1 The paper's example: the circle `K = ∂Δ²` (§5.1)

The paper validates the pipeline on the circle, where `K×K` is a 9-vertex
triangulation of the torus with **18 maximal facets**. The correct answer is
**2 domains** (TC(S¹)=2, so SC_strict = 1), and 1 domain is impossible.

To run it, edit `main.cpp`:

```cpp
int maxInt = 2;              // was 3
int numOfMaxSimplex = 3;     // was 4
komplex.initComplex(numOfMaxSimplex, maxInt + 1);

komplex.K.A[0][0].initSimplex(2, 0, 1);
komplex.K.A[0][1].initSimplex(2, 0, 2);
komplex.K.A[0][2].initSimplex(2, 1, 2);

L.initSubComplexJ(18);       // was 96: circle KxK has 18 maximal facets
```

and drop the six `komplex.graph.addWeight(...)` lines for 4 vertices in
favour of the three circle edges `(0,1)`, `(0,2)`, `(1,2)` (the 1-skeleton is
also rebuilt automatically by `komplex.initAdjMat()`).

`./test_product` (see §5) prints the expected facet counts so you can check
the switch: circle → 18 facets on a 3×3 skeleton, `S²` → 96 facets on 4×4.

### 4.2 `K = ∂Δ³` (S²): the 96-facet search

This is what `main.cpp` does **out of the box**:

* builds the 4 triangular facets `(0,1,2) (0,1,3) (0,2,3) (1,2,3)`,
* forms `K×K` (96 maximal tetrahedra over the 4×4 vertex skeleton),
* runs `OptimizedCovering` with `M=40000, r=0.1, Ns=100`,
* each rank logs its domains to `coreN.txt`.

Expected behaviour on a small machine: the greedy **stalls** — empirically
(`diag_greedy`) only ~13% of AddFacet attempts succeed at `M=1000`, and each
failed attempt costs ~3 s at `M=40000` (see the benchmark below). This is why
the cluster scripts ask for 120 h, and why the exact diagnosis tools of §5
exist: a stall of the randomized search is *not* evidence of impossibility.

```bash
./bench_step     # costs of contiguous() (~70 µs) and d() (~3.8 µs) per walk step
./diag_greedy 1000 20 > diag_greedy.out   # replays AddFacet, logs per-pair success rate
```

---

## 5. Exact search for `K = ∂Δ³` (`diag_sc2`)

`diag_sc2.cpp` is a standalone, MPI-free tool that replaces the randomized
search by an **exact A\*** over the full graph of simplicial maps
`J → K` (moves = changing the image of one vertex; the heuristic
`h = #vertices differing from π₂` is admissible, and single-vertex moves
suffice to refine any contiguity chain). On top of that it runs a
branch-and-bound **minimum-domain-cover search**, so it decides *exactly*
whether a cover with `k` domains exists — and, when found, prints one.

```bash
g++ diag_sc2.cpp SimplexAbstract.cpp SimplexAlpha.cpp LocalSearch.cpp Matrix.cpp -O2 -std=c++17 -o diag_sc2

./diag_sc2 2      # circle:   proves 2 domains suffice, 1 impossible (runs in seconds)
./diag_sc2 3      # S²:       full pipeline, unit-level + simplex-level cover search
./diag_sc2 3 1    # S²:       skip the unit-level cover stage (resume/tune-up mode)
```

Progress goes to stdout; the lines to watch are
`cover with N parts: FOUND` / `no`, `incompatible pairs`, and
`[exact cover] minimum number of JOINTLY-FEASIBLE domains`.

Status of the S² exact search (see `INVES_14.md` for the full write-up):

| Question | Answer |
|---|---|
| Any single facet infeasible? | **No** — all 96 feasible (min chains of 1–8 steps) |
| Incompatible *pairs* of facets? | **0 of 4560** — compatibility graph is complete |
| 1 domain (all of `K×K`)? | **Impossible** (exhaustive proof, 1.9 s) — consistent with TC(S²)=3 |
| **3 domains?** | the decisive open question — search running (`diag_s2.out`) |

So the obstruction to covering `S²×S²` is **joint, not pairwise**: every
facet — and every pair — can live in some planner, but not necessarily in the
same one. That is exactly the situation of the circle (all pairs compatible,
2 domains needed), and it means a failure of the *greedy* dynamics — not a
structural per-facet obstacle — is the most likely explanation for the repo's
stall on `S²`.

---

## 6. Citation

If you use this code, cite the paper it implements:

> C. Ortiz, J. Lara, A. González, G. Borat.
> *A randomized greedy algorithm for motion planning in simplicial complexes.*
> arXiv:2008.13290. https://arxiv.org/abs/2008.13290

Session notes and audits: `~/Downloads/INVES_13.md`, `~/Downloads/INVES_14.md`.