# General-period threshold ranks: minimal norm-4 production backend

Status: production-capable CPU implementation and Python-oracle regression.

The standalone engine is
`src/threshold_rank_integer_period_mc.cpp`.  It extends the cross-channel
bidirectional Newman--Ziff pipeline to every nonsingular 2x2 period matrix
whose entries and quotient order fit the engine's signed 32-bit production
indexing, without changing the frozen cyclic Gaussian engine or its coupling.

## Exact quotient and topology contract

For a row-major integer matrix

```text
P = ((a,b),(c,d))
```

the columns are the two lifted periods and the graph is
`Z^2 / P Z^2`.  A one-time integer column-Hermite reduction constructs

```text
H = P V = ((h11,h12),(0,h22)),  det(V)=+/-1,
0 <= h12 < h11,  h11*h22 = |det(P)|.
```

Every vertex has the unique mixed-radix label

```text
0 <= rx < h11,  0 <= ry < h22,
label = rx + h11*ry.
```

This is a coordinate choice, not an assumption that the quotient is cyclic.
The Smith invariants are recorded separately as

```text
s1 = gcd(a,b,c,d),  s2 = |det(P)|/s1.
```

The engine creates NN forward edges `(1,0),(0,1)` and matching-lattice
forward edges `(1,0),(0,1),(1,1),(1,-1)` from every quotient vertex.  Edge
records retain these lifted displacements.  When a potential-union-find cycle
closes with lifted displacement `Delta`, its period-basis winding is computed
exactly as

```text
w = adj(P) Delta / det(P).
```

The component is cross-wrapping exactly when two rationally independent
winding generators have appeared.  Thus the produced ranks retain the frozen
contract:

- `K_plus`: first black NN rank with one rank-2 component;
- `K_minus`: `N-r+1`, where `r` is the first reverse white matching rank with
  one rank-2 component.

## First production designs

The built-in `--n` values are the nonprimitive norm-4 children from the
maximin/isogeny analysis:

| lineage | first period matrix | second period matrix | Smith group |
|---|---|---|---|
| `65 -> 260` | `((16,-2),(2,16))` | `((14,-8),(8,14))` | `Z/2 x Z/130` |
| `85 -> 340` | `((18,-4),(4,18))` | `((14,-12),(12,14))` | `Z/2 x Z/170` |

The CSV schema is intentionally identical to
`threshold_rank_orientation_mc.cpp`, with Gaussian lineage labels
`(16,2)/(14,8)` and `(18,4)/(14,12)`.  Consequently the existing
`analyze_threshold_rank_orientation.py` full-curve analyzer consumes these
outputs without conversion.  Metadata additionally records both complete
period matrices, both HNFs, and both Smith pairs.

Example:

```bash
clang++ -O3 -DNDEBUG -std=c++17 \
  src/threshold_rank_integer_period_mc.cpp \
  -o build/threshold_rank_integer_period_mc

./build/threshold_rank_integer_period_mc \
  --samples 1000000 --batches 100 --n 260 --threads 1 \
  --seed 20260828 --replica-offset 0 \
  --git-commit "$(git rev-parse HEAD)" \
  --output-prefix results/norm4/n260
```

On Linux, compile with `g++ -fopenmp`; batch scheduling and counter-derived
permutations make `.hist.csv` and `.moments.csv` independent of thread count.
Arbitrary equal-order period pairs are available through
`--first-matrix A B C D --second-matrix A B C D`, with optional lineage labels
`--first-rep A B --second-rep A B`.

## Regression evidence

`tests/test_threshold_rank_integer_period_mc.py` compiles the engine and
checks:

1. exact N=5 all-permutation thresholds;
2. arbitrary-matrix quotient and winding arithmetic;
3. exhaustive threshold-rank invariance under a unimodular period-basis
   change;
4. Smith `(2,130)` and `(2,170)` metadata for N260/N340;
5. N260 integer histograms against the independent general-period Python
   reference, using an explicit HNF-label-to-reference-coordinate map;
6. identical integer outputs with one and two requested threads;
7. the custom arbitrary-matrix CLI.

The prior integer-period reference and cyclic production regressions remain
unchanged and pass alongside the new tests.

## Local benchmark

Apple M4, 10 cores, Apple clang 21.0.0, `-O3 -DNDEBUG`, one thread, 100,000
paired permutations per row:

| backend/design | elapsed in engine metadata | paired replicas/s | paired site updates/s |
|---|---:|---:|---:|
| general period N260 | 2.449 s | 40,827 | 21.23 M |
| general period N340 | 3.291 s | 30,387 | 20.66 M |
| existing cyclic N265 control | 2.241 s | 44,630 | 23.65 M |

At nearly equal size, general-period support costs about 10% in this
single-thread comparison.  Quotient reduction itself is not in the sample
loop; the remaining difference is the generic period-coordinate arithmetic
and the different graphs/crossing-rank distributions.

## Interpretation boundary

The shared HNF label permutation is a legal, deterministic common-random-
number coupling.  Its covariance is not a physical observable and need not
match the cyclic parent/child coupling.  Scientific scoring must use measured
batch covariance.  The first use of N260/N340 is therefore the frozen
`T4=T2^2` radial-curvature and quotient-universality experiment, not a claim
that Smith content is already known to be irrelevant.
