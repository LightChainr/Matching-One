# Norm-5 graph-cover order-statistic coupling pilot

## Decision

The construction is exact, but this pilot gives a **negative promotion
decision** for the frozen norm-5 H4 residual.  On a held-out 2,000-replica
evaluation split of the actual `N=65 -> 325` two-orientation genealogy, the
equal average of all five order-statistic parent fields had variance ratio

```text
0.998199 versus an independent-parent baseline
```

or only `1.0018x` apparent gain.  The paired variance-difference score was
`z=-1.07`.  Weights trained on the first 2,000 replicas made the held-out
variance worse (`ratio=1.00548`).  This is nowhere near Issue #67's `2x`
promotion gate, even before charging for four extra parent evaluations.

This is a method result, not a score of H4, H12, or any percolation mechanism.

## Exact construction and proof

For each parent vertex `j`, let its degree-`Q` child fiber carry iid priorities

```text
U[j,1], ..., U[j,Q] ~ iid Uniform(0,1).
```

Write `U[j,(k)]` for the `k`-th order statistic and define

```text
V[j,k] = F_Beta(k,Q+1-k)(U[j,(k)]).
```

The order-statistic law is `Beta(k,Q+1-k)`.  Its CDF is continuous, so the
probability-integral transform makes `V[j,k]` exactly Uniform.  Fibers for
different `j` use disjoint child variables and are independent.  Therefore,
for every fixed `k`, `{V[j,k]}_j` is an exact iid parent priority field; sorting
it gives an exact uniform parent Newman--Ziff permutation.

The five fields indexed by `k` are correlated with one another.  The result is
five exact-marginal parent couplings, **not** five mutually independent parent
samples.  Nevertheless, for any weights summing to one,

```text
child - r * sum_k w_k parent_k
```

is an unbiased estimator of the same frozen residual, because every parent
term has the correct marginal expectation.

The implementation evaluates the integer-parameter Beta CDF through its
binomial-tail polynomial and consumes the existing exact cover/fiber map.

## Validation

Focused tests cover:

- exact norm-2/norm-5 fiber, edge, and homology-map contracts;
- Uniform mean/variance for every transformed order field at `Q=5`;
- all six parent permutations on a tiny three-site label set;
- unbiased `K_plus` means against exact enumeration on the Gaussian `N=5`
  torus;
- bijective child-priority assignment through an exact norm-5 cover.

All ten focused tests pass.

## Actual-lineage pilot

The pilot uses the frozen raw H4 multiplier

```text
r = (-14/25) * 5^(-13/8) = -0.0409601718395...
```

at `p_ref=0.592746050790`, with the actual two-orientation genealogy:

```text
parent N=65: (8,1) / (7,4)
child N=325: (17,6) / (18,1)
```

Each permutation is converted to the fixed-p matching value by conditioning
on its exact `(K_minus,K_plus)` thresholds.  Replica counters `0..1999` train
the constrained sum-to-one weights; counters `2000..3999` evaluate them.

| parent coupling | residual variance / independent | apparent gain | paired z |
|---|---:|---:|---:|
| independent parent | 1.000000 | 1.0000x | -- |
| equal average, five orders | 0.998199 | 1.0018x | -1.07 |
| trained sum-to-one weights | 1.005480 | 0.9945x | +1.36 |
| best individual order (`k=2`) | 0.996123 | 1.0039x | -1.79 |

The equal-average child/parent correlation is only `0.00183`.  The strongest
individual correlation is `-0.0425` at `k=2`, but its held-out variance gain is
only `0.39%` and is not stable evidence for promotion.

Computing all five parent permutations requires five parent topology
evaluations instead of one.  Even granting the tiny variance change, the
variance-per-wall-time direction is therefore strictly unfavorable in this
reference implementation.

## Consequence

Keep the order-statistic identity as an exact reference and possible building
block for a different observable.  Do not insert it into norm-5 production for
the current frozen H4 semigroup residual.  The direct domain-separated parent
baseline remains the appropriate path unless a materially different coupling
creates substantially larger negative parent/child covariance.

Machine-readable estimates and all replica-level values are in `pilot.json`
and `replicas.csv`.
