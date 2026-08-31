# Trigger pairs cluster around common endpoints: a fixed-count graph baseline

## New empirical result

The real rank-one checkpoints have **4.65–5.39 times** the cooperative
continuation variance predicted by a uniform trigger graph with the **same
number of safe sites and the same number of triggering pairs at every
checkpoint**. Merely counting triggering pairs does not explain how strongly
they overlap. Shared-endpoint organization is the resolved missing structure.

This is zero-new-sample analysis of the N325/N425 20k production blocks at
`e81dd59`. It is not another independent confirmation of those blocks, and
the graph null is a retrospective mechanism comparator, not a pre-data gate.

| N / orientation | observed mean Delta | fixed-edge mean Delta0 | ratio ± cluster SE | excess z |
|---|---|---|---|---|
| 325 first | 1.188143e-4 | 2.552493e-5 | **4.6548 ± 0.0475** | 47.97 |
| 325 second | 1.224303e-4 | 2.556302e-5 | **4.7894 ± 0.0534** | 44.50 |
| 425 first | 7.697726e-5 | 1.428635e-5 | **5.3882 ± 0.0598** | 45.58 |
| 425 second | 7.544237e-5 | 1.411078e-5 | **5.3464 ± 0.0586** | 46.17 |

Here the ratio is `E[Delta]/E[Delta0]`, not the average of noisy per-checkpoint
ratios. Expectations are conditional on the checkpoint being rank one.
The two-coordinate, paired-orientation quadratic forms are 4338.95 and
4206.45 (nominal 2 df); these are asymptotic cluster-covariance summaries,
not exact finite-sample p-values. The mean effects and their SE are the
more useful magnitude statement.

Between **98.49% and 99.21%** of actual at-risk checkpoints individually exceed
the conditional fixed-edge expectation. The average excess numbers of
two-stars are approximately **105.66, 109.71, 160.43, 156.95**, respectively.
This is widespread shared-endpoint concentration, not a signal supported by
one exceptional counterexample.

## Exact dictionary

At a fixed occupied rank-one configuration C, let d be the number of vacant
sites. Let a be the number of vacant sites safe for one insertion. Make a graph
on these a sites:

- an edge is a pair that is individually safe but jointly triggers rank two;
- m is the number of these minimal triggering pairs;
- t_v is a vertex's degree in this trigger graph;
- c_v is its number of safe second insertions.

The archive field `checkpoint_b2_safe_pairs` counts safe pairs, here denoted
`b2_safe`. This is **not** #403's b2 notation for minimal trigger pairs.
The translation is

```text
m = choose(a,2) - b2_safe,
c_v = (a-1)-t_v,
sum_v c_v = 2*b2_safe,
Var_v(c_v) = Var_v(t_v).
```

The already measured same-checkpoint cooperative excess is exactly

```text
Delta(C) = a/[d(d-1)^2] * Var_v(t_v).
```

Thus the new experiment was already measuring a graph degree variance. No
new hidden-state construction or topology classifier is required to read it.

## Fixed-edge null without null simulation

Condition on each observed a,m. Distribute m edges uniformly over the
`M=choose(a,2)` possible pairs: the simple graph model G(a,m).
This keeps both the one-site ceiling and two-site survival count exactly.
It only removes the physical organization of the trigger pairs.

For a fixed vertex, the degree has a hypergeometric distribution with a-1
incident slots among M possible edges. Put p=m/M. Because every graph has
the same mean degree 2m/a,

```text
E_null Var_v(t_v)
  = Var_null(t_fixed_vertex)
  = (a-1) p(1-p) [M-(a-1)]/(M-1)
  = (a-1)^2 p(1-p)/(a+1),   a>=3.
```

Cases a<=2 have zero variance and are treated directly. Multiplying by
`a/[d(d-1)^2]` gives Delta0(C). No graph Monte Carlo is used; this is an exact
finite-population expectation, not an independent-edge approximation.

An equivalent, especially readable observable is the number of two-stars,
`W2=sum_v choose(t_v,2)`: unordered edge pairs that share a center. A triangle
contains three such two-stars. Their exact null expectation is

```text
E_null W2 = a choose(a-1,2) m(m-1)/[M(M-1)].
```

The extra cooperative contribution has the simple exact identity

```text
Delta(C)-Delta0(C)
  = 2 [W2(C)-E_null W2]/[d(d-1)^2].
```

The comparison is therefore a direct overlap measurement, not a fit of an
uninterpreted residual. It needs only the stored integers a, b2_safe and
`sum_v c_v^2`. The script checks this identity on every archived row.

## Covariance, source and limits

One original base permutation is one cluster; both orientations are retained
jointly. Non-risk orientations contribute zero to the influence function of
the conditional means, rather than being treated as independent missing rows.
Each size keeps its full 20-by-20 measurement covariance, including observed
variance, null variance, their difference, stars and graph counts. The ratio
uses their joint covariance. N325 and N425 use independent RNG domains.
Source CSV SHA-256 values and immutable commit are required by the scorer;
metadata hashes, counters and seeds are in the result.

The measured ratio is larger in N425 by 0.7333 ± 0.0764 and 0.5571 ± 0.0793
for the two orientations. This is a descriptive comparison of these two
sizes/geometries, **not** an asymptotic growth law or a fitted exponent.

The eliminated model is the uniform fixed-edge trigger assignment. The result
does not rule out local geometric explanations, prove path memory, identify
a continuum field, or make the degree moment a complete autonomous state.
In particular, pair-graph overlap cannot replace the genuinely minimal
three-site triggers measured separately on the saved N425 configurations.

The next microscopic model must reproduce both the number of trigger edges
and their shared-endpoint concentration. A model that gets H2 and the pair
count right but assigns partners exchangeably misses most of this cooperative
component. That is a sharper target than adding another scalar age fit.

## One overlap coordinate connects branching and three-step survival

The separate real-checkpoint census `f86ab1b` makes the connection explicit.
Let T3 be the number of triangles in the **trigger** graph and c3 the number
of genuine minimal three-site triggers whose proper subsets are all safe.
Inclusion–exclusion on the a individually safe vertices gives

```text
safe_pair_graph_triangles = choose(a,3)-m(a-2)+W2-T3,
actual_safe_triples = choose(a,3)-m(a-2)+W2-T3-c3.
```

For the two saved N425 second-orientation checkpoints, a=d=173 and m=108
are identical. Their trigger two-star counts are **926 and 1466**, their
trigger triangle counts are both **zero**, and c3 is **583 and 509**.
Consequently, from A to B:

```text
Delta W2 = +540,
Delta T3 = 0,
Delta c3 = -74,
Delta actual_safe_triples = 540+74 = 614.
```

The same +540 overlapping edge pairs increase the double-clone probability
by `2*540/[173*172^2]=135/639754`, while the complete three-step survival
increase is `614/choose(173,3)=307/424023`. Exactly `270/307` of that latter
difference comes from pair overlap and `37/307` from the genuine three-site
layer. These are identities for this real pair, not population percentages.

This resolves an otherwise misleading intuition: at fixed trigger-edge count,
more shared-endpoint clustering can **increase branching heterogeneity and
increase three-step survival simultaneously**. Overlapping triggers cover
fewer distinct triples than dispersed triggers. This is not a claim that
three-site effects always have the same sign: here their independent change
reinforces the pair-overlap contribution. The full three-step predictor still
needs the c3 layer.

## Scientific card

- Mechanism changed: fixed-count exchangeable trigger graphs substantially
  underpredict the within-checkpoint cooperative variance.
- Observer / sector: exact common-update/two-clone survival on occupied-NN
  ambient-rank-one checkpoints; trigger-graph two-star excess.
- Source / geometry: original 20k blocks, N325 k0=193 and N425 k0=252,
  two paired Gaussian orientations each.
- Dependency groups: `p334-cooperative-N325-20260831` and
  `p334-cooperative-N425-20260831`; same data as e81dd59/6147e22, no new replica.
- Lifecycle: retrospective analytic baseline scored on production data;
  zero new random draws, cluster covariance retained.
- Not proved: full state closure, continuum memory, a critical exponent, H4
  operator identity, or a physical stochastic-graph replacement theorem.
- Next observation: a geometry-aware trigger model predicts the two-star
  excess and the separately measured minimal-three-site layer together.

## Reproduce

In a clone containing `e81dd59`, with NumPy:

```sh
python scripts/p334_trigger_overlap_baseline.py
```

The scientific calculation reads the immutable CSV blobs directly. No source
Monte Carlo, broad regression suite or earlier synthetic tests are rerun.
