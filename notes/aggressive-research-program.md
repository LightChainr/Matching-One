# Aggressive research program

## Objective

The target is not merely another decimal estimate of the square-lattice site-percolation threshold. The target is to discover and exploit the *structure* that makes the square lattice nearly solvable:

\[
p_c(G)+p_c(G^*)=1,
\]

where `G` is the nearest-neighbour square lattice and `G*` is its site-matching graph (NN+NNN / king graph).

The program treats a closed form as a downstream possibility, not as the primary method. The first objective is to identify an estimator, spectral identity, or convergent sequence whose leading irrelevant fields are cancelled by matching symmetry.

## Current numerical situation

The 2015 periodic-cylinder sequence has 21 exact/high-precision finite-width values and was extrapolated as

\[
p_c(n)=p_c+A_1n^{-4}+A_2n^{-6}+A_3n^{-8}+\cdots.
\]

A 2024 Comment extended cylindrical and helical calculations to width 24 and challenged the final digits of the 2015 estimate. Therefore the repository must store estimates by method and must not encode one extrapolated value as an exact definition.

A direct linear audit of the published `n <= 21` sequence already exposes substantial model dependence. Using all points `n >= 8`, fixed-power fits give:

| correction powers | fitted intercept |
|---|---:|
| `4,6` | `0.592746049501568462` |
| `4,6,8` | `0.592746052559226704` |
| `4,6,8,10` | `0.592746050975478220` |
| `4,6,8,10,12` | `0.592746050900176725` |

These are **not threshold estimates**. They are a diagnostic showing that in-sample agreement and a plausible exponent sequence are insufficient to justify a 14-digit error bar. Model selection must be based on withheld-tail prediction and cross-estimator consistency.

## Primary hypotheses

Each hypothesis is deliberately strong and falsifiable.

### H1 — matching-odd leading correction

For a correctly paired dimensionless observable `Q` on `G` and `G*`, evaluated at complementary probabilities, the leading correction amplitude is odd under matching:

\[
Q_G(p,L)-Q_{G^*}(1-p,L)
 = a_4(p-p_c)L^{y_t-4}+a_6L^{-6}+\cdots.
\]

At the root of the matching difference, the corresponding pseudo-critical estimator has no generic `L^{-4}` bias. Its first surviving bias is `L^{-6}` or smaller.

**Falsification:** fit matched and unmatched estimators over increasing size windows. Reject H1 if the matched estimator retains a stable nonzero `L^{-4}` amplitude or if its rolling-origin error is not improved.

### H2 — a superconvergent matching estimator exists

There is a linear or mildly nonlinear combination of independently defined finite-size estimators,

\[
\widehat p_L=\sum_j w_j p_L^{(j)},\qquad \sum_jw_j=1,
\]

using at least one matching-paired observable, whose `L^{-4}` and `L^{-6}` amplitudes both vanish. The resulting bias begins at `L^{-8}` or beyond.

The weights must be estimated on small/medium sizes and then frozen before testing on larger withheld sizes. A size-dependent unconstrained fit does not count.

**Falsification:** pre-register weights using sizes up to `L_train`; test on larger exact or high-statistics data. Reject if the improvement disappears out of sample.

### H3 — correction exponents form matching-parity sectors

The observed `4,6,8,...` sequence is not a generic polynomial accident. It is the spectrum of irrelevant contributions allowed by the topological-sector eigenvalue identity. Matching decomposes those contributions into even and odd sectors, and an appropriate symmetrization removes one sector.

**Falsification:** compute amplitudes for several boundary conditions and matching pairs. Reject the simple parity version if signs and zeros do not transform consistently.

### H4 — the 2015 last digits are an extrapolation artefact

The numerical roots in the 2015 table are accurate, but the quoted uncertainty on the infinite-width intercept is too small because it does not cover plausible omitted corrections. The true threshold is expected to lie below `0.59274605079210` and in the broad exploratory band

\[
0.592746050789\ < p_c <\ 0.592746050791.
\]

This band is a **research prior**, not an accepted result. It must be replaced by a model-averaged interval supported by withheld-width tests and independent estimator families.

**Falsification:** an independently reproduced width sequence beyond 24, or a demonstrably predictive extrapolation, converges outside the band.

### H5 — no low-complexity elementary closed form

If an exact representation exists, it is more likely to be spectral or constructive than elementary: a limit of critical polynomials, a zero of a Fredholm/spectral determinant, or a distinguished point on a critical manifold. Searching combinations of `pi`, `e`, radicals, and gamma values before stabilizing the numerical interval is low-information work.

**Falsification:** a parameter-free expression survives interval arithmetic, high-height PSLQ controls, and predicts another lattice or finite-size amplitude.

### H6 — an exactly solvable approximant sequence can be constructed

There exists a sequence of self-dual decorated cells or correlated-bond gadgets `G_k` with exactly computable critical manifolds such that

\[
p_c(G_k)\to p_c(\mathbb Z^2_{\rm site})
\]

and the error decreases faster than ordinary basis-size critical-polynomial convergence. The cell construction should preserve the four-terminal connectivity distribution of an occupied square-lattice site to increasing order.

**Falsification:** exhaustive search over bounded cell complexity shows no sequence with stable monotone or accelerated convergence.

## Experimental strategy

### Phase A — provenance and blind extrapolation

1. Store every published finite-size sequence as immutable decimal data with source, boundary condition, estimator definition, and precision provenance.
2. Reproduce each published extrapolation exactly.
3. Hide the largest three or four widths, select models only by prediction of those points, then reveal them.
4. Repeat across `p_pol`, cylindrical, helical, `p_med`, `p_cell`, wrapping, and matching-difference estimators.
5. Report a model ensemble rather than a single preferred fit.

The minimum acceptable output is a table of out-of-sample errors and intercept drift under changes of `n_min`, exponent set, and rational/polynomial ansatz.

### Phase B — exact finite matching identities

Implement simultaneous connectivity tracking for:

- occupied NN clusters on the square lattice;
- vacant NN+NNN clusters on the matching lattice;
- cluster counts;
- horizontal, vertical, both-direction, and either-direction wrapping events.

On a torus, exact finite-size matching relations become property tests. A single simulation stream can evaluate both members of the matching pair, producing strongly correlated estimates and sharply reducing the variance of their difference.

### Phase C — coupled Monte Carlo and GPU campaign

Generate one uniform random number `u_v` per site. At parameter `p`, the site is occupied on `G` when `u_v < p`; at complementary parameter `1-p`, use the same field to define the matching configuration. This coupling should be compared with independent sampling and antithetic variants.

GPU work is justified only after exact small-size tests pass. The useful GPU kernel is a bit-packed, many-replica wrapping/connectivity computation, not a naive single-lattice union-find port. Batch thousands of independent replicas and scan a narrow probability grid around the current threshold band.

Pre-register:

- lattice sizes and aspect ratios;
- number of independent replicas;
- RNG and counter layout;
- root-finding observable;
- covariance estimator;
- finite-size model and withheld sizes.

### Phase D — transfer matrix / spectral extension

The exact transfer-matrix problem is dominated by state count, irregular sparse transitions, hashing, and memory bandwidth. A single consumer GPU is unlikely to be the decisive resource. Priorities are:

1. canonical compact encoding of non-crossing connectivities;
2. symmetry quotienting;
3. matrix-free dominant-eigenvalue iteration;
4. modular or mixed-precision arithmetic where valid;
5. distributed state sharding and deterministic checkpointing.

On the stated 8-core/16-GB host, reproduce widths well below the frontier and profile state growth. Extending beyond width 24 will probably require a high-memory machine rather than a faster GPU.

### Phase E — operator and gadget searches

In parallel:

- derive which irrelevant fields can contribute to the eigenvalue difference in the `c=0` logarithmic CFT;
- test matching parity at the amplitude level;
- enumerate small decorated cells and compare their exact terminal-connectivity distributions with those induced by a square-lattice site;
- search for self-dual critical manifolds that generate convergent exact approximants.

## Statistical and numerical rules

1. Never quote more stable digits than survive all reasonable withheld-tail models.
2. Never convert a deterministic fit residual into a confidence interval.
3. Distinguish arithmetic precision, sampling error, and extrapolation-model error.
4. Freeze model choices before evaluating the largest widths.
5. Publish failed models and negative results.
6. Use interval arithmetic for every claimed exclusion or exact relation.
7. A new decimal is secondary; a new cancellation, identity, exponent rule, rigorous bound, or predictive estimator is primary.

## Hardware map

| Workload | 8 CPU / 16 GB | RTX 5090-class GPU | High-memory CPU node |
|---|---|---|---|
| literature/data audit | sufficient | no | no |
| arbitrary-precision extrapolation | sufficient | no | no |
| exact small torus enumeration | sufficient to moderate sizes | little benefit | useful later |
| coupled Monte Carlo | sufficient for validation | high benefit | optional |
| bit-parallel many-replica scan | slow/moderate | high benefit | optional |
| transfer matrix near current frontier | reproduction only | low/uncertain benefit | likely required |
| self-dual gadget enumeration | sufficient for first search | optional batching | useful for large search |
| rigorous interval bounds | sufficient for prototypes | no | useful for scale |

## First decision gates

- **Gate 1:** the audit script must predict withheld finite widths better than a naive fixed-power fit.
- **Gate 2:** exact matching identities must pass exhaustive tests on small tori.
- **Gate 3:** a matched estimator must show reproducible cancellation or variance reduction.
- **Gate 4:** only then spend GPU time on large-scale Monte Carlo.
- **Gate 5:** only after the numerical interval is stable resume interval-PSLQ searches.
