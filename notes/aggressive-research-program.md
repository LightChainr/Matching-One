# Aggressive research program

## Objective

The target is not merely another decimal estimate of the square-lattice site-percolation threshold. The target is to discover and exploit the structure that makes the square lattice *nearly* solvable:

\[
p_c(G)+p_c(G^*)=1,
\]

where `G` is the nearest-neighbour square lattice and `G*` is its site-matching graph (NN+NNN / king graph).

The program treats a closed form as a downstream possibility, not as the primary method. The first objective is to identify an estimator, spectral identity, scaling function, or convergent exact-approximant sequence whose leading irrelevant fields are cancelled by matching symmetry.

## Current numerical situation

The 2015 periodic-cylinder sequence has 21 exact/high-precision finite-width values and was extrapolated as

\[
p_c(n)=p_c+A_1n^{-4}+A_2n^{-6}+A_3n^{-8}+\cdots.
\]

A 2024 Comment extended cylindrical and helical calculations to width 24 and challenged the final digits of the 2015 estimate. A 2024 Reply retained strong confidence in the leading exponent `4`, weaker confidence in `6`, and substantially weaker confidence in the later correction structure. Therefore this repository must store estimates by method and must not encode one extrapolated value as an exact definition.

A direct linear audit of the published `n <= 21` sequence already exposes substantial model dependence. Using all points `n >= 8`, fixed-power fits give:

| correction powers | fitted intercept |
|---|---:|
| `4,6` | `0.592746049501568462` |
| `4,6,8` | `0.592746052559226704` |
| `4,6,8,10` | `0.592746050975478220` |
| `4,6,8,10,12` | `0.592746050900176725` |

These are **not threshold estimates**. They are a diagnostic showing that in-sample agreement and a plausible exponent sequence are insufficient to justify a 14-digit error bar. Model selection must be based on withheld-tail prediction and cross-estimator consistency.

## Exact finite-size structure to exploit

For an `L x L` torus, define the matching function

\[
M_L(p)=N_L(p)-\widehat N_L(1-p)-L^2\chi_\square(p),
\qquad
\chi_\square(p)=p-2p^2+p^4.
\]

It also has the exact representations

\[
M_L(p)=R_L^x(p)-\widehat R_L^x(1-p),
\qquad x\in\{c,b,e,h\},
\]

in terms of paired wrapping events. In the scaling limit,

\[
M_L(p)=f(z)-f(-z),
\qquad z=b(p-p_c)L^{1/\nu},
\qquad \nu=4/3.
\]

The leading scaling function is therefore odd in `z`; its even Taylor coefficients cancel identically. Corrections to scaling shift the unique finite-size root `M_L(p_L^*)=0`. Empirically this root converges approximately as `L^-4`.

If that exponent is exactly `4`, then

\[
M_L(p_c)\sim aL^{-13/4},
\qquad
M'_L(p_c)\sim sL^{3/4}.
\]

The existing two-size condition

\[
L^{13/4}M_L(p)= (L-1)^{13/4}M_{L-1}(p)
\]

cancels the leading amplitude in `M_L(p_c)` and has shown numerical convergence close to `L^-7`. This is a stronger starting point than an unsupported claim that a simple average of two threshold estimates must cancel the `L^-4` term.

## Primary hypotheses

Each hypothesis is deliberately strong and falsifiable.

### H1 — the useful exact object is the odd matching scaling function

The threshold is only the zero of a richer universal object,

\[
\mathcal M(z)=f(z)-f(-z).
\]

For a fixed torus shape, the normalized odd Taylor invariants

\[
\kappa_{2j+1}
=
\frac{\mathcal M^{(2j+1)}(0)}{\mathcal M'(0)^{2j+1}}
\]

are independent of the nonuniversal metric factor. Reconstructing these invariants across matching pairs and lattice realizations should reveal more exact structure than searching for a closed form of `p_c` alone.

**Falsification:** after controlling finite-size corrections, different matching pairs or microscopic realizations in the same geometry converge to incompatible invariant values.

### H2 — a matching-annihilator hierarchy exists

The known two-size construction removes the leading `L^-13/4` term in `M_L(p_c)` and appears to improve the root bias from approximately `L^-4` to approximately `L^-7`. The strong conjecture is that this can be iterated.

Given measured or theoretically predicted correction exponents `rho_1, rho_2, ...`, choose frozen weights on several sizes so that

\[
\sum_j w_j=1,
\qquad
\sum_j w_j L_j^{-13/4-\rho_r}=0
\]

for the first several correction terms, and solve the resulting weighted matching equation. A three- or four-size annihilator should produce a pseudo-critical sequence beginning at `L^-9`, `L^-10`, or faster.

Weights and exponents must be learned on smaller sizes and frozen before evaluation on larger withheld sizes.

**Falsification:** the apparent `L^-7` improvement does not persist at larger sizes, or higher annihilators amplify noise/model error without out-of-sample gain.

### H3 — a superconvergent multi-observable estimator exists

The annihilator hierarchy can be strengthened by combining independent matching representations of the same exact function:

- cluster-count form;
- horizontal, either, both, and cross-wrapping forms;
- transfer-matrix / critical-polynomial roots;
- open-square `p_med` and `p_cell` estimators.

There exists a frozen combination whose leading correction amplitudes vanish across more than one topology. The requirement of cross-topology prediction is intended to rule out accidental cancellation.

**Falsification:** weights trained on one topology fail on another, or the gain vanishes on withheld sizes.

### H4 — the cubic scaling invariant may be an exact simple constant

Mertens and Ziff found the metric-free ratio

\[
\kappa_3(L)=\frac{M'''_L(p_c)}{M'_L(p_c)^3}
\]

extrapolating numerically to about `-1.67`. A deliberately aggressive candidate is

\[
\boxed{\kappa_3=-\frac53}.
\]

This is not evidence for a closed form of `p_c`; it is a concrete target for the universal scaling function. The candidate must be tested at much larger sizes and with correction-aware fits. Simultaneously estimate `kappa_5`, because a single near-rational number can be accidental whereas a coherent invariant sequence can identify a function class.

**Falsification:** a high-statistics, multi-lattice extrapolation excludes `-5/3` after systematic finite-size uncertainty is included.

### H5 — correction exponents form topology and matching-parity sectors

The observed `4,6,8,...` cylinder sequence and the torus matching-function corrections are not generic polynomial accidents. They arise from irrelevant fields allowed by the selected topological sectors. Matching and topology impose sign/zero rules on their amplitudes.

**Falsification:** amplitudes measured across cylinder, helix, torus, and open square cannot be organized into consistent exponent/sign sectors.

### H6 — the 2015 last digits are an extrapolation artefact

The finite-width roots themselves are accurate; the narrow 2015 uncertainty on the infinite-width intercept is not robust to plausible omitted corrections. The true threshold is expected to lie below `0.59274605079210` and, as a deliberately broad research prior, in

\[
0.592746050789 < p_c < 0.592746050791.
\]

This interval is **not an accepted result**. It is a working range spanning the 2024 disagreement and must be replaced by a model-averaged interval supported by blind prediction and independent estimator families.

**Falsification:** independently reproduced larger-width sequences and a demonstrably predictive extrapolation converge outside this band.

### H7 — no low-complexity elementary closed form

If an exact representation exists, it is more likely to be spectral or constructive than elementary: a limit of critical polynomials, a zero of a Fredholm/spectral determinant, or a distinguished point on a critical manifold. Searching combinations of `pi`, `e`, radicals, and gamma values before stabilizing the numerical interval is low-information work.

**Falsification:** a parameter-free expression survives interval arithmetic, high-height PSLQ controls, and predicts another lattice, correction amplitude, or scaling invariant.

### H8 — finite matching polynomials have no stabilized nontrivial factor

For exactly solvable lattices, every finite-size matching polynomial contains the minimal polynomial of the exact threshold. For square-site percolation, after removal of trivial factors, the relevant factors are conjectured to have increasing degree and no stabilized common divisor:

\[
\gcd\!\left(M_L(p),M_{L'}(p)\right)=1
\]

for generic distinct sufficiently large sizes `L,L'` after normalization.

This does not prove transcendence of the limiting threshold, but it sharply distinguishes an exact finite-cell mechanism from a genuinely limiting one.

**Falsification:** a nontrivial bounded-degree factor persists across sizes and its physical root is compatible with `p_c`.

### H9 — an exactly solvable approximant sequence can be constructed

There exists a sequence of self-dual decorated cells or correlated-bond gadgets `G_k` with exactly computable critical manifolds such that

\[
p_c(G_k)\to p_c(\mathbb Z^2_{\rm site})
\]

and the error decreases faster than ordinary basis-size critical-polynomial convergence. The cell construction should preserve the four-terminal connectivity distribution induced by an occupied square-lattice site to increasing order.

**Falsification:** exhaustive search over bounded cell complexity shows no family with stable monotone or accelerated convergence.

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

On a torus, exact finite-size matching relations become property tests. A single configuration evaluates both members of the matching pair, producing strongly correlated estimates of their difference.

Factor the resulting exact integer polynomials for all feasible small sizes. Record irreducible degrees, coefficient heights, discriminants, root separation, and pairwise gcds.

### Phase C — bidirectional Newman–Ziff campaign

Use one random permutation of the `L^2` sites.

1. Sweep forward, adding black sites and maintaining NN clusters on `G`.
2. Sweep backward through the same permutation, adding white sites and maintaining NN+NNN clusters on `G*`.
3. Pair the microcanonical arrays at complementary occupation counts.
4. Binomially convolve once to obtain `M_L(p)`, its wrapping representations, and analytic derivatives over a dense probability interval.

This bidirectional construction should be superior to independent probability-grid simulations: one permutation contributes to every `p`, preserves the matching covariance, and exposes derivatives needed for Newton roots and universal invariants.

Compare common-random-number, antithetic, and independent couplings. Publish the full covariance matrix of the microcanonical sufficient statistics.

### Phase D — GPU production only after exact validation

GPU work is justified only after exact small-size tests pass. The useful GPU kernel is a bit-packed, many-replica bidirectional sweep, not a naive single-lattice union-find port. Batch thousands of independent permutations and use counter-based random number generation.

Pre-register:

- lattice sizes and aspect ratios;
- number of independent batches and permutations;
- RNG and counter layout;
- root/annihilator observable;
- covariance estimator;
- finite-size model and withheld sizes;
- precision and stopping rules.

### Phase E — transfer matrix / spectral extension

The exact transfer-matrix problem is dominated by state count, irregular sparse transitions, hashing, and memory bandwidth. A single consumer GPU is unlikely to be the decisive resource. Priorities are:

1. canonical compact encoding of non-crossing connectivities;
2. symmetry quotienting;
3. matrix-free dominant-eigenvalue iteration;
4. modular or mixed-precision arithmetic where valid;
5. distributed state sharding and deterministic checkpointing.

On the stated 8-core/16-GB host, reproduce widths well below the frontier and profile state growth. Extending beyond width 24 will probably require a high-memory machine rather than a faster GPU.

### Phase F — operator and gadget searches

In parallel:

- derive which irrelevant fields can contribute to the eigenvalue difference and matching function in the `c=0` logarithmic CFT;
- test topology and matching parity at the amplitude level;
- enumerate small decorated cells and compare their exact terminal-connectivity distributions with those induced by a square-lattice site;
- search for self-dual critical manifolds that generate convergent exact approximants;
- test the factor/gcd conjecture on every exact finite polynomial produced.

## Statistical and numerical rules

1. Never quote more stable digits than survive all reasonable withheld-tail models.
2. Never convert a deterministic fit residual into a confidence interval.
3. Distinguish arithmetic precision, sampling error, covariance, and extrapolation-model error.
4. Freeze model choices before evaluating the largest widths.
5. Publish failed models and negative results.
6. Use interval arithmetic for every claimed exclusion or exact relation.
7. Report condition numbers and noise amplification for every annihilator.
8. A new decimal is secondary; a new cancellation, identity, invariant, exponent rule, rigorous bound, or predictive estimator is primary.

## Hardware map

| Workload | 8 CPU / 16 GB | RTX 5090-class GPU | High-memory CPU node |
|---|---|---|---|
| literature/data audit | sufficient | no | no |
| arbitrary-precision extrapolation | sufficient | no | no |
| exact small torus enumeration | sufficient to moderate sizes | little benefit | useful later |
| bidirectional Newman–Ziff validation | sufficient | optional | no |
| bit-parallel many-replica campaign | slow/moderate | high benefit | optional |
| matching derivatives/invariants | sufficient for validation | high benefit at scale | optional |
| transfer matrix near current frontier | reproduction only | low/uncertain benefit | likely required |
| self-dual gadget enumeration | sufficient for first search | optional batching | useful for large search |
| rigorous interval bounds | sufficient for prototypes | no | useful for scale |

## First decision gates

- **Gate 1:** the audit script must predict withheld finite widths better than a naive fixed-power fit.
- **Gate 2:** exact matching identities must pass exhaustive tests on small tori.
- **Gate 3:** the `L^-4 -> L^-7` annihilator behavior must be reproduced before attempting higher order.
- **Gate 4:** the bidirectional coupling must show a measured variance/throughput advantage.
- **Gate 5:** only then spend GPU time on large-scale Monte Carlo.
- **Gate 6:** only after the numerical interval is stable resume interval-PSLQ searches.

## References driving this program

- J. L. Jacobsen, *Critical points of Potts and O(N) models from eigenvalue identities in periodic Temperley–Lieb algebras*, J. Phys. A 48 (2015) 454003, DOI `10.1088/1751-8113/48/45/454003`.
- S. Mertens and R. M. Ziff, *Percolation in finite matching lattices*, Phys. Rev. E 94 (2016) 062152, DOI `10.1103/PhysRevE.94.062152`.
- S. Mertens, *Exact site-percolation probability on the square lattice*, J. Phys. A 55 (2022) 334002, DOI `10.1088/1751-8121/ac7ed2`.
- Y. Yang and S. Zhou, *Comment on ...*, J. Phys. A 57 (2024) 258001, DOI `10.1088/1751-8121/ad4d2c`.
- J. L. Jacobsen, *Reply to Comment on ...*, J. Phys. A 57 (2024) 258002, DOI `10.1088/1751-8121/ad4d33`.
