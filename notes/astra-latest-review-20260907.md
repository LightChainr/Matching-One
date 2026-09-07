# Latest-delivery audit: restore the measured object before extending the probe

Reviewed 2026-09-07. Repository workspace base `8b5f9d1acb610fd82fc8df6310a96d0fd0ccdf64`;
PR #628 snapshot `a4a1d1f70f82106dd112ee075289e8c0fbaeae4a`. These are
branch deliveries, not merged results. This note does not replace historical artifacts.

## Decision

#612 has returned the bounded N725 forecast score; #613 has returned a
conditional critical-point bridge. Do not order another N725 production.
#617–#625 already contain the new exploratory work. The immediate useful tasks
are two bounded repairs of #628, using existing data, rather than more probes.
The findings below invalidate specific reported controls and interpretations;
they do not settle existence or dimension of the asymptotic shape limit.

## Independently checked failures in PR #628

### 1. The reported bond-duality failure is an implementation failure

In `scripts/probe_invariant_shape/exact_rank_census.py`, the forest path
u→v must be closed with the reverse non-tree edge. The displacement is
`(ax-dx, ay-dy)`, not `(ax+dx, ay+dy)`. The current implementation labels an
elementary contractible plaquette rank 1. Separately, when primal edge i maps
to dual edge pi, dual occupancy is `not occ[i]`, not `not occ[pi]`.

Local enumeration covered all 262,144 L=3 bond configurations. Corrected
forest ranks agree with the existing weighted homology union-find on every
configuration. Correct ranks with geometric dual transport have **zero**
rank-complement failures; pair counts are (0,2):75,460, (1,1):111,224,
(2,0):75,460. The corrected threshold CDF at 1/2 is exactly 1/2.

The submitted edge ordering reproduces its 118,133 failures exactly. Keeping
the bad rank code but sorting the complemented edges changes that to 117,426:
the broken routine is also sensitive to spanning-forest choice. Fixing the
rank but keeping the naive complement still leaves 115,608 failures. Both
repairs are required. This is a finite exact check, not a substitute for the
general topological proof.

### 2. The float CDF does not evaluate the count polynomial

Coefficients already count configurations of cardinality k. Evaluate
`sum(c[k]*p**k*(1-p)**(N-k))`, or first divide c[k] by binomial(N,k) before
using a normalized binomial expectation. `eval_F_float` inserts another
binomial coefficient and leaves its mode scaling unrestored. Site L=3 at
p=1/2 gives **115.9365079365**, versus the exact **43/128**. Bond quantiles
using this evaluator must be regenerated; passing stored-output tests cannot
validate them.

### 3. The same-M/different-F example leaves the probability simplex

For normalized rank-pair masses with rank complement,
`P20+P11+P02=1`, `M=P20-P02`, and `F=P20+P11/2=(1+M)/2`.
Thus equality of the full function M forces equality of F and its quantiles.
It need not determine all three rank-pair probabilities separately.

The proposed `P11 -> 3*P11/2` keeps the other masses unchanged. At L=3,
p=1/2, total mass becomes **593/512**. Its raw M remains -21/64 and raw
F becomes 425/1024, whereas contractual F is 43/128. Normalization gives
M=-168/593 and F=425/1186: M changes. There is no normalized counterexample
here. A valid same-M/different-joint-law example is a different claim and
must have the same F. Equality of only M(1/2), as in #625, is also a different
and much weaker constraint; this audit does not dismiss that question.

### 4. N725 is not the advertised orientation/quantile object

`n725_zflow.cos_four_theta` evaluates cos(2θ). For (26,7), (23,14), true
cos(4θ) values are 0.4958535077 and -0.5780680143. Correct spin-0 weights are
**(0.5382777069, 0.4617222931)**, not (-1.1326530612, 2.1326530612).
The latter weights have spin-4 leakage -1.7944485137.

The script also combines orientation CDFs before inversion. The lineage
object combines separately inverted orientation quantiles. These operations
do not commute. Negative weights additionally remove the general guarantee
that a combined CDF is monotone; this audit does not claim actual
nonmonotonicity in the committed histogram.

Its nine levels replace 0.4,0.6 with 0.25,0.75. The declared nine deciles plus
two anchors require eleven unique levels. A regenerated analysis should use
pooled quantiles, form the nonlinear normalized shape from those quantiles,
and carry aligned delete-one batches through the whole map. The average of
batchwise normalized shapes is a different estimator.

### 5. Covariance rank and tangent-chart angles do not identify a limit model

Full-rank sampling covariance at one N does not reject one-dimensional
variation of the mean shape across N. Counterexample: observations
`y_N = m + f(N)*g + epsilon`, where epsilon is uniformly one of ±e_i in d
dimensions. The mean lies on a line while noise covariance is I/d.

The raw g and normalized-shape displacement ΔZ inhabit different charts.
For `W=Q(b)-Q(a)` and `Z(u)=(Q(u)-Q(a))/W`, use

```
D Z_Q[g](u) = (g(u)-g(a))/W
              - (Q(u)-Q(a))*(g(b)-g(a))/W**2.
```

g currently lives on deciles, so values at 1/4 and 3/4 are unavailable without
an additional interpolation assumption. Either declare that sensitivity or
use anchors present on both grids. Remove actual location/scale directions
in the Q chart, not an unrelated span of 1 and u. A sign-free direction
comparison uses abs(cosine); 159.6 degrees corresponds to 20.4 degrees
between unoriented lines. Tiny exact tori alone cannot establish an
asymptotic mismatch. The existing note already limits the N-regime claim;
the mathematical object must still be corrected before interpreting angles.

## What the recent deliveries can and cannot say

#612 / PR #614 reports N725 amplitude -4.8083463e-4 (SE 2.2941e-6), inside
the predeclared practical 5% tolerance by its interval rule, although displaced
about 5.56 SE from the point forecast. This supports that finite practical
forecast, not a universal exact exponent. Its chart-adjusted curvature is
about 4.25% above prediction; the old 55% discrepancy is not a new-scale
observation. These are reviewed delivery numbers, not recomputed in this
audit; the counterchecks above concern PR #628's separate reanalysis.

The #612 chronology needs a narrow correction: published prediction values
are not exposure to observed outcomes. Document when the estimator/code was
fixed relative to reading any N725 histogram, quantile, score or summary.
Do not upgrade evidence status automatically. New estimands introduced after
that exposure, including this shape exploration, remain exploratory.

#613's bridge remains conditional on the stated geometric and percolation
inputs. The condition shortest period / log N → infinity is sufficient for
the proposed off-critical union-bound argument; it is not an if-and-only-if
necessity statement for quantile convergence. Fixed-width counterexamples
do not prove failure for every O(log N) width. Preserve this distinction
when incorporating the proof delivery.

## Reproduce this audit

Requires Python with NumPy and a separate checkout of the pinned PR #628
snapshot. The reviewed checkout is never modified.

```
python3 scripts/astra628_counterchecks.py --probe /path/to/pr628-checkout --full-bond-census
```

Output: `results/astra628-counterchecks/latest.json`. The full-census flag is
required for the exhaustive assertions. This work creates no production
block, adds no independent statistical evidence vote, and merges no PR.
