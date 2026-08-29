# P205 quotient-character prism

## Result

The next P205 experiment should not be more N325/N425 sampling.  It should be
a three-coordinate mechanism selector:

| N | first / second | Smith first / second | H4 | H8 | H12 |
|---:|:---|:---|:---:|:---:|:---:|
| 25 | `(5,0)` / `(4,3)` | `(5,5)` / `(1,25)` | + | + | + |
| 50 | `(7,1)` / `(5,5)` | `(1,50)` / `(5,10)` | + | - | + |
| 125 | `(11,2)` / `(10,5)` | `(1,125)` / `(5,25)` | + | + | - |

Here each sign is the exact sign of
`cos(s theta_first)-cos(s theta_second)`.  N50 flips only H8 and N125 flips
only H12.  This is a Hadamard-like angular code: H4, H8 and H12 are different
one-dimensional lines in a three-dimensional measurement space, rather than
three similar residuals differentiated mostly by noise.

The exact matrices are Gaussian period matrices
`[[a,-b],[b,a]]`.  In row order they are:

- N25: `[[5,0],[0,5]]` versus `[[4,-3],[3,4]]`;
- N50: `[[7,-1],[1,7]]` versus `[[5,-5],[5,5]]`;
- N125: `[[11,-2],[2,11]]` versus `[[10,-5],[5,10]]`.

Every pair changes Smith type at fixed area.  The design therefore probes the
angular mechanism and quotient arithmetic together, on purpose.

## Why the completed N325/N425 block is not enough

The completed same-parent block was a successful quotient-control experiment,
but not a harmonic selector.  Each N supplied one affine-null residual, the raw
M noise was comparable to the fixed-model residuals, and all three harmonics
survived.  The observed separations were only

- `chi2(H4)-chi2(H8) = 2.2528963`;
- `chi2(H12)-chi2(H8) = 0.2419415`.

Adding samples would shrink the same two coordinates.  The prism changes the
geometry of the question: it adds independent sign flips chosen before target
reveal.

## Maximin objective

The generator enumerates all canonical sums-of-two-squares representations
with `25 <= N <= 2000`, pairs equal-area representations only when their Smith
invariants differ, and constructs for each harmonic

`f_s(N) = N^(-13/8) [cos(s theta_1)-cos(s theta_2)]`.

Each harmonic has one shared unknown amplitude, so it is a line, not a fixed
point.  For every pair of harmonic lines the objective takes both directed
noise-whitened projection distances and then the smaller one.  The campaign
score is the minimum over H4/H8/H12, divided by measured CPU cost.

Noise and runtime come from the completed P205 production rather than a generic
proxy:

- mean 10M pair-contrast variance: `5.314603523220254e-8`;
- mean 10M pair-contrast SE: `2.3053423874167268e-4`;
- CPU seconds per site update: `2.0709907904351133e-7`.

At 10M per pair the selected prism has planning maximin noncentrality
`5.7484551`; its two-pair optimum has only `0.7581049`.  The third coordinate is
the experiment, not redundancy.  At the frozen 12M per pair the planning value
is `6.8981462`, just beyond the df=2 95% reference threshold.  Estimated cost is
497 CPU seconds, or about 39 seconds on two 8-thread lanes.  The reference
amplitude `0.7885` only prices the sample count; it is not a target mean.

## Fingerprints and expected ratios

Relative to the N25 coordinate, including the frozen radial factor, the exact-
geometry planning fingerprints are approximately:

| model | N25 | N50 | N125 |
|:---|---:|---:|---:|
| H4 | +1 | +0.324210 | +0.040960 |
| H8 | +1 | -0.324210 | +0.123349 |
| H12 | +1 | +0.324210 | -0.110035 |

Those signed ratios are the mechanism prediction.  Only one overall amplitude
is fitted per row.  No exponent, quotient offset, or correction coefficient is
available to rescue a failed row.

## Finite-size risk and the bridge designs

N25 is intentionally aggressive.  This campaign asks whether the harmonic
character is already visible as a mechanism; it does not claim that N25 is in
the asymptotic regime.  If the one-amplitude model fails for all three rows, do
not add a fitted correction.  Use one of two already-enumerated bridges:

- information-first medium-N: N338 `(17,7)/(13,13)`, N400
  `(20,0)/(16,12)`, N500 `(22,4)/(20,10)`; it has the same three sign roles
  after reordering, but the P205-calibrated 10M maximin is only `0.0673542`;
- sign-isomorphic scale bridge: N400/N450/N500 with sign code
  `+++ / +-+ / ++-`; its 10M maximin is only `0.0372987`.

The small prism is roughly two orders of magnitude more informative at the
observed noise scale.  The bridge exists to distinguish a small-N breakdown
from a failure of the angular character, not to dilute the first experiment.

## Minimal sufficient statistics

For each pair, preserve synchronized batches and both orientations' marginal
`K_plus` and `K_minus` threshold-rank histograms.  That is enough to reconstruct
the two fixed-p matching values, their common-field contrast, and its covariance.
Independent seeds across N make the final 3x3 contrast covariance diagonal.

Roots, derivatives, joint plus/minus moments, fitted exponents and new full-
curve fits are not needed for the primary score.  The full frozen execution and
score contract is in
`predictions/p205_quotient_character_prism_20260829.yaml`; no production run was
started as part of this design work.
