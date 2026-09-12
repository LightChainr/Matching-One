# N580: the missing covariance recovered, without another experiment

Date: 2026-09-12. Analysis implemented in PR #703.

## Executed, not proposed

GitHub Actions run 34679359201, job 103514945963, completed successfully.
Command: `python scripts/recover_n580_covariance.py`.
Head: `eb78265a7016c99c75f5ecc4a162973c448bf4ee`;
actual tested merge checkout: `dda2a8e2dcb9e5a4a9215cb1f1877a4aa016b189`.
Runner: Python 3.11.16, mpmath 1.4.1. The 14 inference checks and 2 recovery
checks passed in the same job. Full-repository CI is separate from this job.

The reader consumed the EXISTING 100 aligned delete-one rows in
`results/aspect-ladder-n580/shards/rung_r{1,2,4}.json`. It checked replay
seed/offset/sample/batch metadata and recovered the previously stored standard
errors and cov(r1,r2), cov(r1,r4). No Monte Carlo or histogram reconstruction
was run. All historical results and freezes are unchanged.

The result summary is in
`results/research-control-20260912/n580-complete-covariance-summary.json`.
Fields there were extracted from the successful job stdout, then all eight
pure-ray statistics were independently recomputed locally from the recovered
vector/matrix; the largest D difference was zero. That local calculation was
NOT a second shard replay or independent experiment.

## Recovered matrix and decision

The missing cov(r2,r4) is 5.5621485097801135e-9, correlation
0.11001108422637242. The full matrix has rank 3 and condition number 1.96789.
It is not singular; the exact-support bug fixed in #703 does not invalidate
this well-conditioned calculation.

| Pure model, retrospective three-rung test | D, df=2 | nominal Gaussian-reference p | equivalent sigma | declared 3-sigma decision |
|---|---:|---:|---:|---|
| bare aspect ratio | 10.8644599 | 0.00437333 | 2.84990 | not rejected |
| no modulus dependence | 89.8585102 | 3.07236e-20 | 9.21640 | rejected |
| area scaling | 57.5124933 | 3.24579e-13 | 7.28379 | rejected |
| Q4 weight-4 shape | 55.5877105 | 8.49732e-13 | 7.15288 | rejected |
| weight-12 E12 | 132.267324 | 1.89896e-29 | 11.26749 | rejected |
| weight-12 E4 cubed | 132.226117 | 1.93849e-29 | 11.26568 | rejected |
| weight-12 delta | 513.659218 | 2.88614e-112 | not resolved by inverse-erf conversion | rejected |
| weight-8 E8 | 115.804378 | 7.13508e-26 | 10.51803 | rejected |

The previous missing-covariance ambiguity for bare_aspect_ratio is resolved:
it is **not rejected at the declared nominal 3-sigma cutoff**. This does not
confirm that law, make it unique physically, or license a post-hoc 2-sigma
cutoff. The original frozen two-rung comparison remains a different test with
its recorded underpowered verdict; this is the retrospective three-rung test
on the same block, conditional on its pure-amplitude/readout assumptions.

The raw divided difference is

    (m4 - 3*m2 + 2*m1)/6 = -0.00046630613986373587
    SE                    =  0.0001499079748299409
    nominal z             = -3.110615965513005.

This 1-df contrast and the 2-df line test need not cross the same cutoff.
The sign alone cannot reject an arbitrary signed multiple of a convex shape;
any class-level sign claim must additionally declare the amplitude sign.
No continuum operator or percolation threshold is identified.

## What the nuisance sensitivity says, and does not say

Allowing a second column lambda*(-v1,+v2,-v4), lambda=1148/21025, assumes
that spin-8 has the SAME modulus shape and a common ratio to spin-4. Under
that enlarged, unbounded plane the nominal residuals are bare 0.25187/1,
area 6.54806/1, Q4 4.32890/1. The small-ratio bounds matter: at |rho|<=1,
minimum D is 5.72524, 50.81440, 48.88288 respectively. These bounded-cone
minima are not assigned an ordinary chi-square calibration.

This is a sensitivity to a changed model, not a rescue of a rejected pure
ray and not a measurement identifying H8. A bound on |rho| across rungs does
not imply a COMMON rho; the latter is a separate stronger assumption. In
particular the older `required_spin8_ratio` wording that called common rho
an implication of a common bound is not logically valid.

## Decision for the queue and manuscript

Cancel any proposed N580 replay whose sole purpose is this covariance entry.
Use the recovered matrix and clearly separated pure-line/nuisance readings in
P3. Replace 'undetermined because covariance missing' with the narrow decision
above. Keep prospective and retrospective evidence separate. The manuscript
still needs the support/df, amplitude-sign, shared-ratio and prior-art language
corrected; a numerical recovery is not publication readiness.

The scientific next step is not an unbudgeted larger amplitude ladder. First
resolve the angular/Smith-class identifiability in #589 and the within-model
shape question in #622. All p-values here are Gaussian-reference diagnostics
using estimated jackknife covariance, not exact finite-sample coverage.
