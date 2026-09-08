# What remains after using C as the shape clock?

This is post-reveal analysis of the coordinator's N100 three-modulus
production, frozen at `4c1ec50`: 2M common-random replicas per shape, 200
paired batches, fixed `p_ref=.59274605079`. The exact design and competing
three-modulus nulls are `b9e4ea1` and `964d770`. No samples are added.

## A useful gauge, not a field decomposition

For each four-vector Y=(A,E,C,W), set

`D=Y(4i)-Y(2i)`, `U=Y(1/2+i)-Y(2i)`,
`r_C=U_C/D_C`, `R=U-r_C D`.

Thus R_C=0 exactly. After the common offset is removed, the profile pair is
`[0,1,r_C]` (the measured clock) and `[0,0,1]` (the remaining shear response).
The corresponding loadings are D and R. This is saturated on three shapes;
its value is interpretability, not a discovery of two states or fields.

The measured clock is `r_C=-.27798175 +/- .01937689`. Its residual is

| coordinate | R | SE |
|---|---:|---:|
| A | -8.99609e-4 | 3.50155e-4 |
| E | 9.52915e-5 | 2.44801e-4 |
| W | 1.57452e-4 | 6.58762e-5 |

The complete delete-one covariance gives `chi2=17.8234/3`, nominal
`p=.0004783`. This is a mechanism-level localization of the already rejected
common-secant model, not an independent repeated hypothesis test. Choosing
C and a deformation direction after reveal stays explicitly post-reveal.

## The exact thermal interpretation

For threshold births tau1,tau2 on a size-N permutation, integrating each
binomial tail gives `(N+1-tau)/(N+1)`. Therefore for one orientation

`integral A(p) dp=1-2C`, `integral E(p) dp=1-W`.

The constant one cancels in each signed orientation contrast, yielding
`integral P4[A(p)]dp=-2C` and `integral P4[E(p)]dp=-W`. Consequently the
clock quotient has the exact identities

`integral R_A(p)dp=0`, `integral R_E(p)dp=-R_W`.

The first is an imposed zero-area gauge, not an empirical E4 pass. The
second identifies the even-area coordinate with rank-one lifetime deformation.
Its estimate is `-1.57452e-4 +/- 6.58762e-5`, only 2.39 SE; zero even area
still survives at alpha .01. Fixed-p morphology can nevertheless be jointly
resolved through covariance.

The coordinator's separate `7b30648:results/etop-thermal-redistribution/REPORT.md`
(read, not recalculated in this branch; see also its
`notes/etop-clock-redistribution-n100.md`) resolves odd zero-area redistribution:
with z=N^(3/8)(p-p_ref), its first two signed moments are
`-2.93635e-4 +/- 4.6967e-5` and `8.48557e-4 +/- 1.0592e-4`. Its largest lobe
is near p=.32125, not p_ref; hence these finite thermal deformations are not
automatically a critical-field fingerprint. This companion uses the same
histograms and supplies no additional independent evidence.

## Whitening says what is—and is not—identified

The E/W residual correlation is -.76317. Simply reporting their marginal
z-scores hides information. In the frozen E4 residual, the exact Schur
decomposition assigns .21843 chi-square to C and 58.0503 to A/E/W conditional
on an exactly E4-compatible C mean. That stronger premise is different from
estimating the empirical clock; compatibility does not prove it.

For the empirical-clock residual, covariance whitening gives components
`(3.32591,-.23358,2.58982)` in the declared noise eigenbasis. Their squared
norm is the 3-df discrepancy, not three independent discoveries. Source-
selected single-axis nuisance diagnostics give A-only p=.000163, E-only
p=.00873 and W-only p=.03479, without multiplicity correction. Thus a
W-carried deformation still survives at .01; the loading is not uniquely
identified, and E is not established as a separate component.

The covariance-matched training readout is

`Psi=-.0224598 R_A+.222506 R_E+R_W`.

This choice is frozen for a future independent block. Its source SNR is
selected, so it must not be reported with a new one-df p-value.

## Explicit next-block predictions, not a request for production

The least restrictive direction-transfer hypothesis permits an arbitrary
new amplitude but requires two wedges to vanish:

`R_A,100 R_E,next-R_E,100 R_A,next=0`,
`R_A,100 R_W,next-R_W,100 R_A,next=0`.

Keep full source and target uncertainty; do not treat the learned source
direction as exact. A stronger, optional clock-normalized hypothesis predicts

`(R_A,R_E,R_W)/D_C=(.440455,-.0466554,-.0770895)`.

Its complete source covariance is committed. This is a newly proposed
scale/geometry relation, not inferred homogeneity or an implied N400 run.
The next block must declare its lineage/Smith change and compute its own r_C.

## Reproduction and provenance

`results/p267-clock-shape-quotient/input.json` contains the exact original
200x12 batch matrix, mean/covariance, raw hashes, acquisition contract and
SHA256 of the supplied coordinator score. The score was uncommitted when
first read; its exact bytes are now confirmed at source commit
`7b30648be558df0652a7ff22143cc87ed399d042`.

```
python3 -m pip install -r requirements-p267-clock.txt
python3 scripts/p267_clock_shape_quotient.py --source results/p267-clock-shape-quotient/input.json --json results/p267-clock-shape-quotient/score.json --report results/p267-clock-shape-quotient/REPORT.md
python3 -m unittest discover -s tests -p test_p267_clock_shape_quotient.py -v
```

Three new algebra/covariance tests pass. No Monte Carlo or old test suite was
rerun. The report's scientific card and machine predictions retain the common
dependency group and source-selection boundary.
The numeric analysis used NumPy 2.4.6 and SciPy 1.18.0 locally; these optional
dependencies are isolated from the repository's exact-arithmetic runtime.
