# P334: predominantly own-source curvature, with a small local gain roll-off

The prescribed common-batch join of `c48fa360`, `e7473233` and the saved
new64 first-response means is complete. It consumes only these final
artifacts and appends them to `172fbeb1`; the allocation is `fef79403`.

## New paired readout

On the fixed original00 population, both geometries at both sizes have
more negative **own-source C curvature** than other-source pure curvature
or mixed curvature. Differences retain their shared suffix/prefix errors:

| N / receiver | Own − other pure C, ×10^-8 | Own − mixed C, ×10^-8 |
|---|---:|---:|
| 325 / first | −5.241 ± 1.596 | −5.911 ± 1.903 |
| 325 / second | −6.474 ± 1.472 | −3.756 ± 1.489 |
| 425 / first | −7.675 ± 1.320 | −9.326 ± 1.408 |
| 425 / second | −7.322 ± 1.247 | −6.102 ± 1.483 |

The A_ref own-minus-other contrasts are positive:
`(5.267±1.475, 5.945±1.429, 8.439±1.858, 8.290±1.385)×10^-7`
in the same order. E_ref and W do not uniformly reproduce this resolution;
their complete tensors and contrasts are retained rather than selected away.
The mixed tensor is not established to vanish: for example N325 second-C
is `−2.229±0.924 ×10^-8`. The reliable statement is directional predominance
of the own-source curvature in these specified coordinates.

## Local size of the curvature

With H1 and H2 computed on the same new64 original00 block, define
`r=H_oo F/H_o F`. This is the local fractional change of response gain per
unit natural source parameter. The requested ratio of the quadratic Taylor
contribution to the linear one at `t=1/2` is `r/4`.

| N / receiver | A_ref: H2/H1 | C: H2/H1 | A_ref: quadratic/linear at t=1/2 | C: quadratic/linear at t=1/2 |
|---|---:|---:|---:|---:|
| 325 / first | −0.01809 ± 0.00354 | −0.02055 ± 0.00352 | −0.4523 ± 0.0884% | −0.5137 ± 0.0881% |
| 325 / second | −0.01986 ± 0.00393 | −0.02086 ± 0.00339 | −0.4964 ± 0.0981% | −0.5216 ± 0.0847% |
| 425 / first | −0.02849 ± 0.00337 | −0.02774 ± 0.00314 | −0.7123 ± 0.0843% | −0.6935 ± 0.0785% |
| 425 / second | −0.02217 ± 0.00264 | −0.02298 ± 0.00231 | −0.5544 ± 0.0659% | −0.5744 ± 0.0579% |

Thus the resolved curvature opposes the own-source first response for
both named observables: a small **local gain roll-off** in the frozen source
coordinates. This does not extrapolate the finite curve or bound higher
Taylor terms. All eight first-response denominators retain their signs in
the twenty deletions. No ratio was formed for E, W, or cross responses.

## Finite rectangle versus zero-source mixed derivative

The original8 four-corner rectangle has side length one and equals the
integral of H_fs over the source square. Its paired difference from
H_fs(0) is unresolved for all 16 fixed physical/S/D observer readouts,
with magnitude/SE at most 1.41. For physical C the discrepancies are
`(−3.622±37.30, −2.505±30.73, −12.825±43.44, 16.794±33.27)×10^-12`.
This is a precise same-stream comparison, not an exact equality claim.

The new64-minus-old8 H_fs differences retain shared-prefix dependence.
At N325 second-A_ref the difference is `+7.300±2.909 ×10^-7`, and
second-C is `−6.107±2.404 ×10^-8`; the corresponding H4-normalized D
differences are `+1.209±0.455 ×10^-6` and `−1.148±0.437 ×10^-7`.
These correlated conditional-stream changes are kept visible when interpreting
the weak old mixed response. They are not independent population replication
and are not evidence of a changing physical source law.

## Dependence, definitions and reusable output

- C is `(K1+K2)/(2(N+1))`; A_ref is fixed at `p_ref=.59274605079`.
  Own means ff for the first receiver and ss for the second. Mixed fs has
  no extra factor of two. The source coordinates are physical first/second.
- New64 comparisons and ratios use only original00, zero padded to each
  original 1,000-prefix batch. There is no division by the00 prevalence.
  Mixed curvature is identically supported on00, so old8 all/00 agree as
  targets. The finite rectangle keeps the same full-population denominator.
- Source: Hessian `c48fa360`; original8 rectangle `e7473233`; first means
  from `8ad30617:experiments/p334-mechanism-response-20260831/results-extension/`.
  The complete preceding factor is `172fbeb1`. All use the original20 batches.
- Sign alignment: Hessian `factor` is raw-batch centered, so this join uses
  its **negative `LOO_factor`**. Rectangle and previous factors already use
  `sqrt(19/20)*(LOO−mean_LOO)`. Ratios are recomputed in each deletion.
- Observer claim: weak nonlinearity of conditional mean response in a fixed
  commuting exponential source family. It is not a claim about path memory,
  noncommuting actions, a new field, or global state sufficiency.
- No new MC, DP, fork-gzip read, finite weight, fit, determinant, shape test
  or inverse-covariance test was used. The Python research environment only
  supplies the local array aggregation and common-batch uncertainty.

The [full table](../results/p334-mixed-source-curvature-joint/REPORT.md) and
[score](../results/p334-mixed-source-curvature-joint/score.json) retain all
named comparisons. The compressed factors retain every supplied Hessian,
rectangle corner, new64 first-response matrix and predecessor coordinate;
all cross-covariance is recovered as `F.T @ F` (rank at most19 per size).

Reproduce this thin artifact-only join:

```sh
/Users/lc/python-envs/research-py311/bin/python scripts/p334_mixed_source_curvature_joint.py
```
