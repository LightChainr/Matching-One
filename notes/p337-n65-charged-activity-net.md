# The N65 charged contrast is net timing, not common activity

Status: fixed linear re-expression of the same `1714141` N65 20k archive and
the complete covariance published in `83207c7`.  There is no new simulation,
source choice, projective basis, reference point, or dependency group.

## Frozen decomposition

For each A/B1 and D/B2 source, replace `(W,J_birth,J_exit)` by

```text
W,
J_minus = J_birth-J_exit = dW/dp,
J_plus  = J_birth+J_exit.
```

`J_minus` is the net derivative current; `J_plus` is the positive common
activity.  This is an invertible predeclared linear map, so the full three-way
quadratic is unchanged.  Every covariance entry is propagated from the 20
aligned batches rather than reconstructed from marginal errors.

## Result

For the within-batch `7+4i` minus `8+i` contrast:

| channel | coordinate | value | SE | z | marginal quadratic |
|---|---|---:|---:|---:|---:|
| A | `W` | +0.00214027619 | 0.00237873 | +0.900 | 0.810 |
| A | `J_minus` | +0.09701583219 | 0.0286164 | +3.390 | 11.494 |
| A | `J_plus` | +0.02683092280 | 0.0501771 | +0.535 | 0.286 |
| D | `W` | -0.00006569222 | 0.00136975 | -0.048 | 0.002 |
| D | `J_minus` | -0.01059211221 | 0.00903561 | -1.172 | 1.374 |
| D | `J_plus` | -0.01285117100 | 0.0379713 | -0.338 | 0.115 |

The A full score remains `12.153/3 df`.  Net alone carries `11.494`; activity
alone carries `0.286`.  Once `(W,J_minus)` is present, adding `J_plus` increases
the quadratic by only `0.236`.  The A response is therefore localized to net
source-sink timing, not to an overall increase of charged traffic.

D remains unresolved: its full score is `1.509/3 df`, with neither net nor
activity producing a useful orientation discriminator.

## Covariance eigenmodes

Raw covariance eigenvectors are scale dominated because `W` and the currents
have different units.  The certificate therefore retains the raw covariance
but diagnoses eigenmodes in the standardized correlation matrix, ordered as
`(W,J_minus,J_plus)`.

For A:

| eigenvalue | standardized vector | quadratic contribution |
|---:|---|---:|
| 0.5447 | `( +0.717, -0.669, -0.197 )` | 5.477 |
| 0.9675 | `( -0.076, -0.356, +0.931 )` | 0.628 |
| 1.4878 | `( +0.693, +0.653, +0.307 )` | 6.048 |

The A score is split between the two modes carrying a large `J_minus`
component.  The activity-dominated middle mode contributes only `0.628`.
This agrees with the direct marginal and conditional decomposition: net is
the stable localization, not an artifact of ignoring covariance.

For D:

| eigenvalue | standardized vector | quadratic contribution |
|---:|---|---:|
| 0.1574 | `( -0.676, -0.117, +0.727 )` | 0.037 |
| 0.7634 | `( -0.390, +0.895, -0.219 )` | 1.197 |
| 2.0792 | `( +0.625, +0.431, +0.651 )` | 0.275 |

No D eigenmode carries a substantial score.  The small eigenvalue is driven by
strong W/activity correlation (`rho=0.833`) but the observed contrast has
almost no projection on it, so covariance compression does not manufacture a
D signal.

## Crosswalk to #334

This makes the relationship to the #334 **common activity with
orbit-composition counterflow** precise:

```text
#334 common activity       <-> J_plus,
#334 source-sink residual  <-> J_minus.
```

The morphology agrees: two large positive boundary currents share a common
activity, while the discriminating information survives in their much smaller
net residual.  In the present archive, common A activity cancels between the
two microscopic orientations and the A net timing remains.

The observers are not identical.  This A/B1 charge compares the two lines
inside the axis orbit; D/B2 compares the diagonal pair.  The exact `77aa3fe`
N13/N17 analysis compares the axis and diagonal orbit totals in the H4 sector.
It is a mechanism crosswalk, not a statistical replication.

Dependency accounting is explicit:

- `381984d` conditional line sorting and this decomposition reuse the same
  `1714141` N65 block.  Their quadratics cannot be added.
- `77aa3fe` uses exact N13/N17 tables and is a different dependency group.
  Its structural classification is compared qualitatively; no significance or
  quadratic is pooled with the N65 score.

## Boundary and reproduction

The decomposition localizes the already observed A charged orientation score
to a net timing current.  It does not add evidence beyond the 20k archive,
establish large-N persistence, or identify a continuum operator.

```bash
python3 scripts/score_n65_charged_activity_net.py \
  --reveal results/p337-n65-charged-source-reveal/latest.json \
  --json results/p337-n65-charged-activity-net/latest.json \
  --markdown results/p337-n65-charged-activity-net/latest.md

python3 -m unittest discover -s tests \
  -p 'test_n65_charged_activity_net.py'
```
