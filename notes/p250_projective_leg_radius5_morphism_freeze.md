# P250 radius-five morphism acquisition freeze

## Decision target

This is not another rank campaign.  It asks which statement survives after the
two hand-specific degree-two annihilators are extended out of sample:

1. identity plus complex conjugation;
2. Alexander reflection plus complex conjugation, with `R^0..R^3` convention;
3. a general `5+5` direct sum with no parameter-free cross-hand line map.

The support question comes first.  If either old annihilator fails on the fresh
degree-five shifts, none of those morphism labels is interpreted.

## Minimal new geometry

Shifting the quadratic relation by `u=(3,0)` introduces exactly
`(5,0),(4,1),(3,2)`.  Their closure under C4 and Alexander reflection is the
entire 20-point Manhattan radius-five shell.  Therefore no smaller spatial set
can score the source orientation, all four Alexander conventions, and the
identity comparator without changing the observer.

Retain both hands and charges 1/2: the object being extended is the two-charge
hand block, and the exact C4 map acts nontrivially on those two charges.  The
minimal frozen payload is consequently 80 complex values, or 160 real
coordinates, per batch.  All are almost free once the two projective-leg
indices have been built.

## Frozen fresh domain

```text
samples: 1,200,000
batches: 400
workers: 16
p: 0.59274605079
seed: 25050510120261250
counters: [0, 1200000)
target: Huawei-CodeBuddy-XPk2PZ
```

The seed is new.  The old and new covariance contributions are computed
separately and added; batches are never paired across independent seeds.

This commit deliberately has `execution_authorized=false`.  The runner refuses
the manifest until a later explicit authorization changes both that flag and
the status.  No simulation is part of this task.

## Power and identifiability

Use the old ten-real-dimensional line residuals only as an exploratory
alternative.  With `lambda=max(T2-10,0)`, the weakest current Alexander
residual is `R1`, `lambda=0.9221`.  The sample factor is 15 and the row
information factor from `12` to `20` hand-block rows is `20/12`, giving total
factor 25.  At alpha 0.01 the noncentral-chi-square forecast is:

| candidate | forecast power at 1.2M |
|---|---:|
| identity+conjugation, if true | size `0.01` |
| Alexander `R0`, current residual | `>0.999999` |
| Alexander `R1`, current residual | `0.822` |
| Alexander `R2`, current residual | `0.966` |
| Alexander `R3`, current residual | `>0.999999` |

This is a planning approximation, not a frozen scientific effect.  Its purpose
is only to avoid an underpowered repeat of the 80k line comparison.

## Support-first score

First hold `q_plus` and `q_minus` fixed from the old radius-four stream and
score their four degree-three shifts against the fresh shell, separately by
hand.  Both must have `p>=0.01`.

Only then augment each `12 x 6` hand matrix with the eight new charge-by-shift
rows, re-extract both lines, and score the five parameter-free candidates with
the independent old-plus-new covariance.  The independent hand extensions are
the general `5+5` fallback.  Multiple surviving maps mean non-identification;
the scorer never chooses a map by best p-value.
