# P334 current-k0 geometry pilot freeze

## Decision target

The previous production result at `742a8b0` showed that exact saturation index,
birth-local marks, and the future second-birth site class leave essentially all
of the birth-age slope intact.  They do not observe the rank-one configuration
that actually exists at the prediction layer.  This pilot records that current
state at the already frozen intrinsic layers `k0=193` for N325 and `k0=252`
for N425.

The primary question is deliberately narrow: can four cheap, interpretable
current-geometry summaries -- essential occupied mass, essential carrier
count, occupied frontier, and vacant frontier -- absorb the line-controlled
age slope, and does the N325 geometry/hazard map transfer to held-out N425?

## Exact one-step ceiling

For every rank-one state, `H2` counts vacant sites whose insertion makes the
ambient image rank two.  Because the next member of a random permutation is
uniform among the `N-k0` vacant sites,

```text
P(K2=k0+1 | current configuration) = H2/(N-k0).
```

This makes H2 a one-step sufficient statistic almost by definition.  It is
retained as an exact calibration ceiling and decomposed into same-essential-
carrier (`theta`), joined rank-zero/rank-one (`figure8`), and separate-carrier
triggers.  It is not presented as an explanatory discovery.  The informative
comparison is whether the substantially cheaper size/frontier/carrier vector
already removes the age association.

## Frozen production and score

Both sizes use 50,000 paths in 50 batches.  The two orientations of a size
share the exact counter stream; the two sizes use disjoint seeds and replica
domains.  N325 is the source fit and N425 is held out for cross-size transfer.
The full machine-readable contract, matrices, counters, fields, model ladder,
and gates are in `analysis/p334_current_k0_geometry_pilot_freeze.json`.

All regressions remain linear, line fixed-effect partial scores.  No nonlinear
basis, flexible learner, adaptive sample increase, or post-reveal covariate
selection is allowed.  Delete-one batches preserve the common-randomness
orientation covariance.  The N325 and N425 blocks are independent.

The row also retains the realized `K2`.  This was added before any simulation
row existed: both first launchers stopped with exit 127 because `/usr/bin/time`
was absent, before invoking the runner.  After the frozen age score only, K1/K2
may support a descriptive covariance/explained-trace crosswalk to the external
temporal mode-2/mode-3 result at `5a7f2d9`.  That crosswalk cannot change the
primary pilot decision.

## Claim boundary

Persistence after the cheap vector means only that `(k,rank,ell,cheap current
geometry)` is not lumpable for this one-step hazard.  It does not distinguish
intrinsic history from omitted current shape.  Conversely, disappearance after
H2 is expected from the sampling identity and must not be advertised as a
mechanism confirmation.
