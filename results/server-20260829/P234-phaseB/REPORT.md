# P234 Phase B: cross-cutoff top-partner shear

Status: completed production result for cutoff denominators `d=8,12,16` at
`L=64,96,128,192`, with 100,000 configurations and 100 batches per point.

## Chronology

The `d=8` line was already archived as Phase A. Commit `dcfbd8c` froze the
natural-realized-cutoff transform and six-parameter joint score after seeing
all four `d=8` sizes and only the `L=64` previews at `d=12,16`. At that point
the preview was `chi2=6.747/12`, `p=0.874`, with
`kappa_proxy=2.906 +/- 0.545`. The `d=12,16` results at `L=96,128,192` were
therefore held out from the scorer and model.

A separate normalization audit in commit `ceb7c6e`, completed before the last
`d=12, L=192` file revealed, proposed `8/3` as a bold amplitude conjecture in
the declared natural-energy gauge. It did not derive that value from the
thermal exponent and did not change the frozen scorer.

## Frozen primary score

The score uses all 12 JSON inputs and the model

```text
LL = c/L
LD = B + c/(L delta_realized)
DD = D + s log(2 delta_declared) + c/(L delta_realized).
```

The final held-out result is compatible with the frozen model:

```text
chi2 = 22.413980 / 30 df
p    = 0.838565
```

The primary fitted cutoff-shear quantities, using the full same-stream
covariance, are

```text
LD continuum B       =  0.207802 +/- 0.022880
DD log cutoff slope  = -1.048469 +/- 0.071886
kappa_proxy          =  2.522757 +/- 0.321764
```

Here `kappa_proxy=-s/(2B)`. It remains in the connection-probability field
gauge; it is not by itself a gauge-invariant universal coupling.

## Post-reveal diagnostics

These comparisons were calculated after the 12-point score revealed and are
not additional frozen tests. The fitted log slope is `-14.59` standard errors
from zero. Against the separately prerevealed `8/3` amplitude conjecture,

```text
kappa_proxy - 8/3 = -0.143910
normal z           = -0.447
two-sided p        =  0.655
95% normal CI      = [1.892, 3.153].
```

Thus the final estimate is consistent with `8/3`, but this single-lattice,
gauge-fixed agreement does not establish the conjecture. A gauge-invariant
extra statistic and cross-lattice reproduction remain necessary.

## Reproducibility

- Frozen scorer commit: `dcfbd8c`.
- Inputs: 12 JSON/CSV pairs, one per `(d,L)`, with distinct seeds.
- Every input has 100,000 samples split into 100 batches.
- Machine-readable result: `cross_cutoff_shear_score.json`.
- SHA-256 manifests: this directory's `checksums.sha256` for `d=12,16` and
  the Phase A manifest for the reused `d=8` line.
