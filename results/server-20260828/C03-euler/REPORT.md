# C03/P34 Euler and local-motif controls

The previous wrapping-only GLS program is superseded.  The five wrapping
difference channels were configuration-identical, whereas the centered local
Euler quantities provide genuinely distinct zero-mean controls.

## Exact checks

The uncentered identity

```text
N_black - N_white = q_cross + V - E + F0
```

has zero residual for every configuration on axis `L=3` (512 states), diamond
`L=2` (256 states), and the primitive Gaussian `(2,1)` torus (32 states).
For every fixed occupied-site count `K`, the exact configuration averages of
`V`, `E`, `F0`, the declared `(1,1)` diagonal pair count, and the three-site
right-angle motif equal their hypergeometric expectations.  The full rational
audit is in `exact.json`.

## Pilot-frozen Monte Carlo

Each row used 5,000 pilot configurations to freeze coefficients and a disjoint
20,000-configuration evaluation stream.  No ridge was required.  Ratios are
the variance of the best single estimator (`q_cross`) divided by the variance
of the frozen full Euler-plus-motif estimator.

| axis size | N | variance reduction |
|---:|---:|---:|
| 8 | 64 | 2.319x |
| 12 | 144 | 1.852x |
| 16 | 256 | 1.665x |

The method clears the 2x gate at `L=8`, but not at multiple sizes.  Therefore
the current result supports the redesigned control hierarchy but does not yet
authorize a GPU production campaign under the frozen P34 gate.  The local
motifs add only a small improvement over the three Euler controls at these
sizes; further work should target more informative prespecified motifs or a
microcanonical implementation, not duplicate wrapping channels.

Files `L8.json`, `L12.json`, and `L16.json` retain the frozen weights, pilot and
evaluation seeds, covariance matrices, control means, and estimator variances.
