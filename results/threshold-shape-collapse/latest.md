# Threshold-distribution shape-collapse contract

All values below are synthetic contract fixtures; they are not simulation results.

- quantile: left-continuous generalized inverse;
- center: `q50`;
- scale: `q75-q25`;
- frozen grid: `['1/20', '1/10', '1/4', '1/2', '3/4', '9/10', '19/20']`.

| comparison | location shift | scale ratio | shape SSE | max absolute shape residual |
|---|---:|---:|---:|---:|
| positive affine | `11` | `3` | `0` | `0` |
| tail deformation | `0` | `1` | `1/2` | `1/2` |

## Interpretation boundary

The contract proves positive-affine invariance and deformation sensitivity on synthetic fixtures only. It provides no covariance, p-value, conformal map, or universality result.
