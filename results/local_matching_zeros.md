# Local matching-polynomial zeros near the physical root

Source: `results/exact-zero-map-pilot/roots.csv` via `scripts/local_matching_zeros.py`.
Claim level: C1 descriptive. These are matching-polynomial zeros, not Fisher zeros.

Metrics frozen in `predictions/local_matching_zero_metrics_20260829.yaml` before
any future exact L=6 result is inspected. Named diagnostics are not fit targets.

## Local catalogue

| geometry | L | physical root | nearest nonreal | `L^{3/4} Im` | `L^4 |z-p*|` | Re in (0,1) |
|---|---:|---:|---|---:|---:|:---:|
| axis | 1 | 0.500000000000 | all real | — | — | — |
| axis | 2 | 0.541196100146 | all real | — | — | — |
| axis | 3 | 0.586511455113 | 1.149886+0.188111i | 0.429 | 48.1 | no |
| axis | 4 | 0.590672112331 | 1.035629+0.258100i | 0.730 | 131.7 | no |
| diamond | 1 | 0.707106781187 | all real | — | — | — |
| diamond | 2 | 0.604563277854 | -0.067168+0.407454i | 0.685 | 12.6 | no |
| diamond | 3 | 0.594252321169 | 0.996090+0.308717i | 0.704 | 41.0 | yes |

## Named-diagnostic stability

```text
L^{3/4} |Im|  values: 0.429, 0.730, 0.685, 0.704
L^4 |z-p*|   values: 48.1, 131.7, 12.6, 41.0
complex-zero scaling route: closed
```

Axis `L^{3/4} Im` jumps from 0.429 (L=3) to 0.730 (L=4). Distance scaled
by `L^4` is 48, 132, 12, 41. Nearest nonreal roots usually have real part
outside `(0,1)`. No stable local power is present at available exact sizes.

## Boundary

Do not invent further cloud statistics. Do not score a future L=6 polynomial
against a power fitted here. A later exact size may reopen the route only by
the frozen metrics above, not by a new ad hoc summary.
