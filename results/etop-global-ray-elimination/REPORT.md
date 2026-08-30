# Cross-geometry E_top global-ray elimination

## Answer

The one-global-ray null is **global_ray_eliminated** against
the four-lineage-ray alternative. Its absolute goodness-of-fit is
`chi2=15.5147` on `7` df
(`p=0.0299392`), which alone survives at alpha=0.01. Relative to the four lineage-specific
rays, the gauge-free nested penalty is `Delta chi2=14.7004`
on `3` df (`p=0.00209144`), so the stronger geometry-sharing
constraint is `global_ray_eliminated`.

The fitted global direction is `-30.619` degrees
from the positive A_top axis (`E/A=-0.591835`).
The largest independent-geometry incompatibility is `['P43', 'P50']`
with `Delta chi2=10.1371` on 1 df
and absolute unit-ray determinant `0.389471`.

## Lineage directions

| dependency group | angle from A (deg) | E/A | lineage chi2 / 1 df |
|---|---:|---:|---:|
| P49 | -23.356 | -0.431822 | 0.107203 |
| P43 | -44.661 | -0.988237 | 0.270415 |
| P50 | -21.739 | -0.398746 | 0.4342 |
| P57 | -51.871 | -1.27401 | 0.00251722 |

## Pairwise geometry contrasts

| pair | Delta chi2 / 1 df | p | abs det(unit rays) | angle difference (deg) |
|---|---:|---:|---:|---:|
| ['P43', 'P50'] | 10.1371 | 0.00145317 | 0.389471 | 22.9216 |
| ['P49', 'P43'] | 8.94202 | 0.00278685 | 0.363337 | 21.3053 |
| ['P50', 'P57'] | 2.65914 | 0.102956 | 0.501983 | 30.1313 |
| ['P49', 'P57'] | 2.34366 | 0.125793 | 0.477389 | 28.515 |
| ['P43', 'P57'] | 0.105024 | 0.745881 | 0.125501 | 7.20971 |
| ['P49', 'P50'] | 0.0740945 | 0.785466 | 0.0282055 | 1.61627 |

## Leave-one-lineage-out sensitivity

| held out | profiled predictive chi2 / 2 df | p | decision |
|---|---:|---:|---|
| P49 | 3.83855 | 0.146713 | survives |
| P43 | 12.3078 | 0.00212515 | eliminated |
| P50 | 5.39353 | 0.0674233 | survives |
| P57 | 1.30263 | 0.521361 | survives |

The held-out statistic profiles the ray against the training likelihood
penalty rather than pretending that the training direction is known exactly.
The fixed-center exact-batch sensitivity gives global-ray
`chi2=15.514` and nested
`Delta chi2=14.6997`;
the decision is unchanged.

## Scientific card

- MECHANISM SPACE: one universal cross-geometry A_top/E_top ray versus four geometry-family rays.
- RESULT: the test changes only the geometry-sharing constraint; every size retains a free scalar amplitude and no field name is assigned.
- NOT PROVED: ray incompatibility does not identify a field, operator, exponent or asymptotic geometry law.
- OBSERVER-SECTOR-SOURCE-GEOMETRY: P4(A_top,E_top) | Alexander odd/even state plane | threshold-rank source | P49/P43/P50/P57 production geometries.
- DEPENDENCY GROUPS: four named blocks remain explicit; the global statistic is their block-diagonal joint profile, not eight independent narrative votes.
- UPWEIGHT OBSERVATION: a new independent geometry should be frozen near the largest pairwise determinant direction, then scored against the source-profiled global ray.

## Reproduction

```bash
python3 scripts/etop_global_ray_elimination.py --format json --output results/etop-global-ray-elimination/latest.json
python3 scripts/etop_global_ray_elimination.py --format markdown --output results/etop-global-ray-elimination/REPORT.md
python3 -m unittest discover -s tests -p 'test_etop_global_ray_elimination.py'
```

## Claim boundary

This is a cross-geometry model elimination on already revealed production
blocks. It neither assigns the incompatible directions to fields nor treats
the four dependency groups as separate discoveries. Intrinsic-center covariance
is first-order; the fixed-center exact-batch row is the declared sensitivity.
