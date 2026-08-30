# Production E_top rank-one elimination

## Answer

The production state plane does **not** eliminate a common radial ray.
The four-lineage joint profile gives `chi2=0.814335`
on `4` df
(`p=0.936516`). In contrast, the
zero-even baseline gives `chi2=196.188`
on `8` df
(`p=4.05921e-38`) and is eliminated.
Thus E_top is resolved, but its observed parent-child evolution is compatible
with the same two-component ray rescaled by a lineage-specific nuisance lambda.

## Lineage profiles

| lineage | determinant | lambda | min chi2 / 1 df | p | rank-one | E_top=0 p / 2 df |
|---|---:|---:|---:|---:|---|---:|
| P49-N130-to-N170 | 4.33473e-09 | 0.862635 | 0.107203 | 0.743351 | survives | 1.7001e-09 |
| P43-N185-to-N265 | -1.8554e-09 | 0.754354 | 0.270415 | 0.603054 | survives | 3.65016e-25 |
| P50-N145-to-N290 | 5.47096e-09 | 0.470125 | 0.4342 | 0.509935 | survives | 3.19029e-08 |
| P57-N325-to-N425 | 1.39146e-10 | 1.04106 | 0.00251722 | 0.959985 | survives | 0.0126386 |

Each minimum is global: the scorer enumerates every real root of the
degree-at-most-six analytic profile derivative and the compact boundary
`|lambda|=infinity`. Parent truth and lambda are nuisance parameters, leaving
one degree of freedom per lineage.

## Covariance and sensitivity

Primary covariance: `covariance_intrinsic_center_first_order_influence`.
It is the full same-batch A_top/E_top covariance with the displayed
first-order influence correction for each fitted matching center.
The fixed-center exact-batch sensitivity gives joint rank-one `chi2=0.814322`
(`p=0.936518`), so the decision is unchanged.
Cross-size blocks are block diagonal by the production independence contract;
the four dependency-group labels remain explicit rather than being treated as
eight unrelated evidence rows.

## Scientific card

- MECHANISM SPACE: test whether radial evolution of the two-component A_top/E_top state needs a rotating second direction, rather than whether E_top exists.
- RESULT: all four common-ray rank-one lineages survive separately and jointly, while the global E_top=0 baseline is eliminated.
- NOT PROVED: survival does not identify an asymptotic field or make the center influence exact; no stored cross-size CRN covariance is claimed.
- OBSERVER-SECTOR-SOURCE-GEOMETRY: P4(A_top,E_top) | Alexander odd/even rank plane | threshold-rank source | four explicit Gaussian parent-child lineages.
- DEPENDENCY GROUPS: P49, P43, P50 and P57 are kept as four named production blocks; within each size A_top/E_top covariance is full.
- UPWEIGHT OBSERVATION: a future child with a frozen, large ray-rotation prediction or a stored cross-size CRN covariance can turn this surviving nuisance-ray test into a sharper discriminator.

## Reproduction

```bash
python3 scripts/etop_rank1_elimination.py --source raw --format json --output results/etop-rank1-elimination/latest.json
python3 scripts/etop_rank1_elimination.py --format markdown --output results/etop-rank1-elimination/REPORT.md
python3 -m unittest discover -s tests -p 'test_etop_rank1_elimination.py'
```

The first command rebuilds the rank-plane values and covariance from the
same-batch raw histograms/moments, checks equality with the committed
crosswalk, and records SHA-256 hashes for all 16 input files.

## Claim boundary

intrinsic-center covariance is first-order influence; fixed-center covariance is exact for batch estimators conditional on the plug-in p0. The rank-one result is a model
survival statement at these four production edges, not an asymptotic field
identification. The E_top=0 rejection establishes a resolved companion
direction, not that the companion evolves independently of A_top.
