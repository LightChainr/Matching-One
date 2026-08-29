# Production threshold-shape collapse score

The frozen q05/q10/q25/q50/q75/q90/q95 contract is applied to the equal mixture of
Kminus/Kplus and both Gaussian orientations. Central standardized coordinates are
identically fixed, so the covariance score uses the four tail coordinates.

| lineage | location shift | scale ratio | shape SSE | max | tail chi-square | log10 p |
|---|---:|---:|---:|---:|---:|---:|
| N65_to_N130_norm2 | 1.698446046e-05 | 0.7744517053 | 0.0008234603669 | 0.02537476288 | 518694/4 | -112628 |
| N85_to_N170_norm2 | -1.837997036e-05 | 0.7739617347 | 0.0005167091184 | 0.02009137703 | 301790/4 | -65527.7 |
| N65_to_N325_norm5 | 1.154390567e-05 | 0.5510517693 | 0.002667981581 | 0.04566204907 | 3.19199e+06/4 | -693125 |
| N85_to_N425_norm5 | -1.407440759e-05 | 0.5502491476 | 0.001962410474 | 0.03913297309 | 1.98114e+06/4 | -430192 |

## Post-reveal shape-flow direction

- same-parent norm-2 versus norm-5 cosine: N65 `0.999994618`, N85 `0.994738715`;
- norm-5/norm-2 projection scale: N65 `1.799980`, N85 `1.938567`;
- effective positive correction power: N65 `0.631527`, N85 `0.429020` (5/8 predicts ratio `1.804095`).

## Interpretation

A small raw standardized residual is not automatically statistical collapse: the high-statistic archives resolve finite-size tail deformation. Compare the two independent norm-5 lineages as one mechanism pattern, not as additive confirmations.

This is retrospective production evidence from existing archives. It is not a
prospective universality test and does not compare distinct microscopic models.
