# P218 production empirical coalescence diagnostic

## Answer

The committed production block does not resolve a stable second response direction at the explicit 3-sigma diagnostic gate.
The source-input determinant has `|z|=0.645`
and the held-out-parent determinant has `|z|=1.458`.
Therefore plug-in eigenvalue gaps and eigenvector geometry below are reported as
the requested diagnostic, but are not promoted to an empirical Jordan certificate
when the response-rank gate fails.
The plug-in sequence also lacks the requested joint signature: the gap grows
from 0.635 to 1.128 on the two source transitions while J2/J1 worsens from
0.305 to 3.049. That is not positive coalescence evidence, but the rank failure
prevents treating it as a powered elimination.

## Covariance-whitened response geometry

| generation | sizes | response angle (deg) | response cond | determinant z |
|---:|---|---:|---:|---:|
| 0 | [65, 85] | 0.9672 | 118.9 | 4.201 |
| 1 | [130, 170] | 0.3412 | 336.8 | 0.6451 |
| 2 | [260, 340] | 2.798 | 41.33 | 1.458 |
| 3 | [520, 680] | 12.15 | 9.401 | 0.7032 |

## Plug-in transfer geometry

| transition | relative eigenvalue gap | eigenvector angle (deg) | eigenbasis cond | J2/J1 |
|---|---:|---:|---:|---:|
| generation_0_to_1 | 0.6347 | 60.44373457763565 | 1.7166646674416493 | 0.305 |
| generation_1_to_2 | 1.128 | 41.41086633581216 | 2.645664450160727 | 3.049 |
| generation_2_to_3 | 1.226 | 58.74072432216383 | 1.7768625696215339 | 2.584 |

## Held-out N520/N680 model table

| source-frozen class | held-out chi-square / rank | p | decision | nearest optimistic separation chi-square |
|---|---:|---:|---|---:|
| normal_diagonalizable | 0.3839 / 4 | 0.9838 | underpowered | 3.572 |
| rank2_Jordan | 1.059 / 4 | 0.9007 | underpowered | 4.471 |
| generic_2x2 | 0.2614 / 4 | 0.9922 | underpowered | 3.572 |

`normal_diagonalizable` is the covariance-whitened stable-eigenbasis
representative; `rank2_Jordan` is `lambda I + s u(Ju)^T`; and
`generic_2x2` is the saturated source map. A generic real 2x2 matrix with
distinct eigenvalues is already diagonalizable, so the first and third rows
test normal/stable versus unrestricted mixing, not disjoint algebraic sets.
The pairwise separation uses target covariance only and is therefore an
optimistic upper bound; values below 9 certify lack of 3-sigma discrimination.

## Scientific card

- MECHANISM SPACE: empirical 2D Jordan coalescence versus stable diagonalizable or generic mixing in the leading-H4/S-prime state.
- RESULT: report plug-in gap, response angle, eigenbasis condition and minimal-polynomial ratio, but gate interpretation on the source response determinant and held-out rival separation.
- NOT PROVED: the empirical curve-transfer matrix is not the microscopic Potts transfer matrix; target reuse is one mechanism analysis, not an independent evidence vote.
- OBSERVER-SECTOR-SOURCE-GEOMETRY: (N^(13/8)P4[D], U) | thermal H4/S-prime | threshold-rank curves | two dyadic Gaussian lineages.
- DEPENDENCY GROUP: PR277/P154 N65 through N680; scalar and jet scorers are correlated views of the same histograms.
- UPWEIGHT OBSERVATION: resolve the source response determinant and a held-out pairwise model separation above chi2=9 using the same semantic two-coordinate block.

## Reproduction

```bash
python3 scripts/analyze_p218_production_coalescence.py --output results/p218-production-coalescence/latest.json --markdown results/p218-production-coalescence/REPORT.md
python3 -m unittest discover -s tests -p 'test_p218_production_coalescence.py'
```

No new Monte Carlo or exact oracle is used. The script reconstructs every
point and covariance from the committed production histograms and records
their hashes.

## Claim boundary

a second response direction is not resolved at the 3-sigma production gate; plug-in eigen geometry is descriptive only. The transfer diagnostics
are empirical response maps in a frozen curve coordinate, not microscopic
transfer-matrix eigenstates. N520/N680 are held out from basis and source-map
selection but are already revealed production blocks.
