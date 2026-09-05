# P3 evidence tables

**Generated** by `scripts/p3_manuscript_evidence_table.py`. Do not edit by hand;
`tests/test_p3_manuscript_evidence_table.py` fails if this drifts from
`results/p3-projective-inference-manuscript/latest.json`.

Channel `P4_S_prime`, site count 580, 100 batches.

## T1 — the three-rung response

| rung `r` | value | standard error | &sigma; from zero |
|---:|---:|---:|---:|
| 1 | 9.016433e-04 | 2.4914e-04 | 3.62 |
| 2 | 2.910977e-03 | 2.5911e-04 | 11.23 |
| 4 | 4.131808e-03 | 1.9513e-04 | 21.17 |

The frozen design nominated `r=1` as the denominator. Its 3&sigma; Fieller
intervals for the two ratios are

| ratio | lower | upper | width |
|---|---:|---:|---:|
| `r2_over_r1` | 1.636 | 18.908 | 17.272 |
| `r4_over_r1` | 2.395 | 27.468 | 25.074 |

## T2 — the two-entry control

Projective statistic restricted to the two rungs the frozen test used (`r4_over_r1`), against Fieller's *z* squared on the same pair.

| competitor | Fieller *z* | *z*&sup2; | projective *D* (2 entries) | relative deviation |
|---|---:|---:|---:|---:|
| `bare_aspect_ratio` | 0.5019 | 0.25187 | 0.25187 | 1.5e-15 |
| `no_modulus_dependence` | 9.4772 | 89.81806 | 89.81806 | 1.6e-16 |
| `plain_area_scaling` | -2.5589 | 6.54806 | 6.54806 | 1.4e-16 |
| `q4_jordan_weight4` | -2.0806 | 4.32890 | 4.32890 | 2.1e-16 |
| `weight12_E12` | -3.6108 | 13.03818 | 13.03818 | 2.7e-16 |
| `weight12_E4_cubed` | -3.6062 | 13.00465 | 13.00465 | 0.0e+00 |
| `weight12_delta` | 21.1745 | 448.36140 | 448.36140 | 1.3e-16 |
| `weight8_E8` | -3.4780 | 12.09621 | 12.09621 | 1.5e-16 |

Largest relative deviation: **1.5e-15**.

## T3 — frozen verdict against projective verdict

| competitor | ray | frozen &sigma; (1 entry) | projective &sigma; (3 rungs) | frozen | projective | changed |
|---|---|---:|---:|---|---|---|
| `bare_aspect_ratio` | (1.0000, 2.0000, 4.0000) | 0.50 | 2.74 | compatible | compatible | no |
| `no_modulus_dependence` | (1.0000, 1.0000, 1.0000) | 9.48 | 9.22 | excluded | excluded | no |
| `plain_area_scaling` | (1.0000, 4.0000, 16.0000) | 2.56 | 7.13 | compatible | excluded | **yes** |
| `q4_jordan_weight4` | (1.0000, 2.7500, 10.9908) | 2.08 | 7.00 | compatible | excluded | **yes** |
| `weight12_E12` | (1.0000, 32.5156, 2080.3072) | 3.61 | 11.25 | excluded | excluded | no |
| `weight12_E4_cubed` | (1.0000, 20.7969, 1327.6635) | 3.61 | 11.25 | excluded | excluded | no |
| `weight12_delta` | (1.0000, 0.1250, 2.7901e-05) | 21.17 | &infin; | excluded | excluded | no |
| `weight8_E8` | (1.0000, 7.5625, 120.7977) | 3.48 | 10.47 | excluded | excluded | no |

Range of the projective &sigma; over every positive-definite value of the one
covariance entry the frozen artifact did not store:

| competitor | &sigma; range | verdict stable |
|---|---|---|
| `bare_aspect_ratio` | [2.14, 4.77] | **no** |
| `no_modulus_dependence` | [9.21, &infin;] | yes |
| `plain_area_scaling` | [6.10, 8.93] | yes |
| `q4_jordan_weight4` | [5.93, 8.86] | yes |
| `weight12_E12` | [11.13, 11.38] | yes |
| `weight12_E4_cubed` | [11.13, 11.38] | yes |
| `weight12_delta` | [&infin;, &infin;] | yes |
| `weight8_E8` | [10.03, 10.97] | yes |

## T4 — the curvature functional

`f[1,2,4] = (m(4) &minus; 3 m(2) + 2 m(1)) / 6`, a linear functional of the
response: exactly 1 on `r`&sup2;, exactly 0 on any line, and with no denominator
and no matrix inverse anywhere in it.

Measured: **-4.6631e-04 &plusmn; 1.5297e-04**, *z* = **-3.05**, and *z* &isin; [-3.79, -2.62] across the admissible
covariance range.

| competitor | curvature predicted |
|---|---:|
| `bare_aspect_ratio` | 0 (exactly) |
| `no_modulus_dependence` | 0 (exactly) |
| `plain_area_scaling` | 1.000 |
| `q4_jordan_weight4` | 0.790 |
| `weight12_E12` | 330.793 |
| `weight12_E4_cubed` | 211.212 |
| `weight12_delta` | 0.271 |
| `weight8_E8` | 16.685 |

## T5 — the spin-8 amplitude each competitor needs

Leakage coefficient 0.054602 (= 1148/21025, exactly equal and opposite between the two families).

The frozen design assumed: *|A8/A4| well below 1*.

| competitor | required \|A&#8328;/A&#8324;\| to reach `r=2` |
|---|---:|
| `bare_aspect_ratio` | 7.7 |
| `no_modulus_dependence` | 0.4 |
| `plain_area_scaling` | 32.0 |
| `q4_jordan_weight4` | 32.1 |
| `weight12_E12` | 784.6 |
| `weight12_E4_cubed` | 782.8 |
| `weight12_delta` | 222.5 |
| `weight8_E8` | 182.2 |

Smallest requirement 0.39, largest 785.
