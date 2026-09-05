# N100: a new three-modulus response experiment

Each pair has two million new shared-counter permutations in 200 aligned batches. All three shapes were frozen before acquisition; the full 12x12 covariance is retained.

| shape | A_top | E_top | C | W |
|---|---:|---:|---:|---:|
| tau_2i | 0.00099357313 +/- 0.000199 | -0.0002932555 +/- 0.000167 | -0.00062486168 +/- 2.37e-05 | 0.00019874229 +/- 4.64e-05 |
| tau_4i | 0.0016171757 +/- 9.48e-05 | -0.00044115367 +/- 0.000101 | -0.0026673169 +/- 2.62e-05 | 0.0013595302 +/- 5.48e-05 |
| tau_half_plus_i | -7.9385928e-05 +/- 0.00025 | -0.000156851 +/- 0.00013 | -5.7096408e-05 +/- 2.42e-05 | 3.3516291e-05 +/- 2.73e-05 |

## Same-area shape comparison: no N50 calibration or area exponent

| model | chi-square / df | p |
|---|---:|---:|
| affine_E4 | 58.26869 / 4 | 6.70156e-12 |
| affine_height_E4 | 24.4233 / 4 | 6.56908e-05 |
| affine_height_squared | 23.89369 / 4 | 8.38909e-05 |
| affine_log_height | 2330.066 / 4 | 0 |
| free_common_secant | 19.21808 / 3 | 0.00024643 |
| no_shape_response | 12082.66 / 8 | 0 |

The first row is the frozen primary. Other fixed rows are declared comparators; free-common-secant is an exploratory one-parameter relaxation. None is a continuum field count or E4 identity proof.

## Source-informed N50-to-N100 shape transfer

Independent gains are profiled for A/E/C/W, with both uncertain source vectors and all same-stream target covariance. This tests the additional cross-area shape-separability hypothesis without dividing by the weak E denominator.

| model | chi-square / df | p |
|---|---:|---:|
| affine_E4 | 64.79027 / 8 | 5.30959e-11 |
| affine_height_E4 | 91.30763 / 8 | 2.523e-16 |
| affine_height_squared | 117.647 / 8 | 1.01418e-21 |
| affine_log_height | 2603.274 / 8 | 0 |

Gaussian covariance/profile comparisons; same-area primary is a fixed linear null. Models/readouts reuse one new common-random block; scores are not independent evidence.

Prediction freeze: 4c1ec50. Exact geometry: b9e4ea1; three-modulus null: 964d770. No old source p-values are added to the new scores.
