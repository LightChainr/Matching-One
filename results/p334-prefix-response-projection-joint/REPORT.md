# Prefix-response projection: original-block covariance join

All entries are estimates ± one original-20-batch jackknife SE. The source
is the original eight-quartet block, not the separately collected new64 block.

## N325

### Own-source signed covariance loading

| Receiver | Strength share (%) | Contact share (%) | Strength residual | Contact residual |
|---|---:|---:|---:|---:|
| first | 51.6816 ± 8.047 | 90.0062 ± 9.292 | 4.33095e-08 ± 1.209e-08 | 8.95777e-09 ± 9.2e-09 |
| second | 49.09 ± 9.247 | 80.2825 ± 8.076 | 4.389e-08 ± 1.115e-08 | 1.69986e-08 ± 7.792e-09 |

### Additional projected response variance after both baseline clocks

| Receiver | Physical source and response | Clock projection | Contact + clock projection | Increment |
|---|---|---:|---:|---:|
| first | source_first.C | 1.5374e-10 ± 6.055e-11 | 1.88084e-10 ± 5.407e-11 | 3.43435e-11 ± 2.213e-11 |
| first | source_first.W | 4.73683e-13 ± 2.731e-12 | 1.93058e-11 ± 2.249e-11 | 1.88321e-11 ± 2.148e-11 |
| first | source_second.C | 5.66226e-13 ± 1.218e-12 | 2.4657e-12 ± 2.371e-12 | 1.89948e-12 ± 2.437e-12 |
| first | source_second.W | 8.9526e-13 ± 2.584e-12 | 4.65304e-12 ± 5.077e-12 | 3.75778e-12 ± 4.725e-12 |
| second | source_first.C | 1.03031e-12 ± 1.386e-12 | 2.16277e-12 ± 2.63e-12 | 1.13246e-12 ± 1.828e-12 |
| second | source_first.W | 1.76816e-12 ± 3.847e-12 | 1.11055e-11 ± 1.078e-11 | 9.33734e-12 ± 7.978e-12 |
| second | source_second.C | 1.47302e-10 ± 4.598e-11 | 1.87441e-10 ± 5.102e-11 | 4.01394e-11 ± 2.676e-11 |
| second | source_second.W | 1.56197e-11 ± 2.355e-11 | 3.29528e-11 ± 2.762e-11 | 1.73331e-11 ± 1.801e-11 |

### Contact-response cross-moments after projecting out both baseline clocks

| Receiver / contact | source first C | source first W | source second C | source second W |
|---|---:|---:|---:|---:|
| first / joint_safe_mass | 9.54781e-08 ± 3.92e-08 | 5.0227e-08 ± 5.154e-08 | -1.06087e-08 ± 1.784e-08 | -2.94829e-08 ± 2.35e-08 |
| first / own_score_energy | 2.1562e-09 ± 8.804e-10 | -3.36349e-10 ± 9.546e-10 | -4.00206e-10 ± 3.728e-10 | -5.59698e-10 ± 4.716e-10 |
| first / own_safe_degree | 3.06498e-07 ± 1.221e-07 | 1.58911e-07 ± 1.71e-07 | -8.37377e-08 ± 6.361e-08 | -1.39515e-07 ± 8.312e-08 |
| first / own_safe_loop | 1.66149e-08 ± 7.278e-08 | 2.85771e-09 ± 9.846e-08 | -2.13891e-08 ± 3.36e-08 | -4.59266e-08 ± 4.882e-08 |
| second / joint_safe_mass | -1.9209e-08 ± 1.719e-08 | -5.62123e-08 ± 2.699e-08 | 8.2904e-08 ± 3.388e-08 | 3.01227e-08 ± 3.997e-08 |
| second / own_score_energy | 3.12535e-11 ± 2.722e-10 | -8.39432e-10 ± 4.429e-10 | 2.1393e-09 ± 1.013e-09 | 8.74235e-10 ± 1.223e-09 |
| second / own_safe_degree | -3.4478e-08 ± 5.773e-08 | -2.05067e-07 ± 8.596e-08 | 1.886e-07 ± 1.02e-07 | -3.02842e-08 ± 1.484e-07 |
| second / own_safe_loop | -1.53893e-08 ± 3.456e-08 | -4.89475e-08 ± 3.535e-08 | 9.75098e-08 ± 6.737e-08 | -7.61301e-08 ± 7.769e-08 |

### Paired contact-minus-strength gain

| Coordinate | Estimate | Shared-batch SE |
|---|---:|---:|
| first.contact_minus_strength.loading | 3.43517014e-08 | 8.66429e-09 |
| first.contact_minus_strength.loading_share | 0.383246586 | 0.0838521 |
| second.contact_minus_strength.loading | 2.68913039e-08 | 8.86868e-09 |
| second.contact_minus_strength.loading_share | 0.311924907 | 0.0900764 |

## N425

### Own-source signed covariance loading

| Receiver | Strength share (%) | Contact share (%) | Strength residual | Contact residual |
|---|---:|---:|---:|---:|
| first | 51.8618 ± 7.629 | 97.0183 ± 6.974 | 3.97827e-08 ± 8.13e-09 | 2.46416e-09 ± 5.915e-09 |
| second | 53.0963 ± 7.599 | 98.9991 ± 10.05 | 3.5752e-08 ± 8.923e-09 | 7.62933e-10 ± 7.743e-09 |

### Additional projected response variance after both baseline clocks

| Receiver | Physical source and response | Clock projection | Contact + clock projection | Increment |
|---|---|---:|---:|---:|
| first | source_first.C | 1.7206e-10 ± 2.794e-11 | 1.98842e-10 ± 2.853e-11 | 2.67822e-11 ± 1.915e-11 |
| first | source_first.W | 8.74349e-11 ± 5.034e-11 | 1.02086e-10 ± 5.671e-11 | 1.46513e-11 ± 1.425e-11 |
| first | source_second.C | 4.54495e-13 ± 1.18e-12 | 2.90864e-12 ± 3.364e-12 | 2.45414e-12 ± 2.691e-12 |
| first | source_second.W | 5.35638e-12 ± 8.35e-12 | 8.97882e-12 ± 1.034e-11 | 3.62244e-12 ± 4.693e-12 |
| second | source_first.C | 2.0258e-12 ± 2.465e-12 | 5.21231e-12 ± 3.563e-12 | 3.18651e-12 ± 2.552e-12 |
| second | source_first.W | 6.43568e-13 ± 2.063e-12 | 3.26661e-12 ± 3.421e-12 | 2.62304e-12 ± 3.217e-12 |
| second | source_second.C | 1.23783e-10 ± 3.948e-11 | 1.71847e-10 ± 4.215e-11 | 4.80649e-11 ± 1.744e-11 |
| second | source_second.W | 3.88536e-11 ± 3.857e-11 | 5.67917e-11 ± 4.718e-11 | 1.79381e-11 ± 1.928e-11 |

### Contact-response cross-moments after projecting out both baseline clocks

| Receiver / contact | source first C | source first W | source second C | source second W |
|---|---:|---:|---:|---:|
| first / joint_safe_mass | 9.01373e-08 ± 3.723e-08 | 2.18221e-08 ± 3.589e-08 | -1.13566e-08 ± 1.637e-08 | 2.86824e-08 ± 2.814e-08 |
| first / own_score_energy | 1.36972e-09 ± 8.062e-10 | 4.93751e-10 ± 7.724e-10 | -1.29262e-10 ± 3.778e-10 | 1.36228e-10 ± 3.2e-10 |
| first / own_safe_degree | 2.28468e-07 ± 1.117e-07 | 2.92361e-08 ± 1.133e-07 | -2.26856e-08 ± 7.334e-08 | 1.10754e-07 ± 7.584e-08 |
| first / own_safe_loop | 5.1532e-08 ± 6.21e-08 | -1.02318e-07 ± 8.145e-08 | 3.53894e-08 ± 2.817e-08 | 5.96239e-08 ± 3.817e-08 |
| second / joint_safe_mass | -8.20742e-09 ± 2.238e-08 | -4.72376e-09 ± 3.167e-08 | 1.07475e-07 ± 2.363e-08 | 6.01814e-08 ± 3.96e-08 |
| second / own_score_energy | 2.71627e-10 ± 4.653e-10 | -4.40343e-10 ± 6.227e-10 | 1.8078e-09 ± 5.64e-10 | 1.36564e-09 ± 8.005e-10 |
| second / own_safe_degree | -7.42307e-10 ± 6.737e-08 | -6.12582e-08 ± 9.412e-08 | 2.38216e-07 ± 9.885e-08 | 1.68164e-07 ± 1.246e-07 |
| second / own_safe_loop | 3.22731e-08 ± 3.703e-08 | -3.95499e-08 ± 3.057e-08 | -2.24849e-09 ± 3.196e-08 | -6.59795e-09 ± 7.213e-08 |

### Paired contact-minus-strength gain

| Coordinate | Estimate | Shared-batch SE |
|---|---:|---:|
| first.contact_minus_strength.loading | 3.73184914e-08 | 7.75902e-09 |
| first.contact_minus_strength.loading_share | 0.45156471 | 0.0931636 |
| second.contact_minus_strength.loading | 3.49891081e-08 | 4.95006e-09 |
| second.contact_minus_strength.loading_share | 0.459027475 | 0.0666044 |

## Interpretation boundary

Exploratory fixed six-predictor projections of receiver-rank0, within-rank-cell covariance, zero padded to each full population. Loading shares are signed response shares, not variance fractions, probabilities, causal attribution, or an out-of-sample closure result. Contact-after-clock projection increments are plug-in estimates with supplied refit LOO errors, not confidence lower bounds or R-squared. The clock-only loading identity is tautological; partial contact cross-moments retain their sampling uncertainty. All physical receivers, physical sources and C/W responses stay jointly correlated with ce20158a. The inherited covariance factor has rank at most19 per size; no high-dimensional inverse or omnibus test is used. The separate new64 block is neither read nor pooled.
