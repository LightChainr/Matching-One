# Rank-preserving next-label first/completion response

## N325

| Sector / mask | canonical Gamma +/- SE | integrated Gamma +/- SE |
|---|---:|---:|
| all / all | 0.001381329391 +/- 4.46227e-05 | 2.264784181e-05 +/- 1.5092e-06 |
| all / neither | 9.418646854e-05 +/- 2.68289e-05 | 2.139314015e-06 +/- 1.29037e-06 |
| all / mixed | 0.00108446236 +/- 3.18291e-05 | 1.811089285e-05 +/- 6.65283e-07 |
| all / both | 0.0002026805627 +/- 1.30382e-05 | 2.397634946e-06 +/- 2.26171e-07 |
| 01+10 / all | 0.0006475659804 +/- 3.35131e-05 | 1.138231491e-05 +/- 9.87569e-07 |
| 01+10 / neither | 5.235125255e-05 +/- 2.7083e-05 | 1.744832905e-06 +/- 8.91017e-07 |
| 01+10 / mixed | 0.0004445845101 +/- 2.32502e-05 | 8.012136086e-06 +/- 4.00563e-07 |
| 01+10 / both | 0.0001506302177 +/- 1.21476e-05 | 1.625345916e-06 +/- 1.84375e-07 |

## N425

| Sector / mask | canonical Gamma +/- SE | integrated Gamma +/- SE |
|---|---:|---:|
| all / all | 0.0009973691605 +/- 4.00065e-05 | 1.199353818e-05 +/- 8.23242e-07 |
| all / neither | 7.798500499e-05 +/- 2.6945e-05 | 1.573230908e-06 +/- 6.59656e-07 |
| all / mixed | 0.0007824119701 +/- 2.57717e-05 | 9.219230964e-06 +/- 3.7502e-07 |
| all / both | 0.0001369721854 +/- 5.03118e-06 | 1.201076306e-06 +/- 8.46813e-08 |
| 01+10 / all | 0.0004651817403 +/- 3.19245e-05 | 5.653950779e-06 +/- 5.76796e-07 |
| 01+10 / neither | 3.643684147e-05 +/- 2.26316e-05 | 5.474056076e-07 +/- 4.94174e-07 |
| 01+10 / mixed | 0.0003311502468 +/- 1.9566e-05 | 4.27682432e-06 +/- 2.48958e-07 |
| 01+10 / both | 9.759465207e-05 +/- 5.72811e-06 | 8.297208514e-07 +/- 8.50134e-08 |

For each prefix Z, neither/safe-safe contributes pi_safe(Z)^2 Cov_U(m(U)|Z,safe), averaged over prefixes; no division by a pooled safety rate.

Three masks partition original quartets by next-rank changes in either orientation; raw orientation blocks use this same paired mask. Signed Gamma allocations are not mean fractions or complete within-safe variance fractions. All original cells/prefix denominators and cross-provider covariance remain common; no PSD clipping, test suite, covariance inversion or additional production.
