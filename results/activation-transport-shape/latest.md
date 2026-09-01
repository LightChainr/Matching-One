# K1/K2 translation versus zero-area deformation

This reuses the ten frozen threshold-rank archives. It generates no samples and fits no exponent.

For each activation, `D=A fbar+R`: `A fbar` is the pooled-density translation tangent with the same exact area as `D`, while `R` has exactly zero integral.

| N | K1 deformation energy share | K2 deformation energy share | K2 D(p_bar) | K2 translation | K2 deformation | reading |
|---:|---:|---:|---:|---:|---:|:---|
| 65 | 0.531 | 0.529 | 2.756718e-04 | 1.169467e-03 | -8.937957e-04 | translation and deformation partially cancel |
| 85 | 0.546 | 0.612 | 1.443149e-04 | 8.648846e-04 | -7.205697e-04 | translation and deformation partially cancel |
| 130 | 0.575 | 0.658 | 7.812570e-05 | 5.347473e-04 | -4.566216e-04 | translation and deformation partially cancel |
| 145 | 0.573 | 0.644 | 7.886284e-05 | 4.768925e-04 | -3.980296e-04 | translation and deformation partially cancel |
| 170 | 0.571 | 0.610 | 7.851370e-05 | 3.980805e-04 | -3.195668e-04 | translation and deformation partially cancel |
| 185 | 0.639 | 0.789 | 6.455317e-06 | 3.090669e-04 | -3.026116e-04 | translation and deformation partially cancel |
| 265 | 0.637 | 0.833 | -3.885980e-06 | 1.991529e-04 | -2.030388e-04 | shape residual overturns positive translation |
| 290 | 0.579 | 0.596 | 5.164857e-05 | 2.233044e-04 | -1.716559e-04 | translation and deformation partially cancel |
| 325 | 0.752 | 0.894 | -5.808382e-06 | 1.262642e-04 | -1.320726e-04 | shape residual overturns positive translation |
| 425 | 0.659 | 0.870 | -4.431895e-06 | 1.032796e-04 | -1.077115e-04 | shape residual overturns positive translation |

## Decision reading

The same-area translation term is positive at every scored K2 root point. The deformation is negative there and accounts for the local negative values at N=265,325,425. Thus those lobes are activation-distribution reshaping, not a reversal of the integrated translation direction.

The JSON retains every Bernstein coefficient, the exact area/closure audits, within-N covariance, and the full cross-N covariance under aligned dependency-group deletion. Derived translation and deformation coordinates from one archive remain one evidence block.

This finite-archive split does not decide whether the reshaping is carried by rank-one lifetime, plateau line, landing pivotal, or another typed field; it identifies the next mechanism-sensitive information axis.
