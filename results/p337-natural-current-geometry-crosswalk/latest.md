# Geometry-aware crosswalk of the natural A current

The model direction is selected from N65/N85 only. N145 is a held-out two-component diagnostic.

All Gaussian quotients have exact `tau=i`, so `E4(i)` is common. The varying exact descriptor is `Re z_axis=cos(4 theta)`. The charged source obeys `q_A^2=(u+H_F3)/2`, giving one additional projective scalar `u/2=1/2`.

| model | training chi2 / df | N145 predictive chi2 / 2 | predicted N145 pair | pair residual / SE |
|---|---:|---:|---:|---:|
| pure_N_law | 23.1338 / 3 | 42.5927 | 0.00000000 | 6.476 |
| one_H4_geometry_covector | 2.5752 / 3 | 1.5191 | 0.01524124 | 1.119 |
| H4_geometry_plus_A_projective_scalar | 2.4544 / 2 | 1.0293 | 0.01560046 | 1.014 |

The exact H4 model predicts N145 components `-0.00751213, +0.00772911` against observed `-0.01072378, +0.00958769`. Its full predictive score is `1.519/2`; the pair residual is only `1.119` predictive SE.

Adding the charged scalar improves training chi-square by only `0.121` and gives held-out `1.029/2`. It is not selected by these data.

In central-value accounting, exact angle rotation contributes `22.8%` of the apparent rebound relative to the N85-anchored radial target; the remaining `77.2%` is a scale-curvature remainder. But that remainder is not resolved: the joint H4 geometry model already gives an acceptable held-out score. Geometry is necessary; extra curvature is suggested centrally but not required statistically.

The next clean geometry is the N170 angle-flip child `(11+7i,13+i)` of N85 under common multiplier `1+i`: H4 flips exactly, the charged scalar stays fixed, and the area is close to N145. No matching projective-birth archive exists, so it is selected but not run.
