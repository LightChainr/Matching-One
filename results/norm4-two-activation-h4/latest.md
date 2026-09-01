# K1/K2 two-activation H4 decomposition

This is a retrospective reanalysis of archived threshold-rank histograms. It generates no Monte Carlo samples and fits no exponent.

Exact convention: `K1=K_minus`, `K2=K_plus`, `M(p)=-1+F1(p)+F2(p)`.

| N | status | K1 share | interaction | angular delta p1 (z) | angular delta p2 (z) | angular root gap | closure residual |
|---:|:---|---:|:---|---:|---:|---:|---:|
| 65 | scoreable | 0.668 | reinforcing | -6.865100e-05 (-91.79) | -3.414953e-05 (-45.51) | -1.028006e-04 | -1.080e-10 |
| 85 | scoreable | 0.698 | reinforcing | -4.176365e-05 (-75.55) | -1.809705e-05 (-30.77) | -5.986075e-05 | -3.947e-11 |
| 130 | scoreable | 0.742 | reinforcing | -1.924726e-05 (-37.85) | -6.702057e-06 (-10.34) | -2.594933e-05 | -4.049e-12 |
| 170 | scoreable | 0.795 | reinforcing | -1.229058e-05 (-21.86) | -3.169899e-06 (-6.33) | -1.546048e-05 | -1.621e-12 |
| 260 | scoreable | 0.979 | reinforcing | -5.469189e-06 (-7.75) | -1.187904e-07 (-0.20) | -5.587979e-06 | -1.139e-13 |
| 340 | scoreable | 0.852 | reinforcing | -4.607924e-06 (-8.17) | -7.989227e-07 (-1.70) | -5.406847e-06 | -1.046e-13 |

## Descriptive result

Across the 6 scoreable sizes, magnitude classification is K1-dominant at 6, K2-dominant at 0, and shared at 0; 0 sizes show opposite-sign component cancellation.
Using the within-size delete-one standard errors only, K1 has |z|>=2 at 6 sizes and K2 at 4. Unresolved component point estimates remain next-target clues rather than confirmed transitions.
These counts map the decomposition; they are not independent-evidence votes. The dependency groups and their member sizes are recorded explicitly in JSON.

The `nonlinear closure residual` is the observed orientation root gap minus `delta_p1+delta_p2`.  Its smallness diagnoses the local linearization only; it does not establish a continuum exponent or operator identity.

## Joint rank coordinates

Every scoreable archive also reports means, variances and covariance for `C=(K1+K2)/2` and `G=K2-K1` separately for both orientations and their pooled mixture.  These quantities use the archived paired integer moments, not a product reconstructed from marginal histograms.

## Covariance and provenance

All nonlinear coordinates are recomputed after deleting the same batch from both orientations.  The JSON contains each size's full covariance over the declared estimate vector and a cross-size covariance over the decision metrics. Entries from distinct counter intervals are zero by design; shared streams use aligned delete-one covariance.  Input SHA256 values and source commits are stored under each size.

This is a canonical Phase-D state-coordinate decomposition. It does not construct the Phase-E `J_top` versus `J_bulk` comparison and is not an outward-rounded interval or SOS certificate.
