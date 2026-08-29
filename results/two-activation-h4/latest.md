# K1/K2 two-activation H4 decomposition

This is a retrospective reanalysis of archived threshold-rank histograms. It generates no Monte Carlo samples and fits no exponent.

Exact convention: `K1=K_minus`, `K2=K_plus`, `M(p)=-1+F1(p)+F2(p)`.

| N | status | K1 share | interaction | angular delta p1 (z) | angular delta p2 (z) | angular root gap | closure residual |
|---:|:---|---:|:---|---:|---:|---:|---:|
| 65 | scoreable | 0.669 | reinforcing | -6.658850e-05 (-21.51) | -3.289970e-05 (-10.65) | -9.948831e-05 | -1.041e-10 |
| 85 | scoreable | 0.715 | reinforcing | -3.907400e-05 (-16.17) | -1.559013e-05 (-6.20) | -5.466416e-05 | -3.428e-11 |
| 130 | scoreable | 0.730 | reinforcing | -1.952877e-05 (-6.46) | -7.205511e-06 (-2.46) | -2.673428e-05 | -3.968e-12 |
| 145 | scoreable | 0.713 | reinforcing | -1.734329e-05 (-8.82) | -6.983942e-06 (-3.48) | -2.432724e-05 | -6.093e-12 |
| 170 | scoreable | 0.702 | reinforcing | -1.544424e-05 (-6.94) | -6.552739e-06 (-3.12) | -2.199698e-05 | -2.735e-12 |
| 185 | scoreable | 0.944 | reinforcing | -8.869319e-06 (-6.72) | -5.219741e-07 (-0.41) | -9.391294e-06 | -3.354e-13 |
| 265 | scoreable | 0.952 | cancelling | -5.492486e-06 (-7.06) | 2.747598e-07 (0.33) | -5.217726e-06 | -1.726e-13 |
| 290 | scoreable | 0.644 | reinforcing | -6.396891e-06 (-4.19) | -3.530834e-06 (-2.34) | -9.927725e-06 | -4.947e-13 |
| 325 | scoreable | 0.872 | cancelling | -2.595696e-06 (-1.67) | 3.805462e-07 (0.24) | -2.215150e-06 | -6.100e-15 |
| 425 | scoreable | 0.905 | cancelling | -2.501064e-06 (-1.99) | 2.626158e-07 (0.22) | -2.238448e-06 | -7.164e-15 |

## Descriptive result

Across the 10 scoreable sizes, magnitude classification is K1-dominant at 9, K2-dominant at 0, and shared at 1; 3 sizes show opposite-sign component cancellation.
Using the within-size delete-one standard errors only, K1 has |z|>=2 at 8 sizes and K2 at 6. The opposite-sign K2 point estimates at N=265,325,425 are individually unresolved, so the apparent large-size sign change is a next-target clue rather than a confirmed transition.
This count is a map of the decomposition, not an independent-evidence vote: N=65,85,130,170 share one counter stream and remain one dependency group.

The `nonlinear closure residual` is the observed orientation root gap minus `delta_p1+delta_p2`.  Its smallness diagnoses the local linearization only; it does not establish a continuum exponent or operator identity.

## Joint rank coordinates

Every scoreable archive also reports means, variances and covariance for `C=(K1+K2)/2` and `G=K2-K1` separately for both orientations and their pooled mixture.  These quantities use the archived paired integer moments, not a product reconstructed from marginal histograms.

## Covariance and provenance

All nonlinear coordinates are recomputed after deleting the same batch from both orientations.  The JSON contains each size's full covariance over the declared estimate vector and a cross-size covariance over the decision metrics. Entries from distinct counter intervals are zero by design; shared streams use aligned delete-one covariance.  Input SHA256 values and source commits are stored under each size.
