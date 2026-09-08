# Does one finite thermal coordinate explain both odd and even profiles?

A single sign-definite profile can always be matched to another equal-area profile by a cumulative-quantile map. The substantive question is whether the SAME map explains A and E, even after allowing each its own amplitude. The cumulative cross-moments remove that map exactly.

Candidate: U_j(p)=r_j phi'(p)D_j(phi(p)), j=A,E, with phi increasing and fixing0/1. Define J_m(D)=integral(F_A/A0)^m D_E. Then J_m(U)-r_E J_m(D)=0 for every m, without a small-warp or polynomial approximation.

Measured area amplitudes: r_A=-0.277981748, r_E=-0.142339516. The even target area is not used as a denominator.

| cumulative power m | target-source remainder | paired jackknife SE |
|---:|---:|---:|
| 1 | -0.000140630969 | 5.16821e-05 |
| 2 | -0.000140943841 | 5.17667e-05 |
| 3 | -0.000123223509 | 4.53739e-05 |
| 4 | -0.000105332442 | 3.89431e-05 |
| 5 | -9.03841508e-05 | 3.3562e-05 |
| 6 | -7.84458213e-05 | 2.92437e-05 |

Joint: chi-square=53.914357/6, nominal p=7.67641e-10.

The single oriented-area remainder Omega(U)-r_A*r_E Omega(D) is 3.19381521e-07 +/- 1.23524e-07; chi-square=6.685216/1, nominal p=0.00972158.

The six finite-polynomial integrals use 354-node Gauss-Legendre (degree at most706), with all source and target uncertainty propagated by deleting each of200 aligned batches. No new sampling is used.

## Why the odd curve alone is not enough

For the actual empirical histograms, exact integer Bernstein subdivision proves D_A>0 and U_A<0 for every0<p<1. After dividing by their signed areas, both are positive normalized profiles. Thus the odd data alone admit a unique cumulative-quantile map, despite rejecting low-degree tangent models. The finite common-map failure comes from the E channel, not from trying an insufficient velocity polynomial. This existence statement allows degenerate endpoint derivatives and does not prove positivity of the underlying population curves.

## Interpretation boundary

These are necessary conditions for one common signed-profile transport, not for every imaginable redefinition of two observables. A separate map per observer, ordinary observable relabelling without its Jacobian, or a genuine extra response component remain different mechanisms. The oriented area and six moments reuse one data block and must not be counted as separate confirmations.
