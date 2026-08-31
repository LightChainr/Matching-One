# Does one finite thermal coordinate explain both odd and even profiles?

A single sign-definite profile can always be matched to another equal-area profile by a cumulative-quantile map. The substantive question is whether the SAME map explains A and E, even after allowing each its own amplitude. The cumulative cross-moments remove that map exactly.

Candidate: U_j(p)=r_j phi'(p)D_j(phi(p)), j=A,E, with phi increasing and fixing0/1. Define J_m(D)=integral(F_A/A0)^m D_E. Then J_m(U)-r_E J_m(D)=0 for every m, without a small-warp or polynomial approximation.

Measured area amplitudes: r_A=-0.210955162, r_E=-0.401563476. The even target area is not used as a denominator.

| cumulative power m | target-source remainder | paired jackknife SE |
|---:|---:|---:|
| 1 | 2.45776279e-05 | 4.2108e-05 |
| 2 | 2.45777343e-05 | 4.22212e-05 |
| 3 | 2.15657551e-05 | 3.68483e-05 |
| 4 | 1.85523915e-05 | 3.14202e-05 |
| 5 | 1.60308789e-05 | 2.68954e-05 |
| 6 | 1.4001913e-05 | 2.32943e-05 |

Joint: chi-square=3.9008588/6, nominal p=0.690091.

The single oriented-area remainder Omega(U)-r_A*r_E Omega(D) is -5.45946384e-09 +/- 9.2802e-09; chi-square=0.3460867/1, nominal p=0.556337.

The six finite-polynomial integrals use 1404-node Gauss-Legendre (degree at most 2806), with all source and target uncertainty propagated by deleting each of 400 aligned batches. This scorer performs no additional sampling.

## Why the odd curve alone is not enough

One exact Bernstein subdivision does not certify the required signs for both empirical A curves at this scale. This is not a proof of a sign change. The cumulative-moment null remains a necessary common-transport condition without this positivity certificate; a quantile-map existence claim is not made for this target.

## Interpretation boundary

These are necessary conditions for one common signed-profile transport, not for every imaginable redefinition of two observables. A separate map per observer, ordinary observable relabelling without its Jacobian, or a genuine extra response component remain different mechanisms. The oriented area and six moments reuse one data block and must not be counted as separate confirmations.
