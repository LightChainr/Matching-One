# Exact profile-distance certificate

This utility compares normalized rational polynomial densities without floating point.

For densities `rho_1` and `rho_2`, it computes

```text
density L2 squared = integral (rho_1-rho_2)^2 dp,
CDF CvM squared    = integral (F_1-F_2)^2 dp
```

exactly as `Fraction` values. Both inputs must integrate exactly to one. Symmetry and identity-zero are checked directly.

For the frozen N=4 synthetic density relative to the uniform profile,

```text
density L2 squared = 1/20,
CDF CvM squared    = 1/840.
```

The control Gram calculation reveals an additional exact identity:

```text
rho_frozen = (rho_uniform + rho_Beta(2,2))/2.
```

Consequently the two centered difference vectors are exactly collinear and their 2x2 Gram determinant is zero, not merely numerically small. This explains the matrix `[[1/20,1/10],[1/10,1/5]]` without assigning continuum significance to the synthetic fixture.

## Boundary

These are exact synthetic polynomial distances. No production histogram, bootstrap uncertainty, cross-model universality test, or tail conclusion is included. Issue #28 remains open.

