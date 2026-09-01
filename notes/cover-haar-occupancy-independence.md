# Haar cover occupancy independence

Issue #67 declares additive cover couplings

```text
H0: V = frac(sum_m U_m),
H1: 1-V.
```

For cover degree `Q>=2`, remove any one coordinate `U_i`.  The modulo-one sum
of the remaining independent uniforms is Haar uniform on the circle.  Adding
the fixed value of `U_i` only translates that uniform variable, so `V` is
independent of `U_i`.

It follows at every `0<p<1` that both H0 and H1 parent occupancies are
independent of each individual child occupancy `1[U_i<p]`.  By linearity,

```text
Cov(parent occupancy, Q^-1 sum_i 1[U_i<p]) = 0.
```

The checked script independently evaluates the relevant clipped-simplex
volumes with exact rational arithmetic.  It freezes `Q=2,5` at `p=2/5,3/5`;
in all four cases the joint probability is exactly `p^2` for H0 and H1.

## Interpretation

This is a useful negative control.  H0/H1 use every uniform in the fiber, but
their same-threshold single-site occupancy covariance is exactly zero.  Any
variance reduction for a nonlinear wrapping or homology observable must come
from higher-order spatial dependence, not from linear covariance with the
fiber occupancy mean.

## Reproduction

```text
python3 scripts/cover_haar_occupancy_independence.py
python3 -m unittest tests/test_cover_haar_occupancy_independence.py -v
```

## Boundary

The result does not determine nonlinear topological covariance, cover
threshold-rank permutations, wall time, residual-variance gain, or a
production coupling.  Issue #67 remains open.
