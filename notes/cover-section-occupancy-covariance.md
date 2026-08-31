# Occupancy covariance of cover sections

Issue #67 proposes exact-marginal parent fields built from the `Q` uniforms in
each cover fiber.  Before measuring any topological event, the section and
antithetic-section choices already have an exact Bernoulli occupancy
covariance with the child fiber mean.  This certificate freezes that narrow
calculation.

Let

```text
X_m = 1[U_m < p],       Xbar = Q^-1 sum_m X_m.
```

## Direct section

For `Y=X_0`, only the shared zeroth child contributes to the covariance:

```text
Cov(Y,Xbar) = p(1-p)/Q,
Corr(Y,Xbar)^2 = 1/Q.
```

This is the squared-correlation ceiling within the declared pair consisting
of the direct and antithetic sections; it is not a theorem about every
measure-preserving function of an entire fiber.

## Antithetic section

For `Y_anti=1[1-U_0<p]`, the shared success intervals overlap by
`max(0,2p-1)`.  Exact subtraction of `p^2` gives

```text
Cov(Y_anti,Xbar) = -min(p,1-p)^2/Q,
Corr(Y_anti,Xbar)^2 = min(p,1-p)^4 / (Q p^2(1-p)^2).
```

The sign is negative, as desired for the negative-multiplier residual in the
parent issue, but its magnitude is strictly smaller than the direct section
away from `p=1/2`.  At `p=2/5` or `3/5`, the squared correlations are `2/9`
for `Q=2` and `4/45` for `Q=5`, versus `1/2` and `1/5` for the direct section.

## Reproduction

```text
python3 scripts/cover_section_occupancy_covariance.py
python3 -m unittest tests/test_cover_section_occupancy_covariance.py -v
```

## Boundary

These are single-site occupancy identities.  They do not cover H0/H1
additive couplings, determine covariance of a nonlinear topological
observable, include wall time, prove a residual-variance gain, or recommend a
production coupling.  Those gates in Issue #67 remain open.
