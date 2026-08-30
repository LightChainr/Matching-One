# Exact threshold-profile symmetry certificate

This slice certifies reflection symmetry and the center of the frozen N=4 synthetic threshold profile using rational arithmetic only.

The rank-mixture weights are

```text
[1/8, 3/8, 3/8, 1/8].
```

They satisfy `w_k=w_{N+1-k}`. Direct power-basis substitution then verifies

```text
rho(1-p)=rho(p),
F(1-p)=1-F(p),
F(1/2)=1/2.
```

For the frozen density

```text
rho(p)=1/2+3p-3p^2,
rho'(p)=3-6p.
```

The derivative is strictly positive below `1/2` and strictly negative above `1/2`. The endpoint densities are both positive. Therefore the CDF is strictly increasing, the median is uniquely `1/2`, and the density has the unique mode `1/2`.

The mode gate intentionally accepts only a nonconstant quadratic density with an exactly linear derivative. Constant/plateau, wrong-sign, off-center, and higher-degree cases fail closed; higher-degree production profiles will require exact root isolation rather than extrapolation from this fixture.

## Boundary

This is an exact synthetic low-degree certificate. It contains no production histogram, empirical median or mode, bootstrap/profile distance, tail-window fit, or universality claim. Issue #28 remains open.

