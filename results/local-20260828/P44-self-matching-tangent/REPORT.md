# P44 self-matching tangent exact result

The explicit family

```text
p_even = 1/2 + t + lambda
p_odd  = 1/2 + t - lambda
```

has legal domain `|t+lambda|<=1/2`, `|t-lambda|<=1/2` and exact matching
exchange `(t,lambda)->(-t,-lambda)`.

All `2^10=1024` configurations of the `(3,1)` quotient were enumerated. For
each of the five wrapping channels the exact center response matrix, with rows
`(Rplus,Rminus)` and columns `(t,lambda)`, is

```text
[[0,    0],
 [15/8, 5/4]].
```

At `t=0` every channel has

```text
Rminus(lambda) = (5/4) lambda - 4 lambda^5.
```

The only root in the probability interval `[-1/2,1/2]` is `lambda*=0`; the
other real-root magnitude is `(5/16)^(1/4)=0.747674...`. Thus this family has
no nonzero exchange-odd improved point on the minimum exact quotient.

The nontrivial improved-action target is instead the matching-even H4
amplitude. Its frozen first design is N=130, `(11,3)` versus `(9,7)`, at
nonnegative `lambda = 0, 1/8, 1/4, 3/8`; fit the first three in `z=lambda^2`
and retain `3/8` as a no-refit lack-of-fit point.

Evidence boundary: exact finite-quotient lattice algebra, not a continuum-field
identification and not a stochastic N=130 amplitude result.
