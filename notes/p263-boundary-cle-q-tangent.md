# Exact boundary CLE Q-tangent at percolation

## Result

Gefei Cai's general-`kappa` boundary four-point ODE gives an exact,
function-valued calibration of the `Q -> 1` tangent program.  Let

```text
sqrt(Q) = -2 cos(4 pi/kappa),
h = 8/kappa - 1,
```

on the branch with `kappa=6` at `Q=1`.  Direct differentiation gives

```text
d kappa/dQ |1 = -3 sqrt(3)/(2 pi),
d h/dQ     |1 =  sqrt(3)/(3 pi).
```

Write Cai's cross-ratio equation as `L_kappa[U_kappa]=0` and divide by
the same fixed factor 12 used in the paper's `kappa=6` specialization.
For

```text
T(lambda) = partial_Q U_kappa(lambda) |Q=1,
```

the exact tangent equation is

```text
L_6[T] = (3 sqrt(3)/(2 pi)) (partial_kappa L)_6[U_0].
```

Thus the logarithmic/tangent observable is not an arbitrary extra
homogeneous solution: its forcing term is fixed by the ordinary
percolation four-point function.

## Explicit operators

With polynomial coefficients written in `lambda`,

```text
L_6 =
  9 lambda^2(1-lambda)^2 D^3
  + 6 lambda(1-lambda)(1-2lambda) D^2
  + (8lambda(1-lambda)-6) D
  + 4(2lambda-1),

(partial_kappa L)_6 =
  (9/2) lambda^2(1-lambda)^2 D^3
  + 11 lambda(1-lambda)(1-2lambda) D^2
  + (-1-(2/3)lambda(1-lambda)) D
  - 2(2lambda-1).
```

The machine-readable artifact records the same coefficients as exact
rationals, in ascending powers of `lambda`.

## A universal logarithmic coefficient

The high Frobenius branch has

```text
V_{3h+1}(lambda,Q)
  = lambda^r sum_{n>=0} c_n(Q) lambda^n,
r = 3h+1 = 24/kappa-2,
c_0 = 1.
```

At percolation,

```text
r = 2,
partial_Q r |1 = sqrt(3)/pi.
```

Consequently,

```text
partial_Q V_{3h+1} |1
  = (sqrt(3)/pi) V_2(lambda) log(lambda)
    + lambda^2 sum_n (partial_Q c_n)|1 lambda^n.
```

The coefficient `sqrt(3)/pi` is exact and survives every analytic
`Q`-dependent amplitude rescaling.  It is therefore a sharper lattice
target than a fitted effective logarithmic coefficient.  The first
ordinary coefficients are

```text
V_2(lambda) = lambda^2(
    1 + lambda/3 + 37 lambda^2/198 + 112 lambda^3/891 + ...).
```

The exact dual-number recurrence in
`scripts/p263_boundary_tangent_ode.py` also returns the `Q` derivatives
of these coefficients.

## Full Green function versus cross-ratio function

Cai's full four-point function is

```text
G = K(x_1,x_2,x_3,x_4)^(2h) U(lambda),
```

where `K` is the conformal prefactor in the paper.  Therefore

```text
partial_Q G |1 = K^(2/3) [
    T(lambda)
    + (2 sqrt(3)/(3 pi)) log(K) U_0(lambda)
].
```

The second term is compulsory.  A lattice comparison that differentiates
only the cross-ratio function silently drops a known piece of the tangent.

## Normalization gauge and the score to use

The differential equation alone leaves `T` ambiguous up to a homogeneous
solution whenever the `Q` dependence of the continuum normalization is
not fixed.  This is an amplitude gauge, not a failure of the tangent
program.  Eliminate it with an anchored shape derivative,

```text
partial_Q log[U(lambda)/U(lambda_anchor)] |1,
```

or with a universal connectivity ratio.  The high-branch
`(sqrt(3)/pi) log(lambda)` coefficient is already gauge invariant.

For a lattice `Q`-score experiment, the comparison vector should contain

```text
measure-score covariance
+ explicit projector/field derivative
+ conformal-prefactor derivative,
```

and should be projected off the single amplitude-gauge direction before
the cross-ratio chi-square is formed.

## Research consequence

This turns the boundary benchmark into a two-stage exact test:

1. verify the inhomogeneous ODE residual over several cross ratios;
2. verify the universal high-branch logarithmic coefficient
   `sqrt(3)/pi` after amplitude projection.

Failure of the first test diagnoses `Q`-score/projector semantics.  Passing
the ODE but missing the logarithmic coefficient diagnoses a wrong physical
branch or boundary connectivity normalization.  No bulk Potts OPE fit is
needed for either decision.

## Source

- Gefei Cai, *Boundary four-point connectivities of conformal loop
  ensembles*, arXiv:2603.28161v2:
  <https://arxiv.org/abs/2603.28161>

The source paper supplies the general-`kappa` ODE and the `kappa=6`
solutions.  The `Q` differentiation, explicit forcing operator,
amplitude-gauge separation, and `sqrt(3)/pi` tangent coefficient above
are derived here.

