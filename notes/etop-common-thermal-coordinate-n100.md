# A can change coordinates; A and E cannot share that explanation

## Scientific result

The N100 thermal redistribution is not exhausted by an unknown common
area-preserving thermal coordinate. A one-field interpretation is flexible:
**the actual empirical A profiles admit a cumulative-quantile map exactly**.
The informative obstruction appears only when E must use that same coordinate,
even after giving it an independent amplitude. Six coordinate-free cumulative
moments give **chi-square 53.91436/6, nominal p=7.68e-10**.

This advances the earlier finite fixed-p shape splitting into a whole-profile
statement. It is a statement about a specified signed-profile transport class,
not every possible change of scalar observables or a continuum field count.

All calculations reuse `7b30648` / PR484's three N100 shape pairs and 200
aligned batches. They are post-reveal analysis, not a new independent experiment.
There is no new Monte Carlo, engine work, installation or old suite rerun.

## 1. Why failure of low-order temperature shifts is not the final argument

Let `D=Y(4i)-Y(2i)`, `U=Y(1/2+i)-Y(2i)`, and normalize the odd area by
`r_A=U_A,0/D_A,0=-.277981748`. Set `R_A=U_A-r_A D_A` and
`z=s(p-p_ref)`, `s=N^(3/8)`. The tangent model

```text
R_A = (1/s) d_p[v(z) D_A],   v(z)=sum_l theta_l z^l
```

has the exact integrated-moment relation

```text
R_A,j = -j sum_l theta_l D_A,j-1+l.
```

The lowest `degree+1` moments determine theta; higher moments then supply
overidentifying restrictions. Delete-one paired-batch covariance propagates
the fitted clock, source moments and transport coefficients together.

| velocity | theta, in increasing powers of z | remaining odd chi-square / df |
|---|---|---:|
| constant: translation tangent | .07188292 | 481.502 / 5 |
| linear: translation plus normalized dilation | .05033264, -.05105995 | 217.845 / 4 |
| quadratic: nonlinear tangent | .12690859, -.07251449, -.05291981 | 148.243 / 3 |

These restrictions fail strongly, but **a finite map need not have a small
or low-degree tangent velocity**. We therefore remove the entire map rather
than escalating polynomial order. These are profile tangents; ordinary
scalar-observable reparameterization has no Jacobian. Dilation plus an
amplitude correction that fixes the integral produces the model above.

## 2. The relaxed finite candidate

Allow independent signed amplitudes and one shared increasing map:

```text
U_j(p) = r_j phi'(p) D_j(phi(p)),   j=A,E,
phi(0)=0, phi(1)=1.
```

Its area constraints fix `r_j=U_j,0/D_j,0`. The data give

```text
r_A=-.277981748,   r_E=-.142339516.
```

The E amplitude is deliberately not forced to equal the A clock. We divide
only by the well-resolved source E area, not by the weaker target E area.
Thus the result is not merely the earlier nonzero lifetime residual in disguise.
The class allows arbitrary finite interior-regular maps and degenerate
endpoint derivatives. It is broader than a regular endpoint diffeomorphism.

## 3. The empirical A curves themselves have no such obstruction

For a threshold-rank histogram, `F_k(p)=sum_(j>=k) B_(j,N)(p)`.
Therefore its integrated integer threshold counts immediately give the
Bernstein coefficients of the empirical orientation contrast. The raw signed
counts have a common positive scale `625/(1152*2000000)`.

One exact de Casteljau subdivision at `p=1/2` makes all coefficients of
**D_A and -U_A** nonnegative on both half intervals, with a positive midpoint
and nonzero coefficients on each half. Hence both are strictly positive on
the open unit interval. The complete integer certificate is saved in the
finite-invariants JSON. This is exact for the empirical finite polynomials,
not a theorem about their population expectations or unobserved rare tails.

Their normalized cumulative functions

```text
q_D(p)=integral_0^p D_A / D_A,0,
q_U(p)=integral_0^p U_A / U_A,0
```

are consequently strictly increasing. The unique candidate map is
`phi=q_D^(-1) o q_U`. It matches A identically. A alone cannot identify a
mechanism by excluding low-degree warps; an unrestricted map absorbs the
entire one-field response. Endpoint smoothness is a separate condition.

## 4. E in the cumulative A clock is an additional physical readout

Define moments of the E signed measure in the A cumulative coordinate:

```text
J_m(D) = integral_0^1 q_D(p)^m D_E(p) dp.
```

Change variables through the candidate phi. For every m,

```text
J_m(U) = r_E J_m(D).
```

There is no fitted map or velocity order in this relation. When A has fixed
sign, equality of all moments determines the same finite signed measure on
q in [0,1]. The six moments here are necessary conditions, not a finite
completeness claim. Their observed target-minus-source residuals are:

| m | residual | paired jackknife SE |
|---:|---:|---:|
| 1 | -.0001406310 | .0000516821 |
| 2 | -.0001409438 | .0000517667 |
| 3 | -.0001232235 | .0000453739 |
| 4 | -.0001053324 | .0000389431 |
| 5 | -.0000903842 | .0000335620 |
| 6 | -.0000784458 | .0000292437 |

Their full covariance gives `53.91436/6`, rather than treating six similar
marginal z-scores as six confirmations. Shape contrasts remove shared
amplitude uncertainty, which is why the joint discrepancy is stronger than
the marginal values. All uncertainty comes from the same200 aligned batches.

Each integrand has polynomial degree at most `6*(100+1)+100=706`.
The354-node Gauss-Legendre rule therefore integrates these empirical
polynomials exactly in exact arithmetic; the implementation uses double
precision. This is not a numerical sampling approximation to a new process.

The first moment is equivalently an oriented cumulative area. Set

```text
Omega(D) = integral(F_A D_E - F_E D_A).
```

Then the finite candidate requires `Omega(U)=r_A r_E Omega(D)`.
The observed remainder is `3.193815e-7 +/- 1.23524e-7`, with single-coordinate
nominal p=.00972. It is a compact view of the same data, not additional
independent evidence. The whole cumulative shape supplies more discrimination.

## 5. Mechanism update and next observation

**A useful working picture is now two relative response profiles, not merely
one profile with a mysterious temperature scale.** The odd cumulative curve
can itself define a clock; the even measure still changes in that clock.
This supplies a coordinate-free target for subsequent scales or geometries.

What is not established: an exact continuum field, a particular number of
operators, the critical-window origin of the change, failure of arbitrary
independent maps for the observers, or failure of ordinary no-Jacobian
scalar-observable relabelling. Those are different model classes.

The next informative datum is the same `J_m` vector at a homothetic larger
scale, with full source-target covariance and the same observer definitions.
This separates persistence of the relative E shape from a finite thin-geometry
effect. The current measurement stays full-p; it does not automatically become
a critical CFT fingerprint because a convenient coordinate used `N^(3/8)`.

Lifecycle: post-reveal finite-model mechanism analysis / completed /
one reused N100 dependency block / branch-only. Source is the200-batch
seed20260831125401, offset267100000000 dataset in PR484. Theory companion:
`c201c78:notes/p267-thermal-warp-invariants.md`.

Reproduction: `scripts/etop_thermal_transport.py` and
`scripts/etop_finite_transport_invariants.py`; outputs are under matching
`results/etop-thermal-transport/` and `results/etop-finite-transport-invariants/`.
