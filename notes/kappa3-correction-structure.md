# Correction structure of the universal derivative ratio kappa_3

Status: new analytic hypothesis motivated by the Mertens-Ziff finite-size data.

## Known starting point

For the matching function in the square-torus scaling limit,

\[
M_L(p)\to \mathcal M(z)=f(z)-f(-z),
\qquad z\propto(p-p_c)L^{y_t},
\qquad y_t=1/\nu=3/4.
\]

Mertens and Ziff explicitly note that

\[
\kappa_3(L)=\frac{M_L'''(p_c)}{M_L'(p_c)^3}
\]

is independent of the linear metric factor and extrapolates numerically to about `-1.67`; their quoted finite-size approach was approximately `L^-1.38` on the available sizes (Phys. Rev. E 94, 062152 (2016), scaling discussion around Eqs. 36-38).

The candidate

\[
\kappa_3^*=-5/3
\]

is therefore a legitimate exact-value hypothesis, but not yet evidence of rationality.

## New observation: analytic nonlinear scaling fields predict L^-3/2

The microscopic probability is not required to be the exact thermal scaling coordinate. Write the true analytic thermal field as

\[
u(\delta)=b_1\delta+b_2\delta^2+b_3\delta^3+\cdots,
\qquad \delta=p-p_c.
\]

The leading finite-size scaling form is

\[
M_L(p)=\mathcal M\!\left(u(\delta)L^{y_t}\right)+\text{irrelevant corrections}.
\]

Because `mathcal M` is odd,

\[
\mathcal M''(0)=0.
\]

At `delta=0`,

\[
M_L'(p_c)=\mathcal M'(0)b_1L^{y_t},
\]

while

\[
M_L'''(p_c)
=\mathcal M'''(0)b_1^3L^{3y_t}
+6\mathcal M'(0)b_3L^{y_t}
+\cdots.
\]

Hence

\[
\frac{M_L'''(p_c)}{M_L'(p_c)^3}
=\frac{\mathcal M'''(0)}{\mathcal M'(0)^3}
+C_{\rm nl}L^{-2y_t}+\cdots.
\]

Since `y_t=3/4`, the analytic nonlinear-coordinate correction is

\[
\boxed{L^{-3/2}}.
\]

The previously reported apparent exponent `1.38` is close enough that a fixed `3/2` correction model should be tested before introducing a free correction exponent.

## A second route also produces 3/2

The conventional two-dimensional percolation correction-to-scaling literature gives

\[
\Omega=72/91,
\qquad D=91/48,
\qquad \omega=D\Omega=3/2.
\]

See R. M. Ziff, Phys. Rev. E 83, 020107(R) (2011), and the later general Potts derivation in Phys. Rev. E 111, 034108 (2025).

This numerical equality does **not** prove the two mechanisms are the same. It creates an identification problem:

- an `L^-3/2` term in `kappa_3(L)` can arise from nonlinear choice of thermal coordinate;
- a genuine irrelevant/correction field may also generate an `L^-3/2` contribution for some observables.

The project should distinguish these rather than simply fitting `omega=1.5`.

## How to distinguish coordinate contamination from genuine irrelevant corrections

Use a second monotone dimensionless observable `U_L(p)` with a universal square-torus scaling function,

\[
U_L(p)=\mathcal U(z)+\cdots.
\]

Examples include a primal wrapping probability with a fixed topological definition.

Instead of treating `p` as the local coordinate, study the **parametric universal curve**

\[
M_L=\Phi(U_L).
\]

Both observables depend on the same nonlinear thermal field `u(p)`, so eliminating `p` eliminates the analytic reparameterization of the thermal coordinate at leading scaling order.

Prediction:

- `kappa_3` formed from p-derivatives may have a leading `L^-3/2` approach;
- shape coefficients of the parametric curve `M(U)` should lose the nonlinear-coordinate contribution and may converge with the leading genuine irrelevant exponent instead.

If the `L^-3/2` term survives unchanged in the coordinate-free curve, it is evidence for a genuine correction field. If it is strongly suppressed, nonlinear thermal-coordinate contamination was important.

## Exact-threshold controls

Before using disputed square-site `p_c`, measure the same derivative ratio on exact-threshold models.

### Square bond percolation

- square torus;
- self-dual;
- exact `p_c=1/2`;
- same continuum shape as the square-site target.

This is the cleanest universality control. The limiting `kappa_3` should agree with square-site percolation if the ratio is truly lattice-independent for fixed shape.

### Self-matching site models

For a self-matching lattice such as triangular or union-jack site percolation,

\[
p_c=1/2
\]

and the finite matching function satisfies

\[
M_L(1/2)=0
\]

for every `L`. The union-jack model is especially useful because it can be simulated with a square macroscopic boundary, matching the target shape.

## Cheap derivative estimator at p=1/2

For any configuration observable `D(C)` whose expectation is `M_L(p)`, let `K` be the number of occupied Bernoulli variables among `N` and define

\[
x=2K-N.
\]

At `p=1/2`, differentiating the Bernoulli weight gives

\[
M_L'(1/2)=\mathbb E[2xD],
\]

and

\[
M_L'''(1/2)
=\mathbb E\left[8\left(x^3-(3N-2)x\right)D\right].
\]

Thus the exact-threshold controls can estimate `kappa_3` directly from one fixed-`p` ensemble without numerical differencing or root finding. A Newman-Ziff microcanonical implementation remains preferable when the whole curve and higher derivatives are also needed.

## Tiny exact-enumeration smoke tests

Independent exhaustive calculations performed while developing this note give the following finite-size values (not asymptotic estimates):

### Triangular site, natural 60-degree rhombus

```text
L=2: kappa_3 = -8/9               = -0.8888888889
L=3: kappa_3 = -2424832/1975509  = -1.2274466985
L=4: kappa_3 ≈ -1.3641177264
```

### Square bond, square torus

```text
L=2: kappa_3 = -2560/2187 ≈ -1.1705532693
L=3: kappa_3 ≈ -1.4555871991
```

These very small systems are only regression checks. They are qualitatively consistent with drift toward a more negative limit but cannot distinguish `-5/3` from nearby values.

## Pre-registered fit families for server output

When accurate `kappa_3(L)` values arrive, compare at least:

\[
\kappa_3(L)=k_*+aL^{-3/2},
\]

\[
\kappa_3(L)=k_*+aL^{-3/2}+bL^{-2},
\]

\[
\kappa_3(L)=k_*+aL^{-\omega}\quad\text{with omega trained only on smaller sizes}.
\]

Score by withheld-size prediction. Do not select `-5/3` because it looks visually attractive.

A useful hypothesis test is the residual

\[
L^{3/2}\left(\kappa_3(L)+5/3\right).
\]

If `-5/3` is exact and the leading correction is `L^-3/2`, this quantity should approach a finite constant. If it drifts systematically, at least one part of the hypothesis is wrong.

## Higher invariants

Measure

\[
\kappa_5=\frac{M^{(5)}}{(M')^5},\qquad
\kappa_7=\frac{M^{(7)}}{(M')^7}
\]

only after the derivative/noise pipeline is validated. A single rational-looking `kappa_3` is weak evidence. A cross-lattice sequence `(kappa_3,kappa_5,kappa_7,...)` is a meaningful fingerprint of the universal odd scaling function.
