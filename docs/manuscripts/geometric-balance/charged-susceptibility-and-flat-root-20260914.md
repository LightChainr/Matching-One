# Magnetic charged susceptibility and the quantitative flattening of the balance CDF

2026-09-14.  Integration of exact charge coordinates, fixed-width transfer, and the already-merged Pinson--Arguin primitive-sector continuum formula.

The key phenomenon is a product of two opposite effects:

- the **conditional charge odds** cross zero increasingly sharply with width/aspect;
- the total probability of being in either charged endpoint sector becomes exponentially small in aspect.

Their product is the ordinary CDF slope at the balance root.  This gives a quantitative mechanism for “balance without concentration.”

## 1. Exact finite identity

Recall

\[
\chi=P_0+P_2,
\qquad
\theta=\log(P_2/P_0),                                        \tag{1.1}
\]

and

\[
M=\chi\tanh(\theta/2),
\qquad
F=(1+M)/2.                                                    \tag{1.2}
\]

At the matching root `theta=0`,

\[
\boxed{
F'(p_*)=\frac14\chi(p_*)\theta'(p_*).}                       \tag{1.3}
\]

Thus a root can be uniquely selected by a huge conditional odds derivative even while the unconditional CDF is nearly flat, if `chi` is sufficiently small.

## 2. Continuum critical charged susceptibility on a long rectangle

Take the critical continuum torus

\[
\tau=i\rho,
\qquad \rho\to\infty.                                        \tag{2.1}
\]

The merged Pinson--Arguin primitive-sector evaluator gives for the horizontal primitive rank-one sector `{1,0}`

\[
\pi_{i\rho}(\{1,0\})
=\frac{
\theta_3(i\rho/6)-\theta_3(3i\rho/2)-2\theta_2(3i\rho/2)
}{2|\eta(i\rho)|^2}.                                         \tag{2.2}
\]

Use

\[
|\eta(i\rho)|^2
=e^{-\pi\rho/6}[1+O(e^{-2\pi\rho})],                         \tag{2.3}
\]

\[
\theta_3(i\rho/6)
=1+2e^{-\pi\rho/6}+O(e^{-2\pi\rho/3}),                       \tag{2.4}
\]

\[
\theta_3(3i\rho/2)=1+O(e^{-3\pi\rho/2}),                     \tag{2.5}
\]

\[
\theta_2(3i\rho/2)
=2e^{-3\pi\rho/8}[1+O(e^{-3\pi\rho})].                       \tag{2.6}
\]

Then

\[
\boxed{
\pi_{i\rho}(\{1,0\})
=1-2e^{-5\pi\rho/24}+O(e^{-\pi\rho/2}).}                     \tag{2.7}
\]

All nonhorizontal primitive rank-one sectors are `O(poly(rho)e^{-pi rho/2})` or smaller.  Therefore the total rank-one probability satisfies

\[
P_1^{cont}(i\rho)
=1-2e^{-5\pi\rho/24}+O(\operatorname{poly}(\rho)e^{-\pi\rho/2}).\tag{2.8}
\]

At `Q=1`, the trivial and cross topological sector probabilities are equal.  Hence

\[
P_0^{cont}=P_2^{cont}=\frac{1-P_1^{cont}}2,                   \tag{2.9}
\]

and the charged susceptibility is

\[
\boxed{
\chi_{cont}(\rho)
=P_0^{cont}+P_2^{cont}
=2e^{-5\pi\rho/24}[1+o(1)].}                                 \tag{2.10}
\]

The exponent is

\[
\frac{5\pi}{24}
=2\pi\frac5{48}
=2\pi x_m,                                                    \tag{2.11}
\]

exactly the magnetic cylinder gap used in the Jacobsen sector argument.

The leading coefficient `2` is fixed by the continuum wrapping formula; it is not a fitted transfer amplitude.

## 3. Fixed-width transfer version

At fixed `w`, let `p_w^{ch}` be the semi-infinite charge root and

\[
I_w^{ch}
=I^0_{4,w}(p_w^{ch})
=I^0_{8,w}(1-p_w^{ch}).                                      \tag{3.1}
\]

For the periodic topological trace, simple dominant eigenvalues give at large longitudinal length `m`

\[
P_0\sim e^{-mI_w^{ch}},
\qquad
P_2\sim e^{-mI_w^{ch}},                                      \tag{3.2}
\]

up to subleading trace eigenvalues.  Thus

\[
\boxed{
\chi(p_w^{ch})
\sim2e^{-mI_w^{ch}}.}                                        \tag{3.3}
\]

Similarly

\[
\theta'(p_w^{ch})
\sim m\Theta'_w(p_w^{ch}).                                   \tag{3.4}
\]

Insert these into the exact factorization (1.3):

\[
\boxed{
F'(p_w^{ch})
\sim\frac12\,m\Theta'_w(p_w^{ch})
\,e^{-mI_w^{ch}}.}                                           \tag{3.5}
\]

This is the fixed-width large-aspect quantitative version of balance without concentration.

The trace-amplitude coefficient one in (3.2) is the same periodic-sector issue discussed in `finite-aspect-charge-root-crossover-20260914.md`; if an exact finite rank transfer introduces a different sector normalization, (3.3)--(3.5) should be rechecked at that interface.

## 4. Continuum fixed-aspect scaling

Critical finite-size scaling gives

\[
I_w^{ch}
=\frac{2\pi x_m}{w}+o(w^{-1}),                               \tag{4.1}
\]

and

\[
\Theta'_w(p_w^{ch})
\sim A_\theta w^{-1/4}.                                      \tag{4.2}
\]

Take

\[
m=\rho w                                                     \tag{4.3}
\]

with `rho` fixed while `w->infinity`, and only afterwards let `rho` become large.  Equations (3.5), (4.1), (4.2) give

\[
\boxed{
w^{-3/4}F'(p_*)
\sim\frac{A_\theta}{2}\,
\rho e^{-2\pi x_m\rho}.}                                    \tag{4.4}
\]

For square-site percolation the current safe-transfer controls suggest

\[
A_\theta\approx3.4                                           \tag{4.5}
\]

in the raw Bernoulli `p` coordinate, so the prefactor in (4.4) is roughly `1.7` before final wide-width extrapolation.

In logit coordinate `h`, the corresponding first derivative amplitude is about `0.82` in the present widths.

## 5. Why odds sharpen while the CDF flattens

At large `rho`,

\[
\theta'(p_*)
\asymp \rho w^{3/4},                                         \tag{5.1}
\]

so the CONDITIONAL odds between the two charged sectors rotate through the root ever more steeply.

But

\[
\chi(p_*)
\asymp e^{-2\pi x_m\rho},                                   \tag{5.2}
\]

so almost no configurations lie in either charged sector at all.  Therefore

\[
F'(p_*)
\asymp \rho w^{3/4}e^{-2\pi x_m\rho},                        \tag{5.3}
\]

which decays exponentially in aspect after removing the ordinary thermal `w^(3/4)` factor.

This is the precise sense in which

```text
charge odds: sharp
unconditioned CDF: flat
```

can occur simultaneously.

## 6. Aspect of maximal normalized CDF slope

The large-`rho` envelope

\[
\rho e^{-2\pi x_m\rho}                                       \tag{6.1}
\]

has its maximum at

\[
\boxed{
\rho_{max}=\frac1{2\pi x_m}
=\frac{24}{5\pi}
\approx1.5279.}                                               \tag{6.2}
\]

The asymptotic formula itself is not expected to be quantitatively accurate all the way down to `rho≈1.5`, but (6.2) highlights that increasing aspect beyond a modest value can already reduce the median density even though the conditional charge sign becomes steeper.

## 7. Important nonuniformity for exponential aspect

Do **not** insert an exponentially growing

\[
\rho=e^{dw+o(w)}/w                                            \tag{7.1}
\]

directly into the continuum leading term (4.1).  The correction

\[
I_w^{ch}
=\frac{2\pi x_m}{w}+Cw^{-3}+\cdots                           \tag{7.2}
\]

is tiny per row but is multiplied by `m`; when `rho` is exponential, it changes the exponent by an enormous amount.

For exponential aspect the correct statement is the fixed-width one (3.5):

\[
\boxed{
F'(p_*)\sim\frac12m\Theta'_w e^{-mI_w^{ch}},}                \tag{7.3}
\]

using the actual finite-`w` charge gap.  The continuum formula (4.4) is a fixed-/moderately-growing-aspect scaling statement, not a uniform exponential-aspect approximation.

This nonuniformity is another example of why tiny finite-width free-energy corrections cannot be multiplied by an exponentially long direction without audit.

## 8. Interface to #767

The common-window charge coordinates now have explicit large-aspect boundary behavior:

\[
\chi(0,\rho)\sim2e^{-2\pi x_m\rho},                           \tag{8.1}
\]

while the dual-odd thermal scaling note proposes

\[
\theta(X,\rho)\approx\rho\mathcal F(X).                      \tag{8.2}
\]

Together,

\[
M(X,\rho)
=\chi(X,\rho)\tanh[\theta(X,\rho)/2]                         \tag{8.3}
\]

is a topology-compatible two-variable crossover representation with both its neutral plateau weight and its charge sign mechanism explicitly separated.

This is much more constrained than a generic bivariate Poisson/copula ansatz.

## 9. Claim boundary

The Pinson--Arguin expansion (2.7)--(2.10) is an analytic consequence of the already-merged continuum primitive-sector formula.  The exact finite factorization (1.3) is rigorous.  The square-site simultaneous finite-size scaling in sections 4--5 uses standard CFT/thermal assumptions and the safe-transfer amplitude sequence.  Section 7 explicitly states the failure of uniformity for exponential aspect.
