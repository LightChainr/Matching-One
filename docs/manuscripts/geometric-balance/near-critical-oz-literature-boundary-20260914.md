# What the 2026 near-critical OZ theorem does and does not buy us

2026-09-14.  Literature boundary audit for #740/#758/#762/#767 after reading the current v3 of Lucas D'Alimonte and Ioan Manolescu, *Near-critical Ornstein--Zernike theory for the planar random-cluster model*, arXiv:2510.13648v3 (last revised 2026-06-23).

This note is deliberately strict about model/observable transport.  The paper is highly relevant to the architecture of our conjectures, but it is **not** a theorem about square-site percolation or complete winding-component activities.

## 1. Exact model boundary

D'Alimonte--Manolescu study the planar random-cluster model on the square lattice with configurations on **edges**, for

\[
1\le q<4,
\qquad p<p_c(q).                                              \tag{1.1}
\]

At `q=1` this is Bernoulli **bond** percolation, not square-site percolation.

Their correlation length in direction `e` is

\[
\xi_p(e)
=\left[
\lim_{n\to\infty}-\frac1n
\log\phi_p(0\leftrightarrow\lfloor ne\rfloor)
\right]^{-1}.                                                 \tag{1.2}
\]

Thus in the notation used on this branch,

\[
\boxed{\kappa_p(e)=\xi_p(e)^{-1}.}                            \tag{1.3}
\]

Their `xi` is a length; our `kappa/tau` is an inverse length / correlation-norm rate.

## 2. The uniform near-critical two-point theorem

Their Theorem 1.1 states, uniformly in direction, `p in [epsilon,pc)` and physical distance `r>=xi_p(e)`,

\[
\boxed{
\phi_p(0\leftrightarrow\lfloor re\rfloor)
\asymp
\pi_1(\xi_p(e))^2
\sqrt{\frac{\xi_p(e)}r}
\exp[-r/\xi_p(e)].}                                          \tag{2.1}
\]

Here `pi_1(R)` is the **critical** one-arm probability to scale `R`.

If

\[
s=r/\xi_p(e)=\kappa_p(e)r,                                   \tag{2.2}
\]

then

\[
\boxed{
\phi_p(0\leftrightarrow re)
\asymp
\pi_1(\xi_p(e))^2 s^{-1/2}e^{-s}.}                           \tag{2.3}
\]

The formula separates three mechanisms:

1. `e^{-s}`: correlation-mass cost;
2. `s^{-1/2}`: one-dimensional OZ/local-CLT sewing;
3. `pi_1(xi)^2`: two near-critical endpoint insertions.

This is exactly the conceptual decomposition we want for a crossover theory.

## 3. Important limitation: this is comparability, not an exact amplitude

The theorem uses

\[
\asymp,                                                       \tag{3.1}
\]

with constants uniform in `p`, direction and distance in the stated range.  It does **not** identify a multiplicative constant with relative error `1+o(1)` uniformly as `p->pc`.

Therefore it cannot by itself supply:

- an exact near-critical OZ residue;
- an exact complete-component prefactor `zeta(p)`;
- a logarithmic centering constant for an extreme-value window when `O(1)` amplitude information matters.

It does rigorously show that a naive fixed-`p` smooth amplitude continued to criticality is the wrong structural picture: the critical one-arm factor appears explicitly.

## 4. The fixed-p prefactor should not be analytically continued unchanged

For fixed subcritical `p`, `xi_p` and `pi_1(xi_p)` are constants as `r->infinity`, so (2.1) has the familiar

\[
r^{-1/2}e^{-r/\xi_p}                                         \tag{4.1}
\]

shape.

But when `p=p_r->pc` and `xi_p->infinity`, the quantity

\[
\boxed{\pi_1(\xi_p)^2}                                       \tag{4.2}
\]

moves with the parameter and belongs at leading prefactor level.

Thus a crossover ansatz of the form

\[
A(p)r^{-1/2}e^{-\kappa(p)r}                                  \tag{4.3}
\]

with `A(p)` treated as a harmless smooth fixed-`p` residue is not uniformly justified near criticality.  Even in the exactly treated bond-FK model, the endpoint insertion factor becomes a critical-arm observable.

This strongly supports the repository decision not to use a fixed-p OZ amplitude as the foundation of #767.

## 5. Brownian scale: D is of order the correlation length

The paper constructs a killed Markov renewal process at a coarse scale `L(p)` comparable uniformly with the correlation length.  Its endpoint local CLT has transverse variance parameter `sigma(p,w)` bounded uniformly away from zero and infinity.

For a connection of physical longitudinal length `r`, the number of renewal steps is

\[
n\asymp r/L(p),                                               \tag{5.1}
\]

while Theorem 4.10 scales transverse displacement by

\[
\sqrt n\,L(p).                                                \tag{5.2}
\]

Therefore the transverse variance is of order

\[
nL(p)^2\asymp rL(p).                                         \tag{5.3}
\]

The same paper identifies

\[
\xi_p(e)
=L(p)\times\text{a factor uniformly bounded above/below}.     \tag{5.4}
\]

Consequently, in the bond-FK theorem,

\[
\boxed{
\operatorname{Var}(X_\perp\mid0\leftrightarrow re)
\asymp r\xi_p(e).}                                            \tag{5.5}
\]

If we parameterize a Brownian bridge by physical longitudinal distance,

\[
\operatorname{Var}(X_\perp(t))\sim D_p(e)\,r\,t(1-t),         \tag{5.6}
\]

then the rigorous scaling analogue is

\[
\boxed{D_p(e)\asymp\xi_p(e)=\kappa_p(e)^{-1}.}               \tag{5.7}
\]

Hence

\[
\boxed{D_p(e)\kappa_p(e)\asymp1}                             \tag{5.8}
\]

uniformly up to criticality in this bond-FK model.

This is strong structural support for the repository relation

\[
D^{-1}=\partial_{yy}\tau                                      \tag{5.9}
\]

having the correct **scale** near criticality.  It does not prove the exact equality or its SITE version.

## 6. Strict Wulff geometry survives uniformly near criticality

D'Alimonte--Manolescu prove for every subcritical `p` that the correlation-length unit ball and its Wulff dual are strictly convex with differentiable boundaries.  Their construction is performed at correlation-length scale and is uniform in `p` and direction at the renewal level.

This independently confirms the qualitative geometric assumptions used in our loop/branch analysis:

- a unique tangent/dual direction exists;
- the transverse endpoint law has a nondegenerate Gaussian scale;
- nearby directional deviations have a genuine local large-deviation penalty.

However their theorem, as stated, does not give us a ready-made uniform numerical lower/upper bound on the second directional derivative that can simply be transplanted to square-site percolation.

## 7. Consequence for the intermediate aspect-ratio crossover

Consider a winding length `w` and let

\[
\xi=\xi(p),
\qquad
s=w/\xi=\kappa(p)w.                                          \tag{7.1}
\]

If a rare-event opportunity count is `m` (or an aspect ratio proportional to it), the **two-point** OZ balance suggested by (2.3) is

\[
1\asymp
m\,\pi_1(\xi)^2s^{-1/2}e^{-s}.                               \tag{7.2}
\]

Taking logarithms,

\[
\boxed{
s+\tfrac12\log s-2\log\pi_1(\xi)
=\log m+O(1).}                                                \tag{7.3}
\]

This immediately separates three crossover scales.

### Regime A: deep rare-event/OZ merging

If

\[
\log m\gg |\log\pi_1(\xi)|                                  \tag{7.4}
\]

(and `log m->infinity`), then at leading order

\[
s\sim\log m.                                                  \tag{7.5}
\]

The critical endpoint dressing shifts only the lower-order centering.

### Regime B: arm-dressed logarithmic crossover

If

\[
\log m\asymp |\log\pi_1(\xi)|,                               \tag{7.6}
\]

then the critical one-arm term contributes at the same logarithmic order as the opportunity entropy.  It cannot be absorbed into an `O(1)` OZ amplitude.

Polynomial aspect ratios are a natural place where this can happen.

### Regime C: fixed/bounded aspect

When `log m=O(1)`, the system is no longer an extreme-value gas of many independent long opportunities.  The full near-critical torus scaling theory is the natural object instead of an OZ birth-window continuation.

This regime split is only a **template** for square-site component activity; equation (7.2) is a bond-FK two-point statement.

## 8. Complete component activity has an additional insertion problem

Our torus birth intensity is not

\[
P(0\leftrightarrow we_1).                                    \tag{8.1}
\]

It counts complete winding components with an anchor/unrooting convention and with all branches included.

Even if the same renewal skeleton controls the long core, a component activity can differ from the two-point function through:

- anchor/root removal;
- closure/seam insertion;
- the requirement that the component be complete;
- local branches and boundary weights;
- possibly multiple leading renewal bands.

The matrix-sewing result on this branch shows that finite local memory itself does **not** create an arbitrary residue, but it does not identify the correct insertion vector for the actual SITE complete-component object.

Therefore the literature theorem supports

\[
\text{near-critical renewal skeleton + }w^{-1/2}\text{ Gaussian sewing},\tag{8.2}
\]

not the equality of the two-point and complete-component amplitudes.

## 9. Model-transfer boundary: bond FK is not square-site

At `q=1`, D'Alimonte--Manolescu give a rigorous near-critical theory for Bernoulli **bond** percolation on the square lattice.

The standard near-critical scaling-limit theorem of Garban--Pete--Schramm is for **site percolation on the triangular lattice**.

Neither result is a theorem asserting the same uniform OZ formula for square-site percolation, and neither identifies the square-site complete-component insertion amplitude.

So for the actual Matching-One model, using (2.1) as a theorem would require a new SITE adaptation or a separate universality theorem strong enough to transport the relevant quantitative observable.  We do not have such a theorem in this audit.

## 10. What can safely enter our research programme now

### Safe structural import / proof template

Use the paper as strong evidence that a successful square-site near-critical renewal theorem should have:

1. coarse slices at the correlation-length scale;
2. a killed Markov renewal process with uniform mass gap;
3. nondegenerate endpoint local CLT;
4. Brownian bridge scale `sqrt(w xi)`;
5. strict Wulff geometry;
6. critical arm insertions at the endpoints.

### Not safe to import as a SITE theorem

Do **not** claim from this paper alone:

- square-site uniform near-critical OZ;
- an exact `1+o(1)` amplitude;
- a SITE value of `zeta(p)`;
- a complete-component activity formula;
- exact square-site critical arm exponents.

## 11. A focused theorem target suggested by the literature

The full D'Alimonte--Manolescu machinery is more than we need for #767.  A strategically smaller square-site target would be:

> **SITE near-critical complete-component comparability.**  For `p<pc`, `w>=C xi_site(p)`, prove uniformly in a compact angular set that the complete winding-component activity satisfies
> \[
> \nu_w(p)
> \asymp
> R_{site}(p,w/\xi)
> \pi_{1,site}(\xi)^2
> (w/\xi)^{-1/2}
> e^{-w/\xi},
> \]
> where `R_site` is bounded above and below uniformly and encodes the anchor/closure insertion.

Even this **comparability-level** theorem would already determine which terms can enter the logarithmic crossover equation and would distinguish the arm-dressed and deep-OZ regimes.  Exact prefactor identification could come later.

A still weaker first target is to prove only

\[
D_p\asymp\xi_p                                                \tag{11.1}
\]

for the SITE point-to-point conditioned cluster uniformly near criticality.  Combined with the branch's curvature/sewing diagnostics, that would sharply constrain the possible complete-component mechanism.

## 12. Bottom line for #767

The 2026 near-critical OZ paper **supports the architecture** of our crossover programme but does not close it for square-site Matching-One.

The most important correction to prior intuition is:

\[
\boxed{\text{near critical: OZ Gaussian sewing is arm-dressed, not merely amplitude-renormalized.}}\tag{12.1}
\]

Accordingly, #767 should not analytically continue a fixed-p component prefactor to `p_c`.  It should either:

1. obtain a square-site near-critical renewal/comparability theorem with explicit arm insertions; or
2. state the common-window theory conditionally in the exact charge-neutral coordinates `(chi,theta,H)` and treat the one-arm/complete-component insertion as an independent scaling input.

Primary sources audited:

- L. D'Alimonte and I. Manolescu, *Near-critical Ornstein--Zernike theory for the planar random-cluster model*, arXiv:2510.13648v3, 2026.
- C. Garban, G. Pete and O. Schramm, *The scaling limits of near-critical and dynamical percolation*, JEMS 2018 / arXiv:1305.5526v4, triangular-lattice site percolation.