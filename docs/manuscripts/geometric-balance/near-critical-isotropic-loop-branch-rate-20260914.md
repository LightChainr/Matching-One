# Conditional near-critical isotropic loop/branch rate: an explicit universal target

2026-09-14.  A deliberately bold but sharply falsifiable consequence of the #758 loop/branch candidate under one extra universality hypothesis: restoration of rotational invariance of the normalized inverse-correlation norm near criticality.

This is **not** asserted as a square-site theorem.  Rotational restoration/Wulff-rounding is rigorous for triangular-lattice site percolation through its near-critical scaling limit, but is not certified here for the square-site NN/matching pair.

## 1. Input hypotheses

Let

\[
\kappa_p=\tau_p(1,0).                                        \tag{1.1}
\]

Assume the #758 candidate span rate

\[
I_p(A)
=\min_{0\le r\le A}
\left\{
\tau_p(1,2r)-\kappa_p+\kappa_p(A-r)
\right\}.                                                     \tag{1.2}
\]

Assume additionally the near-critical isotropy statement

\[
\boxed{
\frac{\tau_p(x)}{\kappa_p}\longrightarrow |x|
\quad\text{locally uniformly in }x,
\qquad p\uparrow p_c.}                                       \tag{1.3}
\]

The local uniformity may be weakened to the compact set of directions/ratios used below.

Then

\[
\frac{I_p(A)}{\kappa_p}
\longrightarrow
J(A)
:=\min_{0\le r\le A}
\left[\sqrt{1+4r^2}-1+A-r\right].                            \tag{1.4}
\]

## 2. The optimal bulge is explicit

Let

\[
g(r)=\sqrt{1+4r^2}-1-r.                                     \tag{2.1}
\]

Then

\[
g'(r)=\frac{4r}{\sqrt{1+4r^2}}-1.                           \tag{2.2}
\]

The unique stationary point is

\[
16r^2=1+4r^2,                                                  \tag{2.3}
\]

so

\[
\boxed{
r_*=\frac1{\sqrt{12}}=\frac1{2\sqrt3}
=0.2886751345948129\ldots.}                                  \tag{2.4}
\]

Since `g'` is negative below `r_*` and positive above it, the constrained minimizer is

\[
\boxed{r_{opt}(A)=\min\{A,r_*\}.}                            \tag{2.5}
\]

Thus the near-critical isotropic mechanism predicts a **universal saturation bulge**, independent of `A` once `A>r_*`.

## 3. Explicit piecewise rate function

For `A<=r_*`, the constraint is active and `r=A`.  Therefore

\[
\boxed{
J(A)=\sqrt{1+4A^2}-1,
\qquad
0\le A\le\frac1{2\sqrt3}.}                                  \tag{3.1}
\]

For `A>=r_*`, insert the stationary point.  Since

\[
\sqrt{1+4r_*^2}=\frac2{\sqrt3},                               \tag{3.2}
\]

and

\[
\frac2{\sqrt3}-\frac1{2\sqrt3}=\frac{\sqrt3}{2},             \tag{3.3}
\]

we obtain

\[
\boxed{
J(A)=A+\frac{\sqrt3}{2}-1,
\qquad
A\ge\frac1{2\sqrt3}.}                                       \tag{3.4}
\]

The derivative matches at the transition because

\[
\left.\frac{4A}{\sqrt{1+4A^2}}\right|_{A=r_*}=1.             \tag{3.5}
\]

So `J` is `C^1`, strictly convex before `r_*`, and exactly linear afterwards.

## 4. Small-span Brownian overlap

Expand (3.1):

\[
J(A)=2A^2-2A^4+O(A^6).                                       \tag{4.1}
\]

Hence

\[
\boxed{I_p(A)\sim2\kappa_p A^2}                              \tag{4.2}
\]

under the isotropic hypothesis.

The Brownian-bridge range moderate-deviation exponent used on this branch is

\[
2A^2/D_p.                                                      \tag{4.3}
\]

Matching (4.2) therefore gives the particularly sharp near-critical prediction

\[
\boxed{D_p\kappa_p\to1}                                      \tag{4.4}
\]

in the physical normalization used by the rate ansatz.

The 2026 bond-FK renewal theorem rigorously supports only the weaker scale relation

\[
D_p\kappa_p\asymp1,                                          \tag{4.5}
\]

so (4.4) is a genuine normalization/universality test rather than a consequence of that literature.

## 5. No-parameter targets for #762

The #762 morphology plan singled out conditioned macroscopic span thresholds `A=1` and `A=2`.

Both satisfy

\[
A>r_*                                                         \tag{5.1}
\]

by a wide margin, so the isotropic loop/branch mechanism predicts that both are already in the **linear branch regime**.

The large-deviation costs are

\[
\boxed{
\frac{I_p(1)}{\kappa_p}\longrightarrow\frac{\sqrt3}{2}
=0.8660254037844386\ldots,}                                  \tag{5.2}
\]

\[
\boxed{
\frac{I_p(2)}{\kappa_p}\longrightarrow1+\frac{\sqrt3}{2}
=1.8660254037844386\ldots.}                                  \tag{5.3}
\]

In particular,

\[
\boxed{
\frac{I_p(2)-I_p(1)}{\kappa_p}\longrightarrow1.}             \tag{5.4}
\]

Thus one extra unit of macroscopic span costs asymptotically exactly one axial mass unit after the core bulge has saturated.

## 6. Morphological core prediction

At any fixed `A>r_*`, the minimizing winding core uses the same normalized bulge

\[
\boxed{r_{core}\to1/(2\sqrt3).}                              \tag{6.1}
\]

The remaining span `A-r_*` is carried by branches.

If `r` represents half of the core transverse range in the #758 geometry, the total core range target is

\[
\boxed{2r_* = 1/\sqrt3=0.5773502691896258\ldots}              \tag{6.2}
\]

in units of the winding length.

Therefore the conditioned `A=1` and `A=2` samples should have asymptotically the **same core-width distribution after normalization by w**, while their branch lengths differ macroscopically.

This is much stronger than checking only the tail probability.

## 7. A parameter-free shape of the entire candidate rate

The normalized curve

\[
J(A)=
\begin{cases}
\sqrt{1+4A^2}-1,& A\le1/(2\sqrt3),\\
A+\sqrt3/2-1,& A\ge1/(2\sqrt3)
\end{cases}                                                   \tag{7.1}
\]

has no free fit parameters once `A` and the rate are normalized by `w` and `kappa_p`.

A near-critical width sequence can therefore challenge the mechanism by estimating the ratios

\[
-\frac1{\kappa_pw}\log P(R_w\ge Aw)                           \tag{7.2}
\]

at several `A`.  A systematic limit different from (7.1) rules out at least one of:

1. normalized Wulff isotropy;
2. the loop/branch variational form (1.2);
3. the identification of the measured range with the variational `A`.

## 8. Why square-site isotropy must remain a separate hypothesis

Critical RSW and strict Wulff convexity do **not** imply that

\[
\tau_p/\kappa_p\to|\cdot|.                                   \tag{8.1}
\]

The latter is rotational restoration, a stronger universality statement.

For triangular-lattice critical/near-critical site percolation, rotational invariance of the scaling limit implies Wulff-rounding and this type of Euclidean normalization.  For the actual square-site Matching-One model, no such theorem was identified in the current audit.

Keeping isotropy as one explicit hypothesis is useful: it turns all the remaining morphology predictions into exact numbers and makes the missing universality input impossible to hide inside a fitted `D` or curvature.

## 9. Strongest immediate test order

The cheapest discriminators are:

1. verify the already conditional identity `D^{-1}=partial_yy tau` away from criticality where both can be measured/certified;
2. approach criticality and test whether `D kappa` drifts toward one;
3. test core bulge saturation near `0.288675 w` at `A=1,2`;
4. test rate differences `I(2)-I(1)≈kappa`;
5. only after these pass, fit the full curve (7.1).

This order distinguishes a sewing/curvature failure from a rotational-universality failure.

## 10. Claim boundary

Everything after assumption (1.3) is elementary convex optimization.  The square-site isotropy assumption itself is conjectural here.  The note is intended as a no-parameter challenge target, not as a theorem claim.