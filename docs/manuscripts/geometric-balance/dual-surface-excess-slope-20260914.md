# Dual surface excess equals the inverse-mass slope

2026-09-14.  Exact finite-width identity plus a regular-point asymptotic consequence of the component-activity derivative theorem already on #739.

This note links three quantities that previously appeared separately:

1. the black subcritical complete-component logit score;
2. the exponentially large complementary white component's boundary--volume cancellation;
3. the mass slope `v=-kappa'(p)` that sets the median-centred Gumbel scale.

## 1. Setup

Fix black NN occupation probability

\[
p<p_c(G4),\qquad q=1-p>p_c(G8).
\]

Let

\[
\nu_w(p)=\nu_w^4(p)=\nu_w^8(q)                                \tag{1.1}
\]

be the exact complementary complete-winding-component intensity.

For a black component selected under component Palm, write

\[
n=|C|,\qquad b=|\partial_4 C|                                \tag{1.2}
\]

with **distinct external boundary sites**.

For the complementary white matching component selected under its component Palm at density `q`, write

\[
N=|C^*|,\qquad B=|\partial_8 C^*|.                             \tag{1.3}
\]

The two Palm laws are very different: black components have `O(w)` mean volume at fixed subcritical `p`, whereas the alternating-barrier theorem gives the white component an exponentially large longitudinal scale of order `1/nu_w`.

## 2. Exact black logit score

Use black logit coordinate

\[
z=\log\frac{p}{q},\qquad \partial_z=pq\partial_p.             \tag{2.1}
\]

A complete black component has activity

\[
p^nq^b.
\]

Therefore the exact finite-width component-Palm derivative identity is

\[
\boxed{
\partial_z\log\nu_w^4(p)
=E_{4,p}[q n-p b].}                                            \tag{2.2}
\]

This is the logit version already underlying the semiconvexity proof in `poisson-birth-windows.md`.

## 3. Exact white dual score in the **same black coordinate**

The natural white occupation parameter is `q`, whose own logit is

\[
z_8=\log\frac q p=-z.                                        \tag{3.1}
\]

At white density `q`, the component activity is

\[
q^N p^B.
\]

Its natural white-logit score is

\[
\partial_{z_8}\log\nu_w^8(q)
=E_{8,q}[pN-qB].                                               \tag{3.2}
\]

But `z_8=-z` and the exact complementary intensity identity says

\[
\nu_w^8(q)=\nu_w^4(p).
\]

Thus differentiating with respect to the **black** coordinate `z` gives

\[
\partial_z\log\nu_w^4(p)
=E_{8,q}[qB-pN].                                               \tag{3.3}
\]

Combining (2.2)--(3.3):

\[
\boxed{
E_{4,p}[q n-p b]
=E_{8,1-p}[q B-p N]
=\partial_z\log\nu_w(p).}                                    \tag{3.4}
\]

This is an exact finite-`w` dual surface-excess identity.

## 4. Regular-point limit gives the mass slope

The existing author-level component-intensity theorem gives

\[
\frac1w\log\nu_w(p)\to-\kappa_4(p),                           \tag{4.1}
\]

and at every differentiability point of `kappa_4`, uniformly for moving `p_w->p`,

\[
\frac1w\partial_p\log\nu_w(p_w)
\to-\kappa_4'(p)=:v_4(p)>0.                                   \tag{4.2}
\]

Since `partial_z=pq partial_p`, equation (3.4) yields

\[
\boxed{
\frac1wE_{4,p}[q n-p b]
\longrightarrow pq\,v_4(p),}                                 \tag{4.3}
\]

and, far more strikingly on the giant complementary side,

\[
\boxed{
\frac1wE_{8,q}[q B-p N]
\longrightarrow pq\,v_4(p).}                                 \tag{4.4}
\]

The white variables `N` and `B` are individually exponentially large in `w`, but their leading bulk pieces cancel in `qB-pN`, leaving an `O(w)` surface excess whose coefficient is exactly the black inverse-mass slope.

## 5. Bulk ratio and surface excess are two orders of the same identity

The reciprocal white-span theorem gives `E N >= E L ~ 1/nu_w`, so `E N` is exponentially larger than `w`.  Dividing the exact identity (3.4) by `E N` gives

\[
\frac{E B}{E N}=\frac p q+O(w\nu_w),                           \tag{5.1}
\]

recovering

\[
\boxed{\frac{E B}{E N}\to\frac p q.}                          \tag{5.2}
\]

Equation (4.4) supplies the **next order** after this bulk cancellation:

\[
\boxed{
E B
=\frac p q E N
+p\,v_4(p)w+o(w).}                                             \tag{5.3}
\]

Here we used `pq v w / q = p v w` when solving `q E B-p E N = pq v w+o(w)` for `E B`.

Thus the same score identity predicts both:

- the leading giant-component boundary/volume ratio `p/q`;
- the subleading `+p v_4 w` boundary excess.

No Brownian range assumption is needed.

## 6. An independent route to the Gumbel scale

The median-centred first-birth Gumbel variable on the black side uses

\[
v_4(p)=-\kappa_4'(p).                                        \tag{6.1}
\]

Usually one would estimate this from nearby masses or from the derivative of the log intensity itself.  Equation (4.4) gives a geometrically different estimator:

\[
\boxed{
v_4(p)
=\lim_{w\to\infty}
\frac{E_{white}[qB-pN]}{pq\,w}.}                              \tag{6.2}
\]

Similarly, (4.3) gives the black-side score estimator

\[
v_4(p)
=\lim_{w\to\infty}
\frac{E_{black}[q n-p b]}{pq\,w}.                              \tag{6.3}
\]

The two are exact dual finite-width views of the same derivative.  Agreement between them is a strong implementation/Palm-normalization check; neither requires numerical finite differences in `p`.

## 7. Why this is useful for #762

A joint morphology computation already interested in occupation count and **distinct** external boundary sites can report

\[
q n-p b
\]

for black component Palm essentially for free.

On the complementary white side, the raw `N,B` values are huge and strongly correlated.  The residual

\[
qB-pN                                                        \tag{7.1}
\]

is the meaningful surface quantity: the leading bulk terms are supposed to cancel.  Measuring `B/N` alone checks only the leading order; measuring (7.1)/`w` additionally tests the derivative/Gumbel normalization.

An edge-incidence boundary count cannot be substituted for `B`; the exact component activity uses distinct external sites and the cancellation coefficient would change.

## 8. Second-derivative warning

Differentiating once more gives exact variance/curvature identities, but on the giant white side `N+B` is exponentially large and must cancel against an equally large score variance to leave the scaled semiconvex curvature.  This is a delicate second-order object and should not be inferred from a modest finite covariance table.

The first-derivative surface-excess identity (3.4) is the robust target.
