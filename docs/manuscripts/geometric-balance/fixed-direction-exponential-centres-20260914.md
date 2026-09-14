# Fixed integer directions: exact exponential-aspect birth centres

2026-09-14.  Author-level continuation of #739/#765.  This closes the fixed-primitive-integer-direction case without rotating the microscopic interaction and without importing a bond OZ theorem.

The proof uses two ingredients already established on this continuation branch:

1. `first-exit-torus-winding-upper-20260914.md` for the full-period upper bound;
2. the finite directional seed + Harris + translation-packing lower bound, which is the directional analogue of `exponential-birth-centres.md`.

A previously tempting direct comparison between quotient simple-cycle arcs and full-plane two-point events is **not** used: periodic reuse of Bernoulli variables makes that comparison unsafe.  The first-exit skeleton is the rigorous replacement.

## 1. Directional mass along a fixed integer vector

Fix a primitive integer vector

\[
u=(a,b)\in\mathbb Z^2,
\]

and choose `v in Z^2` with

\[
\det(u,v)=1.                                                   \tag{1.1}
\]

For `G=G4` (NN) or `G8` (matching), let

\[
\pi^G_p(x)=P_p^G(0\leftrightarrow x).                         \tag{1.2}
\]

The subcritical inverse-correlation norm is

\[
\tau_{G,p}(x)
=\lim_{k\to\infty}-\frac1k\log \pi^G_p(\lfloor kx\rfloor).    \tag{1.3}
\]

For the fixed integer direction put

\[
\xi_{G,u}(p)=\tau_{G,p}(u).                                   \tag{1.4}
\]

Equivalently, after normalizing the common occupied endpoint and applying Harris,

\[
\xi_{G,u}(p)
=\lim_{k\to\infty}-\frac1k\log s_p^G(ku),
\qquad
s_p(x)=\pi_p(x)/p.                                            \tag{1.5}
\]

For Euclidean unit direction `e=u/|u|`,

\[
\tau_{G,p}(e)=\xi_{G,u}(p)/|u|.                              \tag{1.6}
\]

## 2. The exponential torus adapted to u

Let

\[
\Lambda_{n,m}=\langle nu,mv\rangle.                           \tag{2.1}
\]

Its index is

\[
N=nm.                                                          \tag{2.2}
\]

Assume

\[
n\to\infty,
\qquad
\frac{\log m}{n}\to d\in(0,\infty).                          \tag{2.3}
\]

Because

\[
|\det(u,\alpha nu+\beta mv)|=|\beta|m,                       \tag{2.4}
\]

every period with `beta!=0` has Euclidean length at least `|beta|m/|u|`.  Hence the only periods on the `O(n)` scale are the multiples of `nu`.

The physical short length and transverse height are

\[
\ell_n=n|u|,
\qquad
h_n=N/\ell_n=m/|u|.                                           \tag{2.5}
\]

Thus if #765 uses

\[
\frac{\log h_n}{\ell_n}\to\delta,                            \tag{2.6}
\]

then

\[
d=\delta|u|.                                                   \tag{2.7}
\]

## 3. Rigorous full-period upper bound from first-exit certificates

Let

\[
f_G(n,m;p)=P_p^G(r_G>0).                                     \tag{3.1}
\]

`first-exit-torus-winding-upper-20260914.md` proves that for every fixed subcritical `p` and every `epsilon in (0,1)`, once the torus is large enough,

\[
\boxed{
f_G(n,m;p)
\le N C_{p,\epsilon}
\sum_{\lambda\in\Lambda_{n,m}\setminus0}
\exp[-(1-\epsilon)\tau_{G,p}(\lambda)].}                     \tag{3.2}
\]

The proof is local and quotient-safe:

- choose a shortest nonzero-homology occupied closed walk; it is a torus-vertex-simple cycle;
- apply the finite first-exit skeleton to its lift;
- successive witness interiors are disjoint actual quotient SITE variables;
- each local translate injects, so its coefficient is the genuine planar first-exit coefficient;
- site BK bounds a skeleton by the product of these local coefficients;
- a compact certified Wulff body gives a pointwise factor `exp[-h_T(lambda)]`;
- large first-exit boxes exhaust `(1-epsilon)K_p`, whose support function is `(1-epsilon)tau_p`.

For the period lattice (2.1), parallel periods have

\[
\tau_{G,p}(\alpha nu)=|\alpha|n\xi_{G,u}(p),                 \tag{3.3}
\]

while nonparallel periods cost at least `c_p|beta|m/|u|`.  Therefore the theta sum in (3.2) satisfies

\[
\sum_{\lambda\ne0}e^{-(1-\epsilon)\tau_{G,p}(\lambda)}
\le
C_{p,\epsilon}e^{-(1-\epsilon)n\xi_{G,u}(p)}
+e^{-c_{p,\epsilon}m}.                                        \tag{3.4}
\]

Since `N=nm`,

\[
\limsup_{n\to\infty}\frac1n\log f_G(n,m;p)
\le
-\max\{(1-\epsilon)\xi_{G,u}(p)-d,0\}.                       \tag{3.5}
\]

Let `epsilon downarrow 0`:

\[
\boxed{
\limsup\frac1n\log f_G(n,m;p)
\le
-\max\{\xi_{G,u}(p)-d,0\}.}                                  \tag{3.6}
\]

The complete directional mass appears; no half-period loss remains.

## 4. Lower bound: repeat one finite directional seed and close the seam

Fix `p<pc(G)` and `epsilon>0`.  Choose a fixed integer `L` and a finite box around the segment from `0` to `Lu` such that the conditional finite-box connection probability `q` obeys

\[
q>
\exp[-(\xi_{G,u}(p)+\epsilon/2)L].                            \tag{4.1}
\]

This is possible by the definition of the directional mass and finite-box exhaustion.

Write

\[
n=kL+r,
\qquad0\le r<L.                                               \tag{4.2}
\]

Translate the seed `k` times by `Lu` and use a fixed finite connector for the remaining displacement `ru`.  Conditional Harris at the shared endpoints gives, for all large `n`,

\[
P_p(\text{one prescribed u-ring})
\ge e^{-(\xi_{G,u}(p)+\epsilon)n}.                            \tag{4.3}
\]

The ring support lies in a fixed-width tube around `[0,nu]`, so

\[
|S-S|\le C_{p,\epsilon,u}n.                                   \tag{4.4}
\]

The finite-group packing lemma supplies at least `c m` pairwise site-disjoint translates.  Hence

\[
\boxed{
1-f_G(n,m;p)
\le
\exp\left[-c_{p,\epsilon,u}m
 e^{-(\xi_{G,u}(p)+\epsilon)n}\right].}                       \tag{4.5}
\]

## 5. Exact directional winding-event rate

Combining (3.6) and (4.5), for every fixed subcritical `p`,

\[
\boxed{
\lim_{n\to\infty}\frac1n\log f_G(n,m_n;p)
=-\max\{\xi_{G,u}(p)-d,0\}.}                                 \tag{5.1}
\]

At equality `xi=d`, only the zero exponential rate is asserted.

In the Euclidean normalization `e=u/|u|`, `delta=d/|u|`,

\[
\boxed{
\lim\frac1{\ell_n}\log f_G
=-\max\{\tau_{G,p}(e)-\delta,0\}.}                            \tag{5.2}
\]

This is the fixed-direction rate statement requested in #765.

## 6. Strict p-monotonicity in every fixed integer direction

The continuity argument from `exponential-birth-centres.md` extends verbatim along `ku`:

- the finite-seed infimum gives left continuity;
- the uniform subcritical cluster-volume tail gives right continuity after truncating the connecting cluster at `O(n)` sites.

To prove strict decrease, suppose

\[
0<p<q<pc(G),
\qquad
\xi_{G,u}(p)=\xi_{G,u}(q)=A>0.                               \tag{6.1}
\]

Choose `d<A` and `zeta>0` so that the Friedgut--Kalai increment satisfies

\[
\rho\frac{A-d+\zeta}{d}<q-p.                                 \tag{6.2}
\]

Take `m=ceil(e^{dn})`.  From (5.1),

\[
f_G(n,m;p)>\epsilon_n=e^{-(A-d+\zeta)n}                       \tag{6.3}
\]

for large `n`.  The event `r_G>0` is increasing and invariant under the transitive torus translation group on `N=nm` sites, so Friedgut--Kalai forces

\[
f_G(n,m;q)>1-\epsilon_n\to1.                                 \tag{6.4}
\]

But (3.6) with `xi(q)=A>d` gives `f_G(n,m;q)->0`, a contradiction.  Hence

\[
\boxed{p\mapsto\xi_{G,u}(p)\text{ is strictly decreasing on }(0,pc(G)).}\tag{6.5}
\]

The small-`p` path count gives `xi->infinity` as `p->0`.  As `p->pc(G)`, square symmetry and the norm triangle inequality give

\[
\xi_{G,u}(p)
\le(|a|+|b|)\kappa_G(p)\to0.                                 \tag{6.6}
\]

Thus `xi_{G,u}` continuously and strictly maps `(0,pc(G))` onto `(infinity,0)`.

## 7. Fixed-direction birth centres

Let `T_1,T_2` be the two NN rank births on `Lambda_{n,m}`.  Since

\[
P(T_1\le p)=f_{G4}(n,m;p),                                    \tag{7.1}
\]

(5.1) and strict invertibility imply

\[
\boxed{
T_1\xrightarrow{P}a_u(d),
\qquad
\xi_{4,u}(a_u(d))=d.}                                        \tag{7.2}
\]

Digital Alexander duality gives

\[
P_p^{G4}(r=2)=P_{1-p}^{G8}(r=0),                              \tag{7.3}
\]

so

\[
\boxed{
T_2\xrightarrow{P}b_u(d),
\qquad
b_u(d)=1-c_u(d),
\qquad
\xi_{8,u}(c_u(d))=d.}                                        \tag{7.4}
\]

Equivalently, with `e=u/|u|` and `delta=d/|u|`,

\[
\boxed{
\tau_{4,a(e,\delta)}(e)=\delta,
\qquad
\tau_{8,1-b(e,\delta)}(e)=\delta.}                           \tag{7.5}
\]

## 8. Strict matching mass gap in every fixed direction

The endpoint-direction-independent enhancement argument gives, on every compact parameter interval inside `(0,pc(G8))`, a `delta_I>0` such that

\[
\boxed{
P_p^{G8}(0\leftrightarrow x)
\ge
P_{p+\delta_I}^{G4}(0\leftrightarrow x)}                     \tag{8.1}
\]

for distant endpoints `x`, uniformly in direction.  Taking the logarithmic rate along `ku` gives

\[
\xi_{8,u}(p)
\le\xi_{4,u}(p+\delta_I).                                    \tag{8.2}
\]

By (6.5),

\[
\xi_{4,u}(p+\delta_I)<\xi_{4,u}(p).                           \tag{8.3}
\]

Therefore

\[
\boxed{
\xi_{8,u}(p)<\xi_{4,u}(p),
\qquad0<p<pc(G8),}                                            \tag{8.4}
\]

for every fixed primitive integer direction.

So the strict complement-odd centre displacement proved axially is not axis-specific.

## 9. What remains of #765

This note closes the fixed primitive integer direction with an exponentially long Bezout-complement period.

The genuinely remaining part of #765 is the **varying-direction** problem `u_n/|u_n|->e`, where the arithmetic of `u_n` changes and one needs uniformity of:

- the directional norm convergence;
- finite first-exit certificates near the moving support point;
- finite-seed lower bounds;
- separation from other comparable period classes.

The projective Poisson hard-core note on this branch gives the correct topology if several direction classes remain competitive.  No new angle scan is required to finish the fixed-direction theorem.