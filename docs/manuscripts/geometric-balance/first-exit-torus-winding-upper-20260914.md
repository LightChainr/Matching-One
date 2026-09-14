# First-exit Wulff certificates give a full-period torus winding upper bound

2026-09-14.  This note is the rigorous bridge from `vector-first-exit-domain-20260914.md` to the directional/exponential-torus questions in #765.  It replaces a tempting but unsafe direct comparison between quotient paths and full-plane two-point events.

The result is general: on any sufficiently large honest integer-period torus, positive homology is bounded by a theta sum over **complete period vectors**, with their full planar correlation-norm cost.

## 1. Setup

Work with a translation-invariant finite-range independent site model `G` on `Z^2` at a fixed subcritical parameter `p`.

Let

\[
\tau_p(x)
\]

be the homogeneous inverse-correlation norm and let

\[
K_p=\{t:t\cdot x\le\tau_p(x)\ \forall x\}                    \tag{1.1}
\]

be its polar/exponential-moment body.  By the first-exit theorem on this branch,

\[
\mathcal D_p=K_p^\circ                                        \tag{1.2}
\]

and every compact `T subset K_p^circ` is contained in one sufficiently large finite first-exit certificate

\[
\mathcal C_S=\{t:B_S(t)<1\}.                                  \tag{1.3}
\]

Let `Lambda` be a rank-two period lattice and

\[
T_\Lambda=\mathbb Z^2/\Lambda,
\qquad
N=[\mathbb Z^2:\Lambda].                                     \tag{1.4}
\]

Assume the shortest period is large enough that every translate of the chosen finite set `S`, enlarged by the finite interaction range, injects into the quotient.

Write

\[
f_\Lambda(p)=P_p(r>0).                                       \tag{1.5}
\]

## 2. Every positive-rank configuration contains a simple nonzero-homology cycle

If `r>0`, choose an occupied closed walk with nonzero ambient homology and minimum edge length.  It cannot repeat a torus vertex other than its start/end.  A repeated vertex would split the walk into two shorter closed walks whose homology classes sum to the original nonzero class; at least one shorter piece would still have nonzero homology.

Hence a minimum witness is a graph-simple occupied cycle.

Lift it from a starting vertex `z` to a self-avoiding path in `Z^2` from a chosen lift of `z` to

\[
z+\lambda,
\qquad
\lambda\in\Lambda\setminus\{0\}.                             \tag{2.1}
\]

The lifted endpoint is a different plane vertex but the same torus vertex.

## 3. Local first-exit skeletons remain BK-disjoint on the quotient

Apply the first-exit decomposition of `vector-first-exit-domain-20260914.md` to the lifted simple cycle path using translates of `S`.

Starting from `y_0=z`, obtain exit displacements

\[
v_1,\ldots,v_k\in\partial_{ext}S                             \tag{3.1}
\]

and a final residual `s in S` such that

\[
\lambda=v_1+\cdots+v_k+s.                                    \tag{3.2}
\]

The `i`th first-exit witness uses only the actual cycle vertices from `y_{i-1}` through the internal neighbour immediately preceding `y_i`; the exit site itself is excluded and becomes the starting site of the next witness.

Because the torus cycle is vertex-simple, these witness sets are disjoint **as quotient site variables**, not merely as plane lifts.  This is the point that fails for a generic pair of quotient connection events but holds for the canonical simple-cycle skeleton.

Each witness is local in one translate of the fixed finite set `S`; injectivity of that translate means its probability is exactly the planar coefficient

\[
b_S(v_i;p).                                                   \tag{3.3}
\]

Repeated site BK on the quotient therefore gives, for any fixed skeleton,

\[
P(\text{that skeleton is witnessed})
\le\prod_i b_S(v_i;p).                                       \tag{3.4}
\]

No event is imposed on the final residual, so the fact that the terminal torus vertex equals the initial one creates no shared witness variable.

## 4. Pointwise renewal bound for a homology vector

Define the first-exit renewal majorant

\[
R_S(x)=
\sum_{k\ge0}
\sum_{v_1,\ldots,v_k}
\sum_{s\in S}
1_{\{x=v_1+\cdots+v_k+s\}}
\prod_i b_S(v_i;p).                                           \tag{4.1}
\]

Section 3 proves that, for a fixed start vertex, the probability that a canonical simple winding witness has lift class `lambda` is at most

\[
R_S(\lambda).                                                  \tag{4.2}
\]

For any `t in C_S`, the same generating-function calculation as in the planar first-exit theorem gives

\[
\sum_xR_S(x)e^{t\cdot x}
=\frac{\sum_{s\in S}e^{t\cdot s}}{1-B_S(t)}.                 \tag{4.3}
\]

Hence pointwise

\[
R_S(\lambda)
\le C_S(t)e^{-t\cdot\lambda},
\qquad
C_S(t)=\frac{\sum_{s\in S}e^{t\cdot s}}{1-B_S(t)}.           \tag{4.4}
\]

Let `T` be any compact subset of `C_S`.  The denominator stays uniformly away from zero, so

\[
C_T:=\sup_{t\in T}C_S(t)<\infty.                              \tag{4.5}
\]

Choosing the best `t in T` for each `lambda`,

\[
\boxed{
R_S(\lambda)
\le C_T e^{-h_T(\lambda)},
\qquad
h_T(\lambda)=\sup_{t\in T}t\cdot\lambda.}                    \tag{4.6}
\]

## 5. General torus winding upper bound

There are `N` possible starting torus vertices.  Union over the nonzero lift class of the canonical simple cycle gives

\[
\boxed{
f_\Lambda(p)
\le N C_T
\sum_{\lambda\in\Lambda\setminus\{0\}}
 e^{-h_T(\lambda)}.}                                          \tag{5.1}
\]

This is already a rigorous finite certificate for any computed first-exit body `T`.

Now fix `epsilon in (0,1)` and take

\[
T_\epsilon=(1-\epsilon)K_p.                                  \tag{5.2}
\]

Because `K_p` is compact with the origin in its interior,

\[
T_\epsilon\Subset K_p^\circ=\mathcal D_p.                    \tag{5.3}
\]

The large-box exhaustion theorem supplies a finite `S_epsilon` with

\[
T_\epsilon\subset C_{S_\epsilon}.                            \tag{5.4}
\]

Its support function is

\[
h_{T_\epsilon}(x)
=(1-\epsilon)h_{K_p}(x)
=(1-\epsilon)\tau_p(x).                                      \tag{5.5}
\]

Therefore, for every sufficiently large honest torus,

\[
\boxed{
f_\Lambda(p)
\le N C_{p,\epsilon}
\sum_{\lambda\in\Lambda\setminus\{0\}}
\exp[-(1-\epsilon)\tau_p(\lambda)].}                          \tag{5.6}
\]

This is the desired full-period upper bound.  No half-period arm and no rotated microscopic interaction appear.

## 6. Fixed-direction exponential torus

Take a primitive integer direction `u` and a Bezout complement `v` with `det(u,v)=1`, and

\[
\Lambda_{n,m}=\langle nu,mv\rangle,
\qquad
\frac{\log m}{n}\to d.                                      \tag{6.1}
\]

Write

\[
\xi_u(p)=\tau_p(u).                                           \tag{6.2}
\]

For parallel periods,

\[
\tau_p(\alpha nu)=|\alpha|n\xi_u(p).                         \tag{6.3}
\]

For periods with nonzero transverse coefficient, norm equivalence and

\[
|\det(u,\alpha nu+\beta mv)|=|\beta|m                        \tag{6.4}
\]

give a cost at least `c_p |beta|m/|u|`.  Standard one-dimensional summation in the parallel coordinate then yields

\[
\sum_{\lambda\in\Lambda_{n,m}\setminus0}
 e^{-(1-\epsilon)\tau_p(\lambda)}
\le
C_{p,\epsilon}e^{-(1-\epsilon)n\xi_u(p)}
+e^{-c_{p,\epsilon}m}.                                        \tag{6.5}
\]

Since `N=nm`, (5.6) gives

\[
\limsup_{n\to\infty}\frac1n\log f_{n,m}(p)
\le
-\max\{(1-\epsilon)\xi_u(p)-d,0\}.                           \tag{6.6}
\]

Letting `epsilon downarrow 0`,

\[
\boxed{
\limsup\frac1n\log f_{n,m}(p)
\le-\max\{\xi_u(p)-d,0\}.}                                  \tag{6.7}
\]

Together with the finite-seed seam-closing lower bound already used axially, this identifies the exact fixed-direction winding exponent.

## 7. General homological free-energy interpretation

Equation (5.6) gives a rigorous version of the previously conjectural energy-versus-opportunity competition.

Define the period theta activity

\[
\Theta_{\Lambda,p}^{(\epsilon)}
=\sum_{\lambda\ne0}
 e^{-(1-\epsilon)\tau_p(\lambda)}.                            \tag{7.1}
\]

Then

\[
\boxed{f_\Lambda(p)\le C_{p,\epsilon}N\Theta_{\Lambda,p}^{(\epsilon)}.}\tag{7.2}
\]

If one period family dominates the theta sum and has `M_Lambda` effective translates/opportunities, the transition is governed at exponential scale by

\[
\tau_p(\lambda_*)-\log M_\Lambda.                             \tag{7.3}
\]

The exact entropy factor depends on the geometry/packing of that family; equation (7.2) supplies the universal upper side without guessing it.

For exponentially separated fixed-direction tori, `M_Lambda` is `m` up to subexponential factors and (7.3) becomes `xi_u(p)-log m/n`.

## 8. Why the direct global two-arc comparison was not used

A quotient simple cycle does split into two internally disjoint arcs, but the union over possible quotient arcs is not automatically distributed like two independent full-plane connection events: different plane lifts of the same quotient site share one Bernoulli variable.

The local first-exit proof avoids that issue completely.  Every factor is supported in one fixed injective translate, and BK is applied to actual disjoint quotient witness variables.  This is the safe route from planar directional mass to a finite torus upper bound.

## 9. Claim boundary

Sections 2--5 are finite product-measure arguments once the first-exit certificate and injectivity are supplied.  The passage to `(1-epsilon)K_p` uses the author-level domain equality/exhaustion theorem in `vector-first-exit-domain-20260914.md`, which in turn uses the standard uniform directional exponential estimate for the subcritical SITE norm.  The fixed-direction lower bound remains the finite-seed/Harris/packing construction, not an OZ prefactor.