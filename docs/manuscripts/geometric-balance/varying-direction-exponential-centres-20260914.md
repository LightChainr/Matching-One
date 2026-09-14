# Varying shortest-period directions: full exponential-aspect centre theorem

2026-09-14.  This note completes the geometric direction step posed in #765 for the genuine exponential-aspect regime.  The short-period direction may vary with the torus and converge only along a subsequence.

The proof combines:

- the full-period first-exit upper bound in `first-exit-torus-winding-upper-20260914.md`;
- a fixed local directional seed that approximates the limiting direction, repeated many times and closed by a small deterministic correction;
- the determinant geometry of a rank-two integer period lattice.

No microscopic rotation, directional OZ prefactor, or angle scan is used.

## 1. Geometry and statement

Let `Lambda_n` be honest rank-two sublattices of `Z^2`.  Let

\[
N_n=[\mathbb Z^2:\Lambda_n]                                  \tag{1.1}
\]

and choose a shortest nonzero period

\[
u_n\in\Lambda_n,
\qquad
\ell_n=|u_n|\to\infty.                                       \tag{1.2}
\]

A shortest period is primitive in `Lambda_n`, although it need not be primitive in ambient `Z^2`.

Assume along the sequence or a chosen subsequence

\[
e_n:=u_n/\ell_n\to e\in S^1,                                \tag{1.3}
\]

and define the transverse height

\[
h_n=N_n/\ell_n.                                               \tag{1.4}
\]

The exponential-aspect hypothesis is

\[
\boxed{
\frac{\log h_n}{\ell_n}\to d\in(0,\infty).}                 \tag{1.5}
\]

For `G=G4` or `G8`, let

\[
f_{G,n}(p)=P_p^G(r_G>0).                                     \tag{1.6}
\]

Let `tau_{G,p}` be the subcritical planar inverse-correlation norm.

**Theorem.**  For every fixed `p<pc(G)`,

\[
\boxed{
\lim_{n\to\infty}\frac1{\ell_n}\log f_{G,n}(p)
=-\max\{\tau_{G,p}(e)-d,0\}.}                                \tag{1.7}
\]

Consequently `p -> tau_{G,p}(e)` is continuous and strictly decreasing from `infinity` to zero on `(0,pc(G))`, and the two NN homology births satisfy

\[
T_{1,n}\xrightarrow P a_e(d),
\qquad
\tau_{4,a_e(d)}(e)=d,                                        \tag{1.8}
\]

\[
T_{2,n}\xrightarrow P b_e(d),
\qquad
b_e(d)=1-c_e(d),
\qquad
\tau_{8,c_e(d)}(e)=d.                                        \tag{1.9}
\]

If only a subsequence of directions converges to `e`, the conclusions hold on that subsequence.  No whole-sequence direction is invented.

## 2. Exponential transverse height automatically removes nonparallel period competition

Complete `u_n` to a lattice basis `(u_n,v_n)` with

\[
|\det(u_n,v_n)|=N_n.                                         \tag{2.1}
\]

Every period is

\[
\lambda=a u_n+b v_n,
\qquad a,b\in\mathbb Z.                                      \tag{2.2}
\]

If `b!=0`, then

\[
|\det(u_n,\lambda)|=|b|N_n.                                 \tag{2.3}
\]

The perpendicular component of `lambda` relative to `u_n` therefore has magnitude

\[
\frac{|\det(u_n,\lambda)|}{|u_n|}
=|b|h_n.                                                       \tag{2.4}
\]

Hence

\[
\boxed{|\lambda|\ge |b|h_n\qquad(b\ne0).}                    \tag{2.5}
\]

At any fixed subcritical parameter, norm equivalence gives

\[
\tau_{G,p}(\lambda)\ge c_{G,p}|b|h_n.                        \tag{2.6}
\]

Under (1.5), `h_n=exp(d ell_n+o(ell_n))`.  Thus every nonparallel homology class has super-`ell_n` cost and is negligible compared with the `O(ell_n)` short-period class.

This also proves a purely geometric uniqueness statement.  If two independent periods both had Euclidean length `O(ell_n)`, their determinant would be `O(ell_n^2)`.  Since the determinant of two independent lattice periods is a nonzero integer multiple of `N_n`, this would force

\[
N_n=O(\ell_n^2),
\qquad
h_n=O(\ell_n),                                                \tag{2.7}
\]

contradicting (1.5).  Therefore exponential elongation itself forbids two nonparallel `O(ell_n)` period classes.

The multi-direction hard-core crossover from `projective-poisson-hardcore-crossover-20260914.md` is relevant to geometries where several direction costs genuinely remain comparable, but **not** to the present shortest-period exponential-aspect regime.

## 3. Upper bound from the full-period theta estimate

Fix `p<pc(G)` and `epsilon in (0,1)`.

The first-exit torus theorem gives

\[
f_{G,n}(p)
\le
N_n C_{p,\epsilon}
\sum_{\lambda\in\Lambda_n\setminus0}
\exp[-(1-\epsilon)\tau_{G,p}(\lambda)].                      \tag{3.1}
\]

The parallel periods are `a u_n`.  Homogeneity gives

\[
\tau_{G,p}(a u_n)=|a|\ell_n\tau_{G,p}(e_n).                  \tag{3.2}
\]

By continuity of the norm in direction,

\[
\tau_{G,p}(e_n)\to\tau_{G,p}(e).                             \tag{3.3}
\]

The nonparallel terms are bounded by (2.6) and their lattice sum is exponentially negligible on the `ell_n` scale.  Hence

\[
\sum_{\lambda\ne0}
 e^{-(1-\epsilon)\tau_{G,p}(\lambda)}
\le
C_{p,\epsilon}
 e^{-(1-\epsilon)\ell_n\tau_{G,p}(e_n)}
+e^{-c_{p,\epsilon}h_n}.                                      \tag{3.4}
\]

Since

\[
\frac{\log N_n}{\ell_n}
=\frac{\log h_n}{\ell_n}+rac{\log\ell_n}{\ell_n}
\to d,                                                        \tag{3.5}
\]

we obtain

\[
\limsup\frac1{\ell_n}\log f_{G,n}(p)
\le
-\max\{(1-\epsilon)\tau_{G,p}(e)-d,0\}.                      \tag{3.6}
\]

Let `epsilon downarrow 0`:

\[
\boxed{
\limsup\frac1{\ell_n}\log f_{G,n}(p)
\le-\max\{\tau_{G,p}(e)-d,0\}.}                              \tag{3.7}
\]

## 4. Lower bound with one fixed approximating seed direction

The short period `u_n` changes with `n`, so one should not demand a new large finite seed for every arithmetic direction.  Instead choose the seed direction **before** the torus limit.

Fix `eta>0`.  By directional continuity of the norm, choose one nonzero integer vector `w` such that

\[
\left|\frac{w}{|w|}-e\right|<\eta                             \tag{4.1}
\]

and

\[
\frac{\tau_{G,p}(w)}{|w|}
<\tau_{G,p}(e)+\eta.                                          \tag{4.2}
\]

By the definition of the mass along the fixed integer direction `w`, choose a fixed multiple `z=Lw` and a fixed finite box seed from `0` to `z` whose conditional connection probability `q` obeys

\[
-\log q
<|z|[\tau_{G,p}(e)+2\eta].                               \tag{4.3}
\]

All of `w,L,z`, and the seed box are fixed before `n->infinity`.

For large `n`, `e_n` is also within `eta` of `e`.  Choose

\[
k_n=\operatorname{round}\left(\frac{u_n\cdot z}{|z|^2}\right).\tag{4.4}
\]

Then the remainder

\[
r_n=u_n-k_nz                                                   \tag{4.5}
\]

satisfies

\[
|r_n|\le C\eta\ell_n+O(|z|).                                 \tag{4.6}
\]

Repeat the fixed seed `k_n` times and join the last seed endpoint to `u_n` by a deterministic NN path of length at most

\[
|r_n|_1\le\sqrt2|r_n|.                                       \tag{4.7}
\]

NN edges are available in both `G4` and `G8`.

Harris positive association applies even if translated seed boxes overlap after projection.  Therefore the prescribed closed ring has probability

\[
P_p(\text{one }u_n\text{-ring})
\ge
q^{k_n}p^{|r_n|_1+O(1)}.                                     \tag{4.8}
\]

Equations (4.3)--(4.7) imply

\[
\boxed{
P_p(\text{one }u_n\text{-ring})
\ge
\exp[-(\tau_{G,p}(e)+C_p\eta)\ell_n]}                        \tag{4.9}
\]

for all large `n`.

The complete prescribed support lies in a Euclidean region of diameter `O(ell_n)`, so its difference set contains only

\[
|S_n-S_n|=O(\ell_n^2)                                        \tag{4.10}
\]

quotient vertices.  The finite-group packing lemma therefore supplies at least

\[
\frac{N_n}{C\ell_n^2}
=\frac{h_n}{C\ell_n}                                         \tag{4.11}
\]

pairwise site-disjoint translated attempts.

Hence

\[
1-f_{G,n}(p)
\le
\exp\left[
-\frac{h_n}{C\ell_n}
 e^{-(\tau_{G,p}(e)+C_p\eta)\ell_n}
\right].                                                      \tag{4.12}
\]

If `tau(e)<d`, the right side tends to zero.  If `tau(e)>d`, use `1-e^{-x}>=x/2` for small `x` to obtain

\[
\liminf\frac1{\ell_n}\log f_{G,n}(p)
\ge d-\tau_{G,p}(e)-C_p\eta.                                 \tag{4.13}
\]

Let `eta downarrow 0`.  Together with (3.7), this proves (1.7).

## 5. Strict p-monotonicity for every Euclidean direction

The rate theorem now supplies the missing input for a direction-free Friedgut--Kalai argument.

Fix any `e in S^1`.  Choose primitive ambient integer approximants `u_n` with

\[
|u_n|\to\infty,
\qquad
u_n/|u_n|\to e.                                               \tag{5.1}
\]

Choose Bezout complements `v_n` and a transverse multiplier `m_n` so that the torus

\[
\Lambda_n=\langle u_n,m_nv_n\rangle                           \tag{5.2}
\]

has

\[
\frac{\log(N_n/|u_n|)}{|u_n|}\to d.                          \tag{5.3}
\]

Suppose for contradiction that for some

\[
0<p<q<pc(G)
\]

we had

\[
\tau_{G,p}(e)=\tau_{G,q}(e)=A>0.                              \tag{5.4}
\]

Choose `d<A` close enough to `A`, and `zeta>0`, that the Friedgut--Kalai sharp-threshold increment for a transitive event is strictly less than `q-p`:

\[
\rho\frac{A-d+\zeta}{d}<q-p.                                 \tag{5.5}
\]

By (1.7), at parameter `p`,

\[
f_{G,n}(p)>e^{-(A-d+\zeta)\ell_n}                            \tag{5.6}
\]

for large `n`.  The event `r_G>0` is increasing and translation-transitive on the `N_n` sites, so Friedgut--Kalai forces

\[
f_{G,n}(q)\to1.                                               \tag{5.7}
\]

But the same rate theorem at `q`, with `A>d`, gives

\[
f_{G,n}(q)\to0,                                               \tag{5.8}
\]

a contradiction.

Therefore

\[
\boxed{
\tau_{G,p}(e)\text{ is strictly decreasing in }p
\text{ for every }e\in S^1.}                                 \tag{5.9}
\]

The continuity in `p` follows from the same site finite-seed/cluster-volume comparison used axially, with the standard uniform directional norm estimate supplying the passage from rational approximants to `e`.  Small-`p` path counting gives divergence, and

\[
\tau_{G,p}(e)\le\kappa_G(p)(|e_x|+|e_y|)                     \tag{5.10}
\]

forces the limit zero as `p up to pc(G)`.

Thus every directional mass has a unique inverse value for each `d>0`.

## 6. The two birth centres

For NN,

\[
P(T_{1,n}\le p)=f_{4,n}(p).                                  \tag{6.1}
\]

The rate theorem and strict invertibility give

\[
\boxed{
T_{1,n}\xrightarrow P a_e(d),
\qquad
\tau_{4,a_e(d)}(e)=d.}                                       \tag{6.2}
\]

Digital Alexander duality gives

\[
P_p^{G4}(r=2)=P_{1-p}^{G8}(r=0),                              \tag{6.3}
\]

so

\[
\boxed{
T_{2,n}\xrightarrow P b_e(d),
\qquad
b_e(d)=1-c_e(d),
\qquad
\tau_{8,c_e(d)}(e)=d.}                                       \tag{6.4}
\]

This proves the centre formula anticipated in #765 for arbitrary convergent shortest-period directions.

## 7. Strict matching directional mass gap on the whole circle

The matching-enhancement pivotal conversion is endpoint-direction independent.  On every compact `I subset (0,pc(G8))` it gives a positive `delta_I` such that

\[
P_p^{G8}(0\leftrightarrow x)
\ge
P_{p+\delta_I}^{G4}(0\leftrightarrow x)                      \tag{7.1}
\]

uniformly over distant endpoint directions and `p in I`.

Taking directional rates,

\[
\tau_{8,p}(e)
\le\tau_{4,p+\delta_I}(e).                                   \tag{7.2}
\]

By the all-direction strict monotonicity (5.9),

\[
\tau_{4,p+\delta_I}(e)<\tau_{4,p}(e).                         \tag{7.3}
\]

Therefore

\[
\boxed{
\tau_{8,p}(e)<\tau_{4,p}(e)
\quad\text{for every }e\in S^1,
\ 0<p<pc(G8).}                                                \tag{7.4}
\]

So strict matching/NN centre asymmetry holds in every limiting direction, not only on the coordinate axes.

## 8. Relation to the homological free-energy picture

The upper bound is the rigorous energy side:

\[
f_\Lambda(p)
\le
C N\sum_{\lambda\ne0}e^{-(1-\epsilon)\tau_p(\lambda)}.       \tag{8.1}
\]

The lower construction supplies the opportunity side for the distinguished shortest period.  Under exponential transverse height, there are `h_n/poly(ell_n)` disjoint translated attempts.  Thus the effective exponent is

\[
\tau_p(e)-d.                                                   \tag{8.2}
\]

In more complicated geometries where several period families remain on the same exponential scale, the theta upper bound remains valid while the lower-side opportunity count must be analysed family by family.  That is the correct setting for the broader homological-free-energy variational conjecture.

## 9. Claim boundary

The new torus upper bound is proved through finite first-exit SITE certificates and BK.  Its exhaustion to the exact support function uses the standard uniform directional exponential estimate already isolated in `vector-first-exit-domain-20260914.md`.

The lower bound uses only one fixed finite seed direction chosen after `p,e,eta` and before the torus limit, plus deterministic connectors, Harris association and finite-group packing.

No claim is made here about a common near-critical multi-direction window when the exponential transverse-height separation (1.5) fails.