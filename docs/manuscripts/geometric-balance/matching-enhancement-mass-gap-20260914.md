# Matching enhancement and a strict inverse-correlation-length gap

2026-09-14.  Targeted proof reduction motivated by Grimmett--Li's site-enhancement argument.  The conclusion is **conditional on one two-terminal finite-volume adaptation** isolated below; it is not yet promoted to the theorem layer.

## 1. Desired statement

Let `G4` be nearest-neighbour square-site percolation and `G8` its matching graph.  Write

\[
\kappa_4(p)=-\lim_{n\to\infty}\frac1n\log P_p^{G4}(0\leftrightarrow ne_1),
\qquad
\kappa_8(p)=-\lim_{n\to\infty}\frac1n\log P_p^{G8}(0\leftrightarrow ne_1).
\]

Graph inclusion gives only `kappa_8(p)<=kappa_4(p)`.  The target is the strict fixed-parameter statement

\[
\boxed{\kappa_8(p)<\kappa_4(p),\qquad 0<p<p_c(G8).}
\]

A strict mass gap would immediately sharpen the two-centre geometry in the exponential-aspect regime.

## 2. The published enhancement input

Grimmett--Li, *Percolation critical probabilities of matching lattice-pairs*, Random Structures & Algorithms 65 (2024), Section 5, introduce `hat G` by adding a facial site inside each nontriangular face and connecting it to all boundary vertices.  Under `P_{p,s}`, original vertices are independently open with probability `p` and facial sites with probability `s`.  Their equation (5.2) gives

\[
\theta(p,0)=\theta(p;G),
\qquad
\theta(p,1)=\theta(p;G_*).
\]

For the finite event `v_0 <-> partial Lambda_n`, Lemma 5.4 maps every pivotal original vertex, by a bounded local modification, to a nearby pivotal facial site.  Summing and using Russo's formula gives their equation (5.4)

\[
\partial_p\theta_n(p,s)
\le g(p,s)\,\partial_s\theta_n(p,s),
\]

with `g` locally bounded on the interior parameter square.  Their equation (5.5) integrates this comparison along characteristics: on every compact parameter rectangle, a fixed amount of enhancement can compensate a strictly positive decrement of the ordinary site density.

The square lattice is transitive and nontriangular, so the local enhancement is essential; their theorem indeed gives `p_c(G8)<p_c(G4)`.  What we need here is a **two-terminal** version of the pivotal comparison, not merely the critical-point conclusion.

## 3. Two-terminal Directional Enhancement Lemma (DEL) — the single missing adaptation

Let

\[
A_n=\{0\leftrightarrow ne_1\text{ in }\widehat G\}
\]

under `P_{p,s}`.  Facial sites are auxiliary; at `s=0` this is the NN two-point event, while at `s=1` connectivity of original vertices is equivalent to connectivity in the matching graph, up to the harmless two-edge replacement through a facial site.

**DEL.**  For every compact `K subset (0,1)^2`, there are `g_K<infinity`, `M<infinity` and `n_0` such that, for all `n>=n_0` and `(p,s) in K`,

\[
\boxed{
\partial_p P_{p,s}(A_n)
\le g_K\,\partial_s P_{p,s}(A_n).}
\]

A version with a multiplicative/subexponential finite-endpoint error is already sufficient for the mass comparison below.

### Why this is genuinely a local adaptation

The printed Grimmett--Li proof of Lemma 5.4 starts from an open non-self-touching path through an original pivotal vertex.  Outside a bounded neighbourhood of that pivotal vertex it retains the two arms of the path; inside, a five-stage planar surgery inserts a translated copy of the essential enhancement pattern so that a nearby facial site becomes pivotal.  The target `partial Lambda_n` supplies only the far endpoint of one retained arm.

For `A_n`, a pivotal original vertex away from the two terminal `O(M)` neighbourhoods again has exactly two retained arms, now ending at `0` and `ne_1`.  The same bounded surgery is therefore available without any change in the local enhancement gadget.  What must be written explicitly is:

1. the replacement of the terminal `partial Lambda_n` arm by the arm ending at `ne_1`;
2. the finitely many local types when the pivotal site lies within `O(M)` of either terminal;
3. a bounded-multiplicity map from original-pivotal configurations to nearby facial-pivotal configurations.

On the square lattice these endpoint cases are finite local configurations rather than a new large-scale geometric problem.  Nevertheless, they should be written before the strict mass claim is promoted.  A citation of the strict critical-point theorem alone is not enough.

## 4. DEL gives a uniform sprinkling comparison for the two-point function

Fix `p in (0,p_c(G8))`.  Choose a compact parameter rectangle containing a path from an interior facial parameter near `0` to one near `1`, and let

\[
\gamma=1/\sup_K g>0.
\]

Along a characteristic with `dp/ds=-gamma`, DEL gives

\[
\frac d{ds}P_{p(s),s}(A_n)
=\partial_sP-\gamma\partial_pP\ge0.
\]

Integrating a fixed positive amount in `s` and using ordinary monotonicity for the small endpoint pieces `s in [0,eta]` and `[1-eta,1]` yields a number `delta=delta(p)>0`, independent of `n`, such that

\[
\boxed{
P_p^{G8}(0\leftrightarrow ne_1)
\ge
P_{p+\delta}^{G4}(0\leftrightarrow ne_1)
}
\]

for all sufficiently large `n` (or the same inequality up to a factor `e^{o(n)}` if one uses the weaker DEL form).  Choose `delta` small enough that `p+delta<p_c(G4)`.

Taking `-n^{-1}log` and passing to the limit gives directly

\[
\boxed{\kappa_8(p)\le\kappa_4(p+\delta).}
\]

No directional union bound, rectangle endpoint count, or OZ prefactor is involved: the finite event is already the defining two-point event for the axial mass.

The #739 branch derives the strict parameter comparison

\[
\kappa_4(p)-\kappa_4(q)
\ge(q-p)\kappa_4(q)/\rho>0,
\qquad p<q<p_c(G4),
\]

from its Friedgut--Kalai birth argument.  Therefore DEL implies

\[
\boxed{\kappa_8(p)<\kappa_4(p).}
\]

If one prefers not to reuse that branch inequality, a precise site reference for strict monotonicity of the Bernoulli inverse correlation length may be substituted.

## 5. Consequence for the two exponential-aspect centres

Let

\[
a(d)=\kappa_4^{-1}(d),
\qquad
c(d)=\kappa_8^{-1}(d),
\qquad
b(d)=1-c(d).
\]

The branch currently proves only `c(d)<=a(d)`, equivalently `a(d)+b(d)>=1`.

If `a(d)<p_c(G8)`, the strict mass gap gives

\[
\kappa_8(a(d))<\kappa_4(a(d))=d.
\]

Since `kappa_8` is strictly decreasing, its solution of `kappa_8(c)=d` satisfies `c(d)<a(d)`.  If instead `a(d)>=p_c(G8)`, then automatically `c(d)<p_c(G8)<=a(d)`.  Hence DEL would give for **every** `d>0`

\[
\boxed{a(d)+b(d)>1.}
\]

In the dual-even/odd birth coordinates this says that the limiting centre displacement

\[
C_\infty(d)=\frac{a(d)+b(d)-1}{2}
\]

is strictly positive throughout the separated-window regime.  Equivalently the limiting two-atom law is not reflection-centred about `1/2`, even though its width is controlled separately by `b-a`.

This strict centre asymmetry is unrelated to the possible nondifferentiability exceptional set in the Gumbel proof.  The two questions must remain separate.

## 6. Quantitative version worth pursuing

The enhancement proof is local and its characteristic slope is uniform on compact `(p,s)` sets.  A quantitative DEL would produce an explicit `delta_I>0` for every compact `I subset (0,p_c(G8))`, and therefore

\[
\kappa_4(p)-\kappa_8(p)
\ge \kappa_4(p)-\kappa_4(p+\delta_I)
\]

on `I`, with a further explicit lower bound from the branch's mass-slope inequality.  Such a certificate would separate directional centre intervals without estimating an OZ amplitude.

## 7. Decision rule

- If the two-terminal DEL is written successfully: promote `kappa_8(p)<kappa_4(p)` and `a(d)+b(d)>1` to theorem-level consequences.
- If the local pivotal surgery fails specifically for a two-terminal event: record the obstruction; the published strict critical-point theorem does not by itself imply a strict mass gap.
- No numerical width scan is needed to decide this question.

Primary source checked: G. Grimmett and Z. Li, *Percolation critical probabilities of matching lattice-pairs*, Random Structures & Algorithms 65 (2024), 832--856, especially equations (5.2)--(5.5) and Lemma 5.4.
