# Matching enhancement and a strict inverse-correlation-length gap

2026-09-14.  Targeted proof reduction motivated by Grimmett--Li's site-enhancement argument.  The conclusion is **conditional on one directional finite-volume adaptation** isolated below; it is not yet promoted to the theorem layer.

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

with `g` locally finite/continuous on the interior parameter square.  Their equation (5.5) integrates this comparison along characteristics: on every compact parameter rectangle, a fixed amount of enhancement can compensate a strictly positive decrement of the ordinary site density.

The square lattice is transitive and nontriangular, so the local enhancement is essential; their theorem indeed gives `p_c(G8)<p_c(G4)`.  What we need here is a finite-direction version below, not merely the critical-point conclusion.

## 3. Directional Enhancement Lemma (DEL) — the single missing adaptation

For `C>1`, let `R_n=[0,n] x [-Cn,Cn]` and let `A_n` be the event that there is an open connection in `hat G` from the left side to the right side of `R_n`, with a fixed `O(1)` thickening if needed to absorb facial sites.

**DEL.**  For every compact `K subset (0,1)^2`, there are `g_K<infinity`, a boundary-thickening constant and `n_0` such that, for all large `n` and `(p,s) in K`,

\[
\partial_p P_{p,s}(A_n)
\le g_K\,\partial_s P_{p,s}(A_n)
+e^{-c_K n},
\]

or, more than sufficiently, the same inequality with an error whose logarithmic effect on `P(A_n)` is `o(n)`.

Why this looks close to the printed proof rather than a new enhancement theorem: Grimmett--Li's Lemma 5.4 is local.  Their five-stage surgery takes an open path through an original pivotal vertex, modifies only a bounded neighbourhood, and inserts a translated copy of the essential enhancement pattern so a nearby facial site becomes pivotal.  In the interior of a long rectangle the two arms from the pivotal vertex terminate on the two crossing sides instead of `v_0` and `partial Lambda_n`; the bounded surgery is otherwise of the same type.  What must still be written carefully is the treatment of pivotal sites within `O(1)` of the two terminal sides and the exact finite-domain convention.

This lemma should be proved directly before the strict mass claim is promoted.  A citation of the critical-point theorem alone is not enough.

## 4. DEL implies a finite-volume sprinkling comparison

Fix `p` in a compact subinterval of `(0,p_c(G8))`.  Work first with facial parameter `s` in `[eta,1-eta]`.  By DEL, take `gamma_K=1/g_K>0` (shrinking it if necessary).  Along a segment with

\[
\frac{dp}{ds}=-\gamma_K,
\]

we have, up to the negligible DEL boundary error,

\[
\frac d{ds}P_{p(s),s}(A_n)
=\partial_sP-\gamma_K\partial_pP\ge0.
\]

Integrating across a fixed positive `s` interval and using monotonicity to move from the interior cutoffs to `s=0,1` gives a number `delta=delta(p)>0`, independent of `n`, such that

\[
\boxed{
P_{p}^{G8}(A_n)
\ge P_{p+\delta}^{G4}(A_n)\,e^{-o(n)}.}
\]

Choose `delta` small enough that `p+delta<p_c(G4)`.  The exact numerical value is irrelevant for strictness.

## 5. The rectangle-crossing exponent is the horizontal mass

For either square-symmetric finite-range graph `G`, let `tau_{G,p}` be its subcritical connection norm.  Reflection symmetry makes `y -> tau(1,y)` even and convex, hence

\[
\tau(1,y)\ge\tau(1,0)=\kappa_G(p).
\]

The left-right crossing probability of `R_n` has logarithmic rate `kappa_G(p)`:

- a point-to-point connection from `(0,0)` to `(n,0)` is a left-right crossing, giving the lower probability bound;
- a crossing contains a pair of boundary vertices whose horizontal displacement is `n+O(1)` and vertical displacement `O(n)`; every such pair costs at least `n kappa_G(p)+O(1)` by the norm inequality above, while the number of endpoint pairs is polynomial in `n`.

Thus

\[
-\frac1n\log P_p^G(A_n)\to\kappa_G(p).
\]

This step uses only the already standard existence/norm property of subcritical inverse correlation length and a polynomial endpoint union bound, not an OZ prefactor.

Taking logarithmic rates in the sprinkling comparison yields

\[
\kappa_8(p)\le\kappa_4(p+\delta).
\]

The #739 branch already derives strict decrease of the NN mass in its occupation parameter at differentiability points, and more quantitatively the comparison

\[
\kappa_4(p)-\kappa_4(q)
\ge(q-p)\kappa_4(q)/\rho>0,
\qquad p<q<p_c(G4),
\]

from the Friedgut--Kalai birth argument.  Therefore

\[
\boxed{\kappa_8(p)<\kappa_4(p).}
\]

If one prefers not to reuse that branch inequality, strict monotonicity of the Bernoulli inverse correlation length in `p` on a fixed finite-range graph is a separate standard input that may be substituted after its precise site reference is fixed.

## 6. Consequence for the two exponential-aspect centres

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

## 7. Quantitative version worth pursuing

The enhancement proof is local and its characteristic slope is uniform on compact `(p,s)` sets.  A quantitative DEL would produce an explicit `delta_I>0` for every compact `I subset (0,p_c(G8))`, and therefore

\[
\kappa_4(p)-\kappa_8(p)
\ge \kappa_4(p)-\kappa_4(p+\delta_I)
\]

on `I`, with a further explicit lower bound from the branch's mass-slope inequality.  Such a certificate would be useful for separating directional centre intervals without estimating an OZ amplitude.

## 8. Decision rule

- If DEL is proved: promote `kappa_8(p)<kappa_4(p)` and `a(d)+b(d)>1` to theorem-level consequences.
- If the local pivotal surgery fails specifically for directional crossing: record the obstruction; the published strict critical-point theorem does not by itself imply a strict mass gap.
- No numerical width scan is needed to decide this question.

Primary source checked: G. Grimmett and Z. Li, *Percolation critical probabilities of matching lattice-pairs*, Random Structures & Algorithms 65 (2024), 832--856, especially equations (5.2)--(5.5) and Lemma 5.4.
