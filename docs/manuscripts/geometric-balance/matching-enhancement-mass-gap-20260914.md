# Matching enhancement gives a strict inverse-correlation-length gap

2026-09-14.  Author-level derivation for the actual square-site NN/matching pair.  This note combines Grimmett--Li's facial-site interpolation with the corrected two-dimensional enhancement rerouting theorem of Balister--Bollobas--Riordan.  It should still receive independent proof review before publication, but the previous `Directional Enhancement Lemma` is closed here rather than left as a conjectural interface.

## 1. Statement

Let `G4` be nearest-neighbour square-site percolation and `G8` its matching graph.  For subcritical parameters write

\[
\kappa_4(p)=-\lim_{n\to\infty}\frac1n\log P_p^{G4}(0\leftrightarrow ne_1),
\qquad
\kappa_8(p)=-\lim_{n\to\infty}\frac1n\log P_p^{G8}(0\leftrightarrow ne_1).
\]

**Theorem (strict matching mass gap).**  For every

\[
0<p<p_c(G8),
\]

one has

\[
\boxed{\kappa_8(p)<\kappa_4(p).}
\]

Consequently, for the exponential-aspect centre functions

\[
a(d)=\kappa_4^{-1}(d),\qquad
c(d)=\kappa_8^{-1}(d),\qquad
b(d)=1-c(d),
\]

one has for every `d>0`

\[
\boxed{c(d)<a(d),\qquad a(d)+b(d)>1.}
\]

This strict centre asymmetry is separate from the possible nondifferentiability exceptional set in the median-centred Gumbel theorem.

## 2. Facial-site interpolation

Following Grimmett--Li, add one facial site to every square face and join it to the four corner vertices; call the resulting bipartite augmentation `hat G`.  Under `P_{p,s}`:

- original square-lattice vertices are independently open with probability `p`;
- facial sites are independently open with probability `s`.

For connectivity between original vertices:

\[
s=0\quad\Longleftrightarrow\quad G4,
\]

and

\[
s=1\quad\Longleftrightarrow\quad G8.
\]

The second equivalence is exact: a passage through an always-open facial site connects any two open corners of one square, which is precisely the clique/matching connectivity of that face; conversely every matching diagonal is represented by the two-edge path through its facial site.

Grimmett--Li use the same interpolation in their equation (5.2).  For the radial finite event `v_0 <-> partial Lambda_n`, their Lemma 5.4 constructs a bounded local map from an original pivotal vertex to a nearby pivotal facial site, yielding by Russo's formula

\[
\partial_p\theta_n(p,s)
\le g(p,s)\partial_s\theta_n(p,s).                              \tag{2.1}
\]

The issue is whether the same uniform comparison holds for the two-point event defining `kappa`.  On the square lattice, the corrected enhancement proof of Balister--Bollobas--Riordan supplies exactly the path-rerouting input needed for this adaptation.

## 3. The two-dimensional rerouting theorem that closes the local geometry

Balister--Bollobas--Riordan, *Essential enhancements revisited* (2014), isolate the flawed general Aizenman--Grimmett path lemma and replace it with a statement they prove for site percolation on `Z^2`.

In their notation, for every fixed inner size `m` there exists a fixed `r=r(m,2)` such that two induced red/green paths that enter `B_{r+2}` from outside, end near the origin, and have no red--green adjacency can be modified **only inside `B_r`** so that they reach prescribed opposite inner terminals while remaining mutually nonadjacent.  Equivalently, an induced path through the origin with two remote endpoints can be rerouted inside a bounded ball while keeping the same remote endpoints and replacing the central segment by a prescribed clean connector.  Their Section 5.1 proves this statement for the square lattice.

They then show how such a bounded path surgery produces a finite-to-one pivotal map: after clearing irrelevant local sites, a pivotal connection consists locally of exactly the two arms of a shortest/induced path.  The rerouting inserts the enhancement gadget between those arms without creating an alternate red--green connection.  Endpoint-near pivotal sites are moved a bounded distance into the bulk by another bounded modification.  All constants are independent of the overall connection length.

We use this **proved square-lattice rerouting statement**, not the false unrestricted lemma identified in that paper.

## 4. Two-terminal pivotal comparison

For integer `n` let

\[
A_n=\{0\leftrightarrow ne_1\text{ in }\widehat G\}.
\]

Call an original vertex `p`-pivotal if flipping that original vertex changes `A_n`, and a facial site `s`-pivotal if flipping the facial site changes `A_n`.

### Lemma 4.1 (bounded pivotal conversion)

Fix a compact parameter rectangle `K subset (0,1)^2`.  There are constants `R<infinity` and `C_K<infinity`, independent of `n`, such that for every original vertex `z`,

\[
P_{p,s}(z\text{ is p-pivotal for }A_n)
\le C_K
P_{p,s}(\text{some facial site in }B_R(z)
          \text{ is s-pivotal for }A_n).                       \tag{4.1}
\]

**Proof.**  Take a configuration in which `z` is original-pivotal.

1. **Deactivate nearby facial sites.**  Close active facial sites in a fixed large ball around `z` one at a time.  If the event is destroyed at one step, the last facial site is already facial-pivotal and we are done.  Otherwise `z` remains original-pivotal with all local facial sites closed.

2. **Remove irrelevant local open sites.**  With `z` open, choose a shortest open path `P` from `0` to `ne_1`; it contains `z`, since `z` is pivotal.  Delete local original sites not needed for a connection, one at a time.  We arrive at a configuration whose intersection with a fixed neighbourhood of `z` is the corresponding portion of an induced path.  Removing `z` splits that path into two mutually nonadjacent arms, one leading to each terminal.  Extra connections outside the modification ball do not defeat pivotality: the usual equivalence-relation/minimality argument of the enhancement proof shows that any such bypass would make one of the retained local path pieces unnecessary before the modification.

3. **Reroute the two arms.**  If `z` is farther than the fixed rerouting radius from both terminals, apply the proved `Z^2` red/green rerouting theorem of Balister--Bollobas--Riordan inside a translate of `B_r(z)`.  It brings the two arms to two prescribed separated inner terminals while retaining their remote endpoints and without a red--green adjacency.

4. **Insert the facial-site gadget.**  Continue the two clean arms, still inside the bounded modification region, to **opposite corners of one square face** and close every other local original/facial site that could connect the two sides.  Let `phi` be the center of that face.  With `phi` closed there is no connection between the two arms; with `phi` open the two opposite corners connect through `phi`.  Hence `phi` is facial-pivotal for the same two-terminal event.

5. **Pivots near the terminals.**  If `z` lies within the fixed rerouting radius of `0` or `ne_1`, use the finite endpoint manoeuvre from the enhancement proof: open a fixed small neighbourhood of the terminal, close all but one exit on a surrounding finite sphere, and thereby move the pivotal bottleneck to a site a fixed distance into the connection path.  The preceding construction then applies.  The target endpoint case is the translate/reflection of the source case.

Every modification changes only a bounded number of original and facial states.  On a compact `K`, the ratio between probabilities of any two such local patterns is uniformly bounded.  The map has bounded multiplicity because its image records a pivotal facial site within a fixed ball and there are only finitely many local patterns.  This proves (4.1).  `square`

### Corollary 4.2 (two-terminal differential inequality)

Summing (4.1) over `z`, using bounded overlap of the balls `B_R(z)`, and applying Russo's formula separately to the original and facial coordinates gives a locally bounded `g_K` such that

\[
\boxed{
\partial_p P_{p,s}(A_n)
\le g_K\partial_sP_{p,s}(A_n)                                  \tag{4.2}
}
\]

for every sufficiently large `n`, uniformly on `K`.

This is the two-terminal version of Grimmett--Li (5.4).  The proof above is the only additional model-specific step.

## 5. A fixed amount of enhancement beats a positive p-sprinkling

Fix `p in (0,p_c(G8))`.  Choose `eta in (0,1/4)` and a compact `p` interval contained in `(0,p_c(G4))` and containing `p`.  Let

\[
g_* = \sup_K g_K<\infty,
\qquad \gamma=1/g_*>0.
\]

Choose a small `L>0` so that

\[
\delta=\gamma L>0,
\qquad p+\delta<p_c(G4),
\qquad \eta+L<1-\eta.
\]

Along

\[
s\in[\eta,\eta+L],
\qquad p(s)=p+\delta-\gamma(s-\eta),
\]

Corollary 4.2 gives

\[
\frac d{ds}P_{p(s),s}(A_n)
=\partial_sP-\gamma\partial_pP\ge0.
\]

Hence

\[
P_{p,\eta+L}(A_n)\ge P_{p+\delta,\eta}(A_n).
\]

Monotonicity in `s` then yields

\[
\boxed{
P_{p,1}(A_n)\ge P_{p+\delta,0}(A_n).                            \tag{5.1}
}
\]

By the endpoint identifications of Section 2,

\[
\boxed{
P_p^{G8}(0\leftrightarrow ne_1)
\ge
P_{p+\delta}^{G4}(0\leftrightarrow ne_1).                      \tag{5.2}
}
\]

The crucial point is that `delta>0` is independent of `n`.

## 6. Strict mass gap

Take `-n^{-1}log` in (5.2) and send `n->infinity`:

\[
\kappa_8(p)\le\kappa_4(p+\delta).                              \tag{6.1}
\]

The #739 branch already proves, for `p<q<p_c(G4)`, the quantitative strict parameter monotonicity

\[
\kappa_4(p)-\kappa_4(q)
\ge(q-p)\kappa_4(q)/\rho>0.                                   \tag{6.2}
\]

Applying (6.2) with `q=p+delta` gives

\[
\kappa_8(p)
\le\kappa_4(p+\delta)
<\kappa_4(p),
\]

which proves the theorem.

Notice that no OZ prefactor, p-analyticity, directional differentiability, or numerical estimate of either mass is required.

## 7. Strict separation of the two exponential-aspect centres

Let

\[
a(d)=\kappa_4^{-1}(d),
\qquad
c(d)=\kappa_8^{-1}(d),
\qquad
b(d)=1-c(d).
\]

If `a(d)<p_c(G8)`, then

\[
\kappa_8(a(d))<\kappa_4(a(d))=d.
\]

Since `kappa_8` is strictly decreasing, the solution `c(d)` of `kappa_8(c)=d` satisfies

\[
c(d)<a(d).
\]

If instead `a(d)>=p_c(G8)`, then automatically

\[
c(d)<p_c(G8)\le a(d).
\]

Thus for every `d>0`,

\[
\boxed{c(d)<a(d),}
\qquad
\boxed{a(d)+b(d)>1.}                                          \tag{7.1}
\]

In the dual-even/odd birth coordinates,

\[
C_\infty(d)=\frac{a(d)+b(d)-1}{2}>0.                           \tag{7.2}
\]

So throughout the separated-window exponential-aspect regime the limiting two-atom law has a strictly positive complement-odd centre displacement.  This statement is independent of whether the mass is differentiable at the inverse images used by the Gumbel scaling theorem.

## 8. Quantitative compact-interval version

The pivotal map changes a bounded number of coordinates.  On a compact `I subset (0,p_c(G8))` and a fixed interior facial interval, its finite-energy constant is uniform.  Therefore the construction gives one `delta_I>0` such that

\[
\kappa_8(p)\le\kappa_4(p+\delta_I),\qquad p\in I,
\]

provided `delta_I` is reduced to keep `p+delta_I<p_c(G4)`.  Combining with (6.2),

\[
\boxed{
\kappa_4(p)-\kappa_8(p)
\ge \delta_I\,\kappa_4(p+\delta_I)/\rho>0,
\qquad p\in I.}                                                \tag{8.1}
\]

This is a qualitative rigorous certificate; the constants inherited from the local surgery are extremely conservative and are not proposed as a numerical mass bound.

## 9. Claim boundary and references

This is an author-level proof assembled from published enhancement machinery plus the #739 mass-monotonicity inequality.  It has not received an independent line-by-line audit.  In particular a reviewer should check the finite-to-one bookkeeping in Lemma 4.1 and the endpoint pivotal relocation, rather than merely accepting the analogy with the radial event.

Primary inputs actually read:

- G. Grimmett and Z. Li, *Percolation critical probabilities of matching lattice-pairs*, Random Structures & Algorithms 65 (2024), 832--856: facial-site interpolation, Lemma 5.4, equations (5.2)--(5.5).
- P. Balister, B. Bollobas, O. Riordan, *Essential enhancements revisited*, arXiv:1402.0834: the corrected enhancement argument, especially Conjecture 7 and its proof for `d=2` in Section 5.1, and the finite-to-one pivotal-map reduction preceding it.

The older unrestricted Aizenman--Grimmett combinatorial lemma is **not** used without the Balister--Bollobas--Riordan correction.
