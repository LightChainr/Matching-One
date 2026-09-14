# Safe site frontiers, no-zero type-B noncrossing partitions, and the dilute periodic TL zero-string module

2026-09-14.  Literature/combinatorics bridge.  This note does not claim that the Bernoulli safe transfer has already been algebraically conjugated to Jacobsen's dilute O(N) transfer matrix.  It identifies a strikingly exact common connectivity module and states the remaining intertwining problem.

## 1. Three descriptions with the same state count

The transparent site safe transfer has, for widths

\[
w=1,2,3,4,5,6,7,8,\ldots,                                   \tag{1.1}
\]

state counts

\[
1,3,7,19,51,141,393,1107,\ldots                              \tag{1.2}
\]

(after the trivial width-one convention is aligned).

These are central trinomial coefficients

\[
T_w=[z^0](1+z+z^{-1})^w,                                     \tag{1.3}
\]

with generating function

\[
\boxed{
\sum_{w\ge0}T_w x^w
=\frac1{\sqrt{(1+x)(1-3x)}}.}                                \tag{1.4}
\]

Independently, Jacobsen's dilute periodic Temperley--Lieb treatment defines the zero-string transfer module `T^(0)` and gives exactly

\[
\boxed{
f_0(x)=\frac1{\sqrt{(1+x)(1-3x)}}
=\sum_{n\ge0}a_nx^n,
\qquad a_n=\dim T^{(0)}.}                                    \tag{1.5}
\]

He explicitly records

\[
\dim T^{(0)}=3\text{ at }n=2,
\qquad
\dim T^{(0)}=7\text{ at }n=3,                                \tag{1.6}
\]

and the asymptotic

\[
a_n\sim\frac12\sqrt{\frac3{\pi n}}3^n.                     \tag{1.7}
\]

This is identical to the safe-frontier count and asymptotic.

A third description is supplied by type-B noncrossing partitions with no zero block.  For a frontier occupancy mask with `r` occupied runs, their count is

\[
|NC_B^{nozero}(r)|
=\binom{2r-1}{r-1},                                          \tag{1.8}
\]

which is exactly the safe multiplicity for that mask.  Summing over cyclic masks gives (1.3).

Thus the evidence points to one combinatorial object seen in three languages:

\[
\boxed{
\text{safe site frontier}
\ \leftrightarrow\ 
NC_B^{nozero}
\ \leftrightarrow\ 
\text{dilute periodic TL zero-string connectivity}.}         \tag{1.9}
\]

## 2. Why “dilute” is the right TL adjective

The ordinary periodic Potts/TL `s=0` reduced sectors have

\[
\frac12\binom{2n}{n}\sim4^n                                 \tag{2.1}
\]

open states and the same number of closed states.  Every loop strand position is present.

The safe SITE frontier, however, has explicit vacant positions in the current row.  Its state growth is `~3^w/sqrt(w)`, not `~4^w/sqrt(w)`.  This is exactly the combinatorics of a **dilute** connectivity basis in which a strand position may be empty.

Jacobsen's dilute O(N) reduced states use the symbol `circle` for an empty site; at `n=2`, the zero-string basis is

```text
circle circle,  (),  )(
```

of dimension three.  This is the same structural trichotomy as the width-two safe SITE frontier: empty, contractible pairing, and across-cut pairing with no winding cycle.

The equality of generating functions therefore has a direct state-semantic explanation, not merely a matching integer sequence.

## 3. Topological meaning of the zero-string condition

In the dilute periodic loop module, `s=0` means there is no string propagating between the two time slices.  Noncontractible loops are controlled separately by the winding-loop weight.

In the safe SITE transfer, a state is retained precisely while no occupied component has acquired nonzero horizontal homology.  In the universal-cover/type-B language this is the no-zero-block condition.

The detailed dictionary should therefore identify

- occupied frontier runs with dilute boundary objects;
- site-component connectivity/deck gains with dilute annular pairings;
- a deck-invariant/zero block with a noncontractible loop/string configuration excluded from the safe module.

The exact local conventions differ from Jacobsen's O(N) model; the claim here is about the connectivity module, not equality of Boltzmann weights.

## 4. Why this may explain the remarkable small state space for square SITE

Jacobsen notes that for some problems, specifically including square-site percolation, the actual number of states needed can be even smaller than the generic dense Potts count `1/2 binom(2n,n)`.

The Bernoulli SITE formulation exposes a natural reason: the appropriate topological frontier includes genuine vacancies and therefore lives on a dilute annular connectivity module whose dimension is central-trinomial rather than central-binomial.

This provides a conceptual explanation for why a direct site/homology transfer can reproduce the same semi-infinite root sequence with only

\[
1,3,7,19,51,\ldots                                           \tag{4.1}
\]

safe states rather than a generic `~4^w` dense TL basis.

## 5. Remaining proof: basis bijection

A full theorem should construct an explicit bijection

\[
\Phi_w:\mathcal S^{safe}_w\to\mathcal B^{dTL}_{w,s=0}        \tag{5.1}
\]

which preserves annular connectivity.

The most promising route factors it through the type-B description:

\[
\mathcal S^{safe}_w
\to NC_B^{nozero}
\to\mathcal B^{dTL}_{w,s=0}.                                 \tag{5.2}
\]

For a fixed occupied mask, contract its runs and use the signed universal-cover partition described in `safe-frontier-central-trinomial-20260914.md`.  Standard dilute-TL diagrams already supply an annular noncrossing representation with vacancies.  The missing step is to align their conventions for the cut/deck gain and prove surjectivity at the discrete square-site level.

## 6. Stronger remaining proof: generator intertwining

Basis equivalence alone does not imply the two transfer matrices are the same.  The deeper target is an intertwiner for local updates.

Let

\[
K^{site}_{4,w}(p),\qquad K^{site}_{8,w}(p)                    \tag{6.1}
\]

be the two substochastic safe Bernoulli kernels.  A dilute-TL representation has local vacancy/connectivity operators.  Seek maps of the form

\[
\Phi K^{site}_{G,w}(p)\Phi^{-1}
=\mathcal T^{dTL}_G(\rho_1(p),\ldots,\rho_9(p))               \tag{6.2}
\]

for some nonintegrable local weights `rho_i(p)`, possibly after a simple diagonal gauge transformation.

If (6.2) holds, several observations become structural rather than empirical:

- central-trinomial state count;
- annular/topological sector semantics;
- Perron magnetic scaling;
- descendant relaxation spectrum;
- the fact that the same connectivity vocabulary supports both NN and matching kernels.

There is no reason to assume the resulting `rho_i(p)` lie on Jacobsen's integrable dilute O(N) manifold; the point is algebraic representation, not exact solvability.

## 7. A possible representation-theoretic payoff

Jacobsen observes in the integrable dilute O(N) model that, for a special noncontractible-loop weight, the spectrum of the one-string sector can embed massively into the zero-string sector.  He explicitly suggests a representation-theoretic explanation.

Our SITE problem exhibits a different but related phenomenon: two microscopic graphs (NN and matching) act on what appears to be the same dilute annular connectivity module, while digital Alexander duality exchanges their topological endpoint sectors.

This suggests that the rapid charge-root convergence may be profitably studied at the module/intertwiner level rather than only through CFT asymptotics.  In particular, the sector-even cancellations responsible for the eigenvalue method could have an algebraic shadow in the common annular module.

This is a research direction, not a conclusion.

## 8. Literature boundary

Jacobsen, arXiv:1507.03027, equations (72)--(73), gives the dilute periodic TL dimensions and generating functions.  His ordinary Potts open/closed reduced sectors instead have `1/2 binom(2n,n)` states each.  The type-B no-zero-block enumeration is classical combinatorics.

The identification of the square-site safe frontier with the dilute zero-string connectivity module, and especially the proposed local-generator intertwining, are new bridge conjectures here.

## 9. Claim boundary

The equality of dimension sequences/generating functions is exact.  The state semantics line up strongly and the width-two/three bases match the expected dilute pattern.  A formal all-width bijection and a transfer-generator conjugacy have not yet been proved.
