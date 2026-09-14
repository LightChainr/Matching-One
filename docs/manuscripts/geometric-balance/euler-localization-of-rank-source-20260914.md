# Euler localization of the torus rank source: a symmetric primal--dual cluster gas

Date: 2026-09-14

Status: exact finite embedded-graph identity for a cellular graph on the torus, plus a conjectural loop/TL interpretation.  This sharpens the Krushkal source dictionary and may provide a direct route toward the modified trace sought in #782.

## 1. Setup

Let `G=(V,E)` be a connected cellular graph embedded on `T^2`.  Let `A subset E` be a spanning FK subgraph and let `A*` denote the complementary dual subgraph: a dual edge is present exactly when the corresponding primal edge is absent.

Write

```text
k(A)  = number of connected components of A,
k(A*) = number of connected components of A*,
r(A)  = rank im[H1(A)->H1(T^2)] in {0,1,2}.
```

The claim below is configurationwise and contains no probabilistic/scaling input.

## 2. Euler derivation

Let `N(A)` be a regular neighbourhood of the primal subgraph.  It retracts to `A`, so

```text
chi(N(A)) = |V|-|A|.
```

Let

```text
g  = genus carried by N(A),
g* = genus carried by the closure of T^2\N(A),
b  = number of common boundary components.
```

Then

```text
|V|-|A| = 2 k(A)  - 2g  - b.                     (2.1)
```

The complementary region retracts to the complementary dual ribbon subgraph `A*`.  Since `G` is cellular on the torus,

```text
|V|-|E|+|F|=0,
```

and therefore

```text
chi(T^2\N(A))
 = |F|-|E\A|
 = |A|-|V|
 = 2 k(A*) - 2g* - b.                              (2.2)
```

Equating the two expressions for `b` gives

```text
k(A)-k(A*)+|A|-|V| = g-g*.                         (2.3)
```

For the torus, Krushkal's topological exponents satisfy

```text
s=2g,
s_perp=2g*,
```

and the rank-source identity from the preceding note gives

```text
r(A)-1 = g-g*.
```

Hence the exact configurationwise formula is

```text
boxed:
r(A)-1 = k(A)-k(A*)+|A|-|V|.                      (2.4)
```

Checks:

- empty primal subgraph: `k(A)=|V|`, `k(A*)=1`, `|A|=0` -> `r-1=-1`;
- full primal subgraph: Euler duality gives `r-1=+1`;
- rank-one states give zero charge.

## 3. The topological source becomes a two-sided cluster fugacity

Start from the ordinary FK weight

```text
W_Q,v(A)=Q^{k(A)} v^{|A|}.
```

Insert the bounded rank source:

```text
W_Q,v,h(A)
 = Q^{k(A)} v^{|A|} exp[h(r(A)-1)].
```

Using (2.4),

```text
W_Q,v,h(A)
 = e^{-h|V|}
   (Q e^h)^{k(A)}
   (e^{-h})^{k(A*)}
   (v e^h)^{|A|}.                                  (3.1)
```

Thus, up to the harmless global factor `e^{-h|V|}`, the rank source is exactly a **primal--dual two-cluster-fugacity random-cluster model**:

```text
Q_primal = Q e^h,
Q_dual   = e^-h,
v_source  = v e^h.
```

A global surface-topology source has been converted into ordinary component and edge counts of the primal/dual pair.

This is stronger than merely knowing that the Krushkal polynomial contains the source: it supplies an explicit finite transfer bookkeeping rule.

## 4. Product invariant

The two cluster fugacities obey

```text
boxed:
Q_primal * Q_dual = Q,
```

independent of `h`.

Therefore varying the topological charge redistributes cluster weight between primal and dual sides while keeping their product fixed.

This strongly suggests that in a medial-loop representation the ordinary bulk loop weight `sqrt(Q)` can remain fixed while `h` is represented by an orientation/face/background-charge bias rather than by changing the loop fugacity itself.

That loop statement is a conjectural representation interpretation; equation (3.1) is exact.

## 5. The graph-polynomial source is the symmetric point

The generic-Q graph-polynomial balance source is already known exactly:

```text
h_Q = -1/2 log Q.
```

Substitute it into the two cluster fugacities:

```text
Q_primal = Q e^{h_Q} = sqrt(Q),
Q_dual   = e^{-h_Q}  = sqrt(Q).
```

Therefore

```text
boxed:
the graph-polynomial topological source is exactly the point where
primal and dual cluster fugacities become equal.
```

The sourced edge fugacity becomes

```text
v_source = v/sqrt(Q).
```

At the square-lattice self-dual FK point `v=sqrt(Q)`, this is simply

```text
v_source=1.
```

Thus the sourced graph-polynomial condition has an unexpectedly symmetric form:

```text
primal cluster fugacity = dual cluster fugacity = sqrt(Q),
rescaled edge fugacity = 1.
```

This is a plausible structural explanation for why the periodic Temperley--Lieb eigenvalue criterion is so natural at the graph-polynomial root.

It is not yet a derivation of Jacobsen's complete transfer formula; seam/through-line and closure normalizations still have to be matched.

## 6. A candidate route to the missing modified trace

Equation (3.1) suggests replacing the abstract question

```text
"how do we insert e^{h(r-1)} into a periodic TL trace?"
```

by the more concrete problem

```text
"how do we assign distinct closure fugacities to primal and dual clusters,
while keeping Q_primal Q_dual=Q?"
```

A connectivity transfer can do this locally in time: whenever a primal or complementary dual component closes, multiply by its declared cluster fugacity.  The only additional care is for components that survive around the periodic direction / noncontractible sectors.

The medial-loop formulation should encode the same asymmetry through the alternating primal/dual faces separated by the loops.  Because the product fugacity is fixed, an oriented-loop/height background charge is a natural candidate language.

This narrows the #782 representation problem substantially:

```text
finite state-sum projector      : solved by (2.4)/(3.1),
local primal-dual transfer       : apparently constructible,
periodic TL/Markov trace closure : still to derive,
massive/TBA continuation         : still open.
```

## 7. Rank-one neutral source remains independent

The charge source `h` changes the relative weight of rank0 and rank2 while leaving rank1 with unit topological factor.

The existing affine-TL seam `alpha` instead acts inside rank one:

```text
z=alpha^2/Q,
```

weighting the number `K` of parallel essential clusters.

Thus the two finite source directions are genuinely complementary:

```text
h : primal--dual cluster fugacity imbalance / topological charge,
alpha : neutral essential-count fugacity.
```

A combined periodic transfer should carry both.

## 8. Relation to the h-odd sector and Matching One

Differentiate (3.1) at `h=0`:

```text
r-1 = k(A)-k(A*)+|A|-|V|.
```

So the matching-odd topological observable is an **Euler imbalance** between primal clusters, dual clusters and occupied edges.

This offers another microscopic interpretation of the first noncommon correction:

> the relevant continuum field must survive after the extensive/local primal--dual Euler pieces cancel down to the bounded torus topology charge.

It also explains why low-order ordinary thermal/geometric corrections can be large in individual sectors yet disappear from the rank-source response.

## 9. Possible relation to the thermal spin-four mechanism

Under the sourced representation, a microscopic anisotropy can affect

```text
primal cluster closure statistics,
dual cluster closure statistics,
edge density.
```

The matching-root differential response is their specific Euler-signed combination.

A sector-even identity-family anisotropy can contribute almost equally to the primal and dual cluster gas and cancel.  The first anisotropic mismatch of the thermal/pivotal structure can then survive as the observed spin-four response.

This is qualitative at present but gives a more concrete cluster-language target than an abstract `x=21/4` label.

## 10. Relation to original graph-polynomial algebra

The previous exact result

```text
Z_2-Q Z_0=0
```

was reinterpreted as a topological charge source `h_Q=-1/2 log Q`.

The new Euler localization shows that this same source is the **primal-dual symmetric cluster-gas point**.  Hence two seemingly different explanations of the graph-polynomial criterion coincide:

```text
topological sector balance
<=>
primal/dual cluster-fugacity symmetry after Euler localization.
```

This equivalence may be useful for proving which transfer eigenvalue sectors are compared by the periodic TL criterion.

## 11. Site-model boundary

Everything above is exact for an edge/FK spanning-subgraph model on a cellular torus.

Square-site Matching One has an analogous digital-Alexander rank charge but does not literally have the edge-complement dual subgraph used in (2.4) without passing to an incidence/decorated cell representation.  Therefore the two-cluster-fugacity localization is presently a Potts/FK continuum/interface tool, not an asserted exact square-site identity.

The distinction is important for #275 and for any direct numerical source implementation.

## 12. Claim boundary

Exact:

```text
r-1=k(A)-k(A*)+|A|-|V|,
W_Q,v,h=e^{-h|V|}(Qe^h)^k(e^{-h})^{k*}(ve^h)^{|A|},
Q_primal Q_dual=Q,
h_Q=-1/2 log Q -> Q_primal=Q_dual=sqrt(Q).
```

Conjectural/interface:

- an oriented-loop/height realization with fixed bulk loop fugacity `sqrt(Q)`;
- direct identification with the needed periodic affine-TL modified trace;
- massive/TBA continuation;
- application as an exact site-percolation transfer identity.

This turns the topological source from a closure-only label into an explicit primal--dual cluster-gas deformation.