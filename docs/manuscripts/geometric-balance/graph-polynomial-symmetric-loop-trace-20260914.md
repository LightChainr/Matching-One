# The graph-polynomial source gives a simple symmetric torus loop trace

Date: 2026-09-14

Status: exact finite weight identity at the self-dual FK point, followed by a representation conjecture for the periodic Temperley--Lieb trace.  This is a specialization of the Euler-localized rank source and is intended as a concrete comparison target for the known torus Potts eigenvalue-amplitude decomposition.

## 1. Boundary loops from the primal--dual pair

Use the notation of `euler-localization-of-rank-source-20260914.md`.

For a spanning FK subgraph `A` on a cellular torus, let

```text
b(A) = number of boundary components of a regular neighbourhood N(A),
g(A),g*(A) = genera carried by N(A) and its complement.
```

The two Euler equations give

```text
k(A)  = [|V|-|A|+b+2g]/2,
k(A*) = [|A|-|V|+b+2g*]/2.
```

Hence exactly

```text
boxed:
k(A)+k(A*) = b(A)+g(A)+g*(A).                    (1.1)
```

The boundary components are the medial interfaces separating primal and dual regions.

## 2. Topological value of `g+g*` on the torus

There are only three ambient homology ranks.

```text
rank 0 : (g,g*)=(0,1),
rank 1 : (g,g*)=(0,0),
rank 2 : (g,g*)=(1,0).
```

Therefore

```text
boxed:
g+g* = 1_(rank != 1).                                (2.1)
```

Combining with (1.1),

```text
k+k* = b + 1_(rank != 1).                              (2.2)
```

This is the precise topological correction to the naive planar statement that cluster counts are controlled only by medial loops.

## 3. Evaluate the graph-polynomial source at the self-dual point

The Euler-localized sourced FK weight is

```text
W_Q,v,h(A)
 = e^{-h|V|}
   (Qe^h)^{k(A)}
   (e^-h)^{k(A*)}
   (ve^h)^{|A|}.
```

At the exact graph-polynomial source

```text
h_Q=-1/2 log Q
```

and the square-lattice self-dual FK point

```text
v=sqrt(Q),
```

we have

```text
Qe^{h_Q}=sqrt(Q),
e^{-h_Q}=sqrt(Q),
ve^{h_Q}=1,
e^{-h_Q|V|}=Q^{|V|/2}.
```

Thus

```text
W(A)
 = Q^{|V|/2} (sqrt(Q))^{k(A)+k(A*)}.                 (3.1)
```

Use (2.2):

```text
boxed:
W(A)
 = Q^{|V|/2}
   (sqrt(Q))^{b(A)}
   (sqrt(Q))^{1_(rank != 1)}.                         (3.2)
```

The first factor is global and can be removed.

## 4. A remarkably simple modified trace rule

After global normalization, the sourced self-dual model has the configuration weight

```text
boxed:
(sqrt Q)^(number of medial boundary loops)
×
{ sqrt Q,  rank 0 or rank 2,
  1,       rank 1. }
```

Equivalently:

> give every medial boundary loop the ordinary Potts loop fugacity `sqrt(Q)`, and give the zero-through-line extreme homology sectors one additional factor `sqrt(Q)` relative to the rank-one through-line sectors.

This is the simplest explicit candidate encountered so far for the topological modified trace underlying the graph-polynomial balance.

It does **not** separately distinguish rank0 from rank2 at this symmetric source; it need not, because the graph-polynomial source has already transformed the original asymmetric FK weights into a primal--dual symmetric ensemble.

Away from `h_Q`, the full two-face/cluster fugacity asymmetry is required.

## 5. Relation to the periodic TL sectors

In a medial periodic transfer:

- rank-one FK states carry noncontractible interfaces / through-line sectors;
- rank0 and rank2 are both zero-through-line in the simple loop count and require an additional topological distinction in the ordinary FK partition function.

Equation (3.2) says that **at the sourced graph-polynomial symmetric point**, the remaining distinction collapses to a single sector factor `sqrt(Q)` for the extreme/zero-through-line topology relative to rank one.

This is exactly the kind of information supplied by a modified Markov trace rather than by the local TL generator.

A concrete algebra task is therefore:

```text
construct/evaluate the periodic TL trace with
contractible/noncontractible boundary-loop weight sqrt(Q)
and extreme-sector bonus sqrt(Q),
then compare its characters/eigenvalue amplitudes with the known graph-polynomial sectors.
```

## 6. Connection to the Richard--Jacobsen torus amplitude decomposition

For the toroidal Potts model the partition function is not a simple trace because noncontractible clusters carry global information.  Richard and Jacobsen decompose it into transfer characters labelled by the number of noncontractible clusters and cyclic-group representations, with nontrivial eigenvalue amplitudes.

That formalism is a natural existing place to test (3.2): instead of inventing a new transfer algebra, reweight their sector amplitudes by the sourced primal--dual rule and ask whether the graph-polynomial eigenvalue identity emerges directly.

This is now a much narrower literature/algebra comparison than “find a twisted massive Potts partition function”.

## 7. Why the symmetric point is special

For general topological source `h`, write

```text
Q_primal = sqrt(Q) e^eta,
Q_dual   = sqrt(Q) e^-eta,
eta=h+1/2 log Q.
```

At `eta=0` the primal and dual faces have equal fugacity; the loop boundaries do not need to remember which side is favoured.  Only the torus topology correction (2.1) survives.

For `eta!=0`, one must retain the signed difference of primal and dual component counts, which is naturally represented by an oriented-loop/height/face-charge variable.

Thus the graph-polynomial source is algebraically special because it sits exactly at the **zero face-charge** point of the two-sided cluster gas.

This suggests the continuous rank source around `h_Q` may be represented as an electric/background-charge perturbation of a fixed-`sqrt(Q)` loop model.

## 8. A possible continuum implication

If `eta` is indeed an electric/height charge in the periodic loop description, then derivatives in the topological source `h` are not arbitrary closure insertions: they probe a well-defined electric/topological sector of the loop theory.

This could connect the repository's rank-source coordinate directly to the charged/map-resolved Potts continuum fields being considered in #585/#782.

This statement is a programme, not an identification of the corresponding Coulomb-gas charge.

## 9. Site-model boundary

The derivation is for an edge/FK model on a cellular torus.  It should not be quoted as an exact square-site transfer identity without an explicit decorated/incidence-graph construction.

Its immediate role is in the generic-Q Potts/TL continuum bridge and in understanding Jacobsen's graph-polynomial eigenvalue criterion.

## 10. Claim boundary

Exact:

```text
k+k*=b+g+g*,
g+g*=1_(rank!=1) on T^2,
self-dual graph-polynomial sourced weight
  ∝ (sqrt Q)^b (sqrt Q)^{1_(rank!=1)}.
```

Conjectural/interface:

- identification of (3.2) with a particular periodic TL/Markov trace normalization;
- electric/height interpretation of `eta` away from the symmetric point;
- massive/TBA continuation.

The next useful work is a direct comparison to existing torus Potts eigenvalue amplitudes, not a new width computation.