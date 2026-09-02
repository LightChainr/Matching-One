# P1 theorem scope

Fill-in date: 2026-09-02. Each statement is only as strong as the proof. Issue and pull-request numbers are remarks, not hypotheses.

## Standing hypotheses

**Embedded graph.** Finite undirected graph honestly embedded in an oriented torus. The random site set is the set of **vertices**. Terminals `L` and `R` that appear after the cut are the two occupied boundary components of that cut; they are deterministic given the cut, not random. Switchable vertices are the vacant sites of the rank-one checkpoint. The residual network after contraction is planar. Cofaciality of `L` and `R` is **not** assumed. There is no L–R edge in the residual network (the cut boundaries are distinct precisely because the occupied graph is rank one). Abstract crossed matching graphs, coincident nearest-neighbour darts or loops, and unsupported degenerate quotients are excluded.

This embedding hypothesis is necessary for the cut construction: without a genuine embedding there is no ambient `H1` image and no essential cycle whose transverse intersection detects the second generator. Removing it leaves the residual connectivity statement unproved; the implementation already rejects those excluded geometries.

**Rank-one checkpoint.** A finite occupied induced subgraph `A` whose ambient first homology image has rank one. Occupied means those vertices are already present and remain present. An essential cycle is a simple occupied cycle whose winding is a primitive nonzero ambient class. Rank one is necessary: rank zero supplies no essential cycle to cut; rank two is already the absorbing event. Admissible cuts are simple occupied essential cycles of `A`.

## Cut construction

Cut the torus along a simple occupied essential cycle `γ`. This produces a cylinder with two occupied boundary copies. Contract every old occupied connected component of the cut graph to a permanent vertex; preserve every vacant site as a distinct switchable vertex together with its incident residual edges. Let `L` and `R` be the two occupied boundary components.

**Already-proved equivalence.** For every future vacant subset `U`,

```text
ambient rank of A ∪ U equals two
    if and only if
L is connected to R in the residual vertex network with precisely U switched on.
```

The residual object is a two-terminal **vertex** network, not an independent-edge reliability model. Fixed-cardinality sampling without replacement is the continuation law. (Remark: this is the cut-network theorem recorded in `notes/p487-cut-network-theorem.md`.)

Removing the contraction step would keep occupied geometry that is already decided and is not part of the future random set. Removing vertex-versus-edge randomness would change the sampling measure; the independent-bond gadget no-go already shows that edge-reliability identities do not import.

## Network state

The continuation state retained after the cut is the **rooted typed residual vertex network**, with the initial cut held fixed as a gauge. It is not a law on future occupations (that law is what the network computes) and it is not a named quotient of the network. Alternate admissible cuts of the same checkpoint induce the same future connecting-subset family but need not be isomorphic as rooted graphs; see the cut-dependence verdict below. The paper treats the network as a proof device and as an update-closed representation in a fixed gauge, not as a physical angular mark.

## Future-rank equivalence

Two residual networks, possibly obtained from different admissible occupied cuts of the same checkpoint, induce the same future rank-two event family if and only if they have the same family of connecting vacant subsets. In the embedded rank-one scope this common family is exactly `{ U : r(A ∪ U) = 2 }`, independently of the cut. Pair-trigger edges and genuine-minimal-triple supports are the length-two and length-three minimal members of that family, hence are likewise cut-independent as event sets.

This does **not** identify residual networks up to rooted graph isomorphism.

## Pair-trigger corollary

Inside the stated embedded nearest-neighbour rank-one scope, pair triggers are bipartite. The two sides are the residual sites that touch exactly one of `L` or `R`; the bipartition is cut geometry, not a colouring of an already measured pair graph. Pair edges are the union of direct L/R vacancy edges and the complete bipartite sets `S_L(C) × S_R(C)` for each neutral old occupied component `C`. Overlapping bicliques are allowed; Ferrers / nested-neighbourhood structure is not implied.

Genuine minimal triples survive as length-three switchable terminal paths with a middle vertex that touches neither `L` nor `R` and with no pair subtrigger.

Removing bipartiteness of the construction (for example by allowing a residual L–R edge) would collapse the rank-one hypothesis: an L–R edge in the residual network is a singleton trigger and is already contracted out of the safe pair graph.

## Update rule

Keep the initial cut fixed. Occupying one more switchable vertex `v` contracts `v` with every adjacent permanent component. If `L` and `R` merge, the state is the rank-two cemetery. Otherwise there remains exactly one variable for each still-vacant original site. Uniform choice among those variables reproduces the original permutation continuation law. Shared-prefix and fork experiments are therefore well-defined on the network.

The representation is update-closed **when the cut is held fixed**. It is not claimed that recutting independently at the successor yields the same isomorphism type. Removing the fixed-cut gauge would require a covariant recut rule that has not been supplied.

## Theorem (survival law does not close branching)

There exist two finite rank-one states on the `4×4` square torus, with eight future vacancies each, that have the same complete unbranched safe-subset polynomial

```text
S(z) = 1 + 7z + 18z² + 20z³ + 8z⁴
```

and therefore the same complete unbranched survival law, but distinct delayed-fork probabilities

```text
B(A) = 95/196,    B(B) = 93/196,    gap = 1/98.
```

The two second moments of the successor-exit counts are `q_A = 29` and `q_B = 25`. Ordinary unbranched continuation is identical; the delayed fork is not. Consequently no function of the complete unbranched survival law determines the delayed-fork observable on any class containing these two states.

(Remark: this is the N16 nucleus of the survival-law no-go, independently certified before the present search. It is imported, not re-proved.)

Removing the fork — keeping only unbranched survival — makes the two states indistinguishable, which is the point of the theorem.

## Theorem (parallel-gadget predictive-class lower bound)

Let `A` and `B` be the two gadgets of the previous theorem, viewed as two-terminal cut-network gadgets. For every integer `k ≥ 1`, place `k` future-vertex-disjoint copies in parallel between the same deterministic terminals, and let `a` of them be type `A`. Parallel safety means every gadget is safe, so the full fixed-cardinality safe-subset polynomial is

```text
S_k(z) = S(z)^k,
```

independent of `a`. The delayed-fork probability is

```text
F_{k,a} = (343 k³ − 182 k² + 25 k + 4a) / (8k (8k−1)²),
```

hence `F_{k,a+1} − F_{k,a} = 1 / (2k (8k−1)²) > 0`. Therefore one complete-survival class splits into at least `k+1` exact branching-predictive classes on `8k` future vertices.

Do not inflate `k`. Do not restate this as a Euclidean dimension bound: one real scalar can encode arbitrarily many discrete labels. The embedded realization is the two-port cylinder-then-glue construction; it does not claim every family member is a nearest-neighbour square-site HNF quotient.

(Remark: this is the parallel-gadget theorem of open PR #549. It is imported, not enlarged.)

Removing disjointness of gadget interiors would allow cross-gadget paths and would invalidate `S_k(z) = S(z)^k`.

## Theorem (r=1 bounded summary is not sufficient)

There exist two connected plane two-terminal vertex-networks `G_A`, `G_B`, each with 7 switchable vertices and with terminals on a common face, such that the complete safe-subset polynomials agree, the singleton and pair trigger counts agree, and the radius-1 terminal-local rooted neighbourhoods are isomorphic as typed graphs, but

```text
P(E2_c2; G_A) = 937/1050 ≠ 313/350 = P(E2_c2; G_B).
```

Incidence:

- A: `[[0,4],[0,6],[1,2],[1,3],[1,6],[2,3],[2,5],[4,5],['L',5],['R',6]]`
- B: `[[0,2],[0,4],[0,6],[1,2],[1,3],[2,5],[3,6],[4,5],['L',5],['R',6]]`

The common polynomial is `S(z) = 1 + 7z + 21z² + 35z³ + 33z⁴ + 15z⁵ + 2z⁶`. Both graphs have `H2 = b2 = 0`. The delayed-fork experiment `E1_c1` equals 1 on both graphs, so the split is not the successor-second-moment observable of the previous two theorems.

**Separating mechanism.** Each graph has exactly two connecting 4-sets (`S_4 = 33`). Their interiors are disjoint on `G_A` (intersection `{5,6}`, the corridor ports) and share the R-adjacent core vertex on `G_B` (intersection `{0,5,6}`). The enumerator `S(z)` cannot see that intersection pattern; `E2_c2` is the mean of `p_surv(remainder, 2)²` after a uniform ordered 2-prefix, and can.

**Not an r=2 witness.** The radius-2 neighbourhoods of this pair differ. Lengthening both corridors to 2 hops equalizes r=1 and r=2 on this S-class but kills the `E2_c2` gap. All nine contracted hop combinations of the exhaustive n≤5 cores are closed at radius 2.

The tuple `(S(z), n, H2, b2, radius-1 neighbourhood)` is therefore not a sufficient statistic for the frozen depth-2 compositional language, already inside the 7-vertex planar two-terminal category. This is a statement about the tested class, not an all-graphs minimality theorem.

Removing the r=1 neighbourhood from the summary still leaves a split (the same pair); the neighbourhood is included because it is the coarsening under test. Removing planarity would leave the pair outside the cut-network image category.

## Cut dependence

```text
CUTS_LAW_EQUIVALENT_NOT_ISOMORPHIC
```

Alternate admissible occupied essential cycles of the same rank-one checkpoint induce the same future connecting-subset family, the same pair-trigger set, and the same genuine-minimal-triple set. They need not be isomorphic as rooted two-terminal networks. Cut reversal on the N16 mask `12463` exchanges the terminals and preserves events. Two occupied-cycle searches on each saved N425 checkpoint return different cycles and identical pair/triple sets. The network is a proof device in a fixed-cut gauge, not a physical angular mark. Details: `notes/p1-cut-dependence-20260902.md`.

## Embedding

```text
TWO_PORT_EMBEDDING_SUFFICES
```

The r=1 non-compression theorem is a theorem about plane two-terminal vertex-networks, the image category of the cut construction. The already-proved two-port rule — place the gadget in a cylinder strip with deterministic cut-boundary terminals, then glue — embeds any such gadget into the finite embedded-graph / rank-one continuation category. No explicit torus occupation of the 7-vertex pair was constructed. The pair is not claimed to be a nearest-neighbour square-site HNF quotient. Details: `notes/p1-n7-torus-embedding-20260902.md`.

## Explicit nonclaims

- no Euclidean latent-dimension lower bound
- no “one real scalar cannot encode the state”
- no continuum / CFT / LCFT / field-count statement
- not an all-graphs theorem
- not a proof that the cut network is a minimal sufficient statistic
- not a population-level #334 effect
- not an r=2 non-compression witness
- not a claim that every nearest-neighbour HNF torus realises the n=7 split
