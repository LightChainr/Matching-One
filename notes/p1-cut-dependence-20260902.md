# Cut dependence of the rank-one residual network

Date: 2026-09-02
Paper track: cut-network / predictive-state
Depends on: already-proved embedded rank-one cut theorem (open PR #491, `notes/p487-cut-network-theorem.md`).

## Verdict

```text
CUTS_LAW_EQUIVALENT_NOT_ISOMORPHIC
```

## Scope

The working category is the one already closed by the cut-network theorem: a genuinely embedded undirected graph on an oriented torus, a rank-one occupied checkpoint, a simple occupied essential cycle, and the residual planar **vertex** network obtained by cutting along that cycle and contracting old occupied components. Random variables are vertices. Sampling is fixed-cardinality without replacement. Independent-edge reliability is excluded.

This note does not reopen that theorem, does not produce a new predictive-class lower bound, and does not treat a change of cut presentation as a new branching witness.

## Precise statement

Let `A` be an occupied rank-one checkpoint and let `γ`, `γ'` be two admissible simple occupied essential cycles. Write `N(A,γ)` and `N(A,γ')` for the rooted two-terminal vertex networks obtained by cutting along each cycle and contracting old occupied components, with terminals `L,R` the two occupied boundary components.

**Event equivalence.** For every future vacant subset `U`,

```text
r(A ∪ U) = 2  iff  L connects to R in N(A,γ) after switching U on
              iff  L connects to R in N(A,γ') after switching U on.
```

Both residual networks therefore induce the **same** future rank-two event family. The pair-trigger graph and the genuine-minimal-triple set are likewise independent of the occupied carrier: they are the length-two and length-three minimal supports of those events. This is the already-proved cut theorem plus the already-recorded cut-independence of the event graph, not a new identity.

**Not a canonical isomorphism.** The rooted graphs `N(A,γ)` and `N(A,γ')` need not be isomorphic as rooted two-terminal networks. Changing the essential carrier can change which old occupied components appear as the terminals versus as neutral blocks, and can change the component-wise biclique presentation. Cut reversal exchanges the two terminals. The fixed cut is a convenient gauge for the update rule; recutting independently at a successor is not claimed to return the same graph-isomorphism representative.

**Equivalence relation used by the paper.** Declare

```text
N ~_law N'  iff  they induce the same family of connecting vacant subsets.
```

Then `N(A,γ) ~_law N(A,γ')` for every pair of admissible occupied essential cycles. The paper should treat a concrete cut-network as a **proof device and an update-closed representation in a fixed gauge**, not as a physical angular mark. A component-incidence or overlap statistic that is to be read as a physical mark still needs an explicit covariant cut rule or a demonstrated cut-invariant definition. That is a presentation constraint, not a failure of the continuation law.

This is not `CUTS_NEED_QUOTIENT`: the unmarked continuation theorem does not require a unique isomorphism type. Any admissible cut already computes the same events. A quotient would be needed only if one insisted on a single rooted graph as an observable.

This is not `CUTS_CANONICALLY_ISOMORPHIC`: reversal and alternate-carrier presentations differ as rooted graphs.

This is not `CUTS_COUNTEREXAMPLE`: no admissible occupied cut of a rank-one checkpoint has been observed, or proved, to change the pair/triple event sets.

## Finite example 1 — cut reversal on the N16 witness

On the `4×4` square torus, take the rank-one occupied mask `12463`. Let `N` be the residual network of the search-selected essential cycle, and let `N_rev` be the residual network of the reversed cycle (same vertices, reversed orientation, negated winding). Then:

- the left-site set of `N` equals the right-site set of `N_rev`, and conversely;
- the neutral-site sets coincide;
- the minimal pair sets coincide;
- the genuine-minimal-triple sets coincide.

The reversed network is therefore anti-isomorphic as a rooted two-terminal graph (terminals exchanged) and is **not** the same rooted isomorphism type, while every continuation event is preserved. This is the existing focused test `test_reversing_cut_exchanges_sides` of PR #491; it is restated here as the cut-dependence example, not re-proved.

## Finite example 2 — two occupied-cycle searches on each saved N425 checkpoint

For each of the two hash-pinned N425 rank-one checkpoints, the default occupied-cycle search and the reversed occupied-cycle search return different simple essential cycles (checkpoint A: original cycle of 29 occupied vertices; alternative cycle of length 31) and **identical** pair and triple event sets. The machine certificate records

```text
alternative_cut_pair_and_triple_sets_equal = true
```

on both checkpoints, with 108 pair edges and 583 / 509 genuine triples unchanged. The event graph is therefore cut-independent on these real graphs; the component presentation (which old blocks are L/R versus neutral, how bicliques are written) is allowed to change with the carrier.

## Paper use

- **Proof device.** Any one admissible cut is enough to prove future rank-two iff residual L–R connection, pair-trigger bipartiteness, genuine triples as length-three switchable terminal paths, and update-closure under a **held-fixed** cut.
- **Not a physical mark.** Do not promote a cut-dependent biclique, shared-vertex, or component-size list to an angular observable without a covariant rule.
- **Relation to the r=1 non-compression lemma.** That lemma lives entirely inside the two-terminal category after a cut has been chosen. Cut non-uniqueness does not weaken it: the n=7 pair is a pair of plane two-terminal vertex-networks, not a pair of torus occupations.

## Explicit nonclaims

- no new predictive-class lower bound from a change of cut;
- no claim that recutting at every successor is isomorphism-invariant;
- no claim that the cut network is a minimal sufficient statistic;
- no Euclidean, continuum, CFT, or field-count statement.
