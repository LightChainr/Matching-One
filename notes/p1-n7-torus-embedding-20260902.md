# Embedding of the n=7 r=1 pair

Date: 2026-09-02
Paper track: cut-network / predictive-state

Supersedes the 2026-09-02 draft of this note that cited parallel-gadget §6 as
if it had already proved surjectivity of the cut construction onto all plane
two-terminal vertex-networks. That citation was incorrect. The general
direction is proved in `notes/p1-plane-tt-realization-lemma-20260902.md`.

## Verdict

```text
GENERAL_REALIZATION_LEMMA
```

The n=7 pair satisfies the hypotheses of that lemma. A finite genuinely
embedded torus graph and a rank-one occupied checkpoint therefore exist for
each of `G_A` and `G_B`. No named square-HNF occupation is produced or claimed.

## Hypotheses, checked on the pair

A plane two-terminal vertex-network is a finite undirected graph with
switchable internal vertices, deterministic terminals `L, R` on a common
face, no edge `{L, R}`, and vertex-activation semantics (terminal connection
= vertex-connectivity after occupying a subset of switchable vertices).
Fixed-cardinality sampling is a measure; it is not used by the realization
lemma.

Both certified graphs meet the definition:

```text
A incidence: [[0,4],[0,6],[1,2],[1,3],[1,6],[2,3],[2,5],[4,5],['L',5],['R',6]]
B incidence: [[0,2],[0,4],[0,6],[1,2],[1,3],[2,5],[3,6],[4,5],['L',5],['R',6]]
```

- finite simple undirected graphs, seven switchable vertices, connected
  carriers, an `L`–`R` path, no `L`–`R` edge;
- `H ∪ {e_{LR}}` plane, equivalently `L, R` on a common face (search-certified);
- `deg(L) = deg(R) = 1`, so no extra rim subdivision is required beyond the
  lemma's uniform `m = max(3, d_L, d_R) = 3`;
- no loops and no parallel edges.

The lemma therefore supplies, for each of `G_A` and `G_B`, a host `G_N` on an
oriented torus, an occupied simple essential cycle `γ` of length 3, and a
rank-one checkpoint `A_N = V(γ)`, whose residual cut-network is rooted-
isomorphic to that graph.

## Why this is not a citation of parallel-gadget §6

Parallel-gadget §6 places *k* copies of two **already realized** N16
cut-network gadgets in disjoint cylinder strips and glues the rims. Its
host graphs exist because those two gadgets were obtained from actual
`4×4` occupations by the cut-network theorem. It does not construct a
torus from an arbitrary plane two-terminal network, and it does not apply
to `G_A`, `G_B` until those graphs are known to be residual cut-networks.

The cut-network theorem is the opposite direction: torus checkpoint to
residual network. Using it alone cannot produce `G_N` from `G_A`.

The missing arrow is the realization lemma, which is now proved independently.

## What is not claimed

- no explicit occupied mask on a named nearest-neighbour square-site HNF
  torus whose residual network is `G_A` or `G_B`;
- not a claim that the E2_c2 split occurs for every HNF torus;
- not an r=2 witness (the radius-2 neighbourhoods of this pair differ);
- not a Euclidean or continuum embedding statement;
- cut-network minimality remains `UNRESOLVED`.

The object of the r=1 non-sufficiency theorem remains the plane two-terminal
category. The realization lemma places that category inside the honestly
embedded rank-one continuation category, without enlarging the theorem to
square lattices.
