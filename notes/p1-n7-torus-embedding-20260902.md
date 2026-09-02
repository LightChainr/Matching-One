# Embedding of the n=7 r=1 non-compression pair

Date: 2026-09-02
Paper track: cut-network / predictive-state
Depends on: cut-network theorem (open PR #491) and the two-port embedded-realization lemma (open PR #549, §6 of `notes/p429-parallel-gadget-lower-bound-20260901.md`).

## Verdict

```text
TWO_PORT_EMBEDDING_SUFFICES
```

## Object of the theorem

The r=1 non-compression statement is a statement about **connected plane two-terminal vertex-networks** with terminals on a common face and no L–R edge. That is exactly the image category of the rank-one cut construction: cut along a simple occupied essential cycle, contract old occupied components, retain vacant sites as switchable vertices. Future ambient rank reaches two iff the residual terminals connect.

The n=7 pair is certified inside that category:

```text
A incidence: [[0,4],[0,6],[1,2],[1,3],[1,6],[2,3],[2,5],[4,5],['L',5],['R',6]]
B incidence: [[0,2],[0,4],[0,6],[1,2],[1,3],[2,5],[3,6],[4,5],['L',5],['R',6]]
```

Both graphs are plane with `G ∪ {L,R}` planar, connected carriers, unique L-neighbour, unique R-neighbour. No explicit torus occupation realizing either 7-vertex network was constructed in the bounded search, and none is constructed here.

## Why two-port embedding is enough

The cut-network theorem already maps every honestly embedded rank-one checkpoint to a planar two-terminal vertex network. The parallel-gadget lower bound then records an embedded-realization lemma that applies to **any** planar two-terminal gadget in that category, not only to the N16 base pair:

1. place the gadget in a strip of a cylinder, with the two deterministic cut-boundary components as the terminals;
2. gadget interiors remain disjoint from any other strip;
3. glue the two cylinder boundaries back together.

The result is an embedded torus graph whose initial occupied subgraph has ambient rank one, and whose future rank reaches two exactly when the two-terminal gadget connects. The lemma does **not** assert that the gadget is a nearest-neighbour square-site HNF quotient, nor that every member of the family is one of the repository's Gaussian HNFs.

That already-proved rule places `G_A` and `G_B` as two-terminal blocks in the same finite embedded-graph / rank-one continuation category as the cut-network theorem. It is therefore unnecessary, for the r=1 non-compression lemma, to exhibit an explicit occupied mask on a named HNF torus whose residual network is isomorphic to `G_A` or `G_B`.

The paper's theorem object remains:

```text
plane two-terminal vertex-networks
```

and not

```text
every nearest-neighbour square-site HNF torus occupation.
```

## What an explicit occupation would have to be, and why it is not claimed

An `EXPLICIT_TORUS_OCCUPATION` certificate would require a finite torus graph, an occupied rank-one mask, and an admissible occupied essential cycle such that the contracted residual network is isomorphic, as a rooted two-terminal vertex-network, to `G_A` (resp. `G_B`). The present pair has seven switchable vertices, each of degree at least two, terminals of degree one in the gadget, and a triangle on the A side. Nothing in the search constructed such a mask. Declaring the pair to be a square-site HNF quotient without that mask would be a category error.

The two-port calculus is the correct embedding statement: the gadgets may be inserted as blocks between deterministic cut boundaries. That is the same rule already used to realize the unbounded parallel-gadget family.

## Explicit nonclaims

- no explicit torus occupation of the n=7 pair;
- not a nearest-neighbour square-site HNF quotient unless independently constructed;
- not a claim that the split occurs for every HNF torus;
- not an r=2 witness (the radius-2 neighbourhoods of this pair differ);
- not a Euclidean or continuum embedding statement.
