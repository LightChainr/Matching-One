# Realization of plane two-terminal vertex-networks as rank-one torus cut-networks

Date: 2026-09-02
Paper track: cut-network / predictive-state
Claim level: C5 for the combinatorial embedding statement below.

This note proves a **new** general realization lemma. It is not a corollary of the
cut-network theorem, and it is not supplied by the parallel-gadget note.

- The cut-network theorem (remark: `notes/p487-cut-network-theorem.md`) is the
  direction *torus rank-one checkpoint → residual two-terminal vertex-network*.
- The parallel-gadget note §6 (remark: `notes/p429-parallel-gadget-lower-bound-20260901.md`)
  constructs an embedded torus for *k parallel copies of two already realized
  N16 gadgets*. It does not prove that every plane two-terminal vertex-network
  arises as a residual cut-network.

Hypotheses of the lemma do not use issue or pull-request numbers.

---

## Definition — plane two-terminal vertex-network

A **plane two-terminal vertex-network** is a triple `N = (H, L, R)` where:

1. `H` is a finite undirected graph (loops and parallel edges permitted).
2. `V(H) = V_sw ∪ {L, R}` with `L ≠ R`. The vertices in `V_sw` are
   **switchable**; `L` and `R` are **deterministic terminals**.
3. `H` is connected, and there is at least one `L`–`R` path.
4. There is **no** edge `{L, R}`. Loops at `L` or at `R`, if present, do not
   affect terminal connection and are discarded before the construction.
5. `H` is plane: it is embedded in the 2-sphere (equivalently, in a closed
   disk) so that `L` and `R` lie on a **common face** `F`. Taking `F` as the
   unbounded face, this is the same as: `H ∪ {e_{LR}}` remains plane, where
   `e_{LR}` is a new edge joining `L` to `R` in `F`.
6. **Activation semantics.** A future configuration is a subset `U ⊆ V_sw`.
   Terminals remain present. The occupied vertex set is `U ∪ {L, R}`. An
   activation **connects the terminals** when `L` and `R` lie in the same
   connected component of the subgraph induced by `U ∪ {L, R}`. Connectivity
   is vertex-connectivity in `H`, not independent-edge reliability.

Loops at switchable vertices, if present, do not affect terminal connection
and are retained. Parallel edges, if present, are permitted in `H`. The
cut-network contraction deletes loops created at the contracted terminals
and duplicate residual edges, so the rooted network obtained after the cut
is the simple graph underlying `N` (switchable loops retained).

**Sampling is irrelevant.** Fixed-cardinality sampling without replacement is
a probability measure on subsets of `V_sw`. The realization lemma is a
statement about graphs, embeddings, and induced residual networks. It does
not mention a measure.

Write `d_L = deg_H(L)` and `d_R = deg_H(R)`, counting parallel edges with
multiplicity and ignoring discarded terminal loops. The path hypothesis and
the absence of `{L, R}` force `d_L ≥ 1` and `d_R ≥ 1`. Every neighbour of
`L` or `R` lies in `V_sw`.

---

## Lemma — every such network is a residual cut-network

**Lemma.** Let `N = (H, L, R)` be a plane two-terminal vertex-network. There
exist

- a finite undirected graph `G_N` honestly embedded in an oriented torus `T`,
- a vertex subset `A_N ⊂ V(G_N)`,
- a simple occupied cycle `γ` of the induced subgraph `G_N[A_N]`,

such that:

1. `A_N` is a rank-one checkpoint: the ambient first-homology image of
   `G_N[A_N]` in `H_1(T; Z)` has rank one, and `γ` is a primitive essential
   representative;
2. cutting `G_N` along `γ`, contracting every old occupied connected
   component, and retaining every vertex of `V_sw` as a distinct switchable
   site produces a rooted two-terminal vertex-network isomorphic to the
   simple graph underlying `N`.

In particular, future terminal connection in `N` coincides with residual
`L`–`R` connection after the cut, and therefore with ambient rank of
`A_N ∪ U` reaching two.

The host `G_N` is a genuinely embedded torus graph. It is **not** claimed to
be a nearest-neighbour square-site graph, nor a column-HNF quotient, nor one
of the repository's Gaussian tori.

---

## Proof

### 1. Strip drawing from the common face

Take a plane embedding of `H` on the 2-sphere and a common face `F` of `L`
and `R`. Place a new Jordan arc `e_{LR}` in the closure of `F`, meeting `H`
only at `L` and `R`. The graph `H ∪ {e_{LR}}` is still plane.

The two faces of `H ∪ {e_{LR}}` incident to `e_{LR}` have boundaries
`W^+ ∪ e_{LR}` and `W^− ∪ e_{LR}`. Deleting `e_{LR}` yields two walks
`W^+`, `W^−` of `H` from `L` to `R`. These walks are simple paths if `H` is
2-connected; they may repeat vertices if `L` or `R` is a cut-vertex. This is
the only use of cofaciality.

Stereographic projection from an interior point of `F` realises `H` in the
plane with `F` unbounded. There is then an orientation-preserving
homeomorphism of the plane sending `L` to `(1/2, 0)`, `R` to `(1/2, 1)`, and
a neighbourhood of infinity (still in `F`) so that the open vertical sides
of the rectangle `Q = [0,1]_x × [0,1]_y` lie entirely in `F`. Restricting to
`Q` yields a drawing such that:

- `L` lies on the bottom side `(0,1) × {0}` as a single vertex;
- `R` lies on the top side `(0,1) × {1}` as a single vertex;
- every switchable vertex lies in the open rectangle `(0,1) × (0,1)`;
- every edge is a simple Jordan arc in `Q`, interiors pairwise disjoint
  except at shared endpoints;
- the open vertical sides `{0} × (0,1)` and `{1} × (0,1)` contain no vertex
  and no interior of an edge;
- vertices of `W^+` other than `L, R` sit in the open strip near `x = 0^+`,
  vertices of `W^−` other than `L, R` sit near `x = 1^−`, and the two walks
  meet `∂Q` only at `L` and `R`.

A cut-vertex belonging to both walks is a single interior vertex. No edge is
forced onto a vertical side. Crossings are not created because the spherical
embedding had none.

This drawing is a topological *st*-strip. It does not require *x*-monotonic
edges.

### 2. Replace terminal vertices by occupied boundary paths

Delete `L` and `R`. Let `v_1, …, v_{d_L}` be the neighbours of `L` in the
cyclic order of the interior rotation at `L` (the unique cyclic order of
`L`-incident darts that does not enter `F`; parallel `L`–`v` edges appear as
consecutive darts). Place `d_L` distinct **deterministic occupied** vertices
`ℓ_1, …, ℓ_{d_L}` on the open bottom side `(0,1) × {0}`, ordered by
increasing *x*, and join each `v_i` to `ℓ_i` by a simple arc that meets the
bottom side only at `ℓ_i`. Do the same on the top side: occupied vertices
`r_1, …, r_{d_R}` and arcs `w_j r_j` for the `R`-neighbours in interior
cyclic order.

These rim vertices are never switchable and are never sampled. They belong
to the checkpoint in every future configuration.

If several `L`-incident edges left `L` in a common angular sector, they still
appear as a consecutive block in this order; cofaciality is what makes that
block well-defined.

**Length.** Let

```text
m = max(3, d_L, d_R).
```

Subdivide the bottom side by additional occupied degree-2 vertices, placed
between or beyond the `ℓ_i` and not incident to any switchable vertex, until
the bottom occupied set `P_L` is a path of exactly `m` vertices. Do the same
on the top, obtaining a path `P_R` of `m` vertices. Add occupied edges along
each side so that `P_L` and `P_R` are induced paths in the host.

The host in `Q` is now: switchable copy of `H − {L, R}`, dangling arcs to
`P_L` and `P_R`, and the two occupied paths. It remains plane. No switchable
vertex lies on `∂Q`. No graph element uses the open vertical sides.

**Cyclic order of half-edges.** Around each bottom occupied vertex the local
rotation is: left along `P_L`, right along `P_L`, and (if the vertex is some
`ℓ_i`) the unique interior dart to `v_i`. The sequence of interior darts
along `P_L` reproduces the interior rotation at `L`. The top side likewise
reproduces the interior rotation at `R`.

**Multiple incidences; terminal degree greater than 1.** Distinct attachment
vertices on `P_L` are used precisely so that the host remains a genuine
embedding: two interior darts are never forced through one boundary point.
No switchable vertex is subdivided. After the later contraction of each
occupied rim to a single terminal, all attachments to `P_L` become
attachments to one vertex `L'`; distinctness of rim vertices is not needed
for the isomorphism type of the residual network, only for a simple
embedding of `G_N`. Terminal degree 1 is the same construction with idle
degree-2 points of the future cycle `γ`.

### 3. Cylinder, then torus

Identify the vertical sides of `Q` by `(0, y) ∼ (1, y)` for `y ∈ [0,1]`.
This is an orientation-preserving homeomorphism of the two sides, which carry
no graph. The quotient is a closed cylinder `C ≅ S^1 × [0,1]`.

- The bottom path `P_L` has its two endpoints identified, hence becomes a
  simple occupied circle `γ_L` of length `m ≥ 3` (the bottom rim).
- The top path `P_R` becomes a simple occupied circle `γ_R` of length `m`.
- Switchable vertices remain in the open cylinder `S^1 × (0,1)`.
- Interior edges do not wrap around the `S^1` factor: they were drawn in the
  open rectangle. The embedding on `C` is still crossing-free.

The two rims are disjoint simple cycles. The occupied graph on `C` is
`γ_L ∪ γ_R`.

Now glue the two rims by an orientation-preserving homeomorphism
`φ: γ_L → γ_R` that matches vertices. Because both rims are discretised
circles of the same length, any such homeomorphism is a **cyclic shift**.
The quotient `T = C / φ` is an oriented torus.

Write `γ` for the image of the rims and `A_N` for `V(γ)`. Write `G_N` for the
image of the cylinder graph. Then `G_N` is a finite graph honestly embedded
in `T`: vertices are distinct, edges are simple arcs meeting only at
endpoints, and `γ` is a simple closed curve of graph edges.

**Geometric identification of edges.** The cyclic shift identifies the edge
`ℓ_i ℓ_{i+1}` with the edge `r_{φ(i)} r_{φ(i+1)}`. The two rim graphs become
a single simple cycle of length `m` on the torus, not a doubled cycle. The
occupied graph is this cycle.

**Local sides of `γ`.** In angular coordinates `(θ, ψ) ∈ R^2 / Z^2` with the
cylinder as `ψ ∈ [0,1]` before the identification `ψ = 0 ∼ ψ = 1`, the curve
`γ` is `{ψ = 0}`. The complement `T \ γ` is the open cylinder
`S^1 × (0,1)`, which is connected, so `γ` is non-separating. Locally `γ` has
two sides: attachments that approached `P_L` do so as `ψ → 0^+`; attachments
that approached `P_R` do so as `ψ → 1^- = 0^-`. Cutting along `γ` will
separate those two local sides into two boundary components.

**Choice of `φ`.** If some switchable vertex `v` is adjacent to both rims
(a singleton trigger in `N`), its two attachment vertices must not be
identified by `φ`; otherwise `G_N` would contain a pair of parallel edges
from `v` to one vertex of `γ`. A cyclic shift of `φ` avoids this whenever
`v` uses at most `m − 1` attachments on at least one rim, because `m ≥ 3`.
If no cyclic shift avoids the collision for every such `v` simultaneously
(possible only if some vertex is adjacent to every rim vertex on both
sides), retain a finite multi-edge embedding: two parallel non-crossing arcs
from `v` to one vertex of `γ`, approaching from opposite local sides. This
is still a genuine embedding of a finite multigraph. After the cut those two
arcs attach to the two distinct boundary copies, hence become a `v`–`L'`
edge and a `v`–`R'` edge rather than a pair of residual parallels.

If `H` itself had parallel `v`–`L` edges, the residual contraction deletes
duplicates and the isomorphism type of the simple graph underlying `N` is
what the lemma claims.

### 4. Occupied homology is exactly rank one

The occupied vertex set is `A_N = V(γ)`. The occupied edge set is `E(γ)`.
The occupied graph is a simple cycle. Its ambient homology class is the class
of the glued rim.

That class is essential and primitive: `γ` is the image of a boundary
component of a cylinder under the identification of the two boundaries, hence
is a non-separating simple closed curve on `T`, and it generates a direct
summand of `H_1(T; Z)`. In the coordinates of §3, `γ = {ψ = 0}` has winding
`(0,1)` up to sign.

There is no other occupied edge, so there is no occupied cycle linearly
independent of `γ`. In particular there is **no occupied transverse cycle**:
a cycle homologous to `(1,0)` would require an occupied path through
`T \ γ`, and `T \ γ` contains no occupied vertex. The ambient rank is
exactly one.

### 5. Cut, contract, recover `N`

Cut `T` along `γ`. This reopens the cylinder `C`, with two occupied boundary
copies: the bottom copy is `γ_L` (the `ψ → 0^+` local side), the top copy is
`γ_R` (the `ψ → 0^−` local side). Each copy is connected.

Contract each occupied connected component to a permanent vertex. There are
exactly two such components. Call them `L'` (image of `γ_L`) and `R'` (image
of `γ_R`). Retain every switchable vertex. Delete loops and duplicate edges
created by contraction (planar graph-minor operations, as in the cut-network
construction).

The residual rooted graph `N'` has:

- switchable vertex set `V_sw`, with all edges of `H − {L, R}` intact;
- an edge from `v ∈ V_sw` to `L'` if and only if `v` was joined to `P_L`,
  if and only if `{v, L}` was an edge of `H`;
- likewise to `R'`;
- no residual `L'`–`R'` edge: the two occupied components remain distinct
  because `γ` is essential of rank one, equivalently because `H` had no
  `{L, R}` edge and the construction added none in the residual.

The obvious bijection `L' ↦ L`, `R' ↦ R`, identity on `V_sw` is a rooted
graph isomorphism `N' ≅` (simple graph underlying `N`).

### 6. Future connection is ambient rank two

Let `U ⊆ V_sw`. Write `A(U) = A_N ∪ U` for the occupied vertex set of the
future configuration, and write `G(U)` for the induced subgraph of `G_N` on
`A(U)`. Terminal connection in `N` is equivalent to residual `L'`–`R'`
connection in `N'` with `U` switched on, by the isomorphism of §5.

It remains to compare residual connection with ambient homology.

`T \ γ` is the open cylinder `S^1 × (0,1)`. The graph `G(U)` consists of the
cycle `γ` together with whatever subgraph of switchable vertices in `U` is
occupied, attached to `γ` along the dangling arcs of §2.

- Suppose `U` connects `L'` to `R'` in `N'`. Then `G(U)` contains a path `π`
  that leaves `γ` through a `P_L`-attachment (`ψ → 0^+`), travels in the
  open cylinder, and returns to `γ` through a `P_R`-attachment (`ψ → 0^−`).
  Closing `π` along an arc of `γ` produces a simple cycle `τ` whose `ψ`-winding
  is `±1`. Thus `[τ]` is independent of `[γ]` in `H_1(T; Z)`, and the ambient
  rank of `A(U)` is two.
- Conversely, suppose `U` does not connect `L'` to `R'`. Then every connected
  component of `G(U) − E(γ)` attaches to at most one local side of `γ`. Any
  cycle of `G(U)` is therefore homologous to an integer multiple of `γ`: a
  path that leaves `γ` and returns on the same local side, union an arc of
  `γ`, is supported in a closed cylinder attached to one boundary component,
  hence has `ψ`-winding `0`. The ambient rank of `A(U)` remains one.

Composing the two equivalences: `U` connects the terminals of `N` if and only
if `r(A_N ∪ U) = 2`.

(Remark: the same comparison, for an arbitrary honestly embedded rank-one
checkpoint rather than for the host constructed here, is the cut-network
theorem. The realization itself — the construction of `(G_N, A_N, γ)` and
the isomorphism `N' ≅ N` — does not invoke it. The local-side argument above
is written out so that the rank-two clause is self-contained for this host.)

### 7. Loops, parallel edges, subdivision

- **Loops at switchable vertices.** Copied into `G_N` and into `N'`. They do
  not meet `γ` and do not affect terminal connection.
- **Loops at `L` or `R`.** Discarded in the definition; they never enter the
  host.
- **Loops on `γ`.** Not created: `γ` is a simple cycle of length `m ≥ 3`.
- **Parallel edges in `H`.** Each copy is attached as a distinct interior
  dart to a distinct rim vertex; after contraction, duplicates are deleted,
  matching the simple residual network used as a two-terminal vertex-network.
- **Parallel edges in `G_N`.** Avoided for singleton triggers by the shift of
  `φ` in §3 when a shift exists; otherwise retained as a genuine multi-edge
  embedding with opposite local sides, which the cut separates. Not otherwise
  created: each switchable–rim incidence uses a distinct rim vertex on its
  own side.
- **Subdivision.** Occupied degree-2 vertices may be added on the rims so
  that each rim has length `m ≥ 3` and so that terminal degree `> 1` uses
  distinct attachment vertices. No switchable vertex is subdivided: that
  would change `V_sw`. Attachment edges are not subdivided.

This completes the proof.

---

## What the lemma does not say

- It does not produce a nearest-neighbour square-site occupation, a column-HNF
  period matrix, or an explicit occupied mask on any named lattice torus.
- It does not claim that every residual cut-network of a square-NN torus is
  plane with cofacial terminals (the cut-network theorem does not assert
  cofaciality).
- It does not identify `G_N` up to isomorphism of torus graphs: the integer
  `m`, the shift of `φ`, and idle occupied vertices are gauge.
- It does not discuss a sampling measure.
- It does not prove that the cut-network representation is a minimal
  sufficient statistic.

---

## Relation to existing notes

| Statement | Direction | Scope actually proved |
|---|---|---|
| Cut-network theorem | torus checkpoint → residual network | honestly embedded rank-one graphs |
| Parallel-gadget §6 | *k* parallel copies of two N16 gadgets → torus | those two gadgets, after they are already cut-networks |
| **This lemma** | plane TT network → torus checkpoint | every finite plane two-terminal vertex-network with cofacial terminals |

The three statements compose. The cut-network theorem plus this lemma say
that, on the class of plane two-terminal vertex-networks with cofacial
terminals and no terminal edge, residual networks and rank-one torus
cut-networks are equivalent as continuation objects. The parallel-gadget
construction remains a special embedding for a product of two particular
gadgets; it is not a surrogate for this lemma.
