# Program D: no bounded-radius summary suffices at any fixed radius

**Status:** THEOREM (exact, family-level) — upgrades #550's r=1 certificate
to **every fixed radius**.  Realises Program D's "no bounded-radius quotient"
target (major outcome D) in the sharpest form the published protocol algebra
supports.

Script: `scripts/radius_insufficiency.py`; data:
`results/radius-insufficiency/latest.json`.

---

## 1. The construction

Take the #549 parallel-gadget bank: `k` future-vertex-disjoint N16 gadgets,
`a` of them of type A.  Attach the bank's future-vertex region to the two
terminals `L, R` by two bridge paths of length `L > r`, every bridge vertex
**already occupied** (so the bridge contributes no future random variable).
For any fixed radius `r`, choose `L > r`:

1. the **complete unbranched survival law** is `S_k(z) = S(z)^k`, independent
   of `a` (the parallel product of the equal safe-subset polynomials);
2. the **radius-r terminal neighbourhood** (graph distance from L and R)
   contains only bridge vertices and no future vertex, hence is identical
   across every class `a` and across the two gadget types;
3. the **single fork** probability `F_{k,a}` is strictly increasing in `a`,
   with gap `F_{k,a+1} - F_{k,a} = 1/[2k(8k-1)^2]`.

Consequently, for **every** fixed `r`, there are `k+1` networks that share the
same unbranched law and the same radius-r terminal neighbourhood yet differ in
branching behaviour.  The branching signal is carried by future vertices
arbitrarily far from the terminals.

## 2. Statement

**Theorem.**  For every radius `r >= 0` and every `k >= 1`, the `k+1` hidden
classes `a = 0..k` of the #549 parallel bank admit a realisation as planar
two-terminal vertex-networks in which the unbranched survival law and the
radius-r terminal neighbourhood are identical, but the single-fork probability
differs between consecutive classes by `1/[2k(8k-1)^2]`.

In particular there is **no** finite-radius terminal-neighbourhood summary that
determines even the depth-1 (one-fork) language.  A fortiori no bounded-radius
quotient of the cut network preserves all branching futures.

## 3. Relation to the repository

* #550 proved this at `r = 1` for a fixed depth-2 composition, with the caveat
  that radius 2 separated its witness.  The construction here removes the
  caveat: by pushing the parallel bank beyond the radius, the failure holds at
  every radius simultaneously, and it holds already for the **single fork**
  (depth 1), not only for a depth-2 composition.
* The gap `1/98` at `k = 1` reproduces #435's exact witness.
* The bridge is declared occupied, so it contributes no factor to `S(z)` and
  no future variable — this is the same "fixed-cut gauge" convention #491/#550
  use (the network is a proof device in a gauge, not a physical mark).

## 4. Boundary

* Exact statement about the declared composition of the published #549 fork
  and the occupied-bridge embedding convention; the full N16 site data is not
  needed because the only moving part is the parallel product of the already
  equal gadgets.
* "Neighbourhood" is graph-distance radius on the vertex network; "future
  vertex" means a vertex in the sampled future region.  No geometric
  (Euclidean) claim about square-site tori is made here.

Files: `scripts/radius_insufficiency.py`,
`results/radius-insufficiency/latest.json`.
