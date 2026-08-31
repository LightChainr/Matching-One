# Minimal pair triggers are bipartite: two white essential cycles

**Result:** the bipartiteness observed in the 22 saved checkpoints is forced
by annular topology. It is not a Ferrers/nested-neighbourhood law. No new
Monte Carlo is used.

## Theorem and exact meaning

Let `B` be a black-NN occupied configuration on a finite square torus, with
ambient homology rank one. Let

```text
S = {v vacant : rank(B+v)=1},
uv in E(T_B) iff u,v in S and rank(B+u+v)=2.
```

Then the minimal-trigger-pair graph `T_B` is bipartite. More concretely there
exist two essential cycles `C0,C1` in a fixed planar incidence realization of
the white matching graph such that

```text
(C0 intersect C1) intersect S = empty,
every trigger pair has exactly one endpoint on each cycle.
```

Consequently all nonisolated trigger vertices lie on the two cycles. Vertices
on neither cycle are isolated; the cycles can contain other isolated vertices.
The resulting partition is not required to have nested neighbourhoods, and
need not be canonical before choosing the two cycles. A connected nontrivial
trigger component has the usual unique bipartition up to exchanging sides.

The theorem uses the configurationwise square-NN/matching rank duality
`r_black+r_white=2`, not graph cyclomatic rank. On honest torus cells this is
the repository's digital-Alexander theorem. The same argument applies to
arbitrary nonsingular integer-period quotients with the unrestricted rank
duality of `73d4960`: loops, parallel spokes and identified face corners are
retained, not collapsed. It needs no integral-index assumption or exponent.

## 1. A fixed white graph with the right deletion semantics

For every unit square retain a face-centre auxiliary vertex. Join it by the
four straight half-diagonal spokes to its white corners. Keep every spoke,
including parallel incidences after quotienting. This is an embedded graph
on the torus, unlike the drawing with both white matching diagonals crossing.

For **every** subset of deleted white sites, this same graph with those site
vertices deleted has the same lifted connectivity and ambient homology image
as the remaining white NN+NNN matching graph:

- a matching edge is replaced by the two spokes in its square;
- a two-spoke passage between remaining corners is replaced by their NN or
  diagonal matching edge;
- each replacement has the same lifted endpoints and displacement.

This is a simultaneous deletion-compatible construction, not a pruning chosen
separately for each pair. Extra contractible cycles in either graph are
irrelevant to ambient image rank.

Choose a simple essential black-NN cycle `gamma`. It exists because the black
ambient rank is one. The white spine is disjoint from `gamma`. Cutting along
a sufficiently thin neighbourhood of `gamma` embeds the white spine in an
annulus. Every essential white cycle separates its two boundary components.

## 2. Capacitated annular cycle lemma — proof, not an assumed Menger lift

Let `G` be any finite graph embedded in the interior of an annulus. Assign a
positive integer capacity `c(v)` to each vertex. The minimum total capacity
of a vertex set meeting every essential cycle equals the maximum number of
essential cycles in a packing in which each vertex is used at most `c(v)`
times. Repeated cycles are permitted when capacity permits.

Here is a direct face-distance proof.

Cap the two annulus boundaries by disks, giving a plane/spherical embedding
with distinguished faces `f_out,f_in`. Its face-incidence network allows a
move from any face incident with `v` to any other such face at cost `c(v)`.
Let `d(f)` be shortest-path distance from `f_out`, and `D=d(f_in)`.

1. A face-incidence path from `f_out` to `f_in` gives a vertex transversal:
   every essential cycle must meet its annulus-crossing route. Conversely,
   deleting a transversal merges `f_out` and `f_in` into one face. The merging
   can be described by a simple face/vertex incidence path using only deleted
   vertices. Thus its cost is at most the transversal cost, and the minimum
   transversal capacity is exactly `D`.
2. For each integer `i=0,...,D-1`, colour a face low when `d(f)<=i`. Take the
   primal edges separating a low from a high face. Around each vertex colour
   changes occur an even number of times, so these edges form an Eulerian
   boundary. This boundary separates `f_out` from `f_in`; at least one of its
   simple cycle components separates these faces and is essential in the
   original annulus. Choose one such cycle `Ci`.
3. Write `d_min(v),d_max(v)` for the extreme incident-face distances. The
   network move through `v` gives `d_max(v)-d_min(v)<=c(v)`. A boundary at
   level `i` can use `v` only when

   ```text
   d_min(v) <= i < d_max(v).
   ```

   Therefore at most `c(v)` of the chosen cycles use it. We have packed `D`
   essential cycles. Every packing has size at most every transversal's
   capacity by counting its intersections with that transversal, proving the
   equality.

This proof accommodates disconnected graphs, non-cellular annular regions,
bridges, loops and parallel edges: after capping, faces are complement regions,
not just separately enumerated boundary walks. At a loop its two incidences
are retained in the even-degree argument. It is **not** a min-max assertion
for arbitrary cycle families in arbitrary graphs.

## 3. Apply capacity one only to safe site insertions

In the white spine assign

```text
c(v)=1   for safe white sites v in S,
c(v)=2   for other white sites and every auxiliary face centre.
```

A capacity-zero transversal is impossible because white rank is one. A
capacity-one transversal would consist of one safe white site, contradicting
its definition and rank duality. Hence `D>=2`. The lemma gives `C0,C1` that
share no capacity-one site; sharing unsafe sites or auxiliary centres is
allowed and harmless.

If `uv` is a trigger edge, deleting `{u,v}` kills all essential white cycles,
so it meets both `C0` and `C1`. Both endpoints are safe and the cycles are
safe-disjoint. Thus one endpoint lies on each cycle. This proves the theorem.
When any trigger edge exists its capacity is two, so `D=2` exactly: the pair
itself is a dual optimum certifying the two-cycle packing.

The annular step is essential. For the family of **all** cycles in abstract
`K4`, every singleton deletion is safe and every pair destroys all cycles;
its pair-trigger graph is `K4`, not bipartite. That is not a lattice witness:
the common two-face/annular essential-cycle geometry is absent.

## 4. Exact realization on existing production configurations

`scripts/p334_trigger_bipartite_cycles.py` reuses only the 22 configurations
already selected at `1b5a9de`. Their committed counters reconstruct the exact
occupied prefix; no new permutations are sampled. It constructs the fixed
white spine and uses an **integer min-cost circulation** to find two essential
cycles with the declared capacities. Every cycle is saved as its actual
site/face spoke walk with deck increments. The saved pairs must cross the
two cycle-site sets. The cycle certificate can be checked without trusting
the circulation optimizer or running new pair enumeration.

All 22 existing configurations have such a certificate. The algorithmic
construction is an illustration of the theorem, not the theorem's proof.
Small new-property enumeration covers every HNF quotient of index 1 through
6, plus axis N9 and Gaussian N10/N13. It checks bipartiteness on rank-one
states with nonempty trigger graphs and saves the first two-cycle certificate
per qualifying quotient. This is not an additional Monte Carlo block or a
repeat of the Ferrers/triple analysis.

The bounded census contains **36 quotients, 3,388 rank-one states and 1,804
nonempty trigger graphs**, with no counterexample. Two focused certificate
tests pass: all 22 archived cycle walks and every retained tiny cycle walk are
checked directly from corner coordinates, adjacency, deck sums and the saved
pair sets, independently of the circulation solver. No full repository suite
was rerun.

## Mechanism change

The zero triangle count is now structural: `T3(trigger)=0`. The already-known
three-step identity therefore becomes, for this model,

```text
safe_triples = C(a,3) - m(a-2) + W2 - c3,
```

where `a` is the safe-site count, `m` the minimal trigger-pair count, `W2`
the trigger two-star count, and `c3` the genuine minimal triple layer.
This removes one term by topology; it does **not** remove `c3` or make the
safe complex flag. The observed side-capacity explanation now has a concrete
two-white-cycle origin, while the rejected Ferrers strengthening stays rejected.

This is a finite exact graph/topology result, not an exponent, an identified
continuum field, temporal memory or an ensemble causal percentage.

Reproduction (managed local research Python; no server):

```bash
/Users/lc/python-envs/research-py311/bin/python scripts/p334_trigger_bipartite_cycles.py
/Users/lc/python-envs/research-py311/bin/python -m pytest -q tests/test_p334_trigger_bipartite_cycles.py
```
