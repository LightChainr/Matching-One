# P429: a dual-cycle blocker certificate for the observed bipartition

## Result first

The two already selected real N425 checkpoints in `6147e22`, subsequently
used for the pair-graph analysis at `c827cd8/1b5a9de`, admit two
vertex-disjoint essential cycles in their white-matching complements.
Every minimal two-site black rank-two trigger must remove one site from
each cycle. Their previously abstract bipartitions therefore have explicit
geometric certificates.

| checkpoint | cycle lengths | nonisolated cycle sides | trigger edges | W2 |
|---|---:|---:|---:|---:|
| 43042514269 | 20, 25 | 12, 14 | 108 | 926 |
| 43042505280 | 19, 43 | 5, 29 | 108 | 1466 |

Both pairs of cycles have windings `(12,-19)` and `(-12,19)` in the
period basis `[[425,268],[0,1]]`. Their complete ordered vertex lists,
occupied masks, source commit and expected summaries are in
`results/p429-dual-cycle-blocker/certificate.json`.

This is a new geometric explanation of two existing configurations, not
new Monte Carlo, a universal bipartiteness theorem, or another population
estimate. The existing failure of componentwise Ferrers remains intact.

## 1. The finite continuation object has an exact dual description

Let S be a black-NN rank-one checkpoint on a torus in the digital-Alexander
identity domain. Let W be its white-matching complement. For U contained
in W, the existing rank identity gives

```
r_B(S union U) + r_W(W minus U) = 2.
```

Let C(W) denote the family of vertex sets of white-matching graph cycles
with nonzero ambient winding. A graph has nonzero ambient H1 image if and
only if at least one of its cycles has nonzero image. Deleting U preserves
such a cycle if and only if at least one member of C(W) avoids U. Thus

```
r_B(S union U)=2
  iff U intersects every member of C(W).
```

The inclusion-minimal completion hypergraph is exactly the blocker
(minimal-transversal hypergraph) of C(W). In particular,

```
K_safe(S) = union_{C in C(W)} simplex(W minus C).
```

This description retains incidence; its face counts alone do not determine
links or branching laws. If v is individually safe, then

```
K_safe(S+v) = link_v K_safe(S).
```

Equivalently, remove v from every minimal trigger set and retain only
inclusion-minimal resulting sets. These are finite identities, not claims
that enumerating the entire hypergraph is efficient or that it is a
finite-dimensional continuum state.

## 2. Sufficient two-cycle theorem

Let B1 be the set of singleton triggers. Suppose two members C0,C1 of
C(W) satisfy

```
C0 intersection C1 is contained in B1.
```

A minimal trigger pair contains only individually safe vertices and must
hit both cycles. It therefore contains one vertex from `C0 minus B1` and
one from `C1 minus B1`, which are disjoint. Any safe vertex outside their
union is isolated in the trigger-pair graph. The trigger graph is
bipartite.

This proves the sufficient theorem. It does not prove that such a cycle
pair exists for every HNF or checkpoint, nor that it is necessary for
bipartiteness. Other essential cycles can forbid some cross-pairs, with
no implication that neighborhoods are nested. Hence the theorem does
not imply complete bipartite or Ferrers structure.

The B1 allowance is necessary for the proposed certificate language.
For the original N16 witness A (row-major mask 12463), cycles
`[4,9,10,11]` and `[8,9,14,15]` have winding `(1,0)` and intersect exactly
at the singleton trigger 9. The corresponding B witness has overlap at
its singleton trigger 8. Demanding completely disjoint cycles would be an
unnecessarily strong search condition in such examples.

## 3. What was actually checked

Source of the real occupied masks:

```
6147e22f53902a94e5f133739f2c1d423691d0b8
results/local-20260831/P334-cooperative-closure/scalar_state_collisions.json
```

The checkpoints were selected by the preceding work, not selected here
for a favorable graph property. All data remain part of the existing
cooperative-production dependency group.

Discovery used a separately written integer-lift BFS. It reconstructed
both complementary ranks, obtained fundamental nonzero-winding cycles,
and searched for a second essential cycle after removing the first.
Both existing N425 states supplied certificates. The discovery calculation
also checked every one of the 14,878 vacant pairs in each checkpoint.

The committed verifier does not invoke that search. It validates the
supplied cycles directly: white membership, simple vertex sequence,
matching-edge adjacency, exact lifted winding, and allowed intersection.
It independently uses potential union-find, not the discovery BFS, to
reconstruct ranks and singleton triggers. `--full-pairs` checks all 29,756
safe pairs across the two checkpoints, black/white rank duality, and every
trigger edge's cycle-side assignment. It reproduces the existing edge,
W2, side-size and isolate counts.

The cycle certificate alone proves bipartiteness conditional on the
existing rank-duality theorem. Full pair enumeration is an independent
finite consistency check and reproduction of the earlier statistics,
not an additional premise or a new statistical significance test.

Twelve focused tests cover both real certificates, full pair regeneration,
wrong windings, repeated/nonwhite vertices, unsafe cycle overlap, invalid
periods, occupation/singleton mismatch, schema/ID errors, and allowed
singleton overlap on the N16 control. No repository-wide suite was run.

## 4. Next theorem-facing task

Try to return these geometric certificates on the already declared
22-checkpoint set, rather than begin another random graph-class census.
Then address the sharp existence question: can two complementary
essential cycles be chosen disjoint outside the common singleton
bottlenecks? For H2=0 the condition is ordinary vertex-disjointness.

A correctly typed annular packing or block-decomposition argument is a
candidate proof route, not imported as a theorem here. Matching diagonals
and the digital cell realization require care. Search failure is not an
odd-cycle counterexample: the certificate is sufficient, not proved
necessary. Return an explicit geometric obstruction and check whether a
more local, componentwise certificate is enough.

The original 84% capacity-baseline decomposition is still a posthoc
comparison of two graphs, not a causal population percentage. The new
result gives those sides actual positions on complementary homology
carriers; it neither removes genuine minimal triples nor makes pair data
close the full continuation process.

## Reproduce

```
python3 scripts/verify_p429_dual_cycle_blocker.py
python3 scripts/verify_p429_dual_cycle_blocker.py --full-pairs
python3 -m unittest discover -s tests -p 'test_p429_dual_cycle_blocker.py' -v
```

Related: #429, #334, #401, #403, #466. The direct finite proof above uses no
arm exponent, continuum identification, or physical Jordan interpretation.
