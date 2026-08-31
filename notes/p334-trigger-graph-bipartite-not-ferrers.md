# Trigger-pair structure: bipartite survives; componentwise Ferrers fails

The two saved real N425 checkpoints both have bipartite minimal-trigger-pair
graphs, but **neither is componentwise chain/Ferrers**. A deterministic additional
check of the first five eligible archived counters in each of four fixed
geometries gives 20/20 bipartite graphs but only 11/20 componentwise chain graphs.
The total bounded census is therefore **22/22 bipartite, 11/22 Ferrers**.

This task tested only those two predeclared graph classes. It did not search
for a favorable subclass, train a classifier, generate new Monte Carlo, or
expand after finding a counterexample.

## Graph and selection semantics

Vertices are individually safe vacant insertions at a rank-one checkpoint.
An edge is an unordered pair of these sites that jointly triggers ambient rank
two. Thus the graph is the complement, on safe vertices, of the safe-pair graph
used in the previous notes. Its wedges are
`W2=sum_v binom(trigger_degree(v),2)`.

The two previously selected N425 counters are `43042514269` / `43042505280`,
with the same seed `20260831430425`, period matrix `[[425,268],[0,1]]`,
`k0=252, age=10, ell=(12,-19), H2=0, b1_safe=173, b2_safe=14770`.

The fixed extension uses, independently within N325/N425 first/second, the
five lowest counters among archived rank-one rows with at least one minimal
trigger pair. Eligibility does not depend on bipartiteness or nesting. Complete
selected rows and trigger adjacency are saved. The production potential
union-find and exact counter permutations are reused; no new sample stream is
run.

## The two existing graphs already reject the Ferrers explanation

| quantity | A | B |
|---|---:|---:|
| safe vertices | 173 | 173 |
| minimal trigger edges | 108 | 108 |
| W2 | 926 | 1,466 |
| triangles | 0 | 0 |
| nonisolated connected components | 1 | 1 |
| component bipartition sizes | 14 × 12 | 5 × 29 |
| isolated vertices | 147 | 139 |
| bipartite | yes | yes |
| componentwise chain/Ferrers | **no** | **no** |

A has same-side sites 1 and 47 with nonnested neighborhoods:
`6 in N(1)\N(47)` and `50 in N(47)\N(1)`. Their connected-component certificate
is the induced trigger path

```
1 — 6 — 159 — 50 — 47.
```

In particular `(1,6)` and `(47,50)` trigger rank two, while `(1,50)` and
`(47,6)` are safe pairs. This is an induced 2K2 **within one connected component**,
not the irrelevant 2K2 obtained by combining different components.

B similarly has same-side sites 8 and 390, exclusive neighbors 15 and 26,
and the induced trigger path

```
8 — 15 — 121 — 26 — 390.
```

These four-site nonnested-neighborhood certificates are minimal in arity for
failure of a chain graph. The five-site paths also explicitly demonstrate that
the vertices lie in the same component. Full neighborhoods are retained in
the machine-readable certificate.

## Fixed five-counter extension

| environment | selected replica counters | bipartite | componentwise Ferrers |
|---|---|---:|---:|
| N325 first | 43032500005, 00012, 00019, 00020, 00023 | 5/5 | 4/5 |
| N325 second | 43032500001, 00003, 00005, 00008, 00013 | 5/5 | 3/5 |
| N425 first | 43042500001, 00004, 00006, 00009, 00013 | 5/5 | 1/5 |
| N425 second | 43042500000, 00003, 00010, 00013, 00015 | 5/5 | 3/5 |

The shortened entries retain the leading `430325` or `430425` counter prefix
of their row. Every failure has an explicit four-site nonnested-neighborhood
certificate and a same-component path. No odd cycle occurs in this bounded
census. **That preserves bipartiteness as a structural conjecture, not a theorem
over all rank-one configurations.**

## Two-graph posthoc conditional baseline: side capacity explains 84% of W2 contrast

Only for the two already selected graphs, condition on each nonisolated
component's bipartition sizes L,R and edge count m. Uniformly selecting m edges
from the L*R possible cross edges gives exactly

```
E W2 = [L*C(R,2) + R*C(L,2)] * m*(m-1) / [L*R*(L*R-1)].
```

Proof: there are `L*C(R,2)+R*C(L,2)` candidate wedges, each requiring two
specified distinct edges; their inclusion probability is
`m*(m-1)/[L*R*(L*R-1)]`. Expectations add across components. The two saved graphs
each have one nonisolated component.

| conditional arithmetic | A: 14×12, m=108 | B: 5×29, m=108 |
|---|---:|---:|
| observed W2 | 926 | 1466 |
| E W2 | 138672/167 = 830.371257 | 1284 |
| observed minus E W2 | 15970/167 = 95.628743 | 182 |

The observed contrast is 540. Its conditional-baseline part is

```
1284 - 138672/167 = 75756/167 = 453.6287425,
(75756/167)/540 = 6313/7515 = 84.0053227%.
```

The residual excess contrast is `14424/167=86.3712575`, or **15.9946773%**.
Thus most of this particular two-graph overlap contrast is already present in
the different bipartition capacities; a Ferrers/nested-neighborhood model is
unnecessary and false here. Both graphs still exceed their conditional
uniform-bipartite expectations. This is **posthoc structural decomposition of
two fixed graphs**, not an independent hypothesis test, causal attribution, or
population-wide percentage.

## Scientific card and reproduction

- Changed mechanism space: componentwise Ferrers is eliminated by real
  same-component witnesses; bipartiteness survives the fixed 22-graph check.
- Positive structural lead: bipartition capacity accounts for 84.0053% of the
  selected A/B W2 difference under the explicit conditional random-graph baseline.
- Source/sector: minimal two-site trigger graphs on occupied-NN ambient rank-one
  checkpoints, fixed N325/N425 Gaussian quotients; no cross-geometry comparison
  is used as an exact state counterexample.
- Dependency: the existing cooperative-production blocks, including the same
  previously selected A/B witnesses; no new data block or independent replication.
- Not established: universal bipartiteness, a dual-annular proof, closure of the
  safe complex at triples, a continuum field or full-state temporal memory.
- Next implication: pursue the weaker two-sided/topological-cut structure rather
  than neighborhood nesting; no further graph search is initiated by this note.

```bash
clang++ -std=c++17 -O2 -Wno-unknown-pragmas \
  src/p334_checkpoint_trigger_graph.cpp -o /tmp/p334-checkpoint-trigger-graph
/Users/lc/python-envs/research-py311/bin/python scripts/p334_trigger_graph_structure.py \
  --binary /tmp/p334-checkpoint-trigger-graph --scope bounded
```

The script reuses committed adjacency when present. Artifacts are
`trigger_graph_raw/`, `trigger_graph_structure_existing.json`, and
`trigger_graph_structure_bounded.json` under the cooperative result directory.
