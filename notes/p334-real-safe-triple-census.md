# Two saved N425 checkpoints: pair-clique closure misses 583 / 509 triple triggers

**The new result is a quantitative real-checkpoint census, not a first discovery
that a safe-insertion complex can be non-flag.** PR #415 / `d09f925` already
contains a genuine N13 minimal-triple example and the three-step overlap
identity. After checking #403 and #429, this calculation instead resolves the
entire third layer for the two large, actually sampled N425 checkpoints saved
by `6147e22`. No new Monte Carlo was run.

Both checkpoints use the same second N425 geometry, exact period matrix
`[[425,268],[0,1]]`, with HNF label v represented by `(v,0)`. They have
`k0=252, d=173, age=10, ell=(12,-19), H2=0, b1_safe=173, b2_safe=14770`.
The seed is `20260831430425`; A/B counters are
`43042514269` / `43042505280`.

## Exact third-layer counts

Let T3 be the number of triangles in the safe-pair graph, b3 the number of
vacant triples whose insertion leaves ambient rank one, and c3 the number of
minimal rank-two-triggering triples: each proper subset is safe but the whole
triple is not. Monotonicity gives exactly `b3=T3-c3`.

| count | A | B | B−A |
|---|---:|---:|---:|
| all vacant triples C(173,3) | 848,046 | 848,046 | 0 |
| safe-pair graph triangles T3 | 830,504 | 831,044 | 540 |
| actual safe triples b3 | 829,921 | 830,535 | 614 |
| minimal nonfaces c3 | 583 | 509 | −74 |

Therefore the exact pair-clique overprediction of three-step survival is

```
A: 583/848046 = 0.0006874627083908185,
B: 509/848046 = 0.0006002032908592222.
```

The true three-step survival difference, B−A, is

```
614/848046 = 307/424023,
```

whereas graph-triangle counts alone give only

```
540/848046 = 90/141341.
```

The remaining `74/848046=37/424023` comes from the differing minimal-triple
layer. This is an additional higher-order contribution even after the pair
overlap/degree-second-moment structure measured in the previous result.

## Explicit minimal nonfaces on the actual lattice

The lexicographically first examples are A sites **{1,25,73}** and B sites
**{8,13,14}** in the declared exact quotient coordinates.

| added subset | A ambient rank | B ambient rank |
|---|---:|---:|
| empty | 1 | 1 |
| first singleton | 1 | 1 |
| second singleton | 1 | 1 |
| third singleton | 1 | 1 |
| first + second | 1 | 1 |
| first + third | 1 | 1 |
| second + third | 1 | 1 |
| all three | 2 | 2 |

Thus each is a size-three minimal nonface, not an unsafe pair with a redundant
third site. The artifact retains the first eight nonfaces per configuration,
the eight subset ranks for the first witness, original occupied prefixes,
seed/counter and exact matrix/line/age.

## What is and is not determined by a pair graph

The safe complex is not equal to the clique complex of its safe-pair graph:
pairwise compatibility does not imply triple safety. Predicting its third
cardinality layer from the graph triangle count requires the missing c3 term.
This is the precise sense of pair-clique nonclosure used here.

We do **not** claim these two checkpoints have identical complete pair graphs
or identical triangle counts; they do not. Nor does this result say that an
arbitrary function of a fully geometrically labelled graph could never recover
the underlying configuration. The established failure is the explicit pairwise
compatibility/clique completion rule. The earlier scalar-state witness remains
the separate proof that matching `(age,ell,H2,b2_safe)` cannot close branching.

Notation is intentional: `b2_safe` counts safe pairs, as in the new production;
#403's original `b2` instead denoted **minimal triggering pairs**. They are not
the same variable.

## Reproduction and scientific card

`src/p334_checkpoint_safe_triples.cpp` includes the unchanged production source
and reuses its exact period geometry, counter permutation and potential
union-find. It copies the two saved checkpoint snapshots and exhausts their
safe-pair triangles; no full production stream or old test suite is rerun.
Snapshot b1/b2/degree-square sums agree with the archived integers, and the
explicit nonface subset ranks are saved. This is reuse of the same backend,
not a claim of independent-oracle validation.

```bash
clang++ -std=c++17 -O2 -Wno-unknown-pragmas \
  src/p334_checkpoint_safe_triples.cpp -o /tmp/p334-checkpoint-safe-triples
/tmp/p334-checkpoint-safe-triples \
  results/local-20260831/P334-cooperative-closure/safe_triple_census.json
```

Local execution took 1.40 s wall / 0.85 CPU s, exit 0.

- Changes: on two real large checkpoints, pair-clique closure loses exactly
  583/509 triple triggers; actual third-step probabilities and their difference
  are now known, not inferred from b1/b2 or graph wedges.
- Observer/sector: three distinct uniform vacant insertions, occupied-NN rank-one
  survival / rank-two absorption, one fixed N425 physical quotient.
- Source/dependency: replays of the two existing
  `p334-cooperative-N425-20260831` checkpoints; no new sample block or independent
  statistical replication.
- Not established: first-ever non-flag example, equal full pair graphs,
  population prevalence, a scale law or temporal memory of the full lattice state.
- Next implication: a third-layer reduced prediction must account for the minimal
  triple layer as well as pair-overlap structure; no new stream is proposed here.
