# Real-lattice scalar-state nonclosure: same (age, ell, H2, b2), different branching

The same two N425 checkpoints now also have an exact
[safe-triangle versus actual-safe-triple census](p334-real-safe-triple-census.md):
their pair graphs miss 583/509 minimal triple triggers under clique completion.

**The archived production contains exact counterexamples to using
`(N, orientation, k0, age, ell, H2, b2)` as a state sufficient for
one-common-update/two-clone survival.** No new Monte Carlo was run. These are
actual N325/N425 lattice checkpoints from result commit `e81dd59`, not abstract
graphs, comparisons across different geometries, or stochastic clone outcomes.

## A concrete N425 witness

Fix the second N425 Gaussian quotient, represented by the exact period matrix
`[[425,268],[0,1]]`, lineage `(19,8)`. Both checkpoints have

```
k0=252, d=173, H2=0, b1=173, b2=14770,
age_steps=10, primitive ell=(12,-19).
seed=20260831430425.
```

| exact quantity | checkpoint A | checkpoint B |
|---|---:|---:|
| archived replica counter | 43042514269 | 43042505280 |
| sum c_v | 29,540 | 29,540 |
| sum c_v^2 | 5,045,796 | 5,046,876 |
| `b1*sum(c_v^2) - (2*b2)^2` | 311,108 | 497,948 |
| unordered 2-stars | 2,508,128 | 2,508,668 |
| exact two-clone survival | 1261449/1279508 | 1261719/1279508 |

Thus the two-clone prediction differs by

```
1080 / [173 * 172^2] = 135/639754 = 0.00021101861027832572.
```

The checkpoint-specific closure baseline `s2^2/s1` is identical, so the
cooperative excess also differs by exactly the same fraction. Their complete
original production rows, seed/counter, occupied prefix labels and bitmask,
period matrix and exact rational outputs are in
`results/local-20260831/P334-cooperative-closure/scalar_state_collisions.json`.
The archived counter is decoded with the frozen runner's uint64 SplitMix/Fisher–Yates
rule; its next site is checked against the saved row. No topology production is
repeated by this extraction.

This proves the precise nonclosure claim: if a function of the displayed state
tuple determined the branching probability for every checkpoint, it would have
to return the same number on A and B. The exact archived outputs are unequal.

## It is present within all four geometries, not a cross-geometry effect

The search groups only within identical N, orientation and k0. The first count
below additionally fixes H2 and b2; the second also fixes age and primitive ell.
Each listed group contains at least two distinct exact degree-square sums.

| environment | H2/b2 collision groups | age/ell/H2/b2 collision groups | selected exact branching difference |
|---|---:|---:|---:|
| N325 first | 1,064 | 663 | 122/566313 |
| N325 second | 1,088 | 665 | 89/377542 |
| N425 first | 1,290 | 458 | 483/2559016 |
| N425 second | 1,313 | 476 | 135/639754 |

One exact witness per environment is saved. Selection maximizes the degree-square
spread within an age/line-matched group; it is post-hoc counterexample selection,
not a hypothesis test or a population effect-size claim. No additional samples
are necessary to establish the finite-state nonclosure.

## What the branching observer reads

Define the **safe-insertion graph** of the current checkpoint: its b1 vertices
are safe vacant-site insertions, and an edge joins v,w when adding both sites
still leaves ambient rank one. It has b2 edges. Monotonicity means an unsafe
singleton cannot belong to a safe pair. The degree of v is the exact number
`c_v` of safe second insertions after v.

Consequently,

```
sum_v c_v = 2*b2,
sum_v c_v^2 = 2*b2 + 2*sum_v binom(c_v,2).
```

The last sum counts unordered 2-stars (wedges with a specified center). Vertex
and edge counts do not determine this overlap count. The two independent clones
read the degree second moment, including the possibility that the clones choose
the same second site. Their exact probability is
`sum_v c_v^2/[d*(d-1)^2]`.

The new observable is therefore detecting **overlap structure among safe
continuations that first- and second-step scalar counts omit**. It need not be
called path memory: the complete current lattice configuration still determines
the process. Nor have these witnesses been shown to share the entire single-chain
survival curve; they prove nonclosure of the explicitly matched state tuple.

## Scientific card

- Changes: exact real-lattice nonclosure of `(age, ell, H2, b2)` for this branching
  prediction, stronger than merely rejecting uniform safe successors.
- Observer/sector: one-common-update/two-one-site-clone survival; black-NN rank-one
  plateau and rank-two absorption; same physical quotient within each witness.
- Source/dependency: zero-sample extraction from the existing
  `p334-cooperative-N325-20260831` / `p334-cooperative-N425-20260831` blocks; not
  additional independent evidence or a new random block.
- Not proved: temporal memory of the full configuration, full trace-equivalence,
  a scale law or field identity.
- Next implication: any proposed reduced state for this continuation language
  must retain or determine this 2-star/degree-second-moment information.

Reproduce with

```bash
/Users/lc/python-envs/research-py311/bin/python scripts/p334_checkpoint_scalar_collision.py
```
