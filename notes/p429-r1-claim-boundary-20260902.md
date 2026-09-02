# Claim boundary — r=1 non-compression witness

**Verdict:** `NO_COMPRESSION_WITNESS_FOUND`

Frozen summary (locked before search): `(S(z), n, H2, b2, radius-1 terminal-local neighbourhood)`.
Successor-hazard moments, including delayed-fork `∑x²`, were **not** in the summary.

Claim level: C5 for the finite two-terminal non-sufficiency statement; C1 for the independent exact-rational verifier.

## Killed

- The tuple `(S(z), n, H2, b2, r=1 neighbourhood)` is **not** a sufficient statistic for the frozen depth-2 compositional language, already on a 7-vertex planar two-terminal pair.
- The split is **not** the delayed-fork / successor-second-moment observable: `E1_c1 = 1` on both graphs.
- Series-parallel graphs are not exempt: 330 r=1 splits exist for TTSP n≤12, smallest at n=8, |E|=12, gap 2/1575.

## Not killed

- Radius-2 neighbourhood of this pair: it **differs**, so this is not an r=2 witness.
- Hidden copies of all exhaustive n≤5 cores, over all nine hop combinations (Lh,Rh ∈ {1,2,3}): no r=2 split (0 splits on 66582 graphs). Hidden L2R2 in particular equalizes r=1 and r=2 on the witness S-class and kills the E2_c2 gap.
- Exhaustive plane two-terminal graphs n≤5: closed at both r=1 and r=2 under the frozen language.
- “The cut network is a minimal sufficient statistic.” Not proved.
- Euclidean latent dimension, scalar-encoding impossibility, continuum memory, CFT/LCFT, field counts.
- Universality for every nearest-neighbour HNF torus.
- Explicit torus occupation of the 7-vertex pair (not constructed). Two-port embedding of a planar two-terminal gadget is already enough for the theorem's category.

## Census actually completed

| family | n reached | r=1 splits |
|---|---|---|
| exhaustive plane two-terminal | 5 | 0 |
| hidden hop grid of those cores | 12 (9 combinations) | 1 (L1R1 only) |
| two-terminal series-parallel | **12** (295654 unique graphs) | **330** |
| Wheatstone+SP, multipath, grids | 12 | 0 |

SP n=11 (59324) and n=12 (212771) unique generation is included. The n=7 witness remains smallest in n.

## Independence / chronology

Exact finite combinatorics. No Monte Carlo, no production block, no post-selection of descriptors. Independent of original-`U` identifiability and of the thermal/contact proof gate. Independent of the parallel-gadget successor-second-moment split: that observable was frozen out of the present summary.

## Paper impact

`NEW_REQUIRED_LEMMA` for the P1 cut-network / predictive-state manuscript: a named bounded summary strictly coarser than the cut network fails a predeclared depth-2 language. Do not promote to an all-graphs minimality theorem.

Cut-dependence verdict `CUTS_LAW_EQUIVALENT_NOT_ISOMORPHIC` and embedding verdict `TWO_PORT_EMBEDDING_SUFFICES` are recorded in the accompanying notes. The cut network may be written as sufficient for unmarked continuation and **not minimal for the declared summary class**.
