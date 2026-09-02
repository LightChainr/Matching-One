# Claim boundary — r=1 bounded-summary insufficiency

**Verdict:** `BOUNDED_SUMMARY_INSUFFICIENT`

Frozen search-protocol token (locked before search, not rewritten):
`NO_COMPRESSION_WITNESS_FOUND`. Same mathematical outcome; manuscript wording
uses the token above.

Frozen summary (locked before search): `(S(z), n, H2, b2, radius-1 terminal-local neighbourhood)`.
Successor-hazard moments, including delayed-fork `∑x²`, were **not** in the summary.

Claim level: C5 for the finite two-terminal non-sufficiency statement; C1 for
the search-independent hard-coded witness verifier (same stdlib primitives as
the search library; **not** a fully independent implementation).

## Killed

- The tuple `(S(z), n, H2, b2, r=1 neighbourhood)` is **not** a sufficient statistic for the frozen depth-2 compositional language, already on a 7-vertex planar two-terminal pair.
- The split is **not** the delayed-fork / successor-second-moment observable: `E1_c1 = 1` on both graphs.
- Series-parallel graphs are not exempt: 330 r=1 splits exist for TTSP n≤12, smallest found there at n=8, |E|=12, gap 2/1575.

## Not killed / UNRESOLVED

- Radius-2 neighbourhood of this pair: it **differs**, so this is not an r=2 witness.
- Hidden copies of all exhaustive n≤5 cores, over all nine hop combinations (Lh,Rh ∈ {1,2,3}): no r=2 split (0 splits on 66582 graphs). Hidden L2R2 in particular equalizes r=1 and r=2 on the witness S-class and kills the E2_c2 gap.
- Exhaustive plane two-terminal graphs n≤5: closed at both r=1 and r=2 under the frozen language.
- Cut-network minimality as a sufficient statistic: **`UNRESOLVED`**. The present lemma shows one named coarser summary fails. It does not prove that the residual network is a minimal sufficient statistic.
- Euclidean latent dimension, scalar-encoding impossibility, continuum memory, CFT/LCFT, field counts.
- Universality for every nearest-neighbour HNF torus.
- Named square-HNF occupation of the 7-vertex pair: not constructed. A general cylinder-then-glue realization lemma now places every plane two-terminal vertex-network, including this pair, in the honestly embedded rank-one category (`GENERAL_REALIZATION_LEMMA`). That is not a lattice occupation.

## Census actually completed

| family | n reached | r=1 splits |
|---|---|---|
| exhaustive plane two-terminal | 5 | 0 |
| hidden hop grid of those cores | 12 (9 combinations) | 1 (L1R1 only) |
| two-terminal series-parallel | **12** (295654 unique graphs) | **330** |
| Wheatstone+SP, multipath, grids | 12 | 0 |

SP n=11 (59324) and n=12 (212771) unique generation is included. The n=7 pair is the **smallest witness found in the declared enumerated families**. This is not a claim of global minimality among all plane two-terminal vertex-networks.

## Independence / chronology

Exact finite combinatorics. No Monte Carlo, no production block, no post-selection of descriptors. Independent of original-`U` identifiability and of the thermal/contact proof gate. Independent of the parallel-gadget successor-second-moment split: that observable was frozen out of the present summary.

The witness script recomputes the frozen experiments from hard-coded incidence. It is search-independent. It is not a second implementation of the enumerator.

## Paper impact

`NEW_REQUIRED_LEMMA` for the P1 cut-network / predictive-state manuscript: a named bounded summary strictly coarser than the cut network fails a predeclared depth-2 language. Do not promote to an all-graphs minimality theorem.

Cut-dependence verdict `CUTS_LAW_EQUIVALENT_NOT_ISOMORPHIC`. Embedding verdict `GENERAL_REALIZATION_LEMMA`. Cut-network minimality `UNRESOLVED`.
