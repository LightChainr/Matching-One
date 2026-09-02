# Bounded summary search certificate (extended to contract limits)

**Verdict:** `NO_COMPRESSION_WITNESS_FOUND`

Frozen before search: summary = (S(z), n, H2, b2, radius-1 terminal-local neighborhood).
Successor-hazard moments were **not** in the summary. Arithmetic is exact `fractions.Fraction`. No Monte Carlo.

## Enumerated class

Connected plane two-terminal vertex-networks with L,R on a common face, no L-R edge. Generators: exhaustive n<=5; path-hidden copies of those cores over the **full contracted hop grid** (Lh,Rh in {1,2,3}, 9 combinations, n<=12); two-terminal series-parallel **generated to n=12** (295654 graphs, so the PROMPT's n=11 and n=12 layers are now covered); Wheatstone+SP compositions n<=12 with SP partners n<=6; multipath n<=12; 2/3/4-row ladders n<=12.

| family | size | r=1 splits | r=2 splits |
|---|---|---|---|
| exhaustive_nle5 | 7398 | 0 (closed) | 0 (closed) |
| hidden_L1R1 | 7398 | 1 (the n=7 witness) | 0 |
| hidden_L1R2 | 7398 | 0 | 0 |
| hidden_L1R3 | 7398 | 0 | 0 |
| hidden_L2R1 | 7398 | 0 | 0 |
| hidden_L2R2 | 7398 | 0 | 0 |
| hidden_L2R3 | 7398 | 0 | 0 |
| hidden_L3R1 | 7398 | 0 | 0 |
| hidden_L3R2 | 7398 | 0 | 0 |
| hidden_L3R3 | 7398 | 0 | 0 |
| series_parallel_nle12 | 295654 | 330 | — |
| wheatstone_sp_nle12 | 577 | 0 | — |
| multipath_nle12 | 271 | 0 | — |
| grids_nle12 | 10 | 0 | — |

H2 and b2 are redundant given (S,n): H2 = n − S_1, b2 = C(n,2) − S_2.

## Witness (radius 1)

- n = 7 switchable vertices, |E(A)| = |E(B)| = 10 (pendant-free; every switchable vertex lies on a simple L–R path)
- S(z) = 1 + 7z + 21z^2 + 35z^3 + 33z^4 + 15z^5 + 2z^6, coefficients `[1, 7, 21, 35, 33, 15, 2, 0]`
- H2 = 0, b2 = 0
- r=1 neighborhoods: **identical** (unique L-neighbor, unique R-neighbor, no L–R edge)
- r=2 neighborhoods: **differ** (this is an r=1 witness only)
- planar (G ∪ {L,R}): A=True, B=True; connected carrier, L–R path exists
- first frozen split: **E2_c2**
  - P(A) = 937/1050
  - P(B) = 313/350
  - gap = 1/525
- all other frozen experiments agree (E0_c1=1/1, E0_c2=1/1, E0_mix=1/1, E1_c1=1/1, E1_c2=1/1, E1_mix=1/1, E2_c1=1/1, E2_mix=33/35)

Smallest-edge variants: replace core A by the 9-edge pendant core `exh_n5_1056` or `exh_n5_413` (both in the same behaviour class) to get 9+10 edges with the same S / r=1 / E2_c2 gap 1/525. Both have a pendant switchable vertex, so the 10-edge pendant-free core `exh_n5_4313` remains the primary A.

### Graph A — parallel 4-paths + a triangle bypass

Core `exh_n5_4313` hidden by 1-hop corridors. Vertices `{0,1,2,3,4,5=ℓ,6=r}`.

```
L — 5 — 2 — 1 — 6 — R
         \ /
          3
L — 5 — 4 — 0 — 6 — R
```

Incidence: `[[0, 4], [0, 6], [1, 2], [1, 3], [1, 6], [2, 3], [2, 5], [4, 5], ['L', 5], ['R', 6]]`

The two 4-mincuts are `{5,2,1,6}` and `{5,4,0,6}`, intersecting only at the corridor ports `{5,6}`.

### Graph B — fan-in: two 4-paths share the R-adjacent core vertex

Core `exh_n5_2451` hidden by the same 1-hop corridors.

```
L — 5 — 2 — 0 — 6 — R
         \     /
          1   4
          |
          3 — 6 — R
```

Incidence: `[[0, 2], [0, 4], [0, 6], [1, 2], [1, 3], [2, 5], [3, 6], [4, 5], ['L', 5], ['R', 6]]`

The two 4-mincuts are `{5,2,0,6}` and `{5,4,0,6}`, intersecting in three vertices `{5,0,6}`.

### Separating mechanism

S_4 = 33, so C(7,4)−33 = 2 connecting 4-sets in each graph. The enumerator S(z) cannot see how those two mincuts intersect: A has disjoint interiors (share only the corridor ports), B additionally shares the R-adjacent core vertex (fan-in). E2_c2 occupies a uniform ordered 2-prefix (always safe, since every 3-set is safe), then each clone independently occupies 2 of the remaining 5 vertices; a clone dies iff the resulting 4-set is one of the two mincuts. The mean of p² depends on the intersection pattern, invisible to S(z) and to the radius-1 ball. Delayed-fork E1_c1 equals 1 on both graphs, so the split is **not** the successor-second-moment observable.

### Torus embedding

The pair lives in the already-closed planar two-terminal cut-network category (PR #491): the cut-network representation maps rank-one torus states onto this category; the two-port calculus can place either gadget as a two-terminal block. No explicit torus occupation realizing these 7-vertex networks was constructed in this search. Consistent with #491, random variables are vertices (not edges) and all sampling is fixed-cardinality without replacement.

## Corroboration inside series-parallel graphs

Two-terminal series-parallel graphs generated to n=12 produce **330** r=1 summary classes that split, all on a depth-2 experiment. The smallest SP split is at n=8, |E|=12, gap 2/1575. The layer-by-layer split counts and the smallest per-layer gaps are recorded in the machine JSON (`enumerated_class.SP_nle12`). The pairs are typically a corridor in series with two different same-S SP cores — the same hiding mechanism already inside the series-parallel subcategory. Wheatstone, multipath and grid families produced no r=1 split.

## Relationship to the repository's other no-go results

- PR #549 (parallel-gadget amplification) splits an identical complete-survival class by **successor-hazard second moments** (q_A=29 vs q_B=25). That observable is **frozen out** of the present summary, and on this witness the delayed-fork E1_c1 agrees (both 1). The present split is therefore independent of #549's mechanism.
- Issue #435/#434 (N16/N17 torus configurations) split on delayed branching inside the full survival law. The present pair lives in the planar two-terminal category and splits on a *different* depth-2 experiment while fixing radius-1 neighborhoods.
- The certified pair and all SP corroborations respect the #491 sampling contract: vertex randomness, fixed-cardinality sampling without replacement, no edge-reliability interpretation.

## Manuscript-ready theorem (only the class actually tested)

> **Theorem (r=1 non-compression on a 7-vertex pair).** There exist two connected plane two-terminal vertex-networks G_A, G_B, each with 7 switchable vertices and with terminals on a common face, such that the complete safe-subset polynomials agree, the singleton and pair trigger counts agree, and the radius-1 terminal-local rooted neighborhoods are isomorphic as typed graphs, but P(E2_c2; G_A) = 937/1050 ≠ 313/350 = P(E2_c2; G_B), where E2_c2 is the frozen experiment 'shared prefix of two distinct uniform vertices, then two independent 2-step continuations, observe terminal disconnection in both clones'. In particular the tuple (S(z), n, H2, b2, radius-1 neighborhood) is **not** a sufficient statistic for the frozen depth-2 compositional language, already inside the 7-vertex planar two-terminal category. This does **not** assert failure of the radius-2 neighborhood, nor a lower bound on Euclidean latent dimension, nor any continuum/CFT statement.

## What this does not show

- The same pair is separated by radius-2 neighborhoods, so it is not an r=2 witness.
- Lengthening both corridors to 2 hops equalizes r=1 and r=2 but kills the E2_c2 gap on this S-class. **All nine** contracted hop combinations of the exhaustive cores are closed at radius 2 (0 splits over 66582 graphs).
- Exhaustive plane two-terminal graphs with n<=5 are closed under the frozen summary at both r=1 and r=2.
- Delayed-fork E1_c1 agrees on the witness (both equal 1). The split is a genuine depth-2 experiment.
- The SP corroborations reach n=12; none is smaller than n=8, so the n=7 exhaustive-core witness remains minimal in n.
- This is not an all-graphs theorem and not a proof that the cut network is a minimal sufficient statistic.

## Reproducibility

- Library: `research/summary_search/bounded_summary_search.py`
- Generators/hunts: `sp_gen12.py` (layered SP to n=12), `sp12_hunt.py` (per-layer hunt; per-layer = global because n is part of the summary key), `hid9_hunt.py` (full 9-combo HID sweep), `run_full_search.py` / `finish_remaining.py` (reference)
- Independent check: `python3 research/summary_search/verify_witness.py` (hardcoded incidence; expected `VERIFY_OK`)
- Machine JSON: `artifacts/bounded_summary_search.json`
