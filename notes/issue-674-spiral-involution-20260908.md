# #674 — Spiral involution vs cardinality-only cancellation in `both-same × both-same`

**Date:** 2026-09-08 · **Machine:** local (Mac) · **Parent:** #650 · **Follows:** PR #668 · **Related:** #665, #651, #635, PRs #662

**Outcome (ticket vocabulary): `cardinality-only`, refined to config-wise cancellation — the rank-1 spiral class is self-coincident across the two colourings, so D = 0 needs no pairing and no counting.** Details below; every integer is exact.

## 0. Provenance

The joint census of PR #668 is **imported, not recomputed** (`results/issue674-spiral-involution/imported-issue668-axis-L{3,4}.json`, provenance in `imported-issue668-provenance.json`, fetched via the PR #668 file API). Integrity: the L=4 question-cell masses re-derive from the imported joints as `[120, 416, 448, 128, 8]` at k = 8..12, total 1120, matching the ticket text exactly.

Fresh compute in this note (Mac-OK, integers only):

- L=3 `2^9` stream dumping the 6 question-cell configs (`scripts/issue674_spiral_involution.py`, dump in `results/issue674-spiral-involution/axis-L3-spiral-dump.json`);
- L=4 `2^16` stream: **2.6 s wall clock** — seconds-scale, so permitted by the ticket gate; actions of 56 candidate torus maps on the 1120-config cell (`results/issue674-spiral-involution/axis-L4-spiral-actions.json`);
- L=5 (`2^25`): **not run — NEED_HUAWEI.**

## 1. The colour-flip kill at L=4 (recorded first, as instructed)

Colour-flip C ↦ C̄ preserves k only at k = N/2. It would force the L=4 cell masses palindromic under k ↔ 16−k:

| k | mass | mirror k′ = 16−k | mirror mass |
|---|------|------------------|-------------|
| 8 | 120 | 8 | 120 |
| 9 | 416 | 7 | 0 |
| 10 | 448 | 6 | 0 |
| 11 | 128 | 5 | 0 |
| 12 | 8 | 4 | 0 |

`[120, 416, 448, 128, 8]` is **not palindromic** → **colour-flip is dead as the D = 0 mechanism at L=4.** The L=3 survival (single mass at k = 6 = 9/2 … at k = 6 = N/2·(9/2 is not integer; the point is k = 6 = N − 3 and the window is the single point k = 6) was a small-L accident of the window position, not a mechanism.

Structural form of the kill (verified, L=3 and L=4): a **black** rank-1 `both-same` spiral on the primal (NN) lattice needs k ≥ 2L sites (min k: 6 at L=3, 8 at L=4 — brute-forced over the full 2^N), so its complement has k′ = N − k ≤ N − 2L < 2L for L = 3 (N = 9 < 4L = 12: always) and for all k > 8 at L=4. Concretely at L=3, the complement of each of the 6 spirals has black side label `none` (rank 0) and white side `both-same` **rank 2** — colour-flip does not even preserve rank-1 × rank-1, at any L tested. The palindromic failure is the observable shadow of this site-count obstruction.

## 2. The 6 rank-1 spirals at L=3 (dump)

All six are k = 6, and in **every one** the black side is a single NN cluster with winding basis a single diagonal vector (rank 1, `both` directions) and the white side is a single matching-lattice cluster, likewise rank-1 diagonal. Occupancy (axis coordinates (x, y), x fastest):

```
mask 011101110: black {(1,0),(2,0),(0,1),(2,1),(0,2),(1,2)}
mask 011110101: black {(0,0),(2,0),(1,1),(2,1),(0,2),(1,2)}
mask 101011110: black {(1,0),(2,0),(0,1),(1,1),(0,2),(2,2)}
mask 101110011: black {(0,0),(1,0),(1,1),(2,1),(0,2),(2,2)}
mask 110011101: black {(0,0),(2,0),(0,1),(1,1),(1,2),(2,2)}
mask 110101011: black {(0,0),(1,0),(0,1),(2,1),(1,2),(2,2)}
```

(each is a diagonal staircase winding once around both axes; the complement is the white side's staircase on the matching lattice).

## 3. Answer to the ticket question

**D = 0 on `both-same × both-same` is per-config, not a pairing and not a mere count.** In all 6 (L=3) and all 1120 (L=4) question-cell configurations, `either = 1` on **both** sides simultaneously: each configuration carries its own +1 and −1, so D(C) = 0 term-by-term with nothing to pair. The two relevant refinements:

1. **The rank-1 `both-same` class is self-coincident across colourings.** At L=3: black-side rank-1 `both-same` = 6, white-side rank-1 `both-same` = 6, and the intersection `both-same × both-same` at rank (1,1) is 6 — the two classes are **the same set of configurations**, not two equal-size sets. Same at L=4: 1120 = 1120 = 1120, with per-k masses `[120, 416, 448, 128, 8]` carried by identical configurations on both sides. So the equal-mass option of the ticket holds, but for the strongest possible reason: set identity. Conjecture (mechanism, not proven here): the complement of a primal diagonal spiral is a matching-lattice diagonal spiral and conversely — a matching-lattice self-duality of the diagonal-winding event; verified as set identity at L=3, 4 only.
2. **Fixed-point-free involutions on the cell exist anyway, but are epiphenomenal.** Testing all 56 candidate torus self-maps (identity; 8 translations; 9 rotations by 90°; 3+3 axis reflections; 9 anti-diagonal reflections at L=3; scaled to 56 at L=4): the axis reflections `x ↦ 2c−x` and `y ↦ 2c−y` (all c) act on the cell as fixed-point-free k-preserving involutions — at L=3 all six of them, with e.g. `refl-x(0)` orbits `{0,2},{1,3},{4,5}` on the dumped configs; at L=4 exactly 10 such maps (4 `refl-x`, 4 `refl-y`, translations `(0,2)`, `(2,0)`), all fixed-point-free and k-preserving on all 1120 configs. Rotations by 90° and general translations act on the cell but with orbits of size 4 / fixed points — they are involutions on the torus (90° turns are order 4) but not pairings. **None of these maps is needed for D = 0**, since cancellation is already config-wise (point 1); they are recorded as cell symmetries, not as the mechanism.

**Verdict: `cardinality-only` — sharpened: the equal cardinality is set identity; D = 0 on the question cell is config-wise.** No involution pairing of distinct configurations is required, and colour-flip (the one candidate that acts between the colourings) is obstructed at L=4 and in fact never preserves rank-1 × rank-1 at L=3 either.

## 4. Q3: does the same mechanism force `both-two = 0`?

**No — `both-two = 0` is independent, and it is a theorem, not a census accident.** On the torus, an x-wrapping cluster contains a curve in homology class (1,0) and a y-wrapping cluster one in class (0,1); these classes have algebraic intersection number 1, so two *disjoint* clusters cannot respectively wind in x and in y: any x-wrapping occupied cluster intersects any y-wrapping occupied cluster, hence merges with it into a single `both-same` cluster. Verified exhaustively:

| L | lattice | configs with some x-wrap and some y-wrap | all `both-same` | `both-two` |
|---|---------|------------------------------------------|-----------------|------------|
| 2 | primal | 5 | yes | 0 |
| 2 | matching | 7 | yes | 0 |
| 3 | primal | 97 | yes | 0 |
| 3 | matching | 265 | yes | 0 |
| 4 | primal | 10165 | yes | 0 |
| 4 | matching | 37679 | yes | 0 |

(All rows from full 2^N streams over the black mask; the same config set re-laboured for the white side is the complement image, so the table covers both colourings. The intersection argument above holds at every L and both lattices, so `both-two ≡ 0` is structural — A-continues is forced — independent of the spiral self-coincidence of §3.)

## 5. Tripwires

- Imported #668 integers re-derived: L=4 cell masses `[120, 416, 448, 128, 8]`, total 1120; L=3 cell size 6 — both match the ticket text bit-for-bit.
- Fresh streams reproduce the imported cell sizes exactly (6 at L=3, 1120 at L=4), cross-validating import and code.
- Integers only throughout; no STATUS edit; no ticket closed; draft PR against main.

## 6. Artifacts

- `scripts/issue674_spiral_involution.py` — L=3/L=4 streamer, cell finder, map-action tester
- `results/issue674-spiral-involution/imported-issue668-axis-L{3,4}.json`, `imported-issue668-provenance.json`
- `results/issue674-spiral-involution/axis-L3-spiral-dump.json` (6 configs + all 33 map actions at L=3)
- `results/issue674-spiral-involution/axis-L4-spiral-actions.json` (1120-config cell masses + all 56 map actions at L=4)
- `tests/test_issue674_spiral_involution.py` — axis L=2,3 (seconds)
