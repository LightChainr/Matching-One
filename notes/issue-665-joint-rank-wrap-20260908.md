# #665 joint dictionary — (k, r_b, r_w) x wrapping-type x 5-names, axis L=3,4

**Ticket:** LightChainr/Matching-One#665 (parent #650). **Date:** 2026-09-08.
**Mode:** new census at axis L=3 (`2^9`) and L=4 (`2^16`) only, seconds-scale. Axis
L=5 (`2^25`) is **NEED_HUAWEI** and was not run on the Mac. Integers only.
No `docs/STATUS.md` edit; no ticket closed or merged.

## What was computed

The two exact views of the same site configurations had never been tabulated
jointly. This census streams every configuration on the axis torus and records

```text
(k, r_black, r_white, wrap_black, wrap_white, label5_black, label5_white)
```

with

- `(r_b, r_w)` — ambient H_1 winding ranks of the black NN (primal) graph and
  the white NN+NNN (matching) graph, as in PR #653's repaired lab
  (`scripts/probe_invariant_shape/exact_rank_census.py`);
- coarse wrap labels `neither / dir0 / dir1 / both` via
  `torus_homology.wrapping_channels` (PR #653's 4×4 classifier, where `both`
  lumps rank-2 cross with rank-1 spirals);
- #646 5-names `none / x / y / both-same / both-two` via
  `scripts/probe635_sector_decomposition.py` (`both-same` = a *single* cluster
  wrapping both axes, rank-1 spiral or rank-2 cross alike; `both-two` = the
  two-distinct-clusters case).

Script: `scripts/issue665_joint_rank_wrap_census.py`. Exact integer joints:
`results/issue665-joint-rank-wrap/axis-L3.json`, `axis-L4.json`.
Test: `tests/test_issue665_joint_rank_wrap.py` (axis L=2,3; seconds).

## Tripwires — both pass at L=3 and L=4

1. **Bernstein.** Collapsing the coarse 4×4 with
   `D = 1{primal either} − 1{matching either}` reproduces the committed
   Bernstein integers bit-for-bit:
   - L=3: `[-1, -9, -36, -78, -90, -36, 36, 36, 9, 1]`;
   - L=4: the committed rung (matching the `collapsed_D` field in the JSON).
2. **`r_b + r_w = 2`** on every configuration (PR #653), zero violations at
   both sizes. Rank-pair totals reproduce PR #653 exactly:
   - L=3: `(0,2)=259, (1,1)=162, (2,0)=91`;
   - L=4: `(0,2)=36559, (1,1)=19932, (2,0)=9045`.

## The question: is `both-same` rank-2, or does `both-same × both-same` contain spirals?

**Answer: `both-same` is NOT the same as rank 2, and `both-same × both-same`
is pure spiral — rank-1 × rank-1, never rank-2.** The exact integers:

- The 5-name `both-same` splits by rank (black side; white side is the
  transpose, as expected from the L=3,4 cell support):

  | L | `both-same`, rank 1 (spiral) | `both-same`, rank 2 (cross) |
  |---|---|---|
  | 3 | 6 | 91 |
  | 4 | 1120 | 9045 |

  So `both-same = spirals ⊎ crosses`: the label is a *cluster-level* statement
  (one cluster, both axes) and a rank-1 spiral cluster satisfies it. The coarse
  `both` label also contains both, so neither label alone is the rank-2 class;
  the rank-2 (cross) class is exactly `both-same` with component rank 2 — the
  `cross` bit of `torus_homology.ComponentHomology`, which neither committed
  table carries as a column.

- The cell `both-same × both-same` (the "spiral pairing" cell of #660/#662)
  is **entirely rank-1 × rank-1**:

  | L | `both-same × both-same` by `(r_b, r_w)` |
  |---|---|
  | 3 | `(1,1): 6` (at k=6 only) |
  | 4 | `(1,1): 1120` (k=8: 120, k=9: 416, k=10: 448, k=11: 128, k=12: 8) |

  The coarse `both × both` cell coincides with it at both sizes — there is no
  rank-2 × rank-2 configuration at L=3,4 (a rank-2 black graph forces
  `r_w = 0`, so rank-2 on both sides is impossible given `r_b + r_w = 2`).
  This is the structural reason `both × both` carries `D = 0`: both sides are
  rank-1 spirals there, `either` is 1 on both sides, and no cross is involved.

- **Rank-2 lives exclusively in the D-carrying exclusive-cross cells.** Every
  configuration with `wrap_black = neither, wrap_white = both` has
  `(r_b, r_w) = (0, 2)` (white cross-wraps: rank 2); every configuration with
  `wrap_black = both, wrap_white = neither` has `(r_b, r_w) = (2, 0)`. L=3
  counts: 259 and 91; L=4: 36559 and 9045. No D-cell ever mixes ranks.

- Corollary for the #635/#660 identity: the both-same amplitude rows of PR
  #646's decomposition are **not** rank-pure. `A_black^{both-same}` contains
  the rank-2 cross mass (D = +1 side) *and* rank-1 spiral mass (D = 0); the
  spiral 6 at axis L=3 k=6 rides inside `both × both` with D = 0 and cancels
  out of `M` — consistent with the Mertens–Ziff pairing (spiraling counts
  match), but the sector label alone does not separate them. Any future
  sector-level transfer-matrix reading (#636, not funded) must split
  `both-same` by rank or it is not measuring one homology class.

- `both-two` has zero mass at every k at L=3,4 (consistent with #651's
  A-continues through axis L=5), and `dir0 × dir0 = dir1 × dir1` at every k
  (L=3: 78 each; L=4: 9406 each), reproducing PR #653's published support.

## Cell support (all configurations, both L)

```text
neither × both    (0,2)   D = −1     259 (L=3) / 36559 (L=4)
both    × neither (2,0)   D = +1      91 (L=3) /  9045 (L=4)
dir0    × dir0    (1,1)   D =  0      78 (L=3) /  9406 (L=4)
dir1    × dir1    (1,1)   D =  0      78 (L=3) /  9406 (L=4)
both    × both    (1,1)   D =  0       6 (L=3) /  1120 (L=4)   ← spirals
```

No other cell of the coarse 4×4 is populated at axis L=3 or L=4; within each
populated coarse cell the 5-name pair is constant (`none×both-same`,
`both-same×none`, `x×x`, `y×y`, `both-same×both-same` respectively). The five
Mertens–Ziff cells are therefore exactly rank-diagonal at these sizes:
D-carrying cells are the rank-exhausting `(0,2)`/`(2,0)` mismatches, and the
D-neutral cells are all rank-1 × rank-1.

## What L=3,4 cannot settle

Whether rank-2 × rank-2 configurations first appear somewhere (they cannot at
any L with `r_b + r_w = 2` on site — that identity forces at most one rank-2
side), whether a sixth coarse cell appears, and whether `both-two` stays
empty at larger L. Axis L=5 (`2^25`) is NEED_HUAWEI; the #657 run already
covers the 5-name and coarse views there, but the joint rank-splitting of
this ticket has not been run at L=5.

**Tripwire standing rule (unchanged):** any claimed decomposition must
reproduce `[-1, -9, -36, -78, -90, -36, 36, 36, 9, 1]` at axis L=3
bit-for-bit, integers / Fraction only, or it is unused.

Related: #635, #640, #646, #650, #651, #653, #657, #660, #662.
