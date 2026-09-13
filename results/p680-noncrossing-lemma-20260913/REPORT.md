# REPORT — p680: NN+NNN strip frontier-partition lemma (noncrossing) — small-width control

**Ticket:** LightChainr/Matching-One#680. **Date:** 2026-09-13. **Outcome:**
`proof` for the noncrossing claim (all w, all R — see
`notes/p680-nn-nnn-noncrossing-lemma-20260913.md`); class-set equality (NN vs
NN+NNN) re-verified at w = 2..6, R = 1..4 in the honest convention; #675's 28/66
counts reconciled as the delayed-closing / full-frontier convention (NN side), where
the two rules' class sets actually differ.

## 1. Honest-convention main table (partitions of occupied frontier sites, all bonds)

Equality column: NN class set == NN+NNN class set. bf = brute-force per-configuration
cross-check (run wherever w*R <= 16).

```
w=2 R=1: NN=   4 NNN=   4 sets_equal=True  bf=True/True
w=2 R=2: NN=   4 NNN=   4 sets_equal=True  bf=True/True
w=2 R=3: NN=   4 NNN=   4 sets_equal=True  bf=True/True
w=2 R=4: NN=   4 NNN=   4 sets_equal=True  bf=True/True
w=3 R=1: NN=   8 NNN=   8 sets_equal=True  bf=True/True
w=3 R=2: NN=   9 NNN=   9 sets_equal=True  bf=True/True
w=3 R=3: NN=   9 NNN=   9 sets_equal=True  bf=True/True
w=3 R=4: NN=   9 NNN=   9 sets_equal=True  bf=True/True
w=4 R=1: NN=  16 NNN=  16 sets_equal=True  bf=True/True
w=4 R=2: NN=  21 NNN=  21 sets_equal=True  bf=True/True
w=4 R=3: NN=  21 NNN=  21 sets_equal=True  bf=True/True
w=4 R=4: NN=  21 NNN=  21 sets_equal=True  bf=True/True
w=5 R=1: NN=  32 NNN=  32 sets_equal=True  bf=True/True
w=5 R=2: NN=  50 NNN=  50 sets_equal=True  bf=True/True
w=5 R=3: NN=  51 NNN=  51 sets_equal=True  bf=True/True
w=5 R=4: NN=  51 NNN=  51 sets_equal=True  bf=None/None   (2^20 configs, BFS only)
w=6 R=1: NN=  64 NNN=  64 sets_equal=True  bf=True/True
w=6 R=2: NN= 120 NNN= 120 sets_equal=True  bf=True/True
w=6 R=3: NN= 127 NNN= 127 sets_equal=True  bf=None/None   (2^18 configs, BFS only)
w=6 R=4: NN= 127 NNN= 127 sets_equal=True  bf=None/None   (2^24 configs, BFS only)
```

Class-count sequences (w = 2..6):

- R = 1: NN [4, 8, 16, 32, 64]  | NN+NNN [4, 8, 16, 32, 64]
- R = 2: NN [4, 9, 21, 50, 120] | NN+NNN [4, 9, 21, 50, 120]
- R = 3: NN [4, 9, 21, 51, 127] | NN+NNN [4, 9, 21, 51, 127]
- R = 4: NN [4, 9, 21, 51, 127] | NN+NNN [4, 9, 21, 51, 127]

All classes noncrossing (linear test) in every cell, both rules. Saturation from
R = 3 at every tested w. Note 51 > C_5 and 127 > C_6: with vacancies the honest
ambient is larger than C_w (vacancy-leap blocks); the Catalan comparison of #675
belongs to the convention below.

## 2. Delayed-closing / full-frontier convention (reconciles #675's 28/66)

Row-by-row build, previous row closed before the new row is added, final row's own
horizontal bonds pending, final row fully occupied (classes = partitions of [w]).
BFS cross-validated against direct per-configuration enumeration at
(3,3),(4,3),(4,4),(5,3),(5,4),(6,3),(6,4) for both rules.

```
w=5 R=3: NN= 28 NNN= 13   w=5 R=4: NN= 28 NNN= 13   (28 = #675's "w=5 over 4 rows")
w=6 R=3: NN= 65 NNN= 24   w=6 R=4: NN= 66 NNN= 24   (66 at R=4; "3 rows" gives 65)
```

Under this convention the NN and NN+NNN class sets are NOT equal (differences in
both directions; e.g. NN-only class 0|1|0 at (3,3), NN+NNN-only class 0011 at (4,4)),
though every class of both rules is noncrossing. The equality claim of Lemma A holds
in the honest convention (section 1), which is what the lemma's wording describes.

## 3. Converse-false witness (cluster level)

Occupied set {(1,0),(0,1),(1,2)} on the 2x3 rectangle: king -> one cluster with
frontier footprint {0,2} (position 1 vacant, diagonal relay); NN -> footprints
{0},{2}. NN cannot make that footprint with the same three sites; the *partition*
{0,2} is still NN-reachable at depth >= 2 via a longer relay, which is why the class
sets coincide while the cluster-level converse stays false.

## 4. What was NOT done

No width above 6 enumerated (per ticket). No weighted counts or realizations-per-state
compared. No Huawei budget used; no #636 started. Full repository CI has not been run
for this commit.
