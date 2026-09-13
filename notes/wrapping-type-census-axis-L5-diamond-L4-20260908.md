# #640 / #651 wrapping-type census: axis L=5 and diamond L=4

CLI job on `XPk2PZ` hit `--max-turns 120` without a PR. The enumerations had already finished. Grok downloaded the artifacts, checked the Bernstein tripwire against PR #649, and records them here. No `docs/STATUS.md`. Tickets left open.

Enumerator: `scripts/wrap_census.cpp`, instrumented from PR #649 `exact_matching_ladder_bruteforce.cpp` (kernels K1/K2). Labels as in `scripts/probe635_sector_decomposition.py`: `none / x / y / both-same / both-two`. Coarse 4×4 folds `both-two` into `both`.

Machine: Huawei DevEnvC `XPk2PZ`, 14 threads, BLAS/OMP 1. Persistent dir `/workspace/wrapping-640-651/`.

---

## Axis L=5 (`N=25`, `2^25 = 33 554 432`)

| | |
|---|---|
| configs visited | 33 554 432 |
| wall | 2.209 s (K1; two kernels agree on collapsed D) |
| tripwire | **PASS** vs `rungs.axis_L5.k1_counts` on PR #649 |
| binomial row sums | \(\binom{25}{k}\) at every k |
| 5×5 reconstructs collapsed D | yes |

Support, all k:

| black × white | count | D |
|---|---:|---:|
| none × both-same | 20 218 851 | −1 |
| both-same × none | 3 685 661 | +1 |
| x × x | 4 482 280 | 0 |
| y × y | 4 482 280 | 0 |
| both-same × both-same | 685 360 | 0 |
| any cell with both-two | **0** | — |

Coarse 4×4 is the same five Mertens–Ziff cells as L=3,4. No sixth cell. `dir0=dir1` at every k.

**#651 verdict: A-continues.** `both-two` has zero mass at every k. `x`/`y` have mass and zero D (they pair). Degeneracy is not a small-L accident on the last independently checkable axis rung.

---

## Diamond L=4 (`N=32`, `2^32 = 4 294 967 296`)

| | K1 | K2 |
|---|---:|---:|
| configs visited | 4 294 967 296 | 4 294 967 296 |
| wall | 323.827 s | 75.363 s |
| tripwire vs PR #649 `diamond_L4` | PASS | PASS |
| 5×5 tables K1 vs K2 | **identical** | **identical** |

Support, all k (both kernels):

| black × white | count | D |
|---|---:|---:|
| none × both-same | 2 739 595 215 | −1 |
| both-same × none | 408 426 325 | +1 |
| x × x | 520 535 376 | 0 |
| y × y | 520 535 376 | 0 |
| both-same × both-same | 105 875 004 | 0 |
| any cell with both-two | **0** | — |

Same five-cell coarse 4×4. MZ-forbidden cells: none.

**#651 on this geometry: A-continues**, same sentence.

---

## What this does not do

- Does not prove degeneracy at every L.
- Does not fund #636 by itself. Sector labels #636 would carry remain the degenerate both-same pair (black NN vs white NN+NNN).
- Does not identify those finite-L sectors with Jacobsen’s transfer-matrix eigenvalues (#637 / remainder of #635).

Artifacts: `results/wrapping-type-census/axis-L5.json`, `diamond-L4-k1.json`, `diamond-L4-k2.json`.
