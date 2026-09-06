# Erratum — `π({1,0})` is not `R_h`

Caught on check of the synthesis (2026-09-05 evening). **Does not enter** `docs/STATUS.md`.

`notes/literature-officer-20260905-issue576-wrapping.md` was already correct. Three later files were not.

## The two numbers

At square aspect, criticality:

| object | value | what it is |
|---|---:|---|
| Pinson / Arguin `π({1,0})` = `P̂((1,0),≥1)` | **0.169 415 435 321** | rank-1 homology class `{1,0}` |
| Newman–Ziff `R_1` | **0.169 415 435** | wraps around one specified axis **but not** the other |
| Newman–Ziff `R_h` | **0.521 058 290** | wraps in a specified direction, **including** simultaneous wrap in the other |
| Akhunzhanov torus polynomials | → `R_h` | specified-direction wrap, `L≤12` except corrupt `L=10` |

`R_h = R_e − R_1`, `R_b = R_e − 2 R_1` hold on the four Newman–Ziff digits. `π({1,0})(i)` matches `R_1` to the 12 digits printed in the #576 note; it does **not** match `R_h`. They are different observables. Scoring `R_h` against the `π({1,0})` table is the mix-up this erratum exists to stop.

Conceptual caveat: `R_1` is “wraps x not y” (any horizontal winding, no vertical). `π({1,0})` is the class `{1,0}` specifically. At `r=1` they agree to 12 digits; do not assume that at `r=2,4` without computing both.

## What was wrong where

- `literature-officer-20260905-synthesis.md` — said the #576 `π({1,0})(r=1)` **is** `R_h`. False. Rewritten in `ff1a23c`.
- `literature-officer-20260905-newman-annulus.md` — same sentence (`π({1,0})(r=1) = 0.521058290` is `R_h`). False. The first erratum commit claimed this rewrite and did not do it. Rewritten on the second check (`b1d251e`).
- `literature-officer-20260905-homology-giant.md` — wrote `π({1,0})` and `R_h` as if they were the same A-specified number `0.521`. False. Same overclaim; rewritten on the second check (`602b91d`).

The #576 written choice stands: matching-odd is a **non-claim** vs Pinson wrapping; `pinson_pi10_ratio` (`2.969` at `r=2`, `5.052` at `r=4`) is the named competitor; keep `R_h` (→ `0.521`) and `π({1,0})` (→ `0.169`) apart.
