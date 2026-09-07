# Verification — #594 Q3 bounded-task-rank no-go pack

Submitted first as [#615](https://github.com/LightChainr/Matching-One/pull/615) from `matching-one-theory-analysis-20260907.zip`. Independent rerun after submit. Does **not** close #594 / #580 / #588 / #609. Does **not** enter `docs/STATUS.md`.

## Cite vs claim

| Pack statement | Status after rerun |
|---|---|
| Theorem 1: two analytic families, identical finite-horizon response, `p_c = 1/2` vs `1/3`, `sup r_ex = 2` | **Claim, verified.** Direct-sum witness. |
| Hidden gap closed form `a+b−2√(ab) cos(π/L)` | **Claim, verified.** max \|num − formula\| = `4.21e-14` on the script grid. |
| `γ_∞ = 4 sinh²((p−p_c)/2)` | **Cite of algebra.** Follows from `a=e^{p−p_c}`, `b=e^{-(p−p_c)}`. Not a separate numerical check. |
| Kalman rank = 2 independent of `L` | **Claim, verified.** `L=4,16,64`. |
| Irreducible slowest-mode mask: generic `~0.484`, `C v=0` / `wᵀB=0` at `~1e-17` | **Claim, verified.** Numbers match the pack JSON; script key names were longer (`masked_readout_Cv0` vs `masked_readout`). JSON on this branch is the rerun. |
| P398 `r_balanced ≈ 3–4` is I/O compressibility, not a critical-state dimension | **Cite of the note.** Follows from Theorem 1 + earlier probe ranks; this pack does not recompute P398. |
| Square-site: CSC necessary; signed-observable route (D) is independent of rank | **Cite of the note.** Argument, not a construction of a square-site affine generator (that remains #594 Q3.1). |
| #609 55% second-difference residual: readout-control design | **Cite of a protocol.** Not a recompute of `A(N)~N^{-0.970}`. Costs no new Monte Carlo. |

## What the scripts actually ran

```
python3 scripts/theory/no_go_theorem.py
python3 scripts/theory/irreducible_masking.py
```

- `p ∈ {0.3, 0.5, 0.7}`: `fam1 == fam2 == vis` at all four sampled times.
- Unbiased gap ~ `1/L²`; biased gap saturates (~0.268 at `L=256`).
- Fine-tuned mask is **not** a robust hidden sector; the pack already says so.

## What this does not do

- Does not answer #594 Q1 (lumping subvariety) or Q2 (memory degree).
- Does not construct the square-site affine generator asked in Q3.1.
- Does not close the no-go for a *generic irreducible* chain except in the weak-coupling / balanced-order sense the pack states.
- Does not replace the N=725 production on #609.
