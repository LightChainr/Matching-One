# Literature officer — verification of the homological-balance root ledger

**Date:** 2026-09-06
**Role:** theory input + independent rerun. Same epistemic level as [#602](https://github.com/LightChainr/Matching-One/pull/602). **Does not enter** `docs/STATUS.md`. **Does not close** #276, #321, #337, #275, #537.
**Source:** owner-supplied pack `matching-one-research.zip`. Submitted as a **separate** PR from #602 / #603 / #604.
**Ledger:** [`notes/homological-balance-root-ledger-20260906.md`](homological-balance-root-ledger-20260906.md)

One-line: the L=3,4 census is real; the F1/F2 split of `L^{-4}` is already the #321 identity; qualitative convergence is a reconstruction of the still-open #276 target, not a close.

---

## Independently rerun

Stack: Python 3 stdlib. Two winding algorithms, no sampling.

| Check | Result |
|---|---|
| L=3, 512 configs, `r_b+r_w=2` | **0 failures.** Rank pairs `(0,2):259, (1,1):162, (2,0):91`. |
| L=4, 65536 configs | **0 failures.** `(0,2):36559, (1,1):19932, (2,0):9045`. |
| Algorithm A vs B (tree-cycle vs BFS lift) | **0 mismatches** at L=3 and L=4. |
| Manual configs (empty, full, point, h-loop, v-loop, diagonal, cross) | all `r_b+r_w=2`, ranks as claimed. |
| `M_L(0)=-1`, `M_L(1)=+1`, sampled monotone | exact. |
| `p_L^H` | L=3 **0.586511455113**, L=4 **0.590672112331**. |
| `M_L(1/2)` | L=3 **`-21/64`**, L=4 **`-13757/32768`**. Square-site is not self-matching. |
| `M_L(0.5927460)` | L=3 **+0.02497**, L=4 **+0.01033**. Sign agrees with `p_L^H < p_c`. Not an `L^{-13/4}` fit (ledger already says L=3,4 are too small). |

JSON: `results/homological-balance-exact-torus/latest.json`.

---

## Cite vs claim

| Ledger statement | Verdict |
|---|---|
| I1 `r_b+r_w=2` on honest tori | **CITE.** Already proved: PR #271, [`notes/digital-alexander-duality-proof.md`](digital-alexander-duality-proof.md). The census is a finite **check**, not a new proof, and does not replace the topological argument for general L. |
| I2–I9 matching function / source algebra / clocks | **CITE** #269 / #337 / #275. Reconstruction of existing exact language. |
| Qualitative `p_L^H\to p_c` | **#276's still-open target.** Architecture (subcritical decay + matching duality + monotonicity) is the one #276 already named, with Duncan–Kahle–Schweinhart arXiv:2011.11903 as template. A reconstructed sketch does **not** close #276. |
| Self-matching `⟹ p_L^H=p_c=1/2` exactly, all L | **Finite identity from I1 + complement involution.** Square-bond analog is already [`notes/square-bond-duality-tiny-torus.md`](square-bond-duality-tiny-torus.md) (`E[D]=0` at `p=1/2`). The site/homological-balance form is the right named lemma for this observable; not a percolation-threshold discovery (Sykes–Essam / Kesten matching already gives `p_c=1/2`). |
| Abstract two-scale theorem E | Taylor + implicit function, **conditional** on (H1,H2,H3). Clean statement of what F1 and F2 must supply. Not a square-site theorem. |
| `L^{-4}=(13/4)+(3/4)` | **Already #321.** `M_L(p_c)\sim L^{-13/4}` and `M_L'(p_c)\sim L^{3/4}` are the empirical/CFT inputs named there. Ledger correctly splits them as F1 (universal 4-arm / `ν=4/3`) and F2 (matching-odd residual). |
| F1 open on square-site | **Correct.** Smirnov/LSW is triangular (and FK-Ising). Do not import `α_4≍L^{-5/4}` onto square-site. |
| F2 = `Q_4ε`, `x=21/4` | **CONJ**, as the ledger says. Deeper than F1. |
| Obstacle J: no model with F1 proved **and** F2 nontrivial | **Landscape observation, not a no-go theorem.** Currently true of the rigorously conformal list (triangular site, FK-Ising: self-dual, F2≡0; square-site: F2 live, CI open; square-bond: self-dual so F2≡0, CI also open). It does **not** prove square-site conformal invariance can never be proved. |

---

## What this does not do

- Not a STATUS edit.
- Not a close of #276 or #321.
- Not a proof of `p_L^H-p_c\sim L^{-4}` on square-site.
- Not an identification of L=3,4 `M_L(p_c)` with the `L^{-13/4}` amplitude.
- Not a competing digital-Alexander proof; #271 remains the proof.
