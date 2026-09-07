# Literature officer — verification of the #599 probe pack

**Date:** 2026-09-06
**Role:** theory input + independent rerun. Same epistemic level as [#602](https://github.com/LightChainr/Matching-One/pull/602). **Does not enter** `docs/STATUS.md`. **Does not close** #599, #598, #600, #601, #593, #584, #582.
**Source:** owner-supplied pack `matching-one-round1-20260906.zip` (bundle `3ac6aae`, branch `research/probe-full-round1-20260906`). Submitted here as a **separate** PR from #602 and from #603.
**Rerun stack:** numpy 2.2.6, scipy 1.15.3. No sampling, no new width, no fit.

One-line: the pack is a #599 atlas, not a close of Phase C. Independently rerun numbers match. Several "theorems" are named facts (Schur / Kalman–Wonham / Callan–Smiley); the durable new objects are the affine-rank lemma on the #549 closed form, the language ladder, and the synthetic no-go.

---

## Independently rerun

| script | result |
|---|---|
| `verify_gate1_independent.py` | **GATE1_VERDICT: PASS.** `r_positive = 10,26,76,232,750`. Unique D0-preserving reflection is the cut. `halves_linked` odd only at w=5,7 (defects 0.306, 0.250), matching #603's teeth. |
| `phase_c_quotient_factor.py` | full vs orbit-weighted quotient rel. diff **2.5e-16 … 4.2e-16** at w=4..8. Plain (unweighted) product **0.61, 0.47, 0.37, 0.40, 0.44**. Odd readout from even sources ~1e-18; from a marked odd source 0.029 (w=5), 0.0063 (w=7). Committed JSON matches this rerun to ~1e-14. |
| `c5_l2pi_fibre.py` | L2(π) full vs quotient **1e-16**. C3 W-fibre vs full **1.3e-15**. |
| `group_selection_demo.py` | triv/triv/W integrand **2.6e-17**; other rows 1e-2. Second-order slope 2.006 vs triv-perturbation control 1.022. |
| `rank_vs_classes_demo.py` | affine family rank = 2 at k=1,4,9,24. `F_{k,a+1}-F_{k,a}` constant; k=1 increment **1.020408e-02 = 1/98**. |
| `c2_root_rank.py` | #549 closed form = successor-sum protocol, k=1..12 **PASS**. Abstract d-root: deg(P)=d; tests d=0..k give rank k+1; `h_A=h_B` degenerates to rank 1. |
| `c7_no_go.py` | even-task data independent of the hidden coupling to **1e-16**; full-chain gap moves. |
| `c10_cheap_tests.py` | T1 dynamic resolutions at w=8: 750 and 750. T2 even-invariance True. |

One numerical caveat, not a lemma failure: Vandermonde at k=24 reports rank 18, not 25. A 25×25 Vandermonde on `linspace(0.5,2.0,25)` is float64-ill-conditioned. The k=1,4,9 rows are exact.

One packaging fix on submit: `phase_c_quotient_factor.py` wrote to `scripts/results/` (`parent.parent`). Aligned to `parents[2]` with the other probe scripts so a repo checkout writes `results/probe-p398-quotient-balanced/`. Committed JSON is the pack's original.

---

## Cite vs claim (against #601 / #603)

| Probe claim | Verdict |
|---|---|
| Phase C: full space ≡ R-orbit quotient to 1e-16 | **Independent confirmation of [#603](https://github.com/LightChainr/Matching-One/pull/603).** Not a close of #598/#600. Spectra are **not** comparable to #603 (counting-measure Gauss–Legendre vs frozen Hankel). Shared claim is the factorization identity. Expected theorem is Kalman / Wonham, as #600 said and as the pack's own C1 snapshot already relabelled. |
| Gate 1 orbit = coarsest D0-admissible lumping | **Instance at five widths**, not a theorem (#601 Q2; D'Angeli–Donno 2013 Prop. 13). Pack already carries this caveat in snapshot v2. |
| Fixed-NC count `C(w, ⌊w/2⌋)` / orbit formula | **CITE** Callan–Smiley Thm 1 + Burnside. OEIS A007123. Pack snapshot v2 already cites this. |
| Finite-group response selection (Thm 1a/1b) | Representation-theoretic vanishing is **Schur / Wigner–Eckart**. Diaconis 1988 Ch. 3E is the group-walk case. Pack C5 correctly claims only the Markov-generator / Duhamel-pointwise form + order-`ℓ` tensor criterion + C2-accident. That is the #601 Q3 language slot, not a new algebraic fact. |
| "even = orbit-constant is a C2 accident" | Standard: `V_triv = Fun(X/K)` for any finite K; other isotypics are fibres. Useful **correction of language** for this repo, not a discovery. |
| Affine `F_{k,a}` ⇒ rank ≤ 2 while classes = k+1 | **Ours as a #549 protocol fact.** The closed form is affine in `a`; increment `1/[2k(8k-1)^2]` reproduces #435's 1/98 at k=1. This is the sharpest new exact object in the pack. |
| C7 even-task no-go | Exact on a **synthetic two-copy C2 family**. Not a percolation-threshold theorem and not a constraint on `p_c`. Keep the claim boundary. |
| `r_linear` ≠ task block-Hankel | Honest negative. Stays blocked on #593's defining matrix. Do not quote 10, 26, 72 as McMillan degrees. |

---

## Position relative to later literature-officer notes

Snapshot v2 absorbed #600/#601/#602/#603 as of 2026-09-06 afternoon. After that, [#602](https://github.com/LightChainr/Matching-One/pull/602) `5f0bf1b` settled the Gate 2 leftover: **do not retarget #584 into `Sp`**. Probe C6 is still correctly deferred (needs archived `Q_N`). If C6 later runs, primary coordinate stays the frozen #582 contract.

---

## What this does not do

- Not a STATUS edit.
- Not a close of #599, #598, #600, #601, #593, #584, #582.
- Not a competing Phase C of #598; #603 remains the frozen-convention record.
- Not a start of the #584 screen or of a w=9 memory run.
- Not an identification of the C7 toy with a site-percolation threshold.
- Not a merge.
