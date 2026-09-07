# Literature officer — Gate 2 must not pick Gate 3's coordinate

**Date:** 2026-09-06
**Role:** theory input. Same epistemic level as [#602](https://github.com/LightChainr/Matching-One/pull/602). **Does not enter** `docs/STATUS.md`. **Does not close** #584, #582, #596, #595, #600. **Does not start** the #584 screen.
**Question:** #596 left a constraint: if Gate 3 residual orientation is computed in the historical channel while Gate 2 says that channel is 16× more expensive than the best available, the two gates will disagree for a reason unrelated to the labels. Settle it before running #584.
**Answer, in one line:** they are different functionals; do not retarget. The 16× is a *purchase* multiplier for `(A4, δ)`, not a reason to move the #582 remainder.

---

## The two objects

Gate 2 (`d5c7dad5`, `notes/p205-projective-channel-design-20260906.md`; script `scripts/p205_projective_channel_design.py` at base `390118b3`) scores histogram channels by the sample multiplier for an α-level Fieller set of `(A4, δ)` to sit inside a declared `|A8_fit/A4_fit|` band. Large-`n` figure of merit:

```text
A4² / var(δ)
```

Ranking, identical at two independent sizes: `combination < Sp < S < M < Dp`. Against historical `M @ p_ref`, best readout saves **16.4× at N=325 and 3.7× at N=425**. Price does not transport (`GATE2_NOT_ESTABLISHED`). Recommendation if a calibration is bought anyway: **`Sp`**, not the combination.

Gate 3 (#584) types the **frozen #582 remainder** — the affine-orthogonal Wasserstein shape residual `r_N(u)` after projecting `span{1, Q_N}` out of a quantile displacement of a declared threshold CDF. #582's default contract is

```text
F_N(p) = [1 + M_N(p)] / 2
```

unless a different typed CDF is declared *before* the flow is computed. #582 already says: do not mix different observable/channel semantics into one flow. #584 already says: freeze `g_N` and `r_N`; do not recompute a PCA basis inside each metadata class.

These are not the same linear functional of a histogram. Gate 2 v1 already taught that whitening by `var(A4)` is the wrong matrix for a `(A4, δ)` decision. Using `Sp` for Gate 3 because it is cheap for `(A4, δ)` is that mistake one level up.

---

## Why the 16× cannot fight the screen

1. **The screen is already paid for.** #584's first pass uses the five committed transitions. No new samples. The 16× is a multiplier on archived 10M/node for a *future* Smith/angular calibration. It does not apply to residual orientation on existing full-law profiles.

2. **Gate 2 refused to "use prospectively".** Ranking transported; price did not. Branch A was rejected because a channel cannot be used at a cost uncertain by two orders of magnitude. An unpriced ranking is not a declared coordinate for a different experiment.

3. **Retargeting after seeing Gate 2 is a post-hoc channel declaration.** The object being typed is the #582 remainder. Changing the CDF after a different gate selected a channel is the non-independence error: selection criterion `C` and subsequent test of labels on the selected coordinate are not independent (Kriegeskorte, Simmons, Bellgowan, Baker, *Nat. Neurosci.* 2009). On this repo that is GOVERNANCE **§2C** (don't misdate a freeze) plus **§2E** (roots, slopes, derivatives, quantiles from the same histograms are one block). Even if the five transitions are a different production block from #205, the *channel* was chosen on #205 for a different functional.

4. **Observer bandwidth is not physical rank** (#419, and #582's own pointer to it). A label that appears only in `M` is an observer statement, not a typed fiber. A label that appears only after moving to `Sp` is a *different remainder*. Neither is a reason to silently swap the contract.

5. **Dictionary-relative realization** (#580/#588/#598). P398 already showed that a structured residual can be missing instantaneous state *or* projected memory, relative to the declared dictionary. Changing the dictionary is a different experiment. #584's existing-data screen does not need to change operationally (#584 comment after #580); interpretation of success already requires A/B/C (instantaneous / memory / descriptive). Adding a Gate-2-selected dictionary would confound that test with observer choice.

---

## Protocol for #584 (declared before the screen)

```text
primary coordinate:
  the frozen #582 contract
  (M / the typed CDF actually used to produce g_N and r_N)

null, first:
  smooth one-parameter Taylor curvature
  (#584 comment 5558832200; #596 comment 5558826403)

then:
  pre-existing labels vs permutation / random / production-block negative control

robustness report, not a selection:
  if the same five transitions admit Sp (and S) CDFs,
  reconstruct r_N in those channels separately
  and report whether the surviving label's alignment transports

do not:
  mix M and Sp residuals into one PCA
  let Gate 2's ranking pick the primary coordinate
  call a label that lives in only one channel a typed fiber
```

How to read the robustness report:

| Pattern | Meaning |
|---|---|
| Label survives `M` and `Sp` | More interesting. Still a `typed context index`, not a state. |
| Label only in `M` | Observer-bandwidth (#419). Do not promote. |
| Label only in `Sp` | A different remainder. Open a *declared* retarget, do not silently swap. |
| Labels equivalent on `M`, separated on `Sp` (or vice versa) | Report the equivalence class. Crossed transition, not a third latent. |

---

## If a crossed transition is later bought

That purchase is #584's "acquisition only if the existing five are insufficient", not Gate 2's N=650 calibration.

- The **label** stays the one declared on the #582 contract.
- The **angular/Smith work**, if any of that block is also used for #589/#583, is bought in `Sp` (Gate 2's ranking, still unpriced).
- Do not require the crossed design to be in `M` because the remainder was typed in `M`, and do not require the remainder to have been typed in `Sp` because future angular work is cheaper there. Those are two contracts on one acquisition, stated separately.

Sitter–Wu 1993 F-optimality remains the named Fieller-as-design paper, and still chooses dose-support of a binary-response experiment, not among threshold-law CDFs (#601 round 2). It does not license this retarget.

---

## What this does not do

- Not a STATUS edit.
- Not a rewrite of Gate 2's ranking or of `GATE2_NOT_ESTABLISHED`.
- Not a close of #584 or a start of the screen.
- Not a claim that `M` is the physically correct full-law channel.
- Not an identification of Gate 2's 16× with Gate 3's residual SNR.
- Not a reason to reimplement Gate 2. The script is at `390118b3`; `d5c7dad5` is the ranking. Job 2 (#600) reuses it.
