# P4 — A matching-odd finite-size signal, and the failure of scalar closure

Draft manuscript notes for **Paper 1** of [issue #644](https://github.com/LightChainr/Matching-One/issues/644).
Portfolio track **P4**, following the `p2`/`p3` manuscript convention.

| File | Role |
|---|---|
| `README.md` (this file) | the draft: claims, evidence, limitations, omissions |

**Status of this document.** A first-pass writing draft only. It promotes nothing,
re-scores nothing, and computes nothing; `docs/STATUS.md` was not touched and the claim
levels quoted below are its current `main` content, cited as-is. Issue #644 stays open.
No venue is selected. Unlike `p2`/`p3` there is **no generated `tables.md` and no
evidence-assembly script yet** — every number below is carried with the artifact path and
JSON field where its full-precision value is committed, which is the substitute contract
until a generation pipeline is written (see §10). All paths are relative to the repository
root.

## 1. Abstract (draft)

We report a real matching-odd orientation signal in finite 2D percolation lattices,
established by a two-block independent synthesis, and we report that every attempt to
close the observed finite-size structure with a *single scalar* — one multiplier, one
width, or one pure exponent — fails. Two preregistered, independent 500M-replica evidence
blocks jointly reject a global zero-effect description
(`chi2 = 31.1857355515/4`, `p = 2.81e-6`), while frozen fixed-H4 predictions remain
compatible (`3.4622795373/4`, `p = .484`). On a frozen norm-5 harmonic discrimination,
the spin-4 alias H4 beats its adversarial alias H12 and its control H8 decisively.
A frozen two-parameter finite-size center-slope correction predicts a held-out doubling
lineage to within `0.67` standard errors where the bare asymptotic multiplier misses by
more than 22; but the same held-out curve rejects one-multiplier transfer of the
three-level thermal-even profile (`9.3520/2`, `p = .0093`), a scalar width-collapse
over the full covariance fails (`24.5004/10`), and the pure exponent law
`P4[S'] ~ N^(-5/4)` is prospectively falsified on new geometry (`52.71634/2`).
The contribution is the combination: the signal is real, a frozen finite-size transfer
works, and *scalar closure in any of its three forms does not*. We state the correlated
block accounting explicitly, and we make no claim about `p_c`.

## 2. Scope and non-goals

- **No threshold claim.** Nothing here addresses `p_c`, its value, or its identification.
  The ticket is explicit that this paper "does not need to explain `p_c` and should not
  try."
- **No unified-theory claim.** The paper does not assert the continuum operator, the state
  dimension, or an LCFT module identification. The live mechanism state is low-rank /
  multicomponent and that is where the paper stops.
- **No new computation.** All numbers are already scored and committed; this draft adds
  arrangement and wording only.
- What the negatives are: prospective falsifications of frozen scalar-closure hypotheses,
  made countable by the frozen-prediction discipline of §5. That discipline is the
  paper's methodological spine, not an appendix.

## 3. Observables and the exact structure that defines them

**The observable.** From the threshold-rank engine's long-form batch histograms, the
intrinsic-center p0 where the direction-average matching function vanishes, and the
primal/matching pair `S = (R_G + R_hat)/2`, `D = (R_G − R_hat)/2` (with
primes = `p`-derivatives at the center, from the exact beta-density identities
`dR_hat(1−p)/dp = −E[beta_pdf(K_minus)]`), the spin-4 orientation contrast is
`P4[X] = (X(θ1) − X(θ2))/Δcos(4θ)` — the first-minus-second orientation pairing is the
matching-odd channel. Definitions and frozen parity powers (`P4[S] ~ N^-1`,
`P4[D] ~ N^-13/8`, `P4[D'] ~ N^-5/8`, `P4[S'] ~ N^-5/4`) are committed in the analyzer
header, `scripts/analyze_matching_parity_derivatives.py`. The primary target statistic of
§5.1 is the matching-odd orientation contrast of the wrapping-channel matching function
(threshold-rank engine, `first − second` in frozen genealogy order).

**Exact channel correction (spine, not a footnote).** For the complementary
primal/matching construction on the torus, cross and either channels obey the exact
topological exchange `S_either = 1 − S_cross`, hence for orientation contrasts

```text
DeltaS_cross = -DeltaS_either        (exact, no fit)
```

recorded with its corrected no-refit score (`0.5700315436/2`, residual z `+0.6672,
-0.1189`) in `results/server-20260828/P43-heldout-fullcurve-500m/analysis/channel_map_corrected_DeltaS.json`
and `notes/issue43-cross-either-channel-map-erratum.md`. The odd-channel companion map has
scale `+1` (`D_either = D_cross`), validated as an executable gate before scoring —
`notes/norm5-harmonic-primary-scorer.md`, gate artifact
`predictions/norm5_harmonic_semantic_gate_20260829.yaml`. Without these exact maps the
apparent sign behaviour of the matching-even sector is misread (the erratum's consequence:
an interpretation fix, not a data problem — #108's two-sector reading reclassified to
prospective *compatibility* after the map).

**Finite Russo identity (spine).** On a finite periodic quotient, the matching slope is
exactly a sum of two nonnegative total pivotal masses:

```text
M'(p) = sum_v P_p(Piv_G(v)) + sum_v P_{1-p}(Piv_hat(v))
```

with `M = R_G(p) − R_hat(1−p)`; this is a chain-rule fact, not a scaling assumption.
Derivation and three independent exact verification paths:
`notes/pivotal-russo-finite-volume.md`, `scripts/exact_pivotal_russo.py`; the committed
bridge analyses use it as a definition —
`results/local-20260829/P100-pivotal-h4-bridge/score.json` (`definitions.Mbar_prime`:
"exact total pivotal mass by Russo") — and the companion marked-pivotal analysis enforces
the control role: `russo_control_role: "mu0 is a regression control, not independent
evidence"` (`results/local-20260829/P100-marked-pivotal-h4/analysis/score.json`).

**Frozen predictions.** Every falsification in §6 was scored against a prediction frozen
before its target was revealed:

```text
predictions/gaussian_norm5_harmonic_discrimination_20260828.yaml   frozen_at 2026-08-28,
    status preregistered_before_child_runs (H4/H12/H8 multipliers exact rational)
predictions/p48_sprime_correction_20260828.yaml                    status frozen_before_target_reveal
predictions/p49_slope_two_sector_145_290_20260828.yaml             status frozen_before_fullcurve_N290_reveal
```

This is what makes a negative result a *prospective falsification* rather than a bad fit
(GOVERNANCE §2C).

## 4. The evidence blocks, and counting one block once

GOVERNANCE §2E: roots, slopes, derivatives, quantiles and score modes taken from the same
histograms are all useful, and all are the *same block*; adding them as independent
evidence inflates the conclusion. The paper's block inventory:

| block | content | primary score |
|---|---|---|
| **I** | `issue43_n185_n265_500m_histograms` — held-out N=185/265 full curves, disjoint counter domains | matching-odd `DeltaM` fixed-H4 / zero-effect |
| **II** | `issue57_norm5_production` — N=325/425 norm-5 harmonic discrimination, 500M/orientation | same, in its own channel |
| **III** | `P50` N=145 → 290 full curve, independent seeds | three-level transfer, slope, channels |

The joint two-block claim (Block I + II) is the Issue #212 synthesis, whose artifact
freezes the independence contract (`raw_data_groups` distinct,
`joint_covariance: block_diagonal`, scope "no other ledger views enter"):
`results/evidence-ledger/issue212-matching-odd-synthesis.json`.

Every other score in this paper is a *view of one block*, and the draft says so at each
use: the H4/H12/H8 discrimination and the width-collapse score both live in Block II (the
thermal-jet report states its secondary scores "reuse the primary raw curves and are not
an additional independent evidence block"); the P48 `S'` pure-law score reads the intrinsic
center analyses of Block I's own runs (`results/server-20260828/P48-new-geometry-score/REPORT.md`
names its targets `analysis/n185.p48.json`, `analysis/n265.p48.json` from
`P43-heldout-fullcurve-500m`); the central-contrast and channel tables of Block III are
diagnostics, not evidence rows. There is no number in this paper that is the product of
four confirmations; the honest headline count is two independent primary blocks plus one
independent third block carrying falsifications.

## 5. Positive results

**5.1 The signal exists (two-block, primary).**

```text
zero effect:   chi2 = 31.18573555150965 / 4,  p = 2.805595267905808e-6
fixed H4:      chi2 = 3.4622795373044295  / 4,  p = 0.48363695393249573
```

Field paths: `joint_scores.zero_effect.chi_square`,
`joint_scores.zero_effect.chi_square_survival_p`, `joint_scores.fixed_H4.*` in
`results/evidence-ledger/issue212-matching-odd-synthesis.json`; prose
`notes/issue212-matching-odd-synthesis.md`. Direct predictive comparison: NLPD
`fixed_H4 − zero_effect = -13.806418789729096` nats (field
`predictive_comparison.delta_nlpd_fixed_H4_minus_zero_effect`); the artifact records
`preferred_lower_nlpd: fixed_H4`. The chi-square *difference* is deliberately not treated
as a likelihood-ratio statistic (frozen covariances need not agree).

Ledger level cited from `main` without promotion: `docs/STATUS.md` line 25, "Square-site
matching-odd orientation signal exists — C3".

**5.2 Frozen H4 transfers; tested aliases are disfavored (Block II).**
Fixed, parameter-free multipliers recomputed from exact rational angular factors and
`5^(-13/8)`:

| model | harmonic | chi2 / df | survival |
|---|---|---:|---:|
| H4 (primary) | `cos_4_theta`, target `-14/25 · 5^(-13/8)` | 0.4163 / 2 | 0.8121 |
| H12 (adversarial alias) | `cos_12_theta`, `23506/15625 · 5^(-13/8)`, opposite child sign | 35.1931 / 2 | 2.280e-8 |
| H8 (even-m control) | `cos_8_theta`, `-1054/625 · 5^(-13/8)` | 16.0120 / 2 | 3.335e-4 |
| zero effect | — | 1.7764 / 2 | 0.4114 |

Source: `results/server-20260829/P57-norm5-500m/primary_score.json` (models table,
full-precision `chi_square` fields `0.41630376401835489631`, `35.193078878051031476`,
`16.012022412848537938`), prose `results/server-20260829/P57-norm5-500m/REPORT.md`,
frozen contract `predictions/gaussian_norm5_harmonic_discrimination_20260828.yaml`.
H4 residuals `+0.645σ, −0.021σ`; H12 residuals `−4.612σ, −3.770σ`.
`docs/STATUS.md` line 27: C3.

**The honesty clause this section must keep:** the two target points alone do **not**
reject a zero child effect (row 4 survives at `0.4114`). §5.2 supports the frozen H4
*alias choice*; it is not an independent discovery of a nonzero H4 amplitude — the
discovery claim is §5.1's, on two blocks. The REPORT says this in its own voice and the
paper repeats it.

**5.3 A frozen finite-size correction predicts a held-out lineage (Block III).**
The two-parameter center-slope correction (scalar + frozen H4), trained only on
`65 → 130` and `85 → 170`, predicts the held-out `145 → 290` slope:

```text
frozen scalar+H4 corrected slope:  signed z = -0.6660797543793089  (p = 0.50536)
raw 2^(3/8) asymptotic baseline:   signed z = -22.690348726453408  (p = 5.58e-114)
```

Fields `slope.frozen_scalar_plus_H4_correction.signed_z` and
`slope.raw_asymptotic_baseline.signed_z` of
`results/server-20260829/P50-n145-n290-fullcurve/analysis/score.json`; table
`results/server-20260829/P50-n145-n290-fullcurve/REPORT.md`; frozen prediction
`predictions/p49_slope_two_sector_145_290_20260828.yaml`.
`docs/STATUS.md` line 30: C3. This is the paper's cleanest positive finite-size transfer:
a prediction with **zero refitted parameters** on the target, missing by 22 standard
errors without the correction and by less than one with it.

## 6. The failure of scalar closure (three forms, prospectively)

**6.1 One multiplier does not close the three-level curve (Block III).**
Joint `DeltaM` at `u = 0, 0.025, 0.05` under the frozen multiplier:
`chi2 = 9.352003684819664 / 2`, `p = 0.009316` (field
`primary_deltaM_transfer.chi_square` / `.chi_square_survival`,
`results/server-20260829/P50-n145-n290-fullcurve/analysis/score.json`).
`docs/STATUS.md` line 29: C3 negative.

Where the failure lives, and why it is a shape failure and not an amplitude failure:
the central contrast residual is only `−1.385` standard errors on its own
(observed ratio `−0.52864747` vs frozen `−0.32420989`); rejection comes from the
strongly-correlated three levels resolving a small deformation. Residual-correlation
eigenvalues `1.23785e-13, 1.76308e-6, 2.99999824` give numerical rank 2 at the frozen
cutoff `1e-10`; the amplitude-like mode contributes `≈ 1.9243` of the chi-square, the
active shape mode `≈ 7.4277` (REPORT.md, "Where the joint DeltaM failure lives").
The result does not overturn §5; it shows a single multiplier does not transport the
local thermal-even curve *shape*.

**6.2 One width does not collapse the thermal jet (Block II).**
Full-covariance norm-5 width-collapse score:
`chi2 = 24.500424925937486902 / 10`, survival `0.006377` (field
`primary_width_collapse.score.chi_square`,
`results/server-20260829/P57-norm5-500m/thermal_jet_score.json`; table
`results/server-20260829/P57-norm5-500m/THERMAL_JET_SCORE.md`). Ordered follow-ons on
the same block: width-corrected `q=2` analytic (`c=8/5`) `22.2386/10`, `p = .013934`;
width-corrected rank-2 Jordan (`c = log5/log2`) `17.0513/10`, `p = .073237`.
`docs/STATUS.md` line 32: C2 negative (scalar width); line 33 records that Jordan is
*compatible, not identified*.

**6.3 One pure exponent does not close `S'` (Block I).**
The parity rule with `x = 21/4` and `y_t = 3/4` gives `P4[S'] ~ N^(-5/4)`; the
scaled `Y_N = N^(5/4) P4[S']` frozen pure law (`A0 = 1.9434247576878727 ±
0.0766048795577632`, trained on N=65/85/130) scores on the new N=185/265 geometries
`chi2 = 52.71633588357711 / 2` (field `channels.P4_S_prime.chi_square`,
`results/server-20260828/P48-new-geometry-score/score.json`; zero-model
`1278.555`; marginals `+5.04σ` and `+6.61σ`).
`docs/STATUS.md` line 31: C3 negative, "prospectively falsified".
The three sibling channels of the same frozen score survive their pure laws
(`P4[S]` `1.139/2`, `P4[D]` `0.281/2`, `P4[D']` `0.088/2`) — `S'` is the unique
clear pure-law failure, so the falsification is specific to the matching-odd
derivative channel, not to finite-size laws in general.
The exponent *arithmetic* is not the broken part (the parity-rule derivation is in
`notes/p48-sprime-correction-discriminator.md`); the constant-amplitude *closure* is.

## 7. Where the mechanism boundary actually sits

The three scalar failures are not generic. Within Block III, the projected derivative
channels split sharply: `P4[D']` transfers under its frozen multiplier to `−0.009σ`
(`p = 0.99275`) while `P4[S']` misses at `+2.70σ` (`p = 0.00703`) — REPORT.md,
"Correlated P4 diagnostics". One projected direction of the derivative response is
already a near-exact eigenchannel; its companion is not governed by the same scalar law.
That is the evidence for treating the response as a small non-diagonal (Jordan-type)
block, and it is the same ordering the independent norm-5 diagnostics found (Jordan
`17.05/10` vs `q=2` `22.24/10`, improvement `5.19`, explicitly a ranking diagnostic
between non-nested fixed models, not a likelihood-ratio test).

The frozen correction models for the failed `S'` channel were themselves preregistered
for N290 and both survive there (`q=2` even-scalar correction `z = 1.30913`,
rank-2 Jordan log `z = 0.14675`; `results/server-20260829/P50-n145-n290-fullcurve/analysis/sprime_frozen_score.json`,
`analysis/score.json` and REPORT.md table), with q=2 first in the declared order.
The paper's mechanism claim stops exactly here: **low-rank non-scalar finite-size
mixing, with a shape mode carrying the rejection; no module identification.**
`docs/STATUS.md` line 68 records the `x = 21/4`-candidate reading of the leading sector;
this draft quotes that level and does not exceed it.

## 8. Limitations

1. **Correlated-block accounting (§4) governs every "combined" reading of the results.**
   The only multi-block claim is the two-block zero rejection (Block I+II). Blocks II and
   III contributions are single-block statements; P48 shares its random block with P43's
   primary score; the thermal-jet and central-contrast tables reuse primary raw curves.
   GOVERNANCE §2E.
2. **The N=580 aspect ladder does not adjudicate the amplitude law (#577).** Verdict
   committed: `"underpowered: more than one competing prediction survives at r4_over_r1"`
   (`results/aspect-ladder-n580/latest.json`, `verdict`; narrative
   `notes/aspect-ladder-n580-result-20260905.md`). Three survivors — bare aspect ratio
   (`4.00`, z `+0.50`), weight-4 modular shape (`10.99`, z `−2.08`), plain area scaling
   (`16`, z `−2.56`) — one clean exclusion (`no_modulus_dependence`, z `+9.5`), 3σ
   Fieller interval `[2.40, 27.47]`. The shape/modulus assignment of the `S'` partner is
   therefore a limitation, not a finding, of this paper; the paper's exclusion claims must
   not be read as selecting among the three survivors. (Cross-rung covariance measured by
   replay: ρ `−0.1648`; moves no verdict — same artifact.)
3. **Block II alone does not reject a zero child effect** (§5.2 clause). The rejection of
   global zero rests on the synthesis of §5.1; alias selection and signal existence are
   different claims with different blocks.
4. **Prospective discipline is per-prediction, not per-paper.** The P48 pure law was
   trained on the earlier sizes N=65/85/130 (a retrospective source) but frozen before
   the N=185/265 target values existed (`predictions/p48_sprime_correction_20260828.yaml`,
   `status: frozen_before_target_reveal`; P48 REPORT.md "Frozen source"); the falsification
   is prospective against that new geometry, and §6.3 says so in those words. The
   `2^(3/8)` baseline of §5.3 is a frozen asymptotic reference, not a rival fit.
5. **No covariance pooling across blocks beyond the block-diagonal contract.** The §5.1
   joint scores add two frozen covariances; the NLPD comparison is offered precisely
   because the raw chi-square difference is *not* a valid LR statistic
   (`notes/issue212-matching-odd-synthesis.md`).
6. **Semantics are load-bearing.** The cross/either sign map is exact but channel- and
   parity-specific (even: `−1`; odd: `+1`); an earlier interpretation error (#108) shows
   the failure mode. The paper inherits the typed-gate requirement (source/target channel
   fields, exact map) as reporting hygiene — erratum "Governance consequence".
7. **The scalar-width failure is scored on a 10-dof full-covariance collapse**; the
   ordered successors survive differently at their own levels (§6.2), and none is
   promoted by this draft.
8. **The three-level transfer rejection is at N=145→290, one lineage.** A second
   independent doubling curve (new geometry, not a rerun) is the natural referee ask; we
   note it rather than pre-empt it.

## 9. Relation to `docs/STATUS.md` (read-only)

Cited levels (as of `main`, this draft changes nothing): line 25 C3 (signal exists);
line 27 C3 (H4 beats aliases); lines 29, 31 C3 negative (multiplier, pure `S'` law);
line 30 C3 (frozen slope correction); line 32 C2 negative (scalar width); line 36 C3
negative (weight-4 amplitude law across aspects — underpowered); line 68 candidate
reading. Levels change on evidence, not on prose.

## 10. What had to be left out (for the next draft)

1. **Generated tables and a machine-checkable assembly.** A
   `scripts/p4_manuscript_evidence_table.py` + `tables.md` + drift test in the `p2`/`p3`
   style, so no number in the manuscript is typed. Until then §5–§7 carry explicit
   `path → JSON field` citations instead.
2. **A `references.bib` pass with source verification.** The current bib
   (10 keys) has no finite-size-scaling or crossing-probability entries; the paper will
   need at least Cardy/Watts crossing-probability and standard FSS references, each read
   and marked per repo `[LIT]`/`PRIMARY_TEXT_READ` convention (p2's README documents that
   workflow). Deliberately not added under a no-new-work ticket — invented bib entries
   would violate the same discipline that makes the falsifications count.
3. **The Q4/Jordan module bridge and the shape/modulus question**, which the ladder
   limitation (§8.2) keeps open; also the #138 conjugation phase-node diagnostic
   construction (P50 REPORT "no new Monte Carlo block is needed" note), which belongs to
   a mechanism paper, not to this one.
4. **A power analysis for the shape-mode rejection** (what curve count separates frozen
   shape predictions at the current SE scale) — analysis, not computation, but new work.
5. **Venue decision** — explicitly deferred by #644 and non-blocking.

## 11. Evidence index

| # | number | committed artifact (field) |
|---|---|---|
| 1 | `31.18573555150965 / 4`, `2.805595267905808e-6` | `results/evidence-ledger/issue212-matching-odd-synthesis.json` `joint_scores.zero_effect` |
| 2 | `3.4622795373044295 / 4`, `0.48363695393249573` | same, `joint_scores.fixed_H4` |
| 3 | NLPD delta `-13.806418789729096` | same, `predictive_comparison` |
| 4 | `0.41630376401835489631 / 2`; `35.193078878051031476 / 2`; `16.012022412848537938 / 2`; `1.7764 / 2` | `results/server-20260829/P57-norm5-500m/primary_score.json` (models) |
| 5 | H4/H12 frozen multipliers | `predictions/gaussian_norm5_harmonic_discrimination_20260828.yaml` (`exact_harmonic_predictions`) |
| 6 | `9.352003684819664 / 2`, `0.00931618710954876`; eigenvalues; mode split `1.9243/7.4277` | `results/server-20260829/P50-n145-n290-fullcurve/analysis/score.json`; `REPORT.md` |
| 7 | `z = -0.6660797543793089` vs `-22.690348726453408` | same, `slope.*` |
| 8 | `24.500424925937486902 / 10` (`0.006377`); `22.2386 / 10`; `17.0513 / 10` | `results/server-20260829/P57-norm5-500m/thermal_jet_score.json` `primary_width_collapse.score`; `THERMAL_JET_SCORE.md` |
| 9 | `52.71633588357711 / 2`; siblings `1.139/0.281/0.088`; `A0 = 1.9434247576878727` | `results/server-20260828/P48-new-geometry-score/score.json` `channels.*`; `REPORT.md` |
| 10 | `0.5700315435551194 / 2` (`+0.6672, −0.1189`) | `results/server-20260828/P43-heldout-fullcurve-500m/analysis/channel_map_corrected_DeltaS.json`; `notes/issue43-cross-either-channel-map-erratum.md` |
| 11 | Russo identity + pivotal bridge | `notes/pivotal-russo-finite-volume.md`; `scripts/exact_pivotal_russo.py`; `results/local-20260829/P100-pivotal-h4-bridge/score.json` |
| 12 | ladder verdict; Fieller `[2.40, 27.47]`; survivors | `results/aspect-ladder-n580/latest.json` (`verdict`, compatible set); `notes/aspect-ladder-n580-result-20260905.md` |
| 13 | P48 frozen correction targets (N290): q2 `1.30913`, Jordan `0.14675` | `results/server-20260829/P50-n145-n290-fullcurve/analysis/sprime_frozen_score.json`; `REPORT.md` |
| 14 | D'/S' split `−0.00909` / `+2.69536` | `results/server-20260829/P50-n145-n290-fullcurve/REPORT.md` (Correlated P4 diagnostics) |
| 15 | reproduction: `scripts/score_matching_odd_synthesis.py`, `scripts/score_norm5_harmonic_primary_typed.py` | `notes/issue212-matching-odd-synthesis.md`, `notes/norm5-harmonic-primary-scorer.md` |
