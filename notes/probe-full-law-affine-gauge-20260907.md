# Atlas — probe HIGH: the affine gauge of the full-law correction (#617)

**Probe:** `notes/probes/probe-HIGH-full-law-affine-gauge-20260907.md` (PR #616)
**Branch:** `claude/p617-affine-gauge`, PR against `claude/matching-one-workspace-pwr5pv`
**Frontier at start:** `8b5f9d1a`. **Committed artifacts referenced:** PR #614 (`f0981a98`), #611 (Astra), PR #605 (Gate 3), #582/#584 (frozen g).
**Standing rules honored:** no second exponent, no arXiv, no production, no `docs/STATUS.md` edit, issue #617 left open.

---

## Verdict table (one kill/accept sentence per program)

| P | program | result | kill/accept sentence |
|---|---------|--------|----------------------|
| P0 | reproduction gate | **PASS** | identity `abs err` 1.42e-13 / 2.91e-13 (tol 1e-12); chart shares 97.49% / 94.89% match to <1e-5; 1.55 turned ON at λ=1 (1.524/1.520) and OFF toward λ=0 (1.403/1.396). The probe may continue. |
| P1 | groupoid as matrices | **LEMMA ACCEPTED** | the g-covector is exhibited in R^9 (biorthogonal to {1, A}, `⟨ell_C, A⟩=0` by construction); transport check: predicted amplitude change −9.7274e-4 vs actual −9.7274e-4; the committed identity is the λ=1 member of the exact attachment family; collapse to the committed /k form at 1.42e-13. "Keeping g fixed" is meaningless until the transported covector is named. |
| P2 | attachment path (21 λ × 2 weightings) | **PASS (with a remainder)** | the fake excess is SMOOTH in λ, peaks at the historical middle chart (1.537→measured/pred 1.5368 spin0 g13), falls monotonically toward the consistent attachment (1.4028 at λ=0). The identity is the whole story of the 55% **as far as the sweep goes** — but λ=0 lands at 1.40, not 1.00: attachment explains ~¼–½ of the fake excess; the rest is the /k contraction plus a genuine residual (P3). |
| P3 | holonomy both triangles | **NO CLOSED-LOOP HOLONOMY; RESIDUAL NON-ZERO** | D1: attachment moves the amplitude 7.97% / 8.20% (NOT <4.25%, and NOT closed-loop holonomy — the covector map is a recompute, reversible by construction). D2: curvature residual statistic 2200.2 / 870.0 on 6 df — NOT consistent with zero: a perpendicular shape component exists in the declared chart. p50 correctly refuses (two sizes, no curvature emitted). |
| P4 | weighting fibre | **1-D FIBRE, NO ORTHOGONAL DIRECTION** | the two stored orientations span the fibre; `w0+w1=1` fixes the affine gauge; no orthogonal weighting exists in the stored data (recorded, not assumed). Along the segment: first amplitudes move up to −3.2%, the curvature −4.99%…+6.02%, the N=725 amplitude 0.02% (≈5σ by its SE). The 12%/25σ historical claim is confirmed as a directional derivative, not two anecdotes; no weighting-invariant sentence exists in the two-weight span. |
| P5 | reconstruction C→C′ | **MIXED (readout-bound)** | gate: production readout reproduces the committed p50 object bit-for-bit (err 0.0). Same-rank grid shift (u∈{0.12..0.92}): ratio moves +0.89σ — STABLE. Coarse rank-5: −57.8σ under the grid-rebuilt direction, −1.41σ under the carried direction — the two direction-transport laws themselves disagree; the 4.25% is a property of the DECLARED readout, not a readout invariant. Reconstruction bias does not manufacture the 4%. |
| P6 | candidate invariants | **1 SURVIVOR (I4); I1/I2/I3/I5 killed** | I1 killed (residual-transport pairings move 68σ/58σ across λ — r_N is chart-bound even after transport); I2 killed (Z-vector single-power χ² 2.6e5–3.7e5, catastrophic at all anchors; Astra's scalar-clean/vector-catastrophic pattern reproduced); I3 killed (anchors agree with each other, but every Z-vector is catastrophic like I2); I4 SURVIVES (predeclared denominator, no sign flip, projective ρ span 15.1% < 50%); I5 killed (independent Aff(1)-quotient disagrees with the pipeline's by 5.8e4σ — the quotient is base-dependent). |
| P7 | declared-chart freeze + the 4% | **4% DIES under the label-null** | chart frozen (spin0 AND equal, always together; reconstruction named). Scalar excess: 4.251% (spin0) / 3.879% (equal). Sign-flip null on batch labels (10⁴ draws, seed 6170425): observed residual statistic 42.8 sits at p=0.647 inside the null world (null median 47.8, max 107.5). The perpendicular share of the residual norm is 24.4% but statistically indistinguishable from label noise. Per the brief's stop rule: the death is written and no replacement is hunted. |
| P8 | N=725 forecast | **SUPPORTED (finite)** | a_hat(290→725): −4.8083e-4 (spin0, se 2.294e-6) / −4.8093e-4 (equal); the 3-se interval lies INSIDE the pre-registered ±5% band around the frozen forecast −4.6808e-4 in BOTH weightings. Supporting a finite forecast is not the model being exact: the displacement is −5.56σ of the measurement and the residual statistic is 1447.5/6 df (spin0), 1495.9/6 (equal). I4's ρ: 651.1 / 651.0 (weighting-stable to 0.02%). |
| P9 | synthetic 1.55, noise-free | **PASSED** | a pure Aff(1)+one-frozen-ray family (rates calibrated to the committed β0/k) MANUFACTURES a naive ratio 1.6455 at the middle-vs-base mismatch (brief's band [1.4, 1.7]); the consistent attachment pulls it to 1.5013; the #612 transported formula closes on the middle-chart curvature to 1.84e-14 relative. 1.55 is coordinate geometry, certified without statistics. |
| P10 | Gate-2 projective vs λ | **PASS (projective kills the fake excess)** | the projective ρ moves only 15.3% across the full λ path (threshold: 50%); no denominator sign flip; Fieller sets bounded at every λ. The projective design is genuinely attachment-stable — #602's warning does NOT need strengthening on this axis. |
| P11 | the allowed #582 sentence | **WRITTEN** | see below. |

---

## The pipeline's groupoid (P1 deliverable, in one paragraph)

The observer's state is a 9-vector Q(N) = F_N^{-1}(u), u ∈ {0.1..0.9}. Five generators act on it:

1. **location** `Q ↦ Q + c·1`;
2. **width** `Q ↦ αQ + β·1` (the fitted step contraction, `k = 1 + h0·β0`);
3. **attachment** `Q_λ = (1−λ)Q_base + λQ_target` — the curvature chart's middle basis vector moves along the lower transition's segment;
4. **orientation weighting** `Q = w0·Q^(0) + w1·Q^(1)`, `w0+w1=1` — a ONE-dimensional fibre (P4);
5. **the g-ray** — the consensus direction is *estimated*, so it is not a gauge generator of the same kind; its covector transports as a linear functional (P1's explicit 9-vector `ell_C`).

The amplitude of g is a coordinate on a bundle whose fibre depends on the estimator: the GLS covector `ell` of the curvature fit satisfies `⟨ell, 1⟩ = ⟨ell, A⟩ = 0` (biorthogonality) and the amplitude is `⟨ell, y⟩` for the second-difference observation `y`. Change the attachment and you change `ell`, hence the coordinate — even with g itself frozen. The #612 identity is the parallel transport of this covector, exact at every λ (the covector expansion, verified to 1.6e-13), and the committed /k form is its λ=1 member.

---

## What remains of the full-law residual

- The 55% second-difference excess was chart geometry: ~97% (g13) / ~95% (g17) of it is the attachment/k term (P0), and the sweep P2 shows the fake multiplier as a smooth function of λ.
- The leftover ~4% (p50 triangle, `measured/chart_corrected` = 1.04251):
  - survives as a number only in the declared readout (P5: stable under same-rank grid shift, not defined at coarse rank);
  - is NOT weighting-stable as a *number* (P4: the curvature moves −5%…+6% along the weighting segment), though both weightings see an excess of the same sign;
  - DIES under the sign-flip label null (P7: p=0.647): the residual statistic is ordinary label noise in the declared chart.
- **`FULL_LAW_SHAPE_IS_CHART_GAUGE` does not follow** — because one candidate survived (I4) and because P7's death is a null result on THIS chart, not a proof of invariance. What follows is the weaker, honest sentence:

> The declared chart carries no evidence of a full-law shape beyond label noise (P7), the fake 1.55 is coordinate geometry certified noise-free (P9), the only surviving invariant is the projective ρ with a predeclared denominator (P6-I4, confirmed over the full attachment path in P10), and the finite N=725 forecast holds in both weightings (P8) — while every other candidate invariant is chart- or estimator-bound (P6).

---

## P11 — the sentence #582 may still say

```text
There is a dominant transferable non-affine direction g in the declared
chart (span{1, Q_base, g_frozen}), with a finite-size forecast that is
chart- and weighting-dependent at the few-percent level, plus a small
full-law residual (~4% on the p50 triangle) that the declared chart's
sign-flip null does not distinguish from label noise and whose remaining
interpretation is a chart/weighting/readout choice, not a typed shape.
```

Anything stronger needs a survivor that is simultaneously attachment-, weighting-, and readout-stable. I4 is attachment-stable (P10) and weighting-stable at 0.02% (P8) but was not tested on a coarse readout, and the 4%-residual that motivated the search died under the null. No second exponent is claimed anywhere.

---

## Holes (punch list, not conclusions)

1. **N-normalisation axis not swept** (P5): the binomial-profile CDF is frozen with the artifacts; the production path exposes no toggle. Recorded as a named hole.
2. **The weighting fibre is one-dimensional** (P4): no orthogonal weighting exists in the stored data; the "orthogonal complement" test is vacuous, recorded as such.
3. **p50 has two sizes** (P3): no curvature, no holonomy, no 4%-analogue; any three-size construction must refuse on p50, and P3's does.
4. **The direction-transport law across readout ranks is undetermined** (P5): rebuilding the consensus construction on a grid and carrying the frozen g by u-interpolation disagree badly; the frozen g is grid-bound.
5. **P7's perpendicular share (24.4%) is not interpretable** beyond the null: the residual covariance's rank structure after the chart projection was not resolved further (the projection is covariance-metric, but the perpendicular noise scale was not independently calibrated).
6. **Vendoring**: `p612_chart_identity.py`, `p612_n725_score.py` (commit `f0981a98` of PR #614) and the two committed JSONs are vendored verbatim into `scripts/probe/_vendored/` and `results/probe-affine-gauge/_reference/` because the frontier has not merged #614. If #614 merges, the vendored copies are pinned references, not forks.

---

## Provenance

- Frozen g: `results/p582-amplitude-law/latest.json` (frontier), never re-fit.
- Committed forecast & p50 anchor: `results/probe-affine-gauge/_reference/p612-n725-score-latest.json` (= PR #614's `results/p612-n725-score/latest.json`).
- Astra negative control: `results/probe-affine-gauge/_reference/astra610-intrinsic-shape-pilot.json` (= PR #611's `results/astra610-independent/intrinsic-shape-pilot.json`).
- N=725 raw: `results/server-20260907/P612-n725-fullcurve/raw/` (vendored from PR #614's branch; required by `load_all(with_725=True)`).
- All programs are deterministic, run from committed blocks only, and cache the production load in `results/probe-affine-gauge/_cache-*.pkl` (delete to rebuild).
