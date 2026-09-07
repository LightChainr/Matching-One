# Probe HIGH — Affine gauge of the full-law correction

**Assign to:** one Agent with high reasoning quality **and** enough CPU to rescore every committed histogram block (numpy / mpmath; a full day of cores is in budget, a new Monte Carlo production is not).

**Not for:** a retrieval-only agent, a “write a note and stop” agent, or anyone who will start a new N-block.

**Literature is out of scope.** Do not search arXiv. Do not open a literature-officer note. If a named paper is already quoted in #613 / #611 / #602, you may cite that quote. Anything not already named is someone else’s job.

---

## Standing

This is a long-horizon *computational geometry* probe of the observer, not a short execution ticket and not a request for one predetermined theorem.

- Does **not** enter `docs/STATUS.md`.
- Does **not** close #582, #584, #609, #610, #612, #275, #579, #595, #596.
- Does **not** start a production. N=725 already exists on PR #614.
- A negative answer is a result. Mark conjecture vs proof. State hypotheses.
- Do not re-derive #612’s identity from scratch; reproduce it, then build on it.
- Do not fit `A1 N^{-ω1} + A2 N^{-ω2}`. Do not name CFT operators. Do not reopen #275.
- Do not prove `Q_N(u) → p_c` rates (rate probe). Do not re-enumerate L=3,4 (exact-controls probe).

Frontier, not `main`:

```text
claude/matching-one-workspace-pwr5pv  @  8b5f9d1a
PR #614   claude/p612-n725-chart
PR #611   analysis/astra-610-independent-20260907
PR #605   Gate 3, merged on frontier
```

Packs against `main` (#602, #606, #607, #615) are citable. Do not rebase them here.

You own **the geometry of the observer**: the groupoid of affine charts, transport of the `g`-amplitude covector, orientation weighting, reconstruction-grid changes, and the leftover ~4% as a coordinate residual. Two parallel mid probes own rates and exact algebraic toys. If a sentence needs a rate or an L=3 fraction, write “see that probe” and keep going.

---

## Why this exists

A 55% second-difference excess was treated as a second physical scale. It was not.

First-difference amplitudes of #582’s direction `g` lived in `span{1, Q_base, g}`; the second divided difference lived in `span{1, Q_middle, g}`. Keeping `g` fixed does **not** keep the amplitude covector fixed.

```text
a_curv = 2/(h0+h1) [ a1 - a0/k + ℓ_C(r1) - ℓ_C(r0)/k ],
k = 1 + h0 β0.
```

| lineage (spin0) | chart share of the 55% | identity |a_err| |
|---|---:|---:|
| gaussian_13 | 97.49% | 1.4e-13 |
| gaussian_17 | 94.89% | 2.9e-13 |

Direct first-amplitude ratios: `1.008 / 1.024`, not `1.54 / 1.56`. **1.55 is retired.** Fitting a second exponent is looking at a residual and then giving it a mode (GOVERNANCE §2C, §2E; #602).

What remains after transport:

```text
P'_chart  = 1.0304774e-3
measured  = 1.0742792e-3
excess    ≈ 4.25%   ≈ 7.3 σ
residual-transport terms ≈ 0.12% of curvature
full-law residual still 1447 / 6 df
```

Weighting: spin0 at N=725 is a 5.56 σ miss inside a ±5% band (~±10.2 SE wide); equal has 3 σ wholly outside. Amplitude shift ~12%, ~25 σ. So `A(N) ~ N^{-0.970}` is a chart- and weighting-dependent finite-size compression, not an observer-independent law.

Astra’s negative control stands: a scalar can look clean while a 7-vector single-power fit is catastrophic (χ² in the thousands). You will keep that control as a *required* output of every candidate invariant, not as a footnote.

---

## North-star

What is the **groupoid of charts and transport maps** the pipeline actually uses, and what (if anything) is invariant?

```text
F_N(p) = [1 + M_N(p)] / 2,
Q_N(u) = F_N^{-1}(u),   u ∈ {0.1, …, 0.9}   (9-vector in an affine space A ≅ R^9).
```

Working suspicion, to prove, falsify, or replace, *with code*:

```text
Every "amplitude of g" is a coordinate on a bundle whose fibre
depends on the estimator and the chart attachment.
Changing the affine chart after seeing a residual is not a new
experiment; it is a change of coordinates on the same histograms.
A scientific claim about full-law flow must be a statement about
an object transported into one declared chart, or about an
object invariant under the groupoid.
If no nontrivial invariant exists, write FULL_LAW_SHAPE_IS_CHART_GAUGE.
```

---

## Settled facts you must not re-derive

Cite the PR/JSON. If your reproduction disagrees past the stated tolerance, **stop the whole probe** and report the break.

1. Chart identity — PR #611, PR #614 (`scripts/p612_chart_identity.py`, `results/p612-chart-identity/latest.json`).
2. Gate 3 — PR #605: no pre-existing discrete label indexes `r_N`. Consensus `g_N` carries 99.5–99.9% of each transition χ².
3. N=725 — PR #614, `(26,7)/(23,14)`, **1e8 paired / orientation** (not 100×100M).
4. Matching involution is between two models; `M_G(1/2) = -21/64` at L=3 is the exact-controls probe’s number. Fitted `g` is not an irrep of `G ⊕ Ĝ`.
5. P398 selection is representation parity; uniform `p` is even; #615 no-go.
6. #608 wrap/`X` algebra is the exact-controls probe.
7. #579 asked for projective inference; #595/#596 Gate 2 replaced A4-SNR by a projective design. You will *use* those coordinates, not rediscover Fieller.
8. #275 stays closed except under its own rule.

---

## Data you must actually load

Do not work from remembered tables. Load:

- the eight committed full-curve histogram blocks used by #582/#612 (three lineages: `gaussian_13` 65-130-325, `gaussian_17` 85-170-425, `p50` 145-290);
- N=725 raw + metadata on PR #614 (`results/server-20260907/P612-n725-fullcurve/`);
- `scripts/p612_chart_identity.py`, `scripts/p612_n725_score.py`, Gate 3 scorer from PR #605;
- delete-one jackknife batch structure (100 batches; correlated views of one block count once — GOVERNANCE §2E).

If a file is missing, stop and say which path. Do not reconstruct quantiles from moments as a silent substitute. Astra’s beta-mixture/Brent reconstruction (agreeing to `3e-14`) is an allowed *independent check*, not a replacement for the production path.

---

## Programs (all of them; none is optional colour)

A program is done only when it has (i) a script, (ii) JSON numbers, (iii) a kill/accept sentence in the atlas. “Looks reasonable” is not done.

### P0 — Bit-for-bit reproduction (gate)

Reproduce #612 identity on the production path, published `g` fixed, middle-chart curvature covector included, mpmath path.

```text
accept iff identity abs err < 1e-12 on both primary lineages
accept iff chart shares match 97.49% and 94.89% to 0.05%
```

If this fails, stop the probe. Everything else is undefined.

Also reproduce, as a *negative* control you will later require of every invariant:

- first-amplitude ratios without transport ≈ `1.008 / 1.024`;
- the “naïve 1.55” obtained by comparing un-transported first amplitudes to middle-chart curvature.

You must be able to **turn 1.55 on and off** by switching attachment point. If you cannot, you do not understand the identity yet.

### P1 — The groupoid, as matrices, not prose

Treat `Q ∈ R^9`. Write explicit 9×9 (or 9×k) matrices for every operation the pipeline uses:

| generator | matrix / map | parameters |
|---|---|---|
| location | `Q ↦ Q + c 1` | `c ∈ R` |
| width | `Q ↦ α Q + β 1` | `α>0`, or the fitted step-contraction `β` |
| attachment | chart at `Q_λ = (1-λ) Q_base + λ Q_target` | `λ ∈ [0,1]`, with `λ=1/2` the “middle” used in #610 |
| orientation weight | `Q = w_0 Q^{(0)} + w_1 Q^{(1)}`, `w_0+w_1=1` | spin0, equal, and the segment between them |
| extra ray | projection onto / away from `g` in the covariance metric | `g` is **estimated**, hence not a gauge generator of the same kind |

Deliverable: a section “the pipeline’s groupoid” with every map named, with the composition law, and with a proof (linear algebra, not rhetoric) that **the amplitude of `g` is a coordinate on a bundle whose fibre depends on the estimator**. Write the #612 identity as parallel transport of the covector dual to `g`. If you cannot exhibit the covector in `R^9`, P1 is not done.

Lemma to accept or kill: “keeping `g` fixed” is meaningless until you say *which* covector is transported. The #612 formula is the unique transport compatible with the three-size divided-difference chart actually used. If you find a different transport that also recovers 97% of 1.55, write both and kill one.

### P2 — Continuous attachment path (this is why you have CPU)

For each three-size lineage, sweep `λ ∈ {0, 0.1, …, 1}` (at least 11 points; 21 is better). At each `λ` extract curvature amplitude in `span{1, Q_λ, g}` and the transported prediction from the #612 identity written at that attachment.

Outputs, per lineage, both weightings:

- `a_curv(λ)`, `a_pred(λ)`, ratio `a_curv/a_pred`;
- the λ that maximises the fake excess (expected near the historical middle-chart);
- the λ that minimises it (expected near a consistent attachment);
- a JSON table, not a screenshot.

Kill: if the fake excess is *not* a smooth function of λ, or if it does not pass through ~1.55 near historical middle and ~1.00 near consistent attachment, the identity is not the whole story and you say so.

### P3 — Holonomy around both triangles

Lineages `65-130-325` and `85-170-425` are triangles in log-N. Parallel-transport the `g`-covector around each triangle (base→mid→target→base, and the opposite order). Report:

- holonomy as a 9×9 (or 3×3 in the `{1,Q,g}` frame) minus identity, plus a scalar angle in the covariance metric;
- comparison to the leftover 4.25%;
- whether reversing the order inverts the holonomy to 1e-12 (connection) or not (path-dependence / estimator-dependence).

Three mutually exclusive diagnostics, each with a killing computation:

| diagnosis of the 4% | kill it by |
|---|---|
| holonomy of `g` | holonomy ≪ 4% in the same units, both triangles |
| shape residual ⟂ `{1,Q,g}` | after transport, the residual’s projection on `{1,Q,g}^⊥` is consistent with 0 in the jackknife metric |
| weighting artefact | the 4% changes sign or dies when weighting runs along the spin0–equal segment |

You may conclude “mixture”. You may not conclude “second exponent”.

`p50` (145-290) has only two sizes: no triangle, no curvature in the same sense. Use it as a **stress test**: any construction that requires three sizes must *refuse* to emit a 4% number on p50, not invent one.

### P4 — Weighting as a covector, not a robustness checkbox

spin0 and equal are two points. Sweep the segment

```text
w(t) = (1-t) w_spin0 + t w_equal,   t ∈ [0,1]
```

and at least one more direction in the two-orientation plane (e.g. the orthogonal complement of spin0, normalised). For each `t`:

- first-difference amplitude of `g` on every transition;
- N=725 forecast residual (if the 725 block can be re-weighted from raw orientations — it can: `(26,7)` and `(23,14)` are stored separately);
- transported curvature residual.

Deliverable: the 12% / 25 σ shift as a *directional derivative* on this plane, not two anecdotes. Either produce a weighting-invariant sentence, or prove that no interesting invariant lives in the two-weight span. **Do not prefer spin0 because it passed.**

If raw per-orientation histograms are missing for an older block, skip that block with a named hole; do not average in spin0 as a default.

### P5 — Reconstruction operator `C → C′` (the leftover of #615 §7)

#615 §7 proposed a readout-control that is now *option 3 of 3* (the first two were “second mode” vs “quantile-reconstruction bias”; the main term was chart transport). You still owe the reconstruction axis, because it is a different groupoid generator from attachment.

Hold process, N, and declared chart fixed. Vary only the inverse-CDF recipe:

- default nine deciles vs a shifted grid `{0.12,…,0.92}` vs a coarser five-point grid;
- any `N`-normalisation the production path uses, toggled if the code path allows.

Report whether the leftover 4% *moves*. If it is readout-stable in the declared chart, reconstruction bias is not the 4%. If it moves, reconstruction is entangled with the 4% and you may not call the 4% “shape”. **No new Monte Carlo** — re-read existing histograms.

### P6 — Candidate invariants, each with a killing test

For every candidate: definition, transformation law under P1’s groupoid, 7-vector χ² (the Astra negative control is **mandatory**), both weightings, both primary lineages, N=725 forecast.

| id | candidate | kill if |
|---|---|---|
| I1 | Gate-3 `r_N` **after** transport into the declared chart | still depends on λ after P2 |
| I2 | `Z_N(u) = (Q(u)-Q(0.5))/(Q(0.8)-Q(0.2))` | 7-vector single-power χ² is “pretty” on one scalar and catastrophic on the vector (then I2 is a cherry-pick, not an invariant) |
| I3 | same with anchors `(0.3,0.7)` and `(0.1,0.9)` | the three Z’s disagree past jackknife error |
| I4 | Fieller / projective amplitude (#579, Gate 2) with a **predeclared** denominator | denominator sign flips on any committed block, or the chart-λ dependence of the projective coordinate is as large as 1.55 |
| I5 | Wasserstein tangent after quotient by `Aff(1)` (what #582 claimed to be) | the pipeline’s number disagrees with an independent `Aff(1)`-quotient implementation by more than jackknife noise |

Accept at most what survives. Zero survivors is an allowed, high-value outcome (`FULL_LAW_SHAPE_IS_CHART_GAUGE`). Do not add I6 after seeing the 4%.

### P7 — Declared-chart freeze, then the 4%

After P1–P6, freeze **one** chart in a one-page protocol card. Default unless you kill it:

```text
declared chart  = span{1, Q_base, g_frozen}
g_frozen        = #582/#584 consensus (name the commit SHA)
transport       = #612 identity, applied before any residual
weightings      = spin0 AND equal, always together
reconstruction  = production inverse-CDF, named
forbidden       = second exponent; chart change after seeing residual
```

If you prefer `Q_middle` or projective, rewrite the identity and say why #612’s formula changes. Two live defaults = P7 failed.

Then, **and only then**, recompute the leftover residual in that chart:

- scalar 4.25% analogue, both lineages, both weightings;
- full 7-vector (or 9-minus-affine) residual with jackknife covariance;
- permutation / sign-flip null on batch labels (exact or 10^4 Monte Carlo *of labels*, not of configs);
- projection onto `{1,Q,g}` vs onto the orthogonal complement.

If the 4% dies, write that and **do not hunt a replacement**. If it survives as a weighting-stable direction orthogonal to `{1,Q,g}`, hand the vector to the owner **unnamed**.

### P8 — N=725 is a forecast, not a sixth transition

In the declared chart, for every surviving invariant of P6:

```text
a_hat(290→725), se, a0, (a_hat-a0)/se,
3se interval vs pre-registered ±5% band,
spin0 AND equal.
```

Equal is not optional. Do not promote “supports this finite forecast” to “the model is exact”. Do not fold 725 into a three-size curvature unless you first declare a new freeze that includes it (and then you have used 725 twice: as freeze and as test — forbidden, GOVERNANCE §2C).

### P9 — Synthetic manufacture of 1.55 (noise-free)

Take a single exact quantile 9-vector (the exact-controls probe’s L=3 `Q_L(u)` is enough; if that PR is not merged, manufacture a monotone 9-vector with `M(1/2) ≠ 0`). Apply two different attachments to a fake two-step family whose *true* motion is pure `Aff(1)` plus a single frozen ray. Show, in code, that a middle-vs-base mismatch manufactures a ratio in `[1.4, 1.7]` while consistent attachment returns `1.00 ± 0.01`.

This is the pedagogical certificate that 1.55 is coordinate geometry. If you cannot manufacture it, P9 failed.

### P10 — Gate 2 projective design vs this groupoid

#595/#596 Gate 2 ranked channels by a projective figure of merit (Fieller / `A4²/var(δ)`), and #602 forbade retargeting Gate 3’s coordinate with that ranking. Compute, on the same eight blocks:

- the Gate-2 projective amplitude of `g` as a function of λ (P2’s path);
- whether projective coordinates *kill* the 1.55 fake excess (they should, if they are attachment-invariant) or merely hide it in a denominator.

If projective coordinates still move by O(50%) along λ, Gate 2 did not buy chart invariance and #602’s warning is stronger than written. Report that. Do not retarget Gate 3.

### P11 — Rewrite #582’s allowed sentence

After P0–P10, the only sentence #582 may still say, unless you have an invariant proof of something stronger:

```text
There is a dominant transferable non-affine direction g
in the declared chart, with a finite-size forecast that
is chart- and weighting-dependent at the few-percent
level, plus a small unresolved full-law residual that
has not been typed.
```

Anything stronger needs a P6 survivor and a P7 residual that is not holonomy. Write the sentence, the surviving invariants table, and the holes (missing raw orientations, missing reconstruction toggle, p50’s two-size limitation) as a punch list, not as a conclusion.

---

## Compute budget (use it)

Expected, on committed data only:

| job | scale |
|---|---|
| P0 reproduction | minutes |
| P2 λ-sweep × 2 lineages × 2 weightings × 11–21 λ | tens of minutes |
| P3 holonomy | minutes once P1 is right |
| P4 weighting segment × N=725 re-weight | tens of minutes |
| P5 reconstruction re-read | tens of minutes |
| P6 five candidates × 7-vector χ² | an hour |
| P8 725 forecasts | minutes (scorer exists) |
| P9 synthetic | minutes |
| P10 projective vs λ | tens of minutes |
| jackknife / 1e4 label permutations | an hour |

If something needs more than ~10 CPU-hours, you are simulating configs. Stop. If something needs a new N, you are off-scope.

---

## Stop rules

- P0 fails → stop the probe.
- No P6 survivor → write `FULL_LAW_SHAPE_IS_CHART_GAUGE` and still do P7–P8 on the declared chart (the residual can exist even if no invariant does).
- 4% dies under P7 → stop hunting.
- 4% survives unnamed → deliver the vector, do not name it, do not fit it.
- No STATUS, no ticket close, no production, no arXiv.

---

## Deliverables

One PR against `claude/matching-one-workspace-pwr5pv`:

```text
notes/probe-full-law-affine-gauge-YYYYMMDD.md     (atlas: P0–P11)
notes/probe-full-law-protocol-card-YYYYMMDD.md    (≤ 40 lines, paste-ready)
scripts/probe/affine_gauge_p0_reproduce.py
scripts/probe/affine_gauge_groupoid.py            (P1 matrices)
scripts/probe/affine_gauge_lambda_sweep.py        (P2)
scripts/probe/affine_gauge_holonomy.py            (P3)
scripts/probe/affine_gauge_weighting.py           (P4)
scripts/probe/affine_gauge_reconstruction.py      (P5)
scripts/probe/affine_gauge_invariants.py          (P6)
scripts/probe/affine_gauge_residual.py            (P7)
scripts/probe/affine_gauge_n725.py                (P8)
scripts/probe/affine_gauge_synthetic155.py        (P9)
scripts/probe/affine_gauge_projective.py          (P10)
results/probe-affine-gauge/*.json                 (one file per program)
```

Empty programs are failures, not “future work”.

## Interface

- Rates, F1/F2, Theorem L → rate probe. Need a rate? Stop that sentence.
- `M(1/2)=-21/64`, wrap/`X`, `Q(u)+Q(1-u)` → exact-controls. Cite their JSON.
- Literature → the literature probe. Do not search.
- #275 → closed.
