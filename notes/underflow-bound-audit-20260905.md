# Underflow-bound audit — what else the ~790-site ceiling capped

Audit of `LightChainr/Matching-One#570`, the follow-up to the binomial-tail
underflow found by a 1M pilot. Premise: the fix lives only in **PR #564**
(commit `844173c`, mode-anchored recurrence in
`analyze_p48_retrospective.tail_and_derivative`), not on `main`. This audit was
run on the PR #564 head; every "clean" below is conditional on that fix landing.

## 1. Was any geometry rejected as "too large" because of the silent zero?

No size decision needs revisiting. Three candidate ceilings were checked:

- **N=1300 ladder (v1/v2 → N=580).** The pilot surfaced the underflow *as a
  visible failure* (the analysis path returned zeros at N=1300; `project_size`
  walked off `(0,1)` and raised — the only reason the bug was found at all).
  But N=1300 was abandoned on **affordability**, not on the bug: even after the
  fix, the frozen yaml states "a decisive ratio there is about three orders of
  magnitude beyond what we can spend". The N=580 choice is driven by "maximum
  shared leverage first, then minimum leakage", all below the ~790 ceiling. The
  silent zero capped a *code path*, not a *decision*.
- **ROADMAP item 4 / #155.** Stops adding replicas to N130/N170 because the
  response matrices "remain nearly rank-one" and "R=8 is non-injective on those
  tori" — an injectivity/geometric reason, not a numerical one.
- **#31.** A high-memory transfer-matrix resource probe — bounded by memory,
  not by underflow.

## 2. Is the same `q**n` / factorial-start pattern anywhere else?

Grep over `scripts/` for `** n`, `**(n`, `math.factorial`, `comb(`, and
`q**n`/`(1-p)**n` start terms. The result is one latent duplicate, everything
else is safe:

- **`score_issue43_full_curve.py::tail`** (float `q**n` recurrence, exactly the
  fixed bug's shape). **Latent, not active**: its `DESIGNS` are N=185 and 265,
  so `0.407**n ≈ 1e-72 … 1e-103`, far above the float64 floor (~1e-308). It
  would underflow only if issue #43 were ever extended past ~790 sites.
- **mpmath tails — safe.** `analyze_threshold_ranks.py::matching_value` /
  `matching_derivative`, `analyze_matching_parity_derivatives_fast.py::tail_*`,
  `threshold_score_modes.py::binomial_weights` all compute `q**n` in `mp.mpf`,
  which has a huge exponent range and does not underflow.
- **Exact rational — safe.** `rigorous_pc_confidence_gate.py` uses `Fraction`;
  `essential_birth_histogram_analysis.py` uses `Decimal`.
- **q-series — safe.** The `q**n` terms in `e4_balanced_pell_estimator.py`,
  `pinson_arguin_kdv.py`, `derive_hexagonal_degree2_hecke.py`,
  `derive_q4_jordan_log_slope_shape.py`, `hexagonal_pell_spin_filter.py`,
  `modulus_shape_discrimination.py` are modular-form series with `|q|<1`,
  summed over a small index — convergent, not a site-count walk.
- **Small-`n` combinatorics — safe.** `comb(...)`/`math.factorial` in
  `observer_bandwidth_*`, `bernoulli_*`, `terminal_reliability_polynomial.py`,
  `logodds_derivative_chain.py` etc. run at path/motif/order scale, not site
  count.

## 3. Above p ≈ 0.9 the old path was already inexact — who keeps that value?

No consumer keeps a high-p value. `project_size` bisects on
`mean_matching(midpoint) < 0.0` — it reads only the **sign** in the upper tail,
and the reported channels are evaluated at the root `p0 ≈ 0.593`. The
`8.9e-7` relative inexactness at `p=0.92` never lands in a kept number.

## 4. Which committed artifacts no longer reproduce bit-for-bit?

`tail_and_derivative` is imported by three scorers besides P49:
`score_axis_pair_annihilator.py`, `score_intrinsic_functional_cocycle.py`,
`score_v14_scalar_root_projector.py`. None of their **committed** outputs is a
tail-derived numerical reveal, so nothing has gone stale to measure:

- `results/server-20260829/P57-norm5-500m/functional_cocycle_score.json` —
  `classification: "frozen pre-target functional-cocycle scorer"` (levels /
  covariance groups / model order). A frozen scorer config, no tail-derived
  numbers.
- `results/v14-scalar-retrospective/fixedp_scalar_projector.csv` — REPORT.md
  says "Status: retrospective design/power analysis only. No claim upgrade."
  Design/power numbers, explicitly not a claim.
- axis_pair has only `predictions/axis_pair_annihilator_semantic_gate_20260830.json`
  (a structural gate) and notes; no scored artifact.

Their `_typed` tests pass precisely because they inject a `mock.Mock()`
calculator and never run the real numerics — confirming the ticket's
observation. The only committed tail-derived reveal in the repository remains
`results/server-20260828/P49-fullcurve-doubling-100m/analysis/score.json`, whose
move is already measured leaf-by-leaf (565/598 leaves moved, worst `6.8e-9`
relative, worst χ² move `2.4e-10`) and recorded under `replay_reproducibility`.

## Conclusion

Nothing was affected. No size decision was capped by the silent zero, no other
active consumer carries the same float underflow, and no committed artifact of
the three sibling scorers has gone stale.

One flag for later (no fix applied — this ticket says "not a rewrite"): if
`score_issue43_full_curve.py` is ever run past ~790 sites, its float `tail`
will underflow the same way. A one-line mode-anchored recurrence, identical to
the P48 fix, removes that latent hazard whenever issue #43 is revisited.
