# P3 — Covariance-aware comparison of finite-size amplitude laws

**Current status (2026-09-12): numerical corrections available; manuscript NOT
submission-ready.** Read [CORRECTIONS-20260912.md](CORRECTIONS-20260912.md)
before using the historical manuscript. Target venue remains a proposal, not
an acceptance/readiness claim. Track: #579.

## Current result

#703 recovered the complete N580 covariance from the existing delete-one shards,
without simulation. The retrospective three-rung pure-line score for bare aspect
is 10.8644599 on 2 df, nominal p=.00437333: not rejected at the declared nominal
3-sigma threshold. Seven other pure lines are rejected under that contract.
The original frozen two-rung underpowered result remains unchanged.

A second completed analysis allows independently varying rung-specific H8/H4
ratios. At a uniform bound B=1, the same seven candidate images remain outside
one nominal 99.73% three-dimensional mean ellipsoid; bare aspect intersects.
This is not a fitted-cone chi-square calibration, a measurement of H8, or a new
independent evidence block.

## Files to use

- `notes/n580-complete-covariance-20260912.md` and
  `results/research-control-20260912/n580-complete-covariance-summary.json`:
  executed covariance recovery and pure-line results.
- `notes/n580-rungwise-leakage-20260912.md` and
  `results/research-control-20260912/n580-rungwise-leakage-summary.json`:
  independently bounded leakage, projection proof, results and limits.
- `scripts/projective_inference.py`: explicit covariance-support policy;
  `scripts/recover_n580_covariance.py`: existing-shard reader;
  `scripts/n580_rungwise_leakage.py`: bounded-image calculation.
- `manuscript.md`, `tables.md`, and the original manuscript-evidence JSON:
  HISTORICAL draft/evidence assembly, retained rather than silently overwritten.
  Their readiness, missing-covariance, common-ratio and sign claims require the
  corrections in the addendum. The old assembler does not generate the new tables.
- `predictions/aspect_ladder_n580_20260905.yaml` and
  `results/aspect-ladder-n580/`: original freeze/results; do not rewrite history.

## Remaining writing work, not more N580 production

Integrate the addendum into the manuscript and a new versioned table assembly.
Separate the frozen two-rung analysis, retrospective pure-line comparison and
rungwise bounded sensitivity. Replace unsupported literature-practice/novelty
claims with the narrow sourced account. The existing #695/#701 retrieval is
input, not proof that all citations were read in full. N650 remains subject to
#589's angular/Smith identifiability question.

No replay is needed for cov(r2,r4). No new size, free exponent or optional
sample top-up follows from this draft.
