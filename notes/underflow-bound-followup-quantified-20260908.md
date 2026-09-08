# Underflow-bound follow-up — the audit re-run on merged main, with one gap quantified

Re-audit of `LightChainr/Matching-One#570`, run 2026-09-08. This note supersedes
the *conditional* framing of `notes/underflow-bound-audit-20260905.md` in one
respect and extends it in another:

1. **The condition resolved.** PR #564 (commit `844173c`, the mode-anchored
   recurrence in `analyze_p48_retrospective.tail_and_derivative`) was merged via
   PR #573 on 2026-09-07 and is now on `main`. The 2026-09-05 audit's "every
   'clean' below is conditional on that fix landing" no longer needs the caveat.
2. **The one gap in that audit is now measured, not just argued.** The 2026-09-05
   note argued from artifact *classification* (frozen config, design analysis,
   structural gate) that no committed output of the three sibling scorers was
   tail-derived and therefore nothing went stale. For one of them that argument
   can be replaced by a measurement — see §4. The verdict does not change, but
   it is now a replay result rather than a category claim.

Everything else in the 2026-09-05 note was re-checked on current `main`
(`8b5bc84`) and stands as written: no geometry was rejected on the silent zero
(N=1300 died on affordability, ROADMAP item 4 / #155 on injectivity, #31 on
memory), the only latent float-recurrence duplicate is
`score_issue43_full_curve.py::tail` (safe at its N=185/265 designs by ~5 and ~7
orders of magnitude of exponent headroom), the high-p inexactness lands only in
sign reads inside bisections, and the mpmath / `Fraction` / `Decimal` / q-series
hits in the grep are all out of underflow's reach.

## 4 (extended). Which committed artifacts moved, measured

### `results/server-20260828/P43-heldout-fullcurve-500m/analysis/primary_score.json` — bit-identical

`score_issue43_full_curve.py` carries its own private copy of the old
float-`q**n` recurrence (`tail`, line ~195), so it is *not* fixed by #564 — but
it is also *not moved* by it, and it does not need fixing at its designs. Two
measurements:

- **Replay.** Running the scorer on the committed N=185/265 raw triples
  reproduces the committed `primary_score.json` **byte-for-byte** (sha256
  `3a545e6c…`, 83/83 leaves identical). Nothing went stale.
- **Numerical distance to exact.** Against a `Fraction`-exact binomial-tail
  evaluation of the same histograms at the frozen `p_ref = 0.592746050790`, the
  worst per-histogram relative error of the float recurrence is
  `1.9e-15` at N=185 and `3.2e-15` at N=265 — pure double-rounding, four orders
  of magnitude below the `6.8e-9` worst leaf move the P49 gate measured for the
  *reassociated* sum. The upper-tail inexactness the ticket cites for the P48
  routine (`8.9e-7` at p=0.92) never arises here: the scorer evaluates only at
  `p_ref`.

The latent hazard stands as flagged on 2026-09-05: extend issue #43 past ~790
sites and this private `tail` underflows exactly like the P48 bug did. One-line
mode-anchored fix, same as #564, when that day comes.

### `results/server-20260829/P50-n145-n290-fullcurve/analysis/score.json` — moved, verdicts intact

This is the gap. `score_p50_fullcurve_n290.py` imports `project_size` and the
P49 helpers (`orientation_values`, `solve_target`), so the committed score **is**
tail-derived and the 2026-09-05 note's classification argument did not cover it.
The ticket predicted exactly this shape of outcome ("every artifact derived from
`tail_and_derivative` moves in its last few digits"), so it was measured rather
than regenerated:

- Replay through the merged (mode-anchored) routine: **72 of 104** float leaves
  moved, worst relative move `1.25e-3` — and that worst leaf is a
  `correlation_eigenvalues[0]` diagnostic, not a decision.
- Every verdict-bearing leaf in the frozen REPORT.md table is intact to
  reporting precision, i.e. the move is far below the printed digits:

  | frozen leaf | committed | replay |
  |---|---:|---:|
  | joint ΔM χ² (2 dof) | 9.35200 | 9.3520036 |
  | raw slope z | −22.6903 | −22.69035 |
  | scalar+H4 corrected slope z | −0.66608 | −0.6660798 |
  | raw root-ratio z | 1.38768 | 1.3876824 |
  | frozen induced root-ratio z | 1.38520 | 1.3851962 |

  χ² `9.35200` against 2 dof sits at p ≈ 0.0093 either way; every sign, every
  σ-distance, and all four "verdict" columns of the report are unchanged.

No regeneration was done and none is needed: the committed file stays as the
dated record, and this note is the measured erratum. If the maintainers want
`replay_reproducibility` blocks like P49's for P50, that is a separate ticket —
this one says only what moved and by how much.

The remaining two siblings stay as classified on 2026-09-05:
`results/server-20260829/P57-norm5-500m/functional_cocycle_score.json` is a
frozen scorer configuration (sizes ≤ 425, far under the ceiling), and
`results/v14-scalar-retrospective/` is an explicitly non-claim design/power
analysis; neither has a tail-derived numerical reveal to move.

## Conclusion

Same verdict as 2026-09-05, now unconditional and with the P50 gap closed by
measurement: no size decision was capped by the silent zero, no other *active*
consumer underflows at any size the project runs or plans, and the only
committed artifact that moved (P50 `score.json`, 72/104 leaves, worst `1.25e-3`
on a non-decision diagnostic) moves none of its verdicts. The two latent items
for whenever they mature: the private `tail` in `score_issue43_full_curve.py`
(past ~790 sites), and an optional P50 `replay_reproducibility` block to match
P49's.
