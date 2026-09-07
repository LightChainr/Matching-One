# P3 — Test the ray, not the ratio

Draft manuscript workspace for portfolio track **P3**, arising from
[issue #579](https://github.com/LightChainr/Matching-One/issues/579).

**Target venue:** *Physical Review E*.

| File | Role |
|---|---|
| [`manuscript.md`](manuscript.md) | the draft, sections 1–8 |
| [`tables.md`](tables.md) | **generated** — every numerical table (T1–T5) |
| `results/p3-projective-inference-manuscript/latest.json` | the machine-readable evidence artifact |
| `results/aspect-ladder-n580/latest.json` | the frozen Fieller scoring (§4.1, §4.3) |
| `results/aspect-ladder-n580-projective/latest.json` | the projective rescoring (§4.4–§4.6) |
| `predictions/aspect_ladder_n580_20260905.yaml` | the frozen design, competitor rays, declared systematic |
| `scripts/projective_inference.py` | the statistic: pseudo-inverse, subspace residual, χ² tail |
| `scripts/aspect_ladder_projective_rescore.py` | the N=580 rescoring |
| `scripts/p3_manuscript_evidence_table.py` | assembles the artifact and renders the tables |
| `tests/test_projective_inference.py` | 8 tests, including the Fieller identity anchor |
| `tests/test_aspect_ladder_projective_rescore.py` | 11 tests on the rescoring |
| `tests/test_p3_manuscript_evidence_table.py` | 11 tests locking the assembly and the tables |
| `notes/aspect-ladder-n580-projective-rescore-20260906.md` | the working note the draft is built from |

## Ground rules honoured here

No Monte Carlo was run for this manuscript and no frozen design was reopened. `p3_manuscript_evidence_table.py`
reads three committed artifacts and derives only statements that follow exactly from them; it rescores nothing
that was not already rescored and committed under #579. The one computation it performs that is not a lookup is
the two-entry control of §4.3 — the projective statistic restricted to the pair of rungs the frozen test used —
and that exists to *remove* a possible explanation for the paper's headline, not to supply one.

No number in the draft is hand-typed into a generated table. `tables.md` is rendered from the artifact and
`test_the_rendered_tables_do_not_drift_from_the_artifact` fails if the two separate. The summary tables inside
`manuscript.md` §4.4–§4.6 restate T3–T5 for readability; the generated tables are authoritative.

## Section readiness

| § | Content | Status |
|---|---|---|
| Abstract | claim and scope | draft complete |
| 1 | the problem: proportions, ratios, conditioning, the design pathology | complete |
| 2.1–2.4 | ray, subspace, residual distance, Fieller as the 2×1 case | complete, proposition proved and tested |
| 2.5 | annihilating functionals | complete |
| 3 | why the ratio framework fails invisibly | complete |
| 4.1–4.2 | the measurement and the eight rays | complete, generated |
| 4.3 | two-entry control | complete, generated, `1.5e-15` |
| 4.4 | verdict table | complete, generated |
| 4.5 | curvature | complete, generated |
| 4.6 | the priced assumption | complete, generated |
| 5 | established / undetermined / not established | complete |
| 6 | the N=650 design consequence | complete; see [#583](https://github.com/LightChainr/Matching-One/issues/583) |
| 7 | reproducibility and preconditions | complete |
| 8 | related work | **one open item**: reference [6] is marked **[LIT]** and needs a primary reading |

## One correction already applied to this draft

The `|A₈/A₄|` column of §4.6 first shipped as the raw gap divided by the leakage coefficient, which is the
leading-order form and drops the `(u + 1)` the exact solution carries. It overstated the requirement by 2.4× for
the bare aspect ratio and by 45× for the weight-12 rays, where it read 783 instead of 17.5. The frozen design
records the leakage *per rung with its sign*, which is what makes the exact solve available at all; the fix is in
`required_spin8_ratio`, and `test_every_survivor_of_the_clean_pair_needs_a_large_spin8` now pins the top of the
column below 20 so the leading-order form cannot come back silently. The dichotomy is unchanged.

## Open items before submission

1. **Reference [6].** The claim that ratio-testing of aspect-ratio amplitudes is standard practice in the
   conformal-invariance finite-size literature is stated from familiarity, not from a primary reading. Either
   read and cite specific instances, or weaken §8 to a statement about this repository's own practice. This is
   the only unresolved item in the draft and it is deliberately marked rather than quietly asserted.
2. **§5.2's undetermined verdict.** One deterministic replay of the N=580 ladder under the current scorer
   recovers `cov(r2, r4)` and settles `bare_aspect_ratio`. The draft is publishable without it — an undetermined
   verdict, reported as undetermined, is a legitimate result — but the replay is cheap and would remove the one
   soft spot.
3. **Figure.** The draft is currently table-only. One figure would carry §4.4–§4.5 well: the three-rung response
   with its jackknife errors, the eight competitor rays normalised through the `r=1` point, and the measured
   concavity as a shaded band. Not required for the argument.
