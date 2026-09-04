# P2 — Certified exclusion of low-degree algebraic threshold relations

Draft manuscript workspace for portfolio track **P2**, requested by
[issue #551](https://github.com/LightChainr/Matching-One/issues/551).

| File | Role |
|---|---|
| [`manuscript.md`](manuscript.md) | the draft, sections 1–9 per the issue outline |
| [`tables.md`](tables.md) | **generated** — every numerical table in the draft |
| `results/p2-algebraic-exclusion-manuscript/latest.json` | the machine-readable evidence artifact |
| `results/pslq-degree4-synthetic-boundary-control/latest.json` | the new sensitivity control (§6.4) |
| `scripts/p2_manuscript_evidence_table.py` | assembles the artifact and renders the tables |
| `scripts/degree4_synthetic_boundary_control.py` | runs the sensitivity control |
| `tests/test_p2_manuscript_evidence_table.py` | locks the assembly against the census artifacts |
| `tests/test_degree4_synthetic_boundary_control.py` | locks the sensitivity control |

## Ground rules honoured here

No census was re-run and no degree or height was expanded.
`scripts/p2_manuscript_evidence_table.py` reads the frozen contract, the provenance manifest and the committed
census artifacts, re-verifies the provenance digest, and derives only statements that follow exactly from those
inputs. Root decisions reuse the repository's existing exact Sturm path
(`scripts/exact_polynomial_root_certificate.py`). No number in the draft is hand-typed: `tables.md` is rendered,
and a regression test fails if it drifts from the artifact.

The one new computation is the §6.4 sensitivity control. It is a synthetic control, not a census extension: it
runs the **unmodified** `degree4_interval_exclusion.run_search` on synthetic intervals, and `GOVERNANCE.md`
places synthetic controls outside the production queue. It covers exactly the intervals whose census result was
a null — read from the census artifacts, not hardcoded — and takes about 24 seconds.

## Section readiness

| § | Section | State | Note |
|---|---|---|---|
| — | Abstract | draft ready | numbers generated; wording is a first pass |
| 1.1 | Exact thresholds in 2D | ready, partly unverified | Sykes–Essam verified; Scullard / Ziff / Suding–Ziff cited from secondary indexing, marked in-text and in `references.bib` |
| 1.2 | Why heuristic searches are insufficient | ready | |
| 1.3 | Contributions | ready | **new section** — foregrounds the four methodological results |
| 1.4 | Scope statement | ready | mirrors the issue's explicit non-claims |
| 2 | Canonical provenance table | ready | Table 1 |
| 2.1 | Interval disjointness | ready | **new result** |
| 2.2 | Provenance limitations | ready | Yang–Zhou tables and Jacobsen 2024 Reply left pending, as in the manifest |
| 3 | Completeness theorem | ready | Prop. 1 and Thm. 2 stated and proved from the committed implementation |
| 4.1 | What the bounds must not be | ready | |
| 4.2 | Complexity of the known exact thresholds | ready | **literature pass done** — Table 6, generated from the certified lattice-native artifact |
| 4.3 | Approach resolution / boundary degree | ready | **new result** |
| 5 | Method | ready | |
| 6.1–6.3 | Results and controls | ready | Tables 3 and 5; Results A–D |
| 6.4 | Sensitivity control | ready | **new result and new computation** — Table 7, Result E |
| 7 | Calibration | ready | one gap closed, two recommendations left open |
| 8 | Discussion and scope | ready | |
| 8.1 | Future work | ready | costed recommendation, deliberately not executed |
| 9 | Reproducibility supplement | ready | artifact list and digests generated |
| — | Target venue | **needs a call** | see below |

## Four results derived while assembling

All follow from committed artifacts and are locked by tests.

**1. The four method intervals are pairwise disjoint.**

```text
mertens-2022-p-cell < mertens-2022-p-med < yang-zhou-2024-corrected < jacobsen-2015-eigenvalue
gaps: 1.03e-10, 5e-13, 2.38e-12
```

At their own quoted precisions the four published estimates are mutually inconsistent, so at least three of the
four intervals do not contain `p_c`. This upgrades the repository's existing "do not pool intervals" policy from
a conservatism to a necessity, and it is what makes the 16 surviving quartics interpretable.

**2. Every surviving quartic survives exactly one interval.**

`max_intervals_per_survivor = 1`, with certified separations of `2.4e-12` (the `p_med` survivor) and
`1.2e-10`–`2.6e-10` (the 15 `p_cell` survivors) from every other interval. The survivors are width artifacts,
not candidate formulas.

**3. Degree 4 is the boundary degree.**

The certified approach resolution of the search class falls `1.53e-4 → 6.98e-8 → 7.10e-10 → <1e-12` across
degrees 1 to 4, while the interval widths span `4e-14` to `1.6e-10`. Degree 4 is the first degree whose closest
member no longer clears one interval width, and it is exactly where survivors appear, on the two widest intervals
only. This turns the ticket's stop rule into a quantitative argument and answers the predictable reviewer
question *why stop at degree 4?*

**4. The nulls are sensitivity-certified.**

Planting a committed quartic root witness inside synthetic intervals of the two widths where the census returned
zero, and re-running the unmodified census path: 8/8 trials pass — every positive trial recovers the planted
quartic, including at the narrowest width (`4e-14`), and no negative trial reports it. The zero-survivor results
on the Jacobsen and Yang–Zhou intervals are therefore certified nulls rather than blind spots. This closes the
calibration gap flagged in §7 and, as far as I can tell, is a control the integer-relation literature does not
routinely perform.

The control deliberately does *not* cover the `p_med` and `p_cell` widths: the census itself found 1 and 15 roots
there, so its sensitivity at those widths is already demonstrated and the question only has force where the
answer was zero. This also keeps the control cheap enough for CI — see below.

## Decisions taken while drafting

- **Ran** the degree-4 sensitivity control. It closes the one calibration gap that directly threatens Results A
  and B, it is a synthetic control rather than a census extension, and it costs 24 seconds.
- **Scoped that control to the null-result widths only.** An earlier version covered all four widths (16 trials,
  88 s). CI runs the full suite on three Python versions under a 20-minute per-job timeout, and the slowest job
  on `main` already uses 770 s of that budget, so a 16-trial control would have left roughly 250 s of headroom.
  Restricting to the two widths where the census reported a null costs nothing scientifically — the census
  demonstrates its own sensitivity on the other two — and is the better-argued design regardless.
- **Did not run** the degree ≤ 6, height ≤ 3 census, even though it is only 409,584 polynomials (a few seconds)
  and would close the last uncovered form in the historical tradition — the `(3,12²)`-type radical. Issue #551
  sequences write-up and review before any degree or height expansion. The hypothesis and its exact cost are
  recorded in the artifact and in §8.1, so the call can be made immediately with the number in hand. **This is
  the one substantive decision I would flag for a second opinion**: the argument for running it now is that a
  referee will ask, and the argument against is that the ticket explicitly said not to.
- **Foregrounded** §2.1 and §4.3, and added §1.3 Contributions. They were not in the requested outline, but they
  generalize beyond this constant and are, in my judgement, the more durable part of the paper.

## Open items

1. **Target venue.** J. Phys. A gives the result its audience; a methods or computational-number-theory venue is
   the better home for Theorem 2, the boundary-degree criterion and the sensitivity control. Since §8's
   methodological contributions are the durable part, I lean methods-oriented, but this needs whoever owns
   literature search.
2. **Primary-source verification** for Scullard 2006, Ziff 2006 and Suding–Ziff 1999. arXiv was unreachable from
   this environment, so they are cited from secondary indexing and flagged as `PENDING` in `references.bib`, in
   line with `data/README.md`'s provenance rules. Sykes–Essam 1964 is verified.
3. **Two remaining calibration recommendations** (§7): an independent second implementation of the exact filter,
   and interval-perturbation sensitivity for the degree-4 near hits. Neither blocks a submission draft.

## Regenerating

```bash
python3 scripts/degree4_synthetic_boundary_control.py \
    --output results/pslq-degree4-synthetic-boundary-control/latest.json
python3 scripts/p2_manuscript_evidence_table.py \
    --output results/p2-algebraic-exclusion-manuscript/latest.json
python3 scripts/p2_manuscript_evidence_table.py --markdown \
    --output docs/manuscripts/p2-algebraic-exclusion/tables.md
python3 -m unittest tests.test_p2_manuscript_evidence_table tests.test_degree4_synthetic_boundary_control
```
