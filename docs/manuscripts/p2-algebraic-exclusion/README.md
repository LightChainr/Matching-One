# P2 — Certified exclusion of low-degree algebraic threshold relations

Draft manuscript workspace for portfolio track **P2**, requested by
[issue #551](https://github.com/LightChainr/Matching-One/issues/551).

| File | Role |
|---|---|
| [`manuscript.md`](manuscript.md) | the draft, sections 1–9 per the issue outline |
| [`tables.md`](tables.md) | **generated** — every numerical table in the draft |
| `results/p2-algebraic-exclusion-manuscript/latest.json` | the machine-readable evidence artifact |
| `scripts/p2_manuscript_evidence_table.py` | assembles the artifact and renders the tables |
| `tests/test_p2_manuscript_evidence_table.py` | locks the assembly against the census artifacts |

## Ground rules honoured here

This is a **writing and assembly** task. No census was re-run and no degree or height was expanded.
`scripts/p2_manuscript_evidence_table.py` reads the frozen contract, the provenance manifest and the ten
committed census artifacts, re-verifies the provenance digest, and derives only statements that follow exactly
from those inputs. Root decisions reuse the repository's existing exact Sturm path
(`scripts/exact_polynomial_root_certificate.py`) at 120-bit isolation. No number in the draft is hand-typed:
`tables.md` is rendered, and a regression test fails if it drifts from the artifact.

## Section readiness

| § | Section | State | Note |
|---|---|---|---|
| — | Abstract | draft ready | numbers generated; wording is a first pass |
| 1.1 | Exact thresholds in 2D | **needs literature pass** | survey of exactly-known planar thresholds and duality/star–triangle mechanisms |
| 1.2 | Why heuristic searches are insufficient | ready | argued from the contract's own design |
| 1.3 | Scope statement | ready | mirrors the issue's explicit non-claims |
| 2 | Canonical provenance table | ready | Table 1 from `data/literature_threshold_sources.json` |
| 2.1 | Interval disjointness | ready | **new** — derived here, see below |
| 2.2 | Provenance limitations | ready | Yang–Zhou tables and Jacobsen 2024 Reply left pending, as in the manifest |
| 3 | Completeness theorem | ready | Prop. 1 and Thm. 2 stated and proved from the committed implementation |
| 4.1 | What the bounds must not be | ready | |
| 4.2 | Historically proposed forms | **needs literature pass** | degrees/heights of previously conjectured planar forms |
| 4.3 | Approach resolution / boundary degree | ready | **new** — derived here, see below |
| 5 | Method | ready | |
| 6 | Results | ready | Tables 3 and 5; Results A–D |
| 7 | Calibration | ready as a plan | the three proposed additions are decisions for the owner, not done work |
| 8 | Discussion and scope | ready | |
| 9 | Reproducibility supplement | ready | artifact list and digests generated |
| — | Target venue | **needs literature pass** | cross-check with whoever owns literature search |

## Two results derived while assembling, not previously stated

Both follow exactly from committed artifacts and are locked by tests. Neither required new census computation.

**1. The four method intervals are pairwise disjoint.**

```text
mertens-2022-p-cell < mertens-2022-p-med < yang-zhou-2024-corrected < jacobsen-2015-eigenvalue
gaps: 1.03e-10, 5e-13, 2.38e-12
```

At their own quoted precisions the four published estimates are mutually inconsistent, so at least three of the
four intervals do not contain `p_c`. This upgrades the repository's existing "do not pool intervals" policy from
a conservatism to a necessity, and it is what makes the 16 surviving quartics interpretable: a survivor's status
is decided by which extrapolation one adopts.

**2. Every surviving quartic survives exactly one interval — and degree 4 is the boundary degree.**

Re-deciding all 16 committed root witnesses against all four intervals gives
`max_intervals_per_survivor = 1`, with certified separations of `2.4e-12` (the `p_med` survivor) and
`1.2e-10`–`2.6e-10` (the 15 `p_cell` survivors) from every other interval.

The certified approach resolution of the search class — how close its best member can get without entering —
falls `1.53e-4 → 6.98e-8 → 7.10e-10 → <1e-12` across degrees 1 to 4, while the interval widths span `4e-14` to
`1.6e-10`. Degree 4 is the first degree whose closest member no longer clears one interval width, and it is
exactly where survivors appear, on the two widest intervals only. A degree-5 census at this height would produce
survivors on every interval by counting alone.

This turns the issue's stop rule ("do not increase degree/height by default") from a scope decision into a
quantitative argument the manuscript can make in its own voice, and it is the strongest available answer to the
predictable reviewer question *why stop at degree 4?*

## Open decisions for the owner

1. **Literature passes** for §1.1, §4.2 and the venue choice — flagged above; nothing else blocks a full draft.
2. **Calibration extensions** (§7). Recommended: a degree-4 synthetic positive control at the boundary degree,
   since the existing synthetic calibration is degree-1 only and the boundary degree is where the paper's
   interpretation is load-bearing. Not recommended: any library, degree or height expansion.
3. **Whether §2.1 and §4.3 should be foregrounded.** They were not part of the requested outline. Both are
   general methodological points that transfer to any integer-relation search against a measured constant, and
   they may be worth more prominence than the census counts themselves.

## Regenerating

```bash
python3 scripts/p2_manuscript_evidence_table.py --output results/p2-algebraic-exclusion-manuscript/latest.json
python3 scripts/p2_manuscript_evidence_table.py --markdown --output docs/manuscripts/p2-algebraic-exclusion/tables.md
python3 -m unittest tests.test_p2_manuscript_evidence_table
```
