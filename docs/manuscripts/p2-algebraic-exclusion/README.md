# P2 — Certified exclusion of low-degree algebraic threshold relations

Draft manuscript workspace for portfolio track **P2**, requested by
[issue #551](https://github.com/LightChainr/Matching-One/issues/551).

| File | Role |
|---|---|
| [`manuscript.md`](manuscript.md) | the draft, sections 1–9 per the issue outline |
| [`tables.md`](tables.md) | **generated** — every numerical table in the draft |
| `results/p2-algebraic-exclusion-manuscript/latest.json` | the machine-readable evidence artifact |
| `results/pslq-degree4-synthetic-boundary-control/latest.json` | quartic sensitivity control (§6.4) |
| `results/pslq-degree6-low-height-*/latest.json` | degree-1..6 height-3 exhaustion, per interval (§6.5) |
| `results/pslq-degree6-low-height-control/latest.json` | planted `(3,12²)` sensitivity control (§6.5) |
| `results/pslq-degree6-low-height-replication-*/latest.json` | the same census by a second implementation (§7) |
| `results/pslq-degree6-implementation-agreement/latest.json` | cell-by-cell comparison of the two (§7) |
| `scripts/p2_manuscript_evidence_table.py` | assembles the artifact and renders the tables |
| `scripts/degree4_synthetic_boundary_control.py` | runs the quartic sensitivity control |
| `scripts/degree6_low_height_exclusion.py` | runs the historical-range exhaustion |
| `scripts/degree6_low_height_control.py` | runs its planted-root control |
| `scripts/degree6_independent_replication.py` | the second implementation, committed as received |
| `scripts/degree6_implementation_agreement.py` | compares the two implementations |
| `tests/test_p2_manuscript_evidence_table.py` | locks the assembly against the census artifacts |
| `tests/test_degree4_synthetic_boundary_control.py` | locks the quartic sensitivity control |
| `tests/test_degree6_low_height_exclusion.py` | locks the historical-range exhaustion and its control |
| `tests/test_degree6_implementation_agreement.py` | locks the two-implementation agreement |

## Ground rules honoured here

The quartic census was not re-run and no degree or height was expanded beyond the frozen classes.
`scripts/p2_manuscript_evidence_table.py` reads the frozen contract, the provenance manifest and the committed
census artifacts, re-verifies the provenance digest, and derives only statements that follow exactly from those
inputs — it computes no census itself. Root decisions reuse the repository's existing exact Sturm path
(`scripts/exact_polynomial_root_certificate.py`). No number in the draft is hand-typed: `tables.md` is rendered,
and a regression test fails if it drifts from the artifact.

Three computations were added, each with its reason for being inside the stop rule:

| Computation | Why it is not a census extension | Cost |
|---|---|---:|
| §6.4 quartic sensitivity control | synthetic control on synthetic intervals, running the **unmodified** `degree4_interval_exclusion.run_search`. Covers exactly the widths whose census result was a null, read from the artifacts rather than hardcoded | ~24 s |
| §6.5 degree ≤ 6 height ≤ 3 exhaustion | a different frozen class, not a widening of the quartic one; it closes the last form in the historical tradition (issue #559) | ~4 s |
| §7 second implementation and its comparison | re-decides the §6.5 class, adding no new class and no new interval; its value is entirely in being written separately | ~50 s to regenerate, 12 s in CI |

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
| 6.4 | Quartic sensitivity control | ready | **new result and new computation** — Table 7, Result E |
| 6.5 | Historical complexity range closed | ready | **new result and new computation** — Table 8, Results F and G |
| 7 | Calibration | ready | two gaps closed — the second partly; one recommendation left open, and it is the more valuable one |
| 8 | Discussion and scope | ready | |
| 8.1 | Future work | ready | the costed recommendation was executed; what remains is genuinely harder |
| 9 | Reproducibility supplement | ready | artifact list and digests generated |
| — | Target venue | decided | J. Phys. A, fallback Experimental Mathematics — reasoning in the manuscript |

## Five results derived while assembling

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

**5. The historical complexity range is closed.**

Every exactly-known planar threshold has degree ≤ 6 and height ≤ 3, and the one form outside `C(≤4, ≤100)` is the
`(3,12²)` site value `x⁶ − 3x⁴ + 1`. Exhausting all `409,584` primitive polynomials of degree ≤ 6 at height ≤ 3
excludes every one of them on all four intervals — the certified screen retains **zero** candidates at any degree
on any interval, and the closest member of the class stays `9.23e-8` away, between `5.8e2` and `2.3e6` interval
widths. Because the screen retains nothing, the Sturm path never runs during the exclusion, so the planted-root
control matters as much as it did at degree 4: planting `x⁶ − 3x⁴ + 1` itself at each frozen width gives 8/8,
with every positive trial retaining exactly one candidate and reporting it, and every negative trial retaining
the same candidate and correctly returning no root.

Together with Results A–C this supports the paper's cleanest sentence: *no algebraic form at the complexity of
any exactly-known planar percolation threshold has a root in any of the four published intervals.*

## Decisions taken while drafting

- **Ran** the degree-4 sensitivity control. It closes the one calibration gap that directly threatens Results A
  and B, it is a synthetic control rather than a census extension, and it costs 24 seconds.
- **Scoped that control to the null-result widths only.** An earlier version covered all four widths (16 trials,
  88 s). CI runs the full suite on three Python versions under a 20-minute per-job timeout, and the Python 3.9
  job is the binding constraint: across three `main` runs it took 996 s, 1000 s and 770 s of that 1200 s budget,
  roughly 40 % slower than 3.11/3.13 on this Fraction-heavy suite. A 16-trial control would have added about
  23 %, putting the typical 3.9 run near 1230 s — over the timeout. The rescope was therefore necessary, not
  merely prudent. It also costs nothing scientifically: the census demonstrates its own sensitivity on the two
  intervals where it found roots, so the question only has force on the other two.

  Measured on this branch: 3.9 1037 s, 3.11 722 s, 3.13 728 s — about +4 % against the `main` cluster on each
  version. **The narrow 3.9 margin is pre-existing on `main`, not introduced here**, and is worth a separate
  look; it is out of scope for this ticket.
- **Ran** the degree ≤ 6, height ≤ 3 exhaustion (#559) after the owner's call. It was initially deferred under
  #551's sequencing rule, then filed as a ticket, then run directly once it was clear the work is 3.6 s of
  single-core arithmetic and needs no remote machine. It closes the last uncovered form in the historical
  tradition, so the paper no longer has to name that gap as a limitation.
- **Did not raise the height** on the null result, per §4.3. Degree ≤ 6 at height ≤ 10 is `890,350,944`
  polynomials per interval and, more importantly, would improve the class's approach resolution toward the
  interval widths; the boundary-degree check has to come first or the null is guaranteed and empty.
- **Imported** the second implementation of the §6.5 census rather than only citing it. It was written
  independently against the same frozen protocol, and it screens at the interval midpoint where the primary
  implementation screens at both endpoints. Committing it, its four artifacts and a cell-by-cell comparison turns
  "a second implementation agrees" from a remark into a checked claim that CI re-derives. The comparison also
  checks the two residuals against the mean value bound, which is what rules out their agreeing by both being
  trivially empty.

  **Scoped the CI cost.** The replication census is about 12 s per interval of pure-Python rational arithmetic.
  Recensusing all four in the test suite would cost roughly 70 s on the Python 3.9 job, which is the binding one
  (see above). The test therefore recensuses one interval — the narrowest, carrying Result A — and covers the
  other three through the agreement check against the primary implementation, which CI already rebuilds in full.
  Measured cost of the new module: 12 s here.
- **Did not overstate what the replication covers.** It reaches the degree ≤ 6 height ≤ 3 census and not the
  quartic census, and both implementations share the Sturm code unchanged. That shared code contributes nothing
  to *this* null — the screens retain zero candidates, so isolation never runs — but it does run in the planted
  controls. §7, Table 9 and the artifact's `claim_boundary` all say so rather than leaving it to be inferred.
- **Foregrounded** §2.1 and §4.3, and added §1.3 Contributions. They were not in the requested outline, but they
  generalize beyond this constant and are, in my judgement, the more durable part of the paper.

## Open items

1. **Primary-source verification** for Scullard 2006, Ziff 2006 and Suding–Ziff 1999. arXiv was unreachable from
   the drafting environment, so they are cited from secondary indexing and flagged `PENDING` in
   `references.bib`, in line with `data/README.md`'s provenance rules. Sykes–Essam 1964 is verified. None of the
   paper's results depend on them; they motivate the search class in §1.1 and §4.2, and §4.2's table is
   generated from the repository's own certified artifact rather than from those citations.
2. **One remaining calibration recommendation** (§7): interval-perturbation sensitivity for the degree-4 near
   hits. The second-implementation recommendation is met for the degree ≤ 6 height ≤ 3 census only (Table 9); the
   quartic census of §6.1–6.3, where Results A–D live, still has a single implementation of both its C++ screen
   and its Sturm decisions, and replicating *that* is the more valuable of the two. Neither blocks submission.

## Regenerating

```bash
python3 scripts/degree4_synthetic_boundary_control.py \
    --output results/pslq-degree4-synthetic-boundary-control/latest.json
python3 scripts/degree6_low_height_exclusion.py --all
python3 scripts/degree6_low_height_control.py \
    --output results/pslq-degree6-low-height-control/latest.json
for interval in jacobsen-2015-eigenvalue mertens-2022-p-med \
                mertens-2022-p-cell yang-zhou-2024-corrected; do
    python3 scripts/degree6_independent_replication.py "$interval" \
        --output "results/pslq-degree6-low-height-replication-$interval/latest.json"
done
python3 scripts/degree6_implementation_agreement.py
python3 scripts/p2_manuscript_evidence_table.py \
    --output results/p2-algebraic-exclusion-manuscript/latest.json
python3 scripts/p2_manuscript_evidence_table.py --markdown \
    --output docs/manuscripts/p2-algebraic-exclusion/tables.md
python3 -m unittest tests.test_p2_manuscript_evidence_table \
    tests.test_degree4_synthetic_boundary_control \
    tests.test_degree6_low_height_exclusion \
    tests.test_degree6_implementation_agreement
```
