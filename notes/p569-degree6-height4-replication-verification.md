# Issue #569 — Height-4 extension of the degree-≤6 second implementation: independent verification

Branch: `replication/p569-degree6-height4-20260913`
Caller: compute1 (an outside caller relative to the author of `scripts/degree6_independent_replication.py`)

## Assignment premise vs. repository state — honesty disclosure

The ticket states the task as: only the artifact field at ~line 274 hardcodes the
module constant `HEIGHT`; thread `height` through to it; run all four intervals at
height 4; extend `scripts/degree6_implementation_agreement.py` to height 4.

On the current `main` (HEAD `d31fa5f`) the following are **already present and
tracked in git** — they were NOT produced by this assignment:

* `scripts/degree6_independent_replication.py:288` already reports
  `"coefficient_height_max": height` (the parameter, not the `HEIGHT` constant).
  `_enumerate_degree`, `class_size`, `_degree_polynomials`, `run_search` and
  `build_result` all already take a `height` argument. **No edit to the script was
  required or made by this caller.**
* `results/pslq-degree6-height4-{id}/latest.json` (primary) and
  `results/pslq-degree6-height4-replication-{id}/latest.json` (second implementation)
  exist for all four frozen intervals:
  `jacobsen-2015-eigenvalue`, `mertens-2022-p-cell`, `mertens-2022-p-med`,
  `yang-zhou-2024-corrected`.
* `scripts/degree6_implementation_agreement.py` already enumerates
  `HEIGHT_DIRS = {3: ..., 4: ...}`; its committed `latest.json` reports
  `heights_compared: [3, 4]`, `cells_compared: 48`, `cells_in_agreement: 48`.

So the ticket's stated premise — that an outside caller must mechanically thread
`height` and that the height-4 extension is absent from `main` — does **not** hold in
this tree.

## Disclosure explicitly required by the ticket

The ticket instructs the caller to declare that this height threading was a mechanical
change made by a party *outside* the implementation's author, weakening the
independent-author property. Because the threading was **already present in `main`**
(the artifact field already reports the parameter, and the height-4 artifacts were
already committed there), that declaration would be false if applied to this caller.

Honest statement: this caller (compute1) did **not** edit
`scripts/degree6_independent_replication.py`. The independent-author property of the
code-writing step is therefore intact — the author's code already parameterized
height, and the height-4 extension was already in `main` before this assignment. This
caller's contribution is an **independent re-run** that confirms the committed
artifacts are byte-for-byte reproducible; it is a verification, not an authorship
change, and introduces no weakening of the independent-author property of the existing
code. The maintainers should reconcile this with the ticket's premise (the height-4
extension appears to have been completed in `main` already).

## Independent re-run performed by this caller

Re-ran `scripts/degree6_independent_replication.py <id> --height 4` for all four
intervals on the managed Python (numpy/scipy/mpmath; exact rational arithmetic, no
binary float in any exclusion claim).

| interval                       | coefficient_height_max | class_size_total | excluded | screen_survivors |
|--------------------------------|------------------------|------------------|----------|------------------|
| jacobsen-2015-eigenvalue       | 4                      | 2,351,328        | True     | 0                |
| mertens-2022-p-cell            | 4                      | 2,351,328        | True     | 0                |
| mertens-2022-p-med             | 4                      | 2,351,328        | True     | 0                |
| yang-zhou-2024-corrected       | 4                      | 2,351,328        | True     | 0                |

Per-degree class sizes (identical across all four intervals, independently recounted
by this caller): `[23, 265, 2639, 24913, 229703, 2093785]`; total **2,351,328**.
This matches the ticket's expected numbers exactly. The certified mean-value screen
retains **0** candidates at every degree on every interval, so root isolation never
runs during the census (consistent with the shared Sturm code contributing nothing to
the null).

Re-ran `scripts/degree6_implementation_agreement.py` (output to a scratch path, not
overwriting the committed artifact). Result: `heights_compared: [3, 4]`,
`cells_compared: 48`, `cells_in_agreement: 48`, `implementations_agree: True`; all
eight `(interval, height)` exclusion verdicts are `True` and in agreement between the
two implementations; and the closest-member residual gap falls within the mean-value
allowance for every interval (the check that rules out the two implementations agreeing
merely because both are trivially empty).

## Secondary figure (not independently recomputed to full precision)

The ticket states the closest member of the whole class sits `8.70e-9` from the
nearest interval, i.e. `54` interval-widths at the narrowest interval (versus `5.8e2`
widths at height 3). This is a min-over-interval absolute value, finer than the
`minimum_absolute_residual_at_midpoint` recorded in the artifact. This caller verified
the primary, checkable claims (class sizes and the all-intervals exclusion with 0
survivors) exactly; the `8.70e-9`/`54`-width distance is recorded in the committed
manuscript/artifacts and was not independently recomputed to min-over-interval
precision in this verification run. No inconsistency with it was observed.

## CI note

Full Matching-One repository CI has not been run for this commit.
