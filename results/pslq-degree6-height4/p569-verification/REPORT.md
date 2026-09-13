# Issue #569 verification — height-4 extension of the degree-≤6 second implementation

See `notes/p569-degree6-height4-replication-verification.md` for the full report and
the honesty disclosure.

## Summary of independent re-run (caller: compute1)

* Re-ran `scripts/degree6_independent_replication.py <id> --height 4` for all four
  frozen intervals.
* Per-degree class sizes (all four intervals): `[23, 265, 2639, 24913, 229703, 2093785]`
  → total **2,351,328** (matches ticket expectation exactly).
* All four intervals: `excluded = True`, `screen_survivors_total = 0` (certified
  mean-value screen retains zero candidates at every degree).
* Re-ran `scripts/degree6_implementation_agreement.py`: `heights_compared = [3, 4]`,
  `cells_compared = 48`, `cells_in_agreement = 48`, `implementations_agree = True`;
  closest-member residual gap within mean-value allowance for every interval.

## Honesty note

The height threading was already present in `main` (`coefficient_height_max` already
reports the `height` parameter; height-4 artifacts already committed). This caller did
**not** edit the replication script. The independent-author property of the code is
intact; this assignment is a verification, not an authorship change. See the notes file
for the explicit disclosure the ticket requested.
