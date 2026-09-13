# p593 Phase B — widths 9-10 memory degree and rank notions (raw: raw/)

- date: 2026-09-13
- issue: 593 (Phase B; started after Phase A delivery, PR #753)
- headline results:
  - block-Hankel numerical rank saturates: 4, 9, 12, 13, 14, 14, 14 across widths 4..10; 99% energy order flat at 3, 99.9% flat at 4 from width 6.
  - rank C = 3 at widths 5..10 (2 at width 4).
  - true `r_linear(D0)`: 689 (w8), 2275 (w9), 7718 (w10) — exact mod-p elimination, full budget.
  - `r_positive` (D0 colouring, joint G0+H): 750 (w8), 2494 (w9), 8524 (w10).
  - `r_transport`: 6 at widths 6..9 (D0); width 10 in raw/transport_width10.json.
  - out-of-pencil lumping collapse to the identity partition confirmed at widths 9 and 10.
- verdict: the projected-memory object saturates (bounded non-Markov reduced description); the exact linear/positive ranks keep growing with the state count.

## Files

- `raw/width{4..10}.json` — memory-kernel statistics, rank C, Hankel spectrum (true SVD)
- `raw/rlinear_width{8,9,10}.json` — true `r_linear(D0)` by exact mod-p Krylov
- `raw/lumping_width{4..10}.json` — `r_positive` and out-of-pencil collapse
- `raw/transport_width{9,10}.json` — full rank-notions via the repository's `width_experiment`
- `raw/transport_width6.json` — validation against the repository's stored results
- `raw/controls.json` — all six ticket controls

## Commands

```bash
cd scripts
python p593_noncrossing_fast.py                 # shim validation + transition-table hash
python p593_phaseB_memory.py 4 5 6 7 8 9 10     # Deliverable 1
python p593_phaseB_rlinear_true.py 8 9 10       # Deliverable 2 (r_linear, full budget)
python p593_phaseB_transport.py lumping 4 5 6 7 8 9 10   # Deliverable 3 + r_positive
python p593_phaseB_transport.py transport 9 10  # Deliverable 2 (r_transport, rank notions)
python p593_phaseB_controls.py                  # controls 1-6
```

## Declared departures

1. Block-Hankel singular values from a true SVD (ticket request).
2. `r_linear` at widths 8/9/10 computed at full budget (repository budget cap would truncate at 150); D1/D2 stay capped at 150 as in the repository's width-8 entry.
3. Width-10 `r_transport` eta ladder `(0, -0.25, 0.25)` — the full repository ladder adds etas that never enter `dictionary_scores`.

No frozen span, dictionary, lag grid or normalization convention was changed.

Full Matching-One repository CI has not been run for this commit.
