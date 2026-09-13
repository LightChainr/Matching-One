# p593 Phase B — P398 widths 9-10: memory degree, rank notions, lumping collapse (2026-09-13)

Issue: #593 (Phase B, started after Phase A was delivered in PR #753).
Everything at `eta = 0` unless stated. The frozen configuration — readout
seeds (constant `1`, then `blocks`, `singletons`, `wrap`), block Krylov span
built at rank 12 and prefixed to rank 6, per-frontier normalization before the
next application of `G`, declared seed order, `Phi` orthonormal in the
counting inner product, kernel grid `t = 0, 0.25, ..., 4.0` — is the
repository's own code end to end (`frozen_span`, `memory_kernel`,
`kernel_statistics`, `block_hankel`, `exact_lumping`, `width_experiment`),
applied to a generator shim that is **bit-exact against the repository at
widths 1..8** (states, targets and rate rows equal; transition-table control
below). Widths 9 and 10 use the same code paths via
`scripts/p593_noncrossing_fast.py`.

Two declared departures, both asked for by the ticket:

1. Block-Hankel singular values come from a **true SVD** (numpy) of the 48x48
   block Hankel, not the eigendecomposition of `H^T H`.
2. `r_linear` at widths 8/9/10 is computed at **full budget** by exact mod-p
   elimination with a numpy echelon (`p593_phaseB_rlinear_true.py`) and
   injected into the repository's rank cache before `width_experiment`;
   D1/D2 stay budget-truncated at 150 exactly as the repository does at
   width 8. Width-10 `r_transport` uses the `eta` ladder `(0, -0.25, 0.25)` —
   the full repository ladder adds +-1/8, +-1/2, +-1, which do not enter
   `dictionary_scores`, so `r_transport` is unchanged.

## Controls (all pass; `raw/controls.json`)

| control | value |
|---|---|
| 1. transition-table SHA-256, widths 1..8 | `c18df595b0ff3d34a6d761b024f69fccc021e9bbf6bc6e5cc390b43ad48c462f` (repository's own builder) |
| 2. row sums / rate positivity / strong connectivity | 0.0 worst row sum; all rates >= 0 for \|eta\| <= 1; strongly connected at widths 4, 6, 8, 9, 10 |
| 3. `exp(tG) 1 = 1` | worst 1.6e-15 (widths 4, 8, 10) |
| 4. uniformization vs dense Taylor expm, width 5 | 1.9e-15 |
| 5. full-rank projection K == 0 | 0.0 exactly (width 4) |
| 5. generator-closed Krylov span K == 0 | 5.75e-30 relative to reference norm 19.05 (width 4; ticket value 6.1e-30) |
| 6. resolvent Schur identity, width 5 | 3.8e-16 (ticket value 7.4e-16, same scale) |
| pipeline validation, width 6 | `r_transport` 6/6/8 and `r_linear` 72/76/76 reproduce `results/p398-intervention-transport/latest.json` digit-for-digit |
| memory widths 4..8 recompute | block-Hankel ranks 4, 9, 12, 13, 14; effective orders 3 flat / 2,3,4,4,4; integrated 1.681, decay 0.136, tail 0.0006 at width 8 — all match the ticket |

## Deliverable 1 — the memory degree

Widths 4..8 recomputed with the same code path (agreement as above); 9 and 10 new.

| width | states | rank C | block-Hankel numerical rank (tol 1e-6) | eff. order 99% | eff. order 99.9% | int. \|K\| | decay time | tail mass > t=2 |
|---|---|---|---|---|---|---|---|---|
| 4 | 14 | 2 | 4 | 2 | 2 | 0.876 | 0.0626 | 1.7e-6 |
| 5 | 42 | 3 | 9 | 3 | 3 | 1.016 | 0.0891 | 3.7e-5 |
| 6 | 132 | 3 | 12 | 3 | 4 | 1.207 | 0.1048 | 1.4e-4 |
| 7 | 429 | 3 | 13 | 3 | 4 | 1.456 | 0.1192 | 3.2e-4 |
| 8 | 1430 | 3 | 14 | 3 | 4 | 1.681 | 0.1363 | 5.7e-4 |
| 9 | 4862 | 3 | **14** | 3 | 4 | 1.873 | 0.1541 | 9.3e-4 |
| 10 | 16796 | 3 | **14** | 3 | 4 | 2.036 | 0.1717 | 1.5e-3 |

Leading normalized singular values of the block Hankel are in
`raw/width{9,10}.json` (top-16 each).

**Answer to the ticket's core question: the numerical order saturates.**
The block-Hankel numerical rank goes 4, 9, 12, 13, 14, **14, 14** across
widths 4..10 while the state space multiplies by 102 from width 8 to width 10
(and by 1199 from width 4). The 99%-energy order is flat at 3 from width 5 on
and the 99.9%-energy order is flat at 4 from width 6 on. Integrated kernel
weight and decay time creep up slowly (1.681 -> 2.036 and 0.136 -> 0.172 from
width 8 to 10) but the pole count does not move. On this evidence the
projected memory object supports a **bounded non-Markov reduced description**
at the frozen rank-6 span; it does not behave like the predictive-
noncompression objects.

## Deliverable 2 — the three rank notions continued

| width | states | `r_linear(D0)` true | `r_positive` (joint `G0`+`H`, D0 colouring) | `r_transport` |
|---|---|---|---|---|
| 4 | 14 | 10 | 10 | 4 |
| 5 | 42 | 26 | 26 | 4 |
| 6 | 132 | 72 | 76 | 6 |
| 7 | 429 | 218 | 232 | 6 |
| 8 | 1430 | **689** (was >= 150) | 750 | 6 |
| 9 | 4862 | **2275** (depth 2152) | **2494** | **6** |
| 10 | 16796 | **7718** (depth 7478) | **8524** | **6** |

(The width 4..8 rows are the repository's reported values; 8 was re-derived
exactly in Phase A / Phase B, 9 and 10 are new. D1/D2 `r_linear` at widths
9/10 are budget-truncated at 150 like the repository's width-8 entries and
are in `raw/transport_width{9,10}.json`.)

- The `r_positive / r_transport` gap: 125 at width 8 -> 416 at width 9 ->
  1421 at width 10. State count x3.4 from width 8 to 9 and x3.5 from 9 to 10
  moved the gap x3.3 and x3.4: **the gap widens linearly in the state count**,
  not slower.
- `r_linear`/states drifts slowly down (0.71, 0.62, 0.55, 0.51, 0.48, 0.47,
  0.46 across widths 4..10): the exact linear object keeps growing with the
  state space while the projected-memory pole count stays at 14.
- `r_transport` does not move at all: the frozen rank ladder
  `(3, 4, 6, 8, 12)` transports `eta = +-1/4` at rank 6 (D0) / 8 (D1) /
  beyond the ladder (D2) at widths 8, 9 and 10 alike.

## Deliverable 3 — the out-of-pencil lumping collapse

Joint lumping against the baseline generator plus the `single_point_join`
tilt (outside the `span{J, D}` pencil):

| width | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|---|---|---|---|---|---|---|---|
| coarsest out-of-pencil lumping (blocks) | 14 | 42 | 132 | 429 | 1430 | 4862 | 16796 |
| collapses to identity | yes | yes | yes | yes | yes | **yes** | **yes** |

**Confirmed at widths 9 and 10**: the coarsest lumping valid for the
out-of-pencil family is the identity partition at every width tested. Not a
small-width accident.

## Boundary

This is a calibration model. Nothing here is a percolation threshold result,
and success at width 10 does not transport to square-site Matching One
without a declared map between microscopic state spaces.

Full Matching-One repository CI has not been run for this commit.
