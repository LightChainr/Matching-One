# probe #625 re-run on main — M(1/2) vs shape Z (2026-09-08)

**Variant:** this is a fresh verification on a clean Mac clone of `main`,
NOT a re-do of the merged #629 probe. Relationship to prior work:

- PR #629 (`probe/625-mhalf-vs-shape` @ `ac0bb20a`) was **merged into main**
  on 2026-09-07; the D1–D3 probe programs and results live on main.
- The original D4 (`scripts/probe/bond_L3_Z.py`) enumerated 2^18 bond
  configs. This re-run is **Mac-budgeted to NOT enumerate 2^18**: D4 here
  imports the repaired bond census from PR #653
  (`repair/632-bond-lab-and-wrapping-l3l4` @ `7109b174`,
  `results/probe-invariant-shape/census-exact.json`; main's copy is the
  pre-#632-repair one with `dual_fail=118133`) and derives the D4 verdict
  from the 3-count census arithmetically. Nothing is enumerated beyond 2^9.
- Game B site L=3 (2^9) and the L=4 spot checks were run locally; all
  numbers below were **reproduced on this machine** and match the values
  committed on main (bit-identical floats in the JSONs).

## D5 verdict

```text
W4 at L=3,4:  KILLED
evidence:     Game B L=3, exact: on the m-pinned line F_lambda(p)=1/2
              (m = E_lambda[r_b]-1 = 0 EXACTLY for every lambda = e^beta),
              Z(0.5; anchors 0.1/0.9) sweeps
              0.73817 (lambda=1/16) -> 0.54242 (lambda=1) -> 0.32192
              (lambda=16); interior-u spread over the 11-point grid =
              0.416 (a1/a9), 0.271 (a2/a8).  L=4 spot (lambda = 1, 5/2,
              2/5): spread 0.178.  Both ~4e7 x the 1e-8 tolerance.
Z_3 vs Z_4:   ||Z_3-Z_4||_inf on {0.1..0.9} = 0.0161 (anchors 0.2/0.8);
              at p_L^H |delta| = 0.0066, at p=1/2 |delta| = 0.0837
              (the same-p row separates because M(1/2) != 0:
              -21/64 at L=3, -13757/32768 at L=4).
bond vs site: INCOMPARABLE — bond X-law at p=1/2 is 3-atom
              (0.28786, 0.42429, 0.28786); inverse-CDF is a staircase;
              Z_bond undefined by the probe's own rule; jumps recorded.
              Census imported from PR #653 (dual_fail=0, M(1/2)=0 exact,
              P(r=0)=P(r=2) exact — duality-oddness checked from counts,
              not assumed).
```

## What was run here (all exact rational arithmetic)

- **D1 precondition** (imported from main's `mhalf_common.stream_joint`):
  Alexander `r_b + r_w = 2` asserted on the stream at L=3, L=4 — holds
  (`dual_fail=0`). Precondition, not a headline.
- **D2** `scripts/probe/mhalf_Z_physical_tables.py` (main's, unmodified):
  reproduced `p_L^H = 0.586511455113` (L=3) / `0.590672112331` (L=4),
  `M(1/2) = -21/64` / `-13757/32768`, `F(1/2) = 43/128` / `19011/65536`,
  `||Z_3 - Z_4||_inf = 0.0161`, self-symmetry rows
  `Q(u)+Q(1-u)-1` (max 0.1712 / 0.1803 vs `2p_L^H - 1` = 0.1730 / 0.1813).
- **D3 Game A** (L=3 fake law, `m` pinned at -21/64 exactly): shape still
  swings, spread 0.667 — the easy game fails too.
- **D3 Game B** (the game that matters; site L=3 = 2^9 configs, Mac-OK,
  full 11-point lambda grid; L=4 spot at lambda = 1, 5/2, 2/5):
  monotone guard asserted before every inversion; `m` pinned at 0 on the
  whole lambda family; `Z` not pinned. **W4 dies for the physical
  interpolation.**
- **D4 (imported, no 2^18 here)** `scripts/probe/bond_L3_Z_imported.py`:
  validates the imported census (configs = 2^18, counts sum, dual_fail=0),
  checks `M_bond(1/2) = P(r=2) - P(r=0) = 0` exactly and `P(r=0) = P(r=2)`
  exactly from the counts, reports the 3-atom jumps and
  `P(r >= 1)(p=1/2) = 0.712142944` (census-derivable). The median-crossing
  `p_wrap = 0.419649` needs the full enumerated `F_bond(p)` and is cited
  from main's committed `bond_L3.json`, not recomputed.
- **Tests:** `python3 tests/test_probe_mhalf_vs_shape.py` — 5 fast OK,
  1 slow skip (Python 3.13.7).

## Deliverables in this PR (all new; 0 modified, 0 deleted)

- `notes/probe-Mhalf-vs-shape-main-rerun-20260908.md` (this file)
- `scripts/probe/bond_L3_Z_imported.py` (D4 from imported census)
- `results/probe-Mhalf-vs-shape/bond_census_imported_653.json`
  (imported #653 census, provenance-recorded)
- `results/probe-Mhalf-vs-shape/bond_L3.json` (regenerated in imported mode)
- `results/probe-Mhalf-vs-shape/physicalZ.json` (reproduced, bit-identical)
- `results/probe-Mhalf-vs-shape/gameAB.json` (reproduced, bit-identical)
- `results/probe-Mhalf-vs-shape/latest.json` (v2 index for this variant)

## Standing

Does not enter `docs/STATUS.md`. Closes nothing (#606, #619, #622, #608,
#625 all stay open). Does not merge. Census numbers are inputs, not
deliverables. No ninth mechanism, no exponent, no `L^-theta` fit.
