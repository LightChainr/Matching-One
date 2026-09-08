# probe664 — what #638's matching-side state count does to #636: a read-only cost-model note (2026-09-08)

**Ticket:** #664 (parent #650).  **Claim level:** C2 — read-only reanalysis of a committed
PR artifact; no new enumeration was run on this machine and none was needed.

**Source artifact:** `results/probe638-matching-boundary-states/latest.json`
(schema `matching-one.probe638.boundary-states.v2`) as delivered by **PR #648**
(`research/probe638-matching-states-20260907`, unmerged), script
`scripts/probe638_matching_boundary_states.py`, notes
`notes/probe638-boundary-states-20260907.md`.  Every number quoted below was
read directly from that JSON on PR #648's head, not from a prose summary.

## 1. The actual counts, quoted from the artifact

Per-width fields `nn_class_count` (site frontier, NN), `nnn_class_count`
(site frontier, NN+NNN = the matching lattice), `nn_catalan` (control), and
the bond anchor fields `bond_planar_class_count` / `bond_cylinder_class_count`
(covered below; identical values):

| w | site NN | site NNN (matching) | Catalan(w) | bond planar | bond cylinder |
|---|---|---|---|---|---|
| 1 | 1  | 1   | 1    | 1    | 1    |
| 2 | 2  | 2   | 2    | 2    | 2    |
| 3 | 5  | 5   | 5    | 5    | 5    |
| 4 | 13 | 12  | 14   | 14   | 14   |
| 5 | 34 | 27  | 42   | 42   | 42   |
| 6 | 88 | 61  | 132  | 132  | 132  |
| 7 | 232 | 142 | 429  | 429  | 429  |
| 8 | 620 | 340 | 1430 | 1430 | 1430 |

Artifact scalars, quoted verbatim:

* `control_ok = true`
* `smallest_crossing_width = None` (explicit negative result — the artifact
  distinguishes this from a missing run)
* `crossing_witness = None`
* `conclusion = "site helical frontier: NN classes 1,2,5,13,34,88,232,620
  (all noncrossing, strict subset of Catalan; missing the nested-pocket
  classes); NNN classes 1,2,5,12,27,61,142,340 = the row-cut class sequence
  (all noncrossing, zero crossing classes through w=8, NNN strictly inside
  NN, and NNN-all-phase == NN-row-cut elementwise per width).  Bond anchor
  with sealed semantics: planar == cylinder == Catalan(w) elementwise for
  w=1..8 (1,2,5,14,42,132,429,1430)."`

Per-width structural fields (`nnn_crossing_class_count = 0`,
`nnn_all_noncrossing = true`, `nnn_subset_nn_strict = true`,
`nnn_equals_nn_phase0_elementwise = true`) hold at every width w = 1..8.

## 2. Tripwire: the Catalan control

#638's tripwire — "NN control = Catalan(w), or the enumerator is unused" —
**passes**.  The verdict in PR #648 is a split verdict, and the split matters:

* **Bond side: Catalan exactly.**  `bond_planar_equals_noncrossing_states_elementwise
  = true` at every width w = 1..8; planar and cylinder counts are identical
  and equal Catalan(w) = 1, 2, 5, 14, 42, 132, 429, 1430 elementwise.  This
  requires the sealed-port semantics (an absent vertical bond seals its old
  block shut); the wrong semantics yields phantom closures 15, 52, 203, …,
  and that dead end is pinned by a test in PR #648.
* **Site side: strictly below Catalan, provably.**  The site helical scan —
  the construction #636 literally specifies — reaches 13 of 14 classes at
  w=4 (missing the nested pocket `{0,3},{1,2}`), and the gap grows: 620 vs
  1430 at w=8.  PR #648 validated the transfer bit-for-bit against an
  independent direct classifier over every occupancy prefix (~4·10^7 prefix
  checks, w ≤ 5, NN and NNN), so the shortfall is structural, not an
  enumerator bug.

So the control as stated is satisfied by the bond transfer and deliberately
*not* satisfied by the site transfer — with the site shortfall independently
verified rather than being evidence of a broken enumerator.

## 3. One sentence on the #636 width ceiling

The matching lattice (NNN) does **not** raise or lower #636's width ceiling —
it leaves it, because the NNN reachable class sequence (1, 2, 5, 12, 27, 61,
142, 340 through w=8) is strictly inside the NN sequence and carries zero
crossing classes through w=8, so the binding constraint on w is whichever
representation #636 adopts: Catalan(w) ≈ 4^w if it uses the bond row
transfer, or the cheaper site sequence (~3.2^w, 620 at w=8) if it uses the
site helical scan it currently specifies.

## 4. Not funding #636

This note does **not** fund #636.  Nothing here builds, weights, or
spectrally analyzes a transfer matrix; #636 owns the weights, the operator
sparsity, and the acceptance gate, and #635 still gates the whole
enterprise.  The only admissible reading of this note is a **recommendation
to the #636 design discussion**, recorded here and nowhere else:

* If a bond representation is chosen, the Catalan(w) cost model stands as
  written in #636.
* If the site representation is chosen, the state-count driver is
  1, 2, 5, 13, 34, 88, 232, 620 (NN) / 1, 2, 5, 12, 27, 61, 142, 340 (NNN) —
  roughly 2.3× fewer states at w=8 — but a state count is not a feasibility
  proof: row-operator sparsity is unmeasured here.
* #636 remains unreleased; no STATUS file was touched; no ticket was closed.
