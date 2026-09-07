# probe638 — boundary connectivity classes: the site frontier, the bond anchor, and what Catalan(w) actually guards (2026-09-07)

Ticket: #638 ("Enumerate the matching-lattice boundary state space — the
single biggest unknown in #636's cost model").  Branch
`research/probe638-matching-states-20260907`, based on
`claude/matching-one-workspace-pwr5pv` @ `3909725`.
Artifact: `results/probe638-matching-boundary-states/latest.json`
(schema `matching-one.probe638.boundary-states.v2`).  Script:
`scripts/probe638_matching_boundary_states.py`.  Tests:
`tests/test_probe638_matching_boundary_states.py` (9 tests, all green,
every docstring names its guarded wrong number).

## What was asked and what it got

#638 asks: which boundary connectivity classes can a width-w transfer
actually reach, under (a) NN adjacency — with the ticket's control that the
count must equal Catalan(w) = 1, 2, 5, 14, 42, 132, 429, 1430 exactly — and
(b) NN+NNN matching-lattice adjacency — with the count, the smallest width
carrying a crossing class, and an explicit witness.  Both parts are now
answered, with one honest split verdict on the control.

## Headline numbers (exact, no floats anywhere in the artifact)

| w | site NN classes | site NNN classes | row-cut classes | Catalan(w) | bond planar | bond cylinder |
|---|---|---|---|---|---|---|
| 1 | 1 | 1 | 1 | 1 | 1 | 1 |
| 2 | 2 | 2 | 2 | 2 | 2 | 2 |
| 3 | 5 | 5 | 5 | 5 | 5 | 5 |
| 4 | 13 | 12 | 12 | 14 | 14 | 14 |
| 5 | 34 | 27 | 27 | 42 | 42 | 42 |
| 6 | 88 | 61 | 61 | 132 | 132 | 132 |
| 7 | 232 | 142 | 142 | 429 | 429 | 429 |
| 8 | 620 | 340 | 340 | 1430 | 1430 | 1430 |

Structural facts (all machine-checked in the artifact):

* Site NNN closure ⊂ site NN closure, every width, strictly from w=4.
* Site NNN all-phase class set == site NN row-cut (phase-0) class set,
  elementwise, every width w=2..8.
* **Zero crossing classes through w=8**, site side, NN or NNN.
* Bond anchor: planar and cylinder closures are elementwise EQUAL to the
  committed `noncrossing_states(w)` for every width 1..8 — Catalan exactly,
  including on the cylinder (the wrap edge sews the two end ports into one
  block but adds no new classes at w<=8).

## The control verdict, stated honestly

The ticket says the NN control "must return Catalan(w) exactly … and if it
does not, the enumerator is broken and nothing else in the report means
anything."  What the probe actually found is sharper:

1. **The control holds — for the bond boundary-state space.**  The classic
   dangling-bond row transfer (w cut ports, one per column) reaches exactly
   the noncrossing partitions, PROVIDED two semantics are pinned down:
   (i) a port whose vertical bond is ABSENT becomes its own isolated
   singleton (the old block it used to carry is sealed shut — it can never
   reconnect from the new side, which is what makes the state space planar);
   (ii) a port with a present vertical bond carries the fresh site's
   cluster.  Get (i) wrong — let an absent port keep its old block — and the
   closure jumps to 15, 52, 203, … at w=4,5,6 with crossing classes such as
   `{0,2},{1,3}` reachable even on the planar strip: phantom growth by
   long-range reconnection that no actual lattice path performs.  Get the
   geometry wrong — no wrap edge — you still get exactly Catalan (planar
   strip); with the wrap edge, still Catalan at w<=8.  This is the correct
   reading of the control, and it passes bit-for-bit, elementwise, w=1..8.

2. **The control fails — for the site frontier, and provably so.**  The
   site helical scan (what #636 literally specifies: "adding sites one at a
   time") reaches only 13 of the 14 noncrossing classes at w=4; the missing
   class is the nested pocket `{0,3},{1,2}` (RGS (0,1,1,0)).  This is not an
   enumerator bug: the transfer was validated against an independent
   direct classifier over ALL occupancy prefixes (2 × 10^7 prefix checks,
   bit-for-bit, w=1..5, NN and NNN).  The reason is structural: on the
   scan frontier, two disjoint multi-site blocks must have all their
   cross-boundary adjacencies landing on the scan's two free seams; the
   nested pocket needs its (0,1) and (2,3) pairs adjacent in the interior
   AND wrapped onto each other, which no site configuration produces —
   whereas {0,1},{2,3} (RGS (0,0,1,1)) IS reachable (one witness: row 0
   sites {0,3} wrap-joined, row 1 sites {1,2}).  Site frontier classes =
   1,2,5,13,34,88,232,620: noncrossing always, but a shrinking fraction of
   Catalan (0.857 → 0.238 by w=8).

3. **Consequence for #636's cost model.**  The Catalan(w) state count is
   the right bound for a BOND-percolation transfer on the width-w cylinder.
   For the SITE transfer that #636 describes, the correct state count
   driver is the site-frontier sequence above (620 states at w=8, ~0.43 ×
   Catalan(8)), and the NNN matching-lattice transfer carries FEWER states
   than NN (340 vs 620 at w=8) — all of them noncrossing.  The matching
   lattice does NOT blow up the state space at the widths reachable here;
   if a crossing class exists at all, it appears at w > 8.

## Independence (GOVERNANCE 2)

* NN/NNN site classes: transfer closure vs. independent direct classifier
  (explicit graph search, no union-find, no shadows) over every occupancy
  prefix of length up to w·max(3,18/w) for w=1..5 — 4718592 prefixes each
  at w<=3, 1048576 at w=4, 491520 at w=5; all agree bit-for-bit.
* Bond anchor: counts compared against `noncrossing_states(w)` from the
  committed codec, which itself cross-checks two internal enumerators
  (RGS generation and block-insertion generation).
* Two transfer bugs were caught and fixed by these checks in this probe:
  the slot-merge bug (merged partner slots kept stale labels) and the
  shadow-merge bug (merged shadows not relabelled).  Both are pinned by
  tests with the exact wrong numbers they produced (13-not-14 claims;
  (0,1,2,3,4)-instead-of-(0,1,2,3,1) at w=5 NNN).

## Crossing-class answer for the NNN side

Smallest width with a crossing boundary class: **none through w=8**.  The
artifact therefore records `smallest_crossing_width: null` and no witness —
a negative result, not a missing run.  If the ticket wants the first
crossing width, the next step is w=9..12 with the same closure machinery
(state counts grow roughly ×1.5 per width; w=9 NN site should be ~1700
states, minutes of work; the bond anchor is the expensive one, ~8-20 min
per variant per width).

## What I could not do / what would finish it

* w=9+ enumeration (both site closures) — trivial compute, deliberately
  left out to keep the first artifact at the ticket's requested widths.
* A proof-level statement of the site-frontier unreachability theorem
  (nested pockets unreachable) — the probe establishes it empirically at
  w=4..8 with an exhaustive prefix-level validation at w<=5 and a
  mechanistic explanation; a written proof would need a seam-counting
  argument.  Worth doing before #636 commits to a site-based state count.
* Whether #636's transfer should be defined on bonds (then Catalan(w) is
  exact and the cost model stands as written) or on sites (then the state
  count is the 1,2,5,13,34,88,232,620 sequence) is a design decision for
  the ticket owner; the probe supplies both enumerations so the choice is
  now data-driven.
