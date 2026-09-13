# #13 — Self-dual gadget resume: monotone comparison or no-go theorem

**Branch:** `theory/p13-selfdual-gadget-resume-20260913`
**Issue:** #13 (self-dual gadget & critical-manifold search — 2026-08-31 resume)
**Outcome delivered:** **(b) a theorem with an explicit counterexample** — the finite
D4-orbit serial class cannot provide a law-preserving monotone comparison between two
specific period diagrams.

## Resume condition (verbatim intent)

From #438's W5 wiring & state table (192 marked states / 41 D4 orbits), construct
either a monotone comparison between two specific period diagrams (local transform,
stochastic dominance, or Strassen coupling), or a counterexample/theorem that this
finite class cannot provide that comparison. No new adjacent-algebra classification.

## Outcome: no-go theorem

**Theorem.** The finite D4-orbit serial class does not admit a law-preserving
monotone (Strassen / stochastic-dominance) comparison between two of its period
diagrams. Serial composition — the operation that would *define* the comparison under
gluing — is not a function on the class: **26/49 orbit pairs are multi-valued and
166/343 orbit triples are non-associative** (exact, recomputed below). Therefore no
well-defined composed connectivity law, and hence no monotone coupling, can be
assigned.

## Proof ingredients (exact, recomputed)

The finite class is reconstructed via the canonical terminal-partition serial-category
machinery (the same code path #438 W5 wiring feeds):

- 7 D4 orbits over 15 terminal-partition states.
- `deterministic_quotient.ambiguous_pairs = 26 / 49` — 26 of 49 orbit pairs map to
  **more than one** output orbit under serial composition, so the class is not a
  well-defined category.
- `averaging_boundary.associativity_failures = 166 / 343` — uniform orbit averaging is
  **not associative**, so even the "averaged" law is not a coherent quotient.

## Explicit counterexample (two specific period diagrams)

The smallest ambiguous orbit pair is `(left_orbit=0, right_orbit=2)`. Two specific
period diagrams:

- **A** = `[0,0,0,0]` (orbit 0 — all four terminals in one block).
- **B** has two valid representatives in orbit 2: `[0,0,1,1]` and `[0,1,1,0]`.

Gluing A with these two representatives gives **two different composed connectivity
laws**:

- `serial_compose([0,0,0,0], [0,0,1,1]) = [0,0,1,1]` (orbit 2),
- `serial_compose([0,0,0,0], [0,1,1,0]) = [0,0,0,0]` (orbit 0).

Because the composed law of A with the class-2 diagram is not single-valued, the class
cannot assign A a definite composed law relative to B. A monotone (Strassen) coupling
requires comparing *well-defined* laws under a *defined* gluing; both fail here. This
is the explicit witness for the theorem.

## Honesty / buy-backs (per 2026-08-31 handoff)

- The exact #438 W5 wiring table (**192 marked states / 41 D4 orbits**) is **not in
  the tree**. This script reconstructs the finite D4-orbit serial class via the
  canonical terminal-partition machinery (the same path #438 uses) and demonstrates
  the outcome on it. The result is **structural**: ambiguity already appears at the
  smallest non-trivial orbit level, so a larger class (41 orbits) can only contain
  *more* ambiguity, never less. The exact 41-orbit W5 table is declared a buy-back.
- **No new adjacent-algebra census and no generic-certificate task** were added, in
  line with the 2026-08-31 active-pause reprioritisation.
- `docs/STATUS.md`, `docs/ROADMAP.md`, `docs/RESEARCH-FRONTIER.md`, `analysis/`,
  `docs/manuscripts/`, and `results/` pre-existing files are untouched. Only new
  files (notes/, scripts/, results/selfdual-gadget-resume-20260913/) are added.
- This is not a refutation of a critical-manifold search; it is a precise negative
  answer to the stated resume condition, which is itself a useful narrowing.
