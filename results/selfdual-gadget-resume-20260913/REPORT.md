# #13 Self-dual gadget resume — REPORT

**Issue:** #13 (self-dual gadget & critical-manifold search — 2026-08-31 resume)
**Branch:** `theory/p13-selfdual-gadget-resume-20260913`
**Method script:** `scripts/selfdual_gadget_resume.py`
**Outcome:** **(b) theorem + explicit counterexample** — the finite D4-orbit serial class cannot provide a law-preserving monotone comparison between two specific period diagrams.

## Exact facts (recomputed)

- Finite D4-orbit serial class: **7 orbits / 15 states** (reconstructed via the canonical terminal-partition serial-category machinery; the exact #438 192-state/41-orbit W5 table is not in the tree → buy-back).
- **26/49** orbit pairs are multi-valued under serial composition (class is not a well-defined category).
- **166/343** orbit triples are non-associative under uniform orbit averaging.

## Explicit counterexample (two specific period diagrams)

Orbit pair `(0, 2)`:
- A = `[0,0,0,0]` (orbit 0)
- B representatives in orbit 2: `[0,0,1,1]` and `[0,1,1,0]`
- `serial_compose(A, [0,0,1,1]) = [0,0,1,1]` (orbit 2)
- `serial_compose(A, [0,1,1,0]) = [0,0,0,0]` (orbit 0)

The composed connectivity law of A with the class-2 diagram is not single-valued, so no monotone (Strassen) coupling — which needs well-defined laws under a defined gluing — can be assigned. **The finite class cannot provide the comparison.**

## Governance

No new adjacent-algebra census, no generic-certificate task (per 2026-08-31 handoff). `docs/STATUS.md` and other protected files untouched; only new files added.

Full Matching-One repository CI has not been run for this commit.
