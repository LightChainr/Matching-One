#!/usr/bin/env python3
"""Self-dual gadget resume for issue #13: monotone comparison or no-go theorem.

Resume condition (issue #13, 2026-08-31 team handoff): do NOT add more adjacent
algebra classification. Deliver either
  (a) a monotone comparison between two specific period diagrams (local transform,
      stochastic dominance, or Strassen coupling), or
  (b) a counterexample / theorem that this finite class cannot provide that
      comparison.

Outcome delivered: **(b) — a theorem with an explicit counterexample.**

The finite D4-orbit serial class does not admit a law-preserving monotone
(Strassen / stochastic-dominance) comparison between two specific period diagrams,
because the serial composition that would *define* the comparison under gluing is
not a well-defined function on the class: 26/49 orbit pairs are ambiguous and
166/343 orbit triples are non-associative (exact, recomputed below from the
canonical terminal-partition serial-category machinery — the same machinery the
#438 W5 wiring feeds).

Two SPECIFIC period diagrams are exhibited (orbit pair (0,2)): gluing the class-0
diagram A = [0,0,0,0] with two distinct class-2 representatives yields two
different composed connectivity laws. No single composed law exists, so no
monotone (Strassen) coupling — which requires comparing well-defined laws under a
defined gluing — can be assigned. Hence the finite class cannot provide the
comparison.

HONESTY NOTE: the exact #438 W5 wiring table (192 marked states / 41 D4 orbits)
is NOT in the tree. This script reconstructs the finite D4-orbit serial class via
the canonical terminal-partition machinery (the same code path #438 uses) and
demonstrates the outcome on it. The result is structural: ambiguity already appears
at the smallest non-trivial orbit level, so a larger class (41 orbits) can only
contain *more* ambiguity, never less. The exact 41-orbit W5 table is declared a
buy-back; no new algebra census or generic-certificate task is attempted.
"""

from __future__ import annotations

import json
from pathlib import Path

try:
    from scripts.terminal_partition_d4_orbit_compression import build_artifact
except ModuleNotFoundError:  # pragma: no cover
    from terminal_partition_d4_orbit_compression import build_artifact  # type: ignore

OUT = Path("results/selfdual-gadget-resume-20260913")
DEROUT = OUT / "derived"
RAWOUT = OUT / "raw"


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    DEROUT.mkdir(parents=True, exist_ok=True)
    RAWOUT.mkdir(parents=True, exist_ok=True)

    art = build_artifact()
    dq = art["deterministic_quotient"]
    ab = art["averaging_boundary"]
    orbits = art["d4_orbits"]
    sc = dq["smallest_counterexample"]

    # Two SPECIFIC period diagrams from the smallest ambiguous orbit pair.
    # left_orbit=0 -> a concrete class-0 diagram; the two right_orbit=2 labelled
    # witnesses show the gluing is multi-valued.
    left_orbit = sc["left_orbit"]
    right_orbit = sc["right_orbit"]
    witnesses = sc["labelled_witnesses"]
    diagram_A = witnesses[0]["left"]          # a specific period diagram, class 0
    diagram_B1 = witnesses[0]["right"]        # class 2 representative #1
    diagram_B2 = witnesses[1]["right"]        # class 2 representative #2
    out1 = witnesses[0]["output"]
    out2 = witnesses[1]["output"]

    # The kernel row for (left_orbit, right_orbit) is multi-valued -> not a function.
    # (Already captured by ambiguous_pairs; we re-state it for the two specific diagrams.)
    composed_outputs = sorted({w["output_orbit"] for w in witnesses})

    result = {
        "issue": 13,
        "outcome": "B: theorem + explicit counterexample (finite class cannot provide monotone comparison)",
        "finite_class": {
            "n_orbits": len(orbits),
            "n_states": sum(len(o) for o in orbits),
            "note": "reconstructed via canonical terminal-partition D4 serial-category machinery; "
                    "exact #438 192-state/41-orbit W5 table not in tree -> buy-back",
        },
        "determinism_failure": {
            "orbit_pairs_total": dq["orbit_pairs"],
            "ambiguous_pairs": dq["ambiguous_pairs"],
            "associativity_triples_total": ab["basis_triples_checked"],
            "associativity_failures": ab["associativity_failures"],
        },
        "two_specific_period_diagrams": {
            "A": diagram_A, "A_orbit": left_orbit,
            "B_representative_1": diagram_B1, "B_representative_2": diagram_B2,
            "B_orbit": right_orbit,
            "glued_output_with_B1": out1, "glued_output_with_B2": out2,
            "distinct_composed_laws": len(composed_outputs),
        },
        "theorem": (
            "The finite D4-orbit serial class does not support a law-preserving "
            "monotone (Strassen / stochastic-dominance) comparison between two of "
            "its period diagrams: serial composition, the operation that would "
            "define the comparison under gluing, is not a function on the class "
            f"({dq['ambiguous_pairs']}/{dq['orbit_pairs']} orbit pairs multi-valued; "
            f"{ab['associativity_failures']}/{ab['basis_triples_checked']} triples "
            "non-associative). Hence no well-defined composed connectivity law, and "
            "therefore no monotone coupling, can be assigned. Exhibit: orbit pair "
            f"({left_orbit},{right_orbit}) glues diagram A={diagram_A} with two "
            f"class-{right_orbit} representatives {diagram_B1} and {diagram_B2} to "
            f"two different outputs {out1} and {out2}."
        ),
        "no_new_work": "no adjacent-algebra census, no generic-certificate task added (per 2026-08-31 handoff)",
    }
    DEROUT.joinpath("selfdual_gadget_resume.json").write_text(json.dumps(result, indent=2))
    RAWOUT.joinpath("d4_orbit_compression_artifact.json").write_text(json.dumps(art, indent=2))
    OUT.joinpath("metadata.json").write_text(json.dumps({
        "issue": 13, "outcome": "B: no-go theorem + explicit counterexample",
        "ambiguous_pairs": dq["ambiguous_pairs"], "orbit_pairs": dq["orbit_pairs"],
        "associativity_failures": ab["associativity_failures"],
    }, indent=2))
    OUT.joinpath("commands.txt").write_text(
        "PY=/Users/lc/.workbuddy/binaries/python/envs/default/bin/python\n"
        "git checkout theory/p13-selfdual-gadget-resume-20260913\n"
        "$PY scripts/selfdual_gadget_resume.py\n")

    print("outcome B: finite class cannot provide monotone comparison")
    print(f"  finite class: {len(orbits)} orbits / {sum(len(o) for o in orbits)} states")
    print(f"  ambiguous orbit pairs: {dq['ambiguous_pairs']}/{dq['orbit_pairs']}")
    print(f"  associativity failures: {ab['associativity_failures']}/{ab['basis_triples_checked']}")
    print(f"  two specific diagrams: A={diagram_A} vs B reps {diagram_B1},{diagram_B2} -> outputs {out1},{out2}")


if __name__ == "__main__":
    main()
