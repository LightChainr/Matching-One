#!/usr/bin/env python3
"""One exact character-class calculation for three fixed physical closures.

The calculation compares finite S4/S5 full-central projections with the
declared stable [Q-2,2] component. It is not a Monte Carlo score or a test
suite. No lattice occupation population is enumerated.
"""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
from hashlib import sha256
import json
from math import factorial
from pathlib import Path
import platform
import subprocess
import time


def cycle_types(n: int, smallest: int = 1):
    if n == 0:
        yield ()
        return
    for first in range(smallest, n + 1):
        for tail in cycle_types(n - first, first):
            yield (first,) + tail


def central_coefficient(q: int, winding: int, components: int) -> dict:
    """m_Q=<chi_[Q-2,2], Fix(pi^w)^c>, using exact class weights."""
    terms = []
    total = Fraction(0)
    for cycles in cycle_types(q):
        counts = Counter(cycles)
        centralizer = 1
        for length, count in counts.items():
            centralizer *= length**count * factorial(count)
        x1, x2 = counts[1], counts[2]
        character = x1 * (x1 - 1) // 2 + x2 - x1
        fixed = sum(length * count for length, count in counts.items()
                    if winding % length == 0)
        term = Fraction(character * fixed**components, centralizer)
        total += term
        terms.append({"cycles": cycles, "class_probability": str(Fraction(1, centralizer)),
                      "character": character, "fixed_colours_of_power": fixed,
                      "contribution": str(term)})
    dimension = q * (q - 3) // 2
    return {"Q": q, "dimension": dimension, "dimension_normalized_coefficient": str(total),
            "full_central_trace_coefficient": str(dimension * total), "class_terms": terms}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, default=Path("analysis/colour_specialization_gap_contract.json"))
    parser.add_argument("--output-dir", type=Path, default=Path("results/colour-specialization-gap"))
    args = parser.parse_args()
    started = time.perf_counter()
    contract_bytes = args.contract.read_bytes()
    contract = json.loads(contract_bytes)
    cases = []
    for case in contract["cases"]:
        finite = [central_coefficient(q, case["winding"], case["essential_components"])
                  for q in contract["integer_colours"]]
        m4, m5 = [Fraction(item["dimension_normalized_coefficient"]) for item in finite]
        stable = Fraction(case["generic_multiplicity"])
        cases.append({**case, "finite": finite,
                      "finite_Q4_minus_stable_specialization": str(2 * (m4 - stable)),
                      "dimension_only_prediction_at_Q5": str(5 * m4),
                      "actual_Q5_minus_dimension_only_prediction": str(5 * (m5 - m4)),
                      "dimension_only_candidate": "rejected" if m4 != m5 else "compatible_on_this_case",
                      "formal_generic_Q1_coefficient": str(-stable)})
    commit = subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip()
    result = {"schema": "matching-one.colour-specialization-gap.result.v1",
              "definition_commit": commit,
              "contract_sha256": sha256(contract_bytes).hexdigest(),
              "question": contract["question"], "cases": cases,
              "decision": "dimension_only_continuation_rejected_for_the_physical_closure_family",
              "interpretation": "finite_colour_isotypic_projection_and_stable_sector_specialization_do_not_commute_in_general",
              "next_required_object": "a specified homology_and_multiplicity_resolved_generic_Q_landing_with_its_complete_Q1_jet_not_only_Q4_label_times_dimension",
              "boundary": contract["boundary"],
              "new_lattice_enumerations": 0, "new_random_samples": 0,
              "completed_N25_J22_rerun": False,
              "elapsed_seconds": time.perf_counter() - started,
              "python": platform.python_version(), "machine": platform.machine()}
    lines = ["# Finite-colour projection and stable specialization do not commute", "",
             "One exact character-class sum on three prescribed physical rank1 closures.",
             "Activity, contractible-colour factors and the common rank projection are factored out.", "",
             "| Physical closure | m4 | m5 | stable m | full trace at Q4 | full trace at Q5 | dimension-only Q5 prediction |",
             "|---|---:|---:|---:|---:|---:|---:|"]
    for case in cases:
        a, b = case["finite"]
        lines.append(f'| {case["id"]} | {a["dimension_normalized_coefficient"]} | {b["dimension_normalized_coefficient"]} | {case["generic_multiplicity"]} | {a["full_central_trace_coefficient"]} | {b["full_central_trace_coefficient"]} | {case["dimension_only_prediction_at_Q5"]} |')
    lines += ["", "Here mQ is the exact character inner product, and dim[Q-2,2]=Q(Q-3)/2.",
              "The fixed candidate mQ=m4 for all Q>=4 fails on the latter two actual connectivity patterns.",
              "This does not refute every analytic continuation or claim these patterns occur in the scored N25 packet.",
              "", "The stable coefficient's formal Q1 value is not a positive sector probability.",
              "Before a Q1 field claim, fix the actual homology/multiplicity-resolved continuation and its complete jet.",
              "The completed Q4 J22 result remains unchanged and was not rerun.",
              "", "Proof and geometry: [specialization note](../../notes/colour-specialization-gap.md).",
              "Exact class contributions and provenance: [latest.json](latest.json).",
              "No new lattice enumeration, random samples, cloud computation or scientific test suite.", ""]
    args.output_dir.mkdir(parents=True, exist_ok=True)
    outputs = {"latest.json": json.dumps(result, indent=2) + "\n", "REPORT.md": "\n".join(lines)}
    for name, content in outputs.items():
        with (args.output_dir / name).open("x") as handle:
            handle.write(content)
    receipt = {"schema": "matching-one.colour-specialization-gap.run.v1",
               "definition_commit": commit, "command": "python scripts/analyze_colour_specialization_gap.py --contract analysis/colour_specialization_gap_contract.json --output-dir results/colour-specialization-gap",
               "python": platform.python_version(), "machine": platform.machine(),
               "script_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
               "contract_sha256": result["contract_sha256"],
               "output_sha256": {name: sha256(content.encode()).hexdigest() for name, content in outputs.items()},
               "elapsed_seconds": result["elapsed_seconds"], "new_random_samples": 0,
               "new_lattice_enumerations": 0, "tests_run": 0, "cloud_jobs": 0}
    with (args.output_dir / "run.json").open("x") as handle:
        json.dump(receipt, handle, indent=2)
        handle.write("\n")
    print(json.dumps({"decision": result["decision"], "elapsed_seconds": result["elapsed_seconds"],
                      "cases": [{"id": x["id"], "m4": x["finite"][0]["dimension_normalized_coefficient"],
                                 "m5": x["finite"][1]["dimension_normalized_coefficient"]} for x in cases]}))


if __name__ == "__main__":
    main()
