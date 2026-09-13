#!/usr/bin/env python3
"""Vertex-subset matching defect polynomial, not a Tutte specialization.

The finite matching polynomial is the Bernoulli generating function of the
wrapping-event difference

    q(ω) = I_wrap(black NN) - I_wrap(white NN+NNN)  ∈ {-1, 0, +1}.

On the implemented square tori the P34 identity says q is configuration-
identical across wrapping channels and

    C_black - C_white = q + V - E + F0,

so the matching function is also M(p) = E[q].  This is a vertex-subset
homology-event generating function.

It is not an edge-subset Tutte / Bollobás–Riordan / Krushkal specialization:
site percolation sums over vertex subsets, and the black/white graphs are
different (NN versus NN+NNN) except on a self-matching triangulation.

On the C4 N=10 triangulation the two graphs coincide, occupation complement
is an involution with q(ω^c) = -q(ω), and the defect polynomial reproduces
the exact Beta(3,3) control 12p^5-30p^4+20p^3-1.

First deliverable of issue #144: exact generating-function identification
plus the Tutte obstruction, verified on axis L=2/3 and C4 N=10.  No deletion-
contraction algorithm is claimed.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List, Sequence

from c4_self_matching_exact import c4_self_matching_torus
from euler_motif_controls import (
    WRAPPING_CHANNELS,
    configuration_identity,
    wrapping_differences,
)
from exact_matching_polynomial import bernstein_to_power, polynomial_string
from integer_period_torus import IntegerTorusGeometry, axis_integer_torus

ROOT = Path(__file__).resolve().parents[1]

# Bernstein a_k = sum_{|ω|=k} q(ω), locked against the committed matching polynomials.
AXIS_BERNSTEIN = {
    2: [-1, -4, -2, 4, 1],
    3: [-1, -9, -36, -78, -90, -36, 36, 36, 9, 1],
}
N10_BERNSTEIN = [-1, -10, -45, -100, -100, 0, 100, 100, 45, 10, 1]
N10_POWER = [-1, 0, 0, 20, -30, 12]


def wrapping_event(geometry: IntegerTorusGeometry, active: Sequence[bool]) -> int:
    wrapping = wrapping_differences(geometry, active)
    if len(set(wrapping.values())) != 1:
        raise AssertionError("wrapping-difference channels are not identical")
    return wrapping["either"]


def bernstein_from_q(geometry: IntegerTorusGeometry) -> List[int]:
    if geometry.n > 16:
        raise ValueError("exhaustive defect enumeration is limited to N<=16")
    counts = [0] * (geometry.n + 1)
    for mask in range(1 << geometry.n):
        active = [bool((mask >> vertex) & 1) for vertex in range(geometry.n)]
        counts[sum(1 for value in active if value)] += wrapping_event(geometry, active)
    return counts


def complement_involution_failures(geometry: IntegerTorusGeometry) -> int:
    failures = 0
    n = geometry.n
    for mask in range(1 << n):
        active = [bool((mask >> vertex) & 1) for vertex in range(n)]
        q = wrapping_event(geometry, active)
        complement = [not value for value in active]
        if wrapping_event(geometry, complement) != -q:
            failures += 1
    return failures


def tutte_obstruction(geometry: IntegerTorusGeometry) -> Dict[str, object]:
    primal = {(min(edge.i, edge.j), max(edge.i, edge.j)) for edge in geometry.primal_edges}
    matching = {(min(edge.i, edge.j), max(edge.i, edge.j)) for edge in geometry.matching_edges}
    return {
        "name": geometry.name,
        "N": geometry.n,
        "vertex_subset_model": True,
        "primal_edges": len(geometry.primal_edges),
        "matching_edges": len(geometry.matching_edges),
        "primal_equals_matching": primal == matching,
        "tutte_specialization": False,
        "reason": (
            "edge_subset_Tutte_does_not_apply_to_site_occupation"
            if primal != matching
            else "self_matching_still_sums_vertex_subsets_by_wrapping_event_not_edge_subsets"
        ),
    }


def axis_defect(L: int) -> Dict[str, object]:
    geometry = axis_integer_torus(L)
    bernstein = bernstein_from_q(geometry)
    # Cross-check P34: q equals C_black-C_white-(V-E+F0) on every configuration.
    identity_failures = 0
    for mask in range(1 << geometry.n):
        active = [bool((mask >> vertex) & 1) for vertex in range(geometry.n)]
        if configuration_identity(geometry, active, mask).residual != 0:
            identity_failures += 1
    power = bernstein_to_power(bernstein)
    return {
        "name": geometry.name,
        "L": L,
        "N": geometry.n,
        "bernstein": bernstein,
        "power": power,
        "power_string": polynomial_string(power),
        "identity_failures": identity_failures,
        "complement_involution_failures": complement_involution_failures(geometry),
        "obstruction": tutte_obstruction(geometry),
        "channels": list(WRAPPING_CHANNELS),
    }


def n10_defect() -> Dict[str, object]:
    geometry = c4_self_matching_torus(3, 1)
    bernstein = bernstein_from_q(geometry)
    power = bernstein_to_power(bernstein)
    return {
        "name": "c4-self-matching-3-1",
        "N": geometry.n,
        "bernstein": bernstein,
        "power": power,
        "power_string": polynomial_string(power),
        "complement_involution_failures": complement_involution_failures(geometry),
        "obstruction": tutte_obstruction(geometry),
        "self_matching": True,
        "beta33": True,
    }


def run_suite() -> Dict[str, object]:
    axis = [axis_defect(2), axis_defect(3)]
    n10 = n10_defect()
    passed = (
        all(row["identity_failures"] == 0 for row in axis)
        and all(row["bernstein"] == AXIS_BERNSTEIN[row["L"]] for row in axis)
        and n10["bernstein"] == N10_BERNSTEIN
        and n10["power"] == N10_POWER
        and n10["complement_involution_failures"] == 0
        and all(not row["obstruction"]["tutte_specialization"] for row in axis + [n10])
    )
    return {
        "schema": "matching defect polynomial v1",
        "identification": "M(p) is the Bernoulli generating function of the wrapping-event q",
        "tutte_specialization": False,
        "passed": passed,
        "axis": axis,
        "n10": n10,
    }


def render_report(payload: Dict[str, object]) -> str:
    lines = [
        "# Matching defect polynomial (vertex-subset, not Tutte)",
        "",
        "Source: `scripts/matching_defect_polynomial.py`.",
        "Claim level: C5 identification `M=E[q]` on the enumerated quotients; C5 obstruction",
        "against an edge-subset Tutte specialization. Issue #144 first deliverable.",
        "",
        "## Identification",
        "",
        "```text",
        "q(ω) = I_wrap(black NN) - I_wrap(white NN+NNN) ∈ {-1,0,+1}",
        "a_k  = sum_{|ω|=k} q(ω)",
        "M(p) = sum_k a_k p^k (1-p)^{N-k} = E[q]",
        "```",
        "",
        "On square tori this is the same `q` as in the P34 identity",
        "`C_black-C_white = q+V-E+F0`.",
        "",
        "## Axis enumerations",
        "",
    ]
    for row in payload["axis"]:
        lines.extend(
            [
                f"### {row['name']} L={row['L']}, N={row['N']}",
                "",
                f"Bernstein `a_k`: `{row['bernstein']}`",
                f"Power basis: `{row['power_string']}`",
                f"P34 identity failures: {row['identity_failures']}",
                f"Complement involution failures: {row['complement_involution_failures']}",
                f"Primal equals matching: {row['obstruction']['primal_equals_matching']}",
                "",
            ]
        )
    n10 = payload["n10"]
    lines.extend(
        [
            "## C4 self-matching N=10",
            "",
            "Primal and matching graphs coincide. Occupation complement is an involution",
            "`q(ω^c)=-q(ω)`, so `M(1-p)=-M(p)`.",
            "",
            f"Bernstein `a_k`: `{n10['bernstein']}`",
            f"Power basis: `{n10['power_string']}`",
            f"Complement involution failures: {n10['complement_involution_failures']}",
            "",
            "## Tutte obstruction",
            "",
            "Site matching is a vertex-subset model. The black and white graphs are NN",
            "versus NN+NNN except on a self-matching triangulation, and even there the",
            "sum is over wrapping events of vertex subsets, not over edge subsets of a",
            "single ribbon graph. No Bollobás–Riordan / Krushkal specialization is used.",
            "",
            "A cheaper deletion-contraction algorithm is **not** established.",
            "",
        ]
    )
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()
    payload = run_suite()
    status = "PASS" if payload["passed"] else "FAIL"
    print("matching defect polynomial " + status)
    for row in payload["axis"]:
        print("{name} L={L}: {power_string}".format(**row))
    print("n10: {power_string}".format(**payload["n10"]))
    if args.report is not None:
        args.report.write_text(render_report(payload), encoding="utf-8")
        print("wrote " + str(args.report))
    return 0 if payload["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
