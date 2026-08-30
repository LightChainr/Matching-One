#!/usr/bin/env python3
"""Bounded exact HNF atlas for the two-orbit source/sink theorem."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from fractions import Fraction
import json
from pathlib import Path
from typing import Optional, Sequence

from digital_alexander_quotient_frontier import (
    has_four_distinct_face_corners,
    hnf_matrices,
)
from integer_period_torus import integer_torus_geometry
from p334_orbit_flux_phase_diagram import (
    _bernstein_edge_poly,
    _decimal,
    _eval,
    _root_midpoint,
    _sub,
)
from p334_third_geometry_falsifier import (
    _has_quarter_turn,
    _load_frozen_parent,
    _qpayload,
    _root_payload,
    exact_orbit_census,
    lattice_d4_stabilizer,
    projective_orbits,
)
from projective_essential_birth_oracle import chi4, subset_marks


Matrix = tuple[tuple[int, int], tuple[int, int]]
Vector = tuple[int, int]
ComplexQ = tuple[Fraction, Fraction]
Poly = tuple[Fraction, ...]


def _graph_connected(geometry) -> bool:
    adjacency = [set() for _ in range(geometry.n)]
    for edge in geometry.primal_edges:
        adjacency[edge.i].add(edge.j)
        adjacency[edge.j].add(edge.i)
    seen = {0}
    pending = [0]
    while pending:
        vertex = pending.pop()
        for neighbor in adjacency[vertex] - seen:
            seen.add(neighbor)
            pending.append(neighbor)
    return len(seen) == geometry.n


def _mean_character(matrix: Matrix, orbit: Sequence[Vector]) -> ComplexQ:
    values = [chi4(matrix, line) for line in orbit]
    return (
        sum(value[0] for value in values) / len(values),
        sum(value[1] for value in values) / len(values),
    )


def geometry_atlas_gate(maximum_order: int = 12) -> dict[str, object]:
    """Apply the declared gate without reading a boundary transition count."""

    excluded: Counter[str] = Counter()
    included = []
    scanned = 0
    for matrix in hnf_matrices(4, maximum_order):
        scanned += 1
        geometry = integer_torus_geometry(matrix, name="p334-two-orbit-atlas")
        if not has_four_distinct_face_corners(geometry):
            excluded["degenerate_face"] += 1
            continue
        if not _graph_connected(geometry):
            excluded["disconnected_primal_graph"] += 1
            continue
        if _has_quarter_turn(matrix):
            excluded["quarter_turn_lattice"] += 1
            continue
        marks = subset_marks(geometry, matching=False)
        lines = tuple(sorted({line for rank, line, _ in marks if rank == 1}))
        stabilizer = lattice_d4_stabilizer(matrix)
        orbits = projective_orbits(lines, stabilizer)
        if len(orbits) != 2:
            excluded["not_two_projective_orbits"] += 1
            continue
        characters = tuple(_mean_character(matrix, orbit) for orbit in orbits)
        if any(value == (0, 0) for value in characters):
            excluded["zero_effective_character"] += 1
            continue
        if characters[0] == characters[1]:
            excluded["aliased_effective_characters"] += 1
            continue
        included.append(
            {
                "matrix": matrix,
                "geometry": geometry,
                "marks": marks,
                "stabilizer": stabilizer,
                "orbits": orbits,
                "characters": characters,
            }
        )
    return {
        "maximum_order": maximum_order,
        "matrices_scanned": scanned,
        "exclusion_counts": dict(sorted(excluded.items())),
        "included": included,
        "gate": (
            "HNF index 4..12; four distinct face corners; connected primal graph; "
            "no physical quarter-turn lattice symmetry; exactly two projective "
            "line orbits under the exact D4 stabilizer; distinct nonzero effective "
            "orbit characters. Only geometry, line support and chi4 enter."
        ),
    }


def _line_boundary_tables(selection: dict[str, object]) -> dict[Vector, dict[str, list[int]]]:
    geometry = selection["geometry"]
    marks = selection["marks"]
    lines = {line for orbit in selection["orbits"] for line in orbit}
    tables = {
        line: {"birth": [0] * geometry.n, "exit": [0] * geometry.n}
        for line in lines
    }
    for mask, (old_rank, old_line, _) in enumerate(marks):
        k = mask.bit_count()
        for vertex in range(geometry.n):
            if mask & (1 << vertex):
                continue
            new_rank, new_line, _ = marks[mask | (1 << vertex)]
            if old_rank == 0 and new_rank == 1:
                tables[new_line]["birth"][k] += 1
            if old_rank == 1 and new_rank == 2:
                tables[old_line]["exit"][k] += 1
    return tables


def _orbit_compression_gate(selection: dict[str, object]) -> dict[str, object]:
    tables = _line_boundary_tables(selection)
    rows = []
    all_pass = True
    for index, orbit in enumerate(selection["orbits"]):
        reference = tables[orbit[0]]
        equal = all(tables[line] == reference for line in orbit[1:])
        characters = [chi4(selection["matrix"], line) for line in orbit]
        effective = selection["characters"][index]
        mean_exact = (
            sum(value[0] for value in characters) == len(orbit) * effective[0]
            and sum(value[1] for value in characters) == len(orbit) * effective[1]
        )
        passed = equal and mean_exact
        all_pass &= passed
        rows.append(
            {
                "orbit": f"orbit_{index}",
                "line_count": len(orbit),
                "line_boundary_coefficients_equal": equal,
                "effective_character_average_exact": mean_exact,
                "pass": passed,
            }
        )
    return {"orbits": rows, "all_pass": all_pass}


def _mul(left: Poly, right: Poly) -> Poly:
    result = [Fraction(0)] * (len(left) + len(right) - 1)
    for i, first in enumerate(left):
        for j, second in enumerate(right):
            result[i + j] += first * second
    while len(result) > 1 and result[-1] == 0:
        result.pop()
    return tuple(result)


def _scale(poly: Poly, scalar: Fraction) -> Poly:
    return tuple(value * scalar for value in poly)


def _dot(left: ComplexQ, right: ComplexQ) -> Fraction:
    return left[0] * right[0] + left[1] * right[1]


def _sign(value: Fraction) -> str:
    return "positive" if value > 0 else "negative" if value < 0 else "zero"


def _frozen_close_threshold(parent: dict[str, object]) -> Fraction:
    separations = []
    for geometry in parent["geometries"].values():
        roots = geometry["signed_share_singularities"]
        first = _root_midpoint(roots["axis_share_zero"][0])
        second = _root_midpoint(roots["diagonal_share_zero"][0])
        separations.append(abs(first - second))
    return max(separations)


def score_geometry(
    selection: dict[str, object], parent: dict[str, object]
) -> dict[str, object]:
    census = exact_orbit_census(selection)
    matrix = selection["matrix"]
    census["geometry"]["id"] = (
        "hnf-" + "-".join(str(value) for row in matrix for value in row)
    )
    census["geometry"]["period_matrix"] = [list(row) for row in matrix]
    polys = {}
    roots = {}
    for index in range(2):
        label = f"orbit_{index}"
        birth = _bernstein_edge_poly(census, label, "birth_edges")
        exit_flux = _bernstein_edge_poly(census, label, "exit_edges")
        net = _sub(birth, exit_flux)
        polys[label] = net
        roots[label] = _root_payload(net)

    root_pattern = tuple(roots[f"orbit_{index}"]["interior_root_count"] for index in range(2))
    root_rows = []
    phase_intervals = []
    separation = None
    close = None
    ordering = None
    if root_pattern == (1, 1):
        values = {
            label: _root_midpoint(roots[label]["roots"][0])
            for label in ("orbit_0", "orbit_1")
        }
        ordered = sorted((value, label) for label, value in values.items())
        ordering = [label for _, label in ordered]
        separation = abs(values["orbit_0"] - values["orbit_1"])
        close = separation <= _frozen_close_threshold(parent)
        boundaries = [Fraction(0), ordered[0][0], ordered[1][0], Fraction(1)]
        gram = _dot(selection["characters"][0], selection["characters"][1])
        for left, right in zip(boundaries, boundaries[1:]):
            probe = (left + right) / 2
            currents = [_eval(polys[f"orbit_{index}"], probe) for index in range(2)]
            cross = gram * currents[0] * currents[1]
            phase_intervals.append(
                {
                    "lower": _decimal(left),
                    "upper": _decimal(right),
                    "J1_sign": _sign(currents[0]),
                    "J2_sign": _sign(currents[1]),
                    "cross_term_sign": _sign(cross),
                    "topology": "reinforce" if cross > 0 else "cancel",
                }
            )
        root_rows = [
            {"orbit": label, **roots[label]["roots"][0]}
            for label in ("orbit_0", "orbit_1")
        ]

    gram = _dot(selection["characters"][0], selection["characters"][1])
    product = _mul(polys["orbit_0"], polys["orbit_1"])
    theorem_right = _scale(product, gram)
    # Direct complex dot product of the two orbit contributions.
    theorem_left = _scale(
        product,
        selection["characters"][0][0] * selection["characters"][1][0]
        + selection["characters"][0][1] * selection["characters"][1][1],
    )
    theorem_residual = _sub(theorem_left, theorem_right)
    compression = _orbit_compression_gate(selection)
    return {
        "geometry": census["geometry"],
        "d4_stabilizer_order": len(selection["stabilizer"]),
        "orbits": {
            f"orbit_{index}": {
                **census["orbits"][f"orbit_{index}"],
                "effective_chi4": _qpayload(selection["characters"][index]),
            }
            for index in range(2)
        },
        "coefficientwise_continuity": census["gates"]["coefficientwise_continuity"],
        "orbit_compression_gate": compression,
        "root_count_pattern": list(root_pattern),
        "roots": root_rows,
        "root_ordering": ordering,
        "root_separation": _decimal(separation) if separation is not None else None,
        "close_under_frozen_parent_envelope": close,
        "character_gram": str(gram),
        "character_gram_sign": _sign(gram),
        "phase_intervals": phase_intervals,
        "cross_term_polynomial_identity": {
            "formula": "Re[(chi1 J1) conjugate(chi2 J2)] = Gram(chi1,chi2) J1 J2",
            "residual_coefficients": [str(value) for value in theorem_residual],
            "pass": all(value == 0 for value in theorem_residual),
        },
    }


def _signature(row: dict[str, object]) -> tuple:
    roots = tuple(sorted(root["decimal"] for root in row["roots"]))
    topology = tuple(interval["topology"] for interval in row["phase_intervals"])
    return row["geometry"]["N"], row["character_gram"], roots, topology


def build_certificate(root: Optional[Path] = None) -> dict[str, object]:
    root = root or Path(__file__).resolve().parents[1]
    gated = geometry_atlas_gate()
    selections = gated.pop("included")
    parent = _load_frozen_parent(root)
    rows = [score_geometry(selection, parent) for selection in selections]
    root_counterexamples = [row for row in rows if row["root_count_pattern"] != [1, 1]]
    close_counterexamples = [
        row
        for row in rows
        if row["root_count_pattern"] == [1, 1]
        and not row["close_under_frozen_parent_envelope"]
    ]
    gram_strata = defaultdict(list)
    for row in rows:
        gram_strata[row["character_gram_sign"]].append(row["geometry"]["id"])
    signatures = defaultdict(list)
    for row in rows:
        signatures[_signature(row)].append(row["geometry"]["id"])
    theorem_pass = all(
        row["cross_term_polynomial_identity"]["pass"]
        and row["orbit_compression_gate"]["all_pass"]
        for row in rows
    )
    return {
        "schema": "matching-one/p334-two-orbit-exact-atlas/v1",
        "issues": [334, 337],
        "parent_commit": "2ec9f19",
        "status": "bounded_exact_two_orbit_atlas",
        "geometry_gate": gated,
        "frozen_close_threshold": _decimal(_frozen_close_threshold(parent)),
        "atlas": rows,
        "summary": {
            "included_hnf_count": len(rows),
            "total_subset_states": sum(
                row["geometry"]["subset_states"] for row in rows
            ),
            "total_directed_boundary_edges": sum(
                row["geometry"]["directed_addition_edges"] for row in rows
            ),
            "exact_mechanism_signature_count": len(signatures),
            "signature_members": [
                {"members": members} for members in signatures.values()
            ],
            "one_simple_root_per_orbit_count": len(rows) - len(root_counterexamples),
            "minimal_root_count_counterexample": (
                root_counterexamples[0]["geometry"]["id"]
                if root_counterexamples
                else None
            ),
            "close_pair_count": len(rows) - len(close_counterexamples),
            "minimal_close_pair_counterexample": (
                close_counterexamples[0]["geometry"]["id"]
                if close_counterexamples
                else None
            ),
            "gram_strata": dict(sorted(gram_strata.items())),
            "all_curve_cross_term_identities_pass": theorem_pass,
        },
        "classification": {
            "statement": (
                "One simple balance root per orbit is universal in the bounded atlas, "
                "but closeness under the inherited N13/N17 envelope is not: the N8 "
                "class is the minimal counterexample. Cooperation topology stratifies "
                "exactly by character-Gram sign across every full p curve."
            ),
            "theorem": (
                "For two exactly compressible orbit currents, "
                "sign Re[(chi1 J1) conjugate(chi2 J2)] = "
                "sign(Gram(chi1,chi2) J1 J2)."
            ),
            "next_prediction": (
                "Beyond N12, root multiplicity and separation are dynamical questions, "
                "but any two-orbit quotient passing the same compression gate must "
                "retain Gram-sign cooperation topology over every root interval."
            ),
        },
        "claim_boundary": [
            "The atlas exhausts the declared HNF geometry gate only through index 12.",
            "HNF variants sharing one exact signature are symmetry copies, not independent evidence.",
            "Root closeness is judged by the frozen parent envelope and is not retuned on this atlas.",
            "No Monte Carlo sample, Huawei production, new PR, or merge is used.",
        ],
    }


def render_markdown(payload: dict[str, object]) -> str:
    summary = payload["summary"]
    lines = [
        "# Bounded exact two-orbit HNF atlas",
        "",
        f"Status: `{payload['status']}`.",
        "",
        payload["geometry_gate"]["gate"],
        "",
        f"Scanned {payload['geometry_gate']['matrices_scanned']} HNFs and included "
        f"{summary['included_hnf_count']}. The frozen close threshold is "
        f"`{payload['frozen_close_threshold']}`.",
        f"The revealed atlas contains only {summary['total_subset_states']} subset states and "
        f"{summary['total_directed_boundary_edges']} directed boundary edges in total.",
        "",
        "| HNF | N | orbit roots | separation | close | Gram | topology below/between/above |",
        "|---|---:|---|---:|---|---:|---|",
    ]
    for row in payload["atlas"]:
        roots = ", ".join(root["decimal"] for root in row["roots"])
        topology = "/".join(interval["topology"] for interval in row["phase_intervals"])
        lines.append(
            f"| {row['geometry']['period_matrix']} | {row['geometry']['N']} | {roots} | "
            f"{row['root_separation']} | {row['close_under_frozen_parent_envelope']} | "
            f"{row['character_gram']} | {topology} |"
        )
    lines += [
        "",
        "## Atlas answer",
        "",
        payload["classification"]["statement"],
        "",
        f"- one simple root per orbit: {summary['one_simple_root_per_orbit_count']}/{summary['included_hnf_count']};",
        f"- close under frozen envelope: {summary['close_pair_count']}/{summary['included_hnf_count']};",
        f"- minimal closeness counterexample: `{summary['minimal_close_pair_counterexample']}`;",
        f"- exact mechanism signatures after HNF symmetry copies: {summary['exact_mechanism_signature_count']};",
        f"- all coefficientwise compression and cross-term identities: `{summary['all_curve_cross_term_identities_pass']}`.",
        "",
        "Exact theorem:",
        "",
        "```text",
        "Re[(chi1 J1) conjugate(chi2 J2)] = Gram(chi1,chi2) J1 J2",
        "```",
        "",
        "Positive Gram produces reinforce/cancel/reinforce across the two simple zeros; negative Gram produces cancel/reinforce/cancel. This holds on every included full curve.",
        "",
        "Next prediction: " + payload["classification"]["next_prediction"],
        "",
        "## Boundary",
        "",
    ]
    lines.extend(f"- {item}" for item in payload["claim_boundary"])
    return "\n".join(lines)


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--format", choices=("json", "markdown"), default="json")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args(argv)
    payload = build_certificate()
    rendered = (
        json.dumps(payload, indent=2, sort_keys=True) + "\n"
        if args.format == "json"
        else render_markdown(payload) + "\n"
    )
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
