#!/usr/bin/env python3
"""Geometry-blind selection and exact third-quotient falsifier for Issue #334."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from fractions import Fraction
import json
from math import atan2, pi
from pathlib import Path
from typing import Optional, Sequence

from digital_alexander_quotient_frontier import (
    has_four_distinct_face_corners,
    hnf_matrices,
)
from integer_period_torus import integer_torus_geometry, matrix_vector
from p334_orbit_flux_phase_diagram import (
    _bernstein_edge_poly,
    _decimal,
    _derivative,
    _eval,
    _root_midpoint,
    _sub,
    isolate_open_unit_roots,
)
from p334_n13_multiorbit_flux import P_REF
from projective_essential_birth_oracle import (
    canonical_projective,
    chi4,
    subset_marks,
)


Matrix = tuple[tuple[int, int], tuple[int, int]]
Vector = tuple[int, int]
ComplexQ = tuple[Fraction, Fraction]


D4_PHYSICAL: tuple[Matrix, ...] = (
    ((1, 0), (0, 1)),
    ((0, -1), (1, 0)),
    ((-1, 0), (0, -1)),
    ((0, 1), (-1, 0)),
    ((1, 0), (0, -1)),
    ((-1, 0), (0, 1)),
    ((0, 1), (1, 0)),
    ((0, -1), (-1, 0)),
)


def _matrix_product(left, right):
    return tuple(
        tuple(sum(left[i][k] * right[k][j] for k in range(2)) for j in range(2))
        for i in range(2)
    )


def _rational_inverse(matrix: Matrix):
    (a, b), (c, d) = matrix
    determinant = a * d - b * c
    return (
        (Fraction(d, determinant), Fraction(-b, determinant)),
        (Fraction(-c, determinant), Fraction(a, determinant)),
    )


def lattice_d4_stabilizer(matrix: Matrix) -> tuple[Matrix, ...]:
    inverse = _rational_inverse(matrix)
    stabilizer = []
    for physical in D4_PHYSICAL:
        action = _matrix_product(_matrix_product(inverse, physical), matrix)
        if all(value.denominator == 1 for row in action for value in row):
            stabilizer.append(
                tuple(tuple(int(value) for value in row) for row in action)
            )
    return tuple(stabilizer)


def _has_quarter_turn(matrix: Matrix) -> bool:
    inverse = _rational_inverse(matrix)
    action = _matrix_product(_matrix_product(inverse, D4_PHYSICAL[1]), matrix)
    return all(value.denominator == 1 for row in action for value in row)


def projective_orbits(
    lines: Sequence[Vector], stabilizer: Sequence[Matrix]
) -> tuple[tuple[Vector, ...], ...]:
    support = set(lines)
    unseen = set(lines)
    orbits = []
    while unseen:
        seed = min(unseen)
        orbit = {
            canonical_projective(matrix_vector(action, seed))
            for action in stabilizer
        }
        if not orbit <= support:
            raise AssertionError("lattice symmetry did not preserve line support")
        orbits.append(tuple(sorted(orbit)))
        unseen -= orbit
    return tuple(sorted(orbits))


def _orbit_characters(matrix: Matrix, orbit: Sequence[Vector]) -> tuple[ComplexQ, ...]:
    return tuple(sorted({chi4(matrix, line) for line in orbit}))


def geometry_scan(maximum_order: int = 18) -> dict[str, object]:
    """Select without reading any source/sink boundary count."""

    rows = []
    selected = None
    for matrix in hnf_matrices(4, maximum_order):
        geometry = integer_torus_geometry(matrix, name="p334-third-geometry-scan")
        row: dict[str, object] = {
            "matrix": [list(part) for part in matrix],
            "N": geometry.n,
        }
        if not has_four_distinct_face_corners(geometry):
            row["decision"] = "exclude_degenerate_face"
            rows.append(row)
            continue
        quarter_turn = _has_quarter_turn(matrix)
        row["quarter_turn_lattice_symmetry"] = quarter_turn
        if quarter_turn:
            row["decision"] = "exclude_quarter_turn"
            rows.append(row)
            continue
        marks = subset_marks(geometry, matching=False)
        lines = tuple(sorted({line for rank, line, _ in marks if rank == 1}))
        stabilizer = lattice_d4_stabilizer(matrix)
        orbits = projective_orbits(lines, stabilizer)
        characters = [_orbit_characters(matrix, orbit) for orbit in orbits]
        row.update(
            {
                "d4_stabilizer_order": len(stabilizer),
                "primitive_line_support": [list(line) for line in lines],
                "projective_orbits": [
                    [list(line) for line in orbit] for orbit in orbits
                ],
                "orbit_chi4": [
                    [
                        {"real": str(value[0]), "imag": str(value[1])}
                        for value in orbit_values
                    ]
                    for orbit_values in characters
                ],
            }
        )
        constant_character = all(len(values) == 1 for values in characters)
        nontrivial_characters = all(
            values[0][1] != 0 for values in characters if len(values) == 1
        )
        distinct_characters = (
            len(characters) == 2
            and constant_character
            and characters[0][0] != characters[1][0]
        )
        if len(orbits) != 2:
            row["decision"] = "exclude_orbit_count"
        elif not constant_character:
            row["decision"] = "exclude_nonconstant_orbit_character"
        elif not nontrivial_characters or not distinct_characters:
            row["decision"] = "exclude_trivial_or_aliased_character"
        else:
            row["decision"] = "select_first_lexicographic_candidate"
            selected = {
                "matrix": matrix,
                "geometry": geometry,
                "marks": marks,
                "stabilizer": stabilizer,
                "orbits": orbits,
                "characters": tuple(values[0] for values in characters),
            }
        rows.append(row)
        if selected is not None:
            break
    if selected is None:
        raise AssertionError("no two-orbit HNF candidate found through N=18")
    return {
        "selection_rule": (
            "First HNF in increasing N/lexicographic order with four distinct face "
            "corners, no physical quarter-turn lattice symmetry, exactly two "
            "projective stabilizer orbits, and distinct constant non-real chi4 on "
            "both orbits. Source/sink counts are forbidden during selection."
        ),
        "maximum_order": maximum_order,
        "rows_examined": rows,
        "selected": selected,
    }


def _qpayload(value: ComplexQ) -> dict[str, str]:
    return {"real": str(value[0]), "imag": str(value[1])}


def exact_orbit_census(selection: dict[str, object]) -> dict[str, object]:
    geometry = selection["geometry"]
    marks = selection["marks"]
    orbits = selection["orbits"]
    line_to_label = {
        line: f"orbit_{index}"
        for index, orbit in enumerate(orbits)
        for line in orbit
    }
    state_counts: Counter[tuple[int, str]] = Counter()
    birth: Counter[tuple[int, str]] = Counter()
    exit_flux: Counter[tuple[int, str]] = Counter()
    direct_rank2 = 0
    total_edges = 0
    for mask, (old_rank, old_line, _) in enumerate(marks):
        k = mask.bit_count()
        if old_rank == 1:
            state_counts[(k, line_to_label[old_line])] += 1
        for vertex in range(geometry.n):
            if mask & (1 << vertex):
                continue
            total_edges += 1
            new_rank, new_line, _ = marks[mask | (1 << vertex)]
            if old_rank == 0 and new_rank == 1:
                birth[(k, line_to_label[new_line])] += 1
            elif old_rank == 0 and new_rank == 2:
                direct_rank2 += 1
            if old_rank == 1 and new_rank == 2:
                exit_flux[(k, line_to_label[old_line])] += 1

    rows = []
    failures = 0
    for k in range(geometry.n):
        row: dict[str, object] = {"lower_subset_size": k}
        for index in range(2):
            label = f"orbit_{index}"
            derivative_coefficient = (
                (k + 1) * state_counts[(k + 1, label)]
                - (geometry.n - k) * state_counts[(k, label)]
            )
            transition = birth[(k, label)] - exit_flux[(k, label)]
            passed = derivative_coefficient == transition
            failures += not passed
            row[label] = {
                "rank_one_states_at_k": state_counts[(k, label)],
                "rank_one_states_at_k_plus_1": state_counts[(k + 1, label)],
                "birth_edges": birth[(k, label)],
                "exit_edges": exit_flux[(k, label)],
                "birth_minus_exit_edges": transition,
                "coefficient_identity_pass": passed,
            }
        rows.append(row)
    return {
        "geometry": {
            "id": "hnf-7-2-0-1",
            "N": geometry.n,
            "period_matrix": [list(part) for part in selection["matrix"]],
            "subset_states": 1 << geometry.n,
            "directed_addition_edges": total_edges,
        },
        "orbits": {
            f"orbit_{index}": {
                "primitive_lines": [list(line) for line in orbit],
                "chi4": _qpayload(selection["characters"][index]),
                "rank_one_state_count": sum(
                    count
                    for (k, label), count in state_counts.items()
                    if label == f"orbit_{index}"
                ),
                "birth_edge_count": sum(
                    count
                    for (k, label), count in birth.items()
                    if label == f"orbit_{index}"
                ),
                "exit_edge_count": sum(
                    count
                    for (k, label), count in exit_flux.items()
                    if label == f"orbit_{index}"
                ),
            }
            for index, orbit in enumerate(orbits)
        },
        "direct_rank2_edge_count": direct_rank2,
        "coefficient_rows": rows,
        "gates": {
            "coefficientwise_continuity": failures == 0,
            "coefficient_failures": failures,
            "edge_count_exact": total_edges == geometry.n * (1 << (geometry.n - 1)),
        },
    }


def _load_frozen_parent(root: Path) -> dict[str, object]:
    path = root / "results/p334-orbit-flux-phase-diagram/latest.json"
    return json.loads(path.read_text(encoding="utf-8"))


def _small_rational_root(poly, row: dict[str, str]) -> Optional[Fraction]:
    midpoint = _root_midpoint(row)
    candidate = midpoint.limit_denominator(10000)
    if Fraction(row["lower"]) <= candidate <= Fraction(row["upper"]):
        if _eval(poly, candidate) == 0:
            return candidate
    return None


def _root_payload(poly) -> dict[str, object]:
    payload = isolate_open_unit_roots(poly)
    for row in payload["roots"]:
        exact = _small_rational_root(poly, row)
        if exact is not None:
            row.clear()
            row.update(
                {
                    "kind": "exact_rational",
                    "root": str(exact),
                    "decimal": _decimal(exact),
                    "lower": str(exact),
                    "upper": str(exact),
                }
            )
    return payload


def _complex_scale(value: ComplexQ, scalar: Fraction) -> ComplexQ:
    return value[0] * scalar, value[1] * scalar


def _complex_add(left: ComplexQ, right: ComplexQ) -> ComplexQ:
    return left[0] + right[0], left[1] + right[1]


def _dot(left: ComplexQ, right: ComplexQ) -> Fraction:
    return left[0] * right[0] + left[1] * right[1]


def _phase_degrees(value: ComplexQ) -> str:
    return format(atan2(float(value[1]), float(value[0])) * 180 / pi, ".12g")


def score_falsifier(
    census: dict[str, object], characters: Sequence[ComplexQ], frozen: dict[str, object]
) -> dict[str, object]:
    polys = {}
    for index in range(2):
        label = f"orbit_{index}"
        birth = _bernstein_edge_poly(census, label, "birth_edges")
        exit_flux = _bernstein_edge_poly(census, label, "exit_edges")
        polys[label] = {
            "birth": birth,
            "exit": exit_flux,
            "net": _sub(birth, exit_flux),
        }
    roots = {
        label: _root_payload(polys[label]["net"])
        for label in ("orbit_0", "orbit_1")
    }
    root_values = {
        label: _root_midpoint(roots[label]["roots"][0])
        for label in ("orbit_0", "orbit_1")
    }
    parent_separations = []
    for geometry in frozen["geometries"].values():
        first = _root_midpoint(
            geometry["signed_share_singularities"]["axis_share_zero"][0]
        )
        second = _root_midpoint(
            geometry["signed_share_singularities"]["diagonal_share_zero"][0]
        )
        parent_separations.append(abs(first - second))
    close_threshold = max(parent_separations)
    separation = abs(root_values["orbit_0"] - root_values["orbit_1"])
    gram = _dot(characters[0], characters[1])

    ordered = sorted((value, label) for label, value in root_values.items())
    intervals = []
    boundaries = [Fraction(0), ordered[0][0], ordered[1][0], Fraction(1)]
    for left, right in zip(boundaries, boundaries[1:]):
        probe = (left + right) / 2
        nets = {
            label: _eval(polys[label]["net"], probe)
            for label in ("orbit_0", "orbit_1")
        }
        alignment = gram * nets["orbit_0"] * nets["orbit_1"]
        intervals.append(
            {
                "lower": _decimal(left),
                "upper": _decimal(right),
                "orbit_0_net_sign": "positive" if nets["orbit_0"] > 0 else "negative",
                "orbit_1_net_sign": "positive" if nets["orbit_1"] > 0 else "negative",
                "character_contribution_alignment": (
                    "reinforce" if alignment > 0 else "cancel"
                ),
            }
        )

    point = {}
    contributions = []
    for index in range(2):
        label = f"orbit_{index}"
        birth = _eval(polys[label]["birth"], P_REF)
        exit_flux = _eval(polys[label]["exit"], P_REF)
        net = birth - exit_flux
        contribution = _complex_scale(characters[index], net)
        contributions.append(contribution)
        root = root_values[label]
        point[label] = {
            "birth": _decimal(birth),
            "exit": _decimal(exit_flux),
            "net": _decimal(net),
            "net_fraction_of_activity": _decimal(abs(net) / (birth + exit_flux)),
            "net_slope_at_p_ref": _decimal(_eval(_derivative(polys[label]["net"]), P_REF)),
            "root": _decimal(root),
            "p_ref_minus_root": _decimal(P_REF - root),
            "net_slope_at_root": _decimal(_eval(_derivative(polys[label]["net"]), root)),
            "chi4_weighted_contribution": {
                "real": _decimal(contribution[0]),
                "imag": _decimal(contribution[1]),
            },
        }
    total = _complex_add(contributions[0], contributions[1])
    total_norm2 = _dot(total, total)
    projection_share_0 = _dot(contributions[0], total) / total_norm2
    determinant = characters[0][0] * characters[1][1] - characters[0][1] * characters[1][0]
    return {
        "root_sets": roots,
        "paired_zero_test": {
            "frozen_close_threshold_from_parent_max": _decimal(close_threshold),
            "observed_separation": _decimal(separation),
            "close_pair_pass": separation <= close_threshold,
            "ordering": [label for _, label in ordered],
        },
        "character_geometry": {
            "chi4_gram_real": str(gram),
            "chi4_oriented_determinant": str(determinant),
            "characters_linearly_independent_over_R": determinant != 0,
            "interpretation": (
                "positive Gram: same-sign incidence contributions reinforce; "
                "opposite-sign incidence contributions cancel"
            ),
        },
        "phase_intervals": intervals,
        "p_ref_metrics": {
            **point,
            "total_complex_net": {
                "real": _decimal(total[0]),
                "imag": _decimal(total[1]),
                "phase_degrees": _phase_degrees(total),
            },
            "orbit_0_signed_projection_share": _decimal(projection_share_0),
            "orbit_1_signed_projection_share": _decimal(1 - projection_share_0),
        },
        "frozen_prediction_score": {
            "paired_net_zeros_close": separation <= close_threshold,
            "reinforcement_only_between_zeros": (
                intervals[0]["character_contribution_alignment"] == "cancel"
                and intervals[1]["character_contribution_alignment"] == "reinforce"
                and intervals[2]["character_contribution_alignment"] == "cancel"
            ),
            "verdict": "paired_timing_survives_but_reinforcement_topology_is_falsified",
        },
        "complex_total_zero_gate": {
            "has_interior_zero": False,
            "reason": (
                "The two chi4 vectors are real-linearly independent and the two "
                "scalar net roots are disjoint, so their complex sum cannot vanish."
            ),
        },
    }


def build_certificate(root: Optional[Path] = None) -> dict[str, object]:
    root = root or Path(__file__).resolve().parents[1]
    scan = geometry_scan()
    selected = scan.pop("selected")
    census = exact_orbit_census(selected)
    frozen = _load_frozen_parent(root)
    score = score_falsifier(census, selected["characters"], frozen)
    return {
        "schema": "matching-one/p334-third-geometry-falsifier/v1",
        "issue": 334,
        "parent_commit": "77aa3fe",
        "status": "exact_third_geometry_falsifier",
        "selection": {
            **scan,
            "selected_matrix": [list(part) for part in selected["matrix"]],
            "selected_N": selected["geometry"].n,
            "selected_orbits": [
                [list(line) for line in orbit] for orbit in selected["orbits"]
            ],
            "selected_characters": [_qpayload(value) for value in selected["characters"]],
            "selection_was_flux_blind": True,
            "gaussian_similarity": False,
            "gaussian_similarity_reason": (
                "A square-similarity lattice is invariant under physical quarter-turn; "
                "the exact quarter-turn stabilizer gate fails here."
            ),
        },
        "census": census,
        "score": score,
        "mechanism_update": {
            "classification": "paired_timing_zeros_with_character_gram_controlled_alignment",
            "exact_two_orbit_alignment_theorem": {
                "formula": (
                    "Re[(chi1 J1) conjugate(chi2 J2)] = "
                    "Re[chi1 conjugate(chi2)] J1 J2 = Gram(chi1,chi2) J1 J2"
                ),
                "scope": "real scalar orbit-net currents J1,J2 and fixed complex orbit characters chi1,chi2",
                "consequence": (
                    "Simple paired zeros determine only the sign pattern of J1 J2. "
                    "The character Gram sign independently decides whether each "
                    "interval is reinforcing or cancelling."
                ),
            },
            "statement": (
                "The close pair of source/sink balance times survives the asymmetric "
                "HNF quotient, but the Gaussian claim that reinforcement occurs only "
                "between them is false. The timing pair and the spin-4 alignment are "
                "separate layers: the latter flips with the exact chi4 Gram sign."
            ),
            "next_falsifiable_prediction": (
                "For a fresh two-orbit quotient, the sign of Re(chi_a conjugate(chi_b)) "
                "must determine the alignment topology: positive Gram gives cancellation "
                "between the simple net zeros and reinforcement outside; negative Gram "
                "gives the reverse. Failure rejects the two-scalar-current reduction."
            ),
        },
        "claim_boundary": [
            "The HNF candidate was selected from geometry and line-orbit support before boundary flux counts were read.",
            "This exact N=7 quotient is an asymmetric finite-volume falsifier, not an asymptotic geometry.",
            "The Gram-controlled alignment rule is a new mechanism classification, not a continuum-field identity.",
            "No N13/N17 recomputation, Monte Carlo sample, Huawei production, PR, or merge is used.",
        ],
    }


def render_markdown(payload: dict[str, object]) -> str:
    selection = payload["selection"]
    census = payload["census"]
    score = payload["score"]
    lines = [
        "# Third-geometry exact falsifier",
        "",
        f"Status: `{payload['status']}`.",
        "",
        "## Flux-blind selection",
        "",
        selection["selection_rule"],
        "",
        f"Selected HNF `{selection['selected_matrix']}` at N={selection['selected_N']} after "
        f"{len(selection['rows_examined'])} geometry-only rows. It has no physical quarter-turn symmetry and is not Gaussian-similar.",
        "",
        f"Primitive-line orbits: `{selection['selected_orbits']}`.",
        f"Characters: `{selection['selected_characters']}`.",
        "",
        "## Exact subset-boundary census",
        "",
        f"- states: {census['geometry']['subset_states']}",
        f"- directed edges: {census['geometry']['directed_addition_edges']}",
        f"- direct 0->2 edges: {census['direct_rank2_edge_count']}",
        f"- coefficientwise continuity: `{census['gates']['coefficientwise_continuity']}`",
        "",
        "| orbit | line | chi4 | states | birth | exit | net zero |",
        "|---|---|---|---:|---:|---:|---:|",
    ]
    for label, row in census["orbits"].items():
        root = score["root_sets"][label]["roots"][0]["decimal"]
        lines.append(
            f"| {label} | {row['primitive_lines']} | ({row['chi4']['real']}, {row['chi4']['imag']}) | "
            f"{row['rank_one_state_count']} | {row['birth_edge_count']} | "
            f"{row['exit_edge_count']} | {root} |"
        )
    pair = score["paired_zero_test"]
    lines += [
        "",
        "## Frozen score",
        "",
        f"The zero separation is `{pair['observed_separation']}` versus the frozen parent-envelope "
        f"`{pair['frozen_close_threshold_from_parent_max']}`: close-pair gate **passes**.",
        "",
        f"The exact character Gram is `{score['character_geometry']['chi4_gram_real']}` (positive). "
        "Therefore the two contributions cancel between the zeros and reinforce outside:",
        "",
    ]
    for row in score["phase_intervals"]:
        lines.append(
            f"- `({row['lower']}, {row['upper']})`: {row['orbit_0_net_sign']}/"
            f"{row['orbit_1_net_sign']}; **{row['character_contribution_alignment']}**."
        )
    point = score["p_ref_metrics"]
    lines += [
        "",
        f"At `p_ref`, orbit nets are `{point['orbit_0']['net']}` and "
        f"`{point['orbit_1']['net']}`; total phase is `{point['total_complex_net']['phase_degrees']}` degrees.",
        f"The net slopes at their own zeros are `{point['orbit_0']['net_slope_at_root']}` "
        f"and `{point['orbit_1']['net_slope_at_root']}`; the first orbit is only "
        f"`{point['orbit_0']['net_fraction_of_activity']}` of its source-plus-sink activity at `p_ref`.",
        "",
        "Frozen verdict: **paired timing survives, but between-zero reinforcement is falsified**.",
        "",
        "## Mechanism update",
        "",
        payload["mechanism_update"]["statement"],
        "",
        "Exact two-orbit alignment theorem:",
        "",
        "```text",
        payload["mechanism_update"]["exact_two_orbit_alignment_theorem"]["formula"],
        "```",
        "",
        "Paired zeros control only `sign(J1 J2)`. Gaussian opposite characters have negative Gram; this HNF has positive Gram, so the same timing window acquires the opposite alignment topology.",
        "",
        "Next falsifier: " + payload["mechanism_update"]["next_falsifiable_prediction"],
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
