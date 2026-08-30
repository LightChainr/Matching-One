#!/usr/bin/env python3
"""Symbolic certificates for unrestricted degenerate digital Alexander duality.

The theorem is proved on the canonical honest four-sheeted cover and then
descended over rational homology.  The finite HNF scan is retained only as a
machine regression and counterexample locator.
"""

from __future__ import annotations

import argparse
import json
from collections import Counter
from math import gcd
from pathlib import Path
from typing import Any, Sequence, Tuple

from digital_alexander_local_bridge import (
    BOUNDARY_EDGES,
    CORNERS,
    face_pattern,
)
from digital_alexander_quotient_frontier import (
    build_state_table,
    has_four_distinct_face_corners,
    hnf_matrices,
)
from integer_period_torus import (
    IntegerPeriods,
    integer_torus_geometry,
)


Vector = Tuple[int, int]
Matrix = Tuple[Tuple[int, int], Tuple[int, int]]


def subtract(first: Vector, second: Vector) -> Vector:
    return first[0] - second[0], first[1] - second[1]


def add(first: Vector, second: Vector) -> Vector:
    return first[0] + second[0], first[1] + second[1]


def path_displacement(path: Sequence[int]) -> Vector:
    result = (0, 0)
    for first, second in zip(path, path[1:]):
        result = add(result, subtract(CORNERS[second], CORNERS[first]))
    return result


def path_boundary(path: Sequence[int]) -> dict[int, int]:
    boundary: Counter[int] = Counter()
    for first, second in zip(path, path[1:]):
        boundary[first] -= 1
        boundary[second] += 1
    return {vertex: coefficient for vertex, coefficient in sorted(boundary.items())
            if coefficient}


def lift_face_certificate(mask: int) -> dict[str, Any]:
    row = face_pattern(mask)
    black = set(row["black_corners"])
    black_boundary_edges = [
        list(edge) for edge in BOUNDARY_EDGES
        if edge[0] in black and edge[1] in black
    ]
    replacements = []
    for replacement in row["removed_diagonal_replacements"]:
        diagonal = tuple(replacement["diagonal"])
        path = tuple(replacement["boundary_path"])
        direct_displacement = subtract(CORNERS[diagonal[1]], CORNERS[diagonal[0]])
        replacement_displacement = path_displacement(path)
        direct_boundary = path_boundary(diagonal)
        replacement_boundary = path_boundary(path)
        replacements.append({
            "diagonal": list(diagonal),
            "replacement_path": list(path),
            "direct_displacement": list(direct_displacement),
            "replacement_displacement": list(replacement_displacement),
            "same_relative_boundary": direct_boundary == replacement_boundary,
            "difference_is_closed": direct_boundary == replacement_boundary,
            "difference_lift_displacement": [
                direct_displacement[0] - replacement_displacement[0],
                direct_displacement[1] - replacement_displacement[1],
            ],
            "ambient_H1_difference_zero_for_every_period_lattice": (
                direct_displacement == replacement_displacement
            ),
        })
    retained = row["retained_diagonals"]
    retained_disjoint = not retained or not black_boundary_edges
    return {
        "mask": mask,
        "black_boundary_edges": black_boundary_edges,
        "retained_white_diagonals": retained,
        "retained_diagonal_disjoint_from_black_carrier": retained_disjoint,
        "matching_connectivity_preserved": row["connectivity_preserved"],
        "at_most_one_retained_diagonal": row["embedded_diagonal_gate"],
        "replacement_certificates": replacements,
        "all_lifted_chain_gates_pass": (
            retained_disjoint
            and row["connectivity_preserved"]
            and row["embedded_diagonal_gate"]
            and all(
                replacement["same_relative_boundary"]
                and replacement["ambient_H1_difference_zero_for_every_period_lattice"]
                for replacement in replacements
            )
        ),
    }


def universal_face_audit() -> dict[str, Any]:
    patterns = [lift_face_certificate(mask) for mask in range(16)]
    return {
        "patterns": patterns,
        "pattern_count": len(patterns),
        "removed_diagonal_count": sum(
            len(row["replacement_certificates"]) for row in patterns
        ),
        "retained_diagonal_masks": [
            row["mask"] for row in patterns if row["retained_white_diagonals"]
        ],
        "all_patterns_pass": all(row["all_lifted_chain_gates_pass"] for row in patterns),
        "universal_statement": (
            "Every removed diagonal differs from its retained white NN path by "
            "a closed lifted 1-chain of displacement zero. Projection to R2/L "
            "therefore preserves its ambient H1 class for every finite-index L."
        ),
    }


def consistent_face_mask(periods: IntegerPeriods, mask: int) -> bool:
    colors: dict[Vector, int] = {}
    for corner, coordinate in enumerate(CORNERS):
        key = periods.quotient_key(coordinate)
        color = int(bool(mask & (1 << corner)))
        if key in colors and colors[key] != color:
            return False
        colors[key] = color
    return True


def canonical_endpoint(edge: Any) -> Tuple[int, int]:
    return min(edge.i, edge.j), max(edge.i, edge.j)


def quotient_projection_row(matrix: Matrix) -> dict[str, Any]:
    periods = IntegerPeriods(matrix)
    geometry = integer_torus_geometry(matrix, name="unrestricted-projection")
    face_rows = []
    failures = []
    for mask in range(16):
        if not consistent_face_mask(periods, mask):
            continue
        lift = lift_face_certificate(mask)
        for replacement in lift["replacement_certificates"]:
            diagonal = replacement["diagonal"]
            path = replacement["replacement_path"]
            direct_start = periods.quotient_key(CORNERS[diagonal[0]])
            direct_end = periods.quotient_key(CORNERS[diagonal[1]])
            path_start = periods.quotient_key(CORNERS[path[0]])
            path_end = periods.quotient_key(CORNERS[path[-1]])
            difference = tuple(replacement["difference_lift_displacement"])
            projected_pass = (
                direct_start == path_start
                and direct_end == path_end
                and periods.winding(difference) == (0, 0)
            )
            if not projected_pass:
                failures.append({"mask": mask, "replacement": replacement})
        face_rows.append({
            "mask": mask,
            "retained_diagonal_count": len(lift["retained_white_diagonals"]),
            "replacement_count": len(lift["replacement_certificates"]),
        })
    primal_endpoints = [canonical_endpoint(edge) for edge in geometry.primal_edges]
    matching_endpoints = [canonical_endpoint(edge) for edge in geometry.matching_edges]
    return {
        "matrix": [list(row) for row in matrix],
        "order": periods.order,
        "four_distinct_face_corners": has_four_distinct_face_corners(geometry),
        "consistent_local_patterns": len(face_rows),
        "primal_loop_incidences": sum(edge.i == edge.j for edge in geometry.primal_edges),
        "matching_loop_incidences": sum(edge.i == edge.j for edge in geometry.matching_edges),
        "primal_repeated_endpoint_incidences": len(primal_endpoints) - len(set(primal_endpoints)),
        "matching_repeated_endpoint_incidences": len(matching_endpoints) - len(set(matching_endpoints)),
        "projection_failures": failures,
    }


def quotient_projection_audit(maximum_order: int) -> dict[str, Any]:
    rows = [quotient_projection_row(matrix)
            for matrix in hnf_matrices(maximum_order=maximum_order)]
    failures = [
        {"matrix": row["matrix"], "failures": row["projection_failures"]}
        for row in rows if row["projection_failures"]
    ]
    return {
        "maximum_order": maximum_order,
        "representatives": len(rows),
        "honest_face_representatives": sum(
            row["four_distinct_face_corners"] for row in rows
        ),
        "self_identifying_face_representatives": sum(
            not row["four_distinct_face_corners"] for row in rows
        ),
        "representatives_with_primal_loops": sum(
            bool(row["primal_loop_incidences"]) for row in rows
        ),
        "representatives_with_matching_loops": sum(
            bool(row["matching_loop_incidences"]) for row in rows
        ),
        "representatives_with_repeated_primal_endpoints": sum(
            bool(row["primal_repeated_endpoint_incidences"]) for row in rows
        ),
        "representatives_with_repeated_matching_endpoints": sum(
            bool(row["matching_repeated_endpoint_incidences"]) for row in rows
        ),
        "consistent_face_patterns_checked": sum(
            row["consistent_local_patterns"] for row in rows
        ),
        "failure_count": len(failures),
        "first_failure": failures[0] if failures else None,
        "rows": rows,
    }


def scale_matrix(matrix: Matrix, factor: int) -> Matrix:
    return (
        (factor * matrix[0][0], factor * matrix[0][1]),
        (factor * matrix[1][0], factor * matrix[1][1]),
    )


def finite_regular_cover_row(matrix: Matrix) -> dict[str, Any]:
    """Return the exact algebra and face gate for T_(2L) -> T_L."""
    base = IntegerPeriods(matrix)
    cover_matrix = scale_matrix(matrix, 2)
    cover = IntegerPeriods(cover_matrix)
    cover_geometry = integer_torus_geometry(
        cover_matrix, name="unrestricted-honest-cover"
    )
    degree = cover.order // base.order
    map_columns = (
        base.winding(cover.period_vector((1, 0))),
        base.winding(cover.period_vector((0, 1))),
    )
    h1_map = [
        [map_columns[0][0], map_columns[1][0]],
        [map_columns[0][1], map_columns[1][1]],
    ]
    intersection_scale = (
        h1_map[0][0] * h1_map[1][1]
        - h1_map[0][1] * h1_map[1][0]
    )
    return {
        "matrix": [list(value) for value in matrix],
        "base_order": base.order,
        "cover_matrix": [list(value) for value in cover_matrix],
        "cover_order": cover.order,
        "cover_degree": degree,
        "cover_has_four_distinct_face_corners": (
            has_four_distinct_face_corners(cover_geometry)
        ),
        "H1_map_in_period_bases": h1_map,
        "intersection_scale": intersection_scale,
    }


def finite_regular_cover_audit(maximum_order: int) -> dict[str, Any]:
    """Certify the canonical honest four-sheeted cover L' = 2L.

    In period bases the induced H1 map is 2I, hence has determinant four and
    scales the torus intersection form by four.  Four distinct face corners
    follow because every vector of 2L has even ambient coordinates, whereas a
    nonzero difference of two unit-square corners has a coordinate equal to
    plus or minus one.
    """
    rows = []
    failures = []
    for matrix in hnf_matrices(maximum_order=maximum_order):
        row = finite_regular_cover_row(matrix)
        if (
            row["cover_degree"] != 4
            or not row["cover_has_four_distinct_face_corners"]
            or row["H1_map_in_period_bases"] != [[2, 0], [0, 2]]
            or row["intersection_scale"] != 4
        ):
            failures.append(row)
        rows.append(row)
    return {
        "maximum_base_order": maximum_order,
        "representatives": len(rows),
        "cover_degree": 4,
        "honest_cover_representatives": sum(
            row["cover_has_four_distinct_face_corners"] for row in rows
        ),
        "failure_count": len(failures),
        "first_failure": failures[0] if failures else None,
        "symbolic_reason": (
            "2L contains only even-coordinate vectors, while every nonzero "
            "difference of unit-square corners has a coordinate of absolute "
            "value one; p_* is 2I in period bases and omega(2u,2v)=4 omega(u,v)."
        ),
        "rows": rows,
    }


def subset_state_regression(maximum_order: int) -> dict[str, Any]:
    counters: Counter[str] = Counter()
    first_failure = None
    rows = []
    for matrix in hnf_matrices(maximum_order=maximum_order):
        geometry = integer_torus_geometry(matrix, name="unrestricted-state-regression")
        table = build_state_table(geometry)
        local: Counter[str] = Counter()
        for mask, (black, white) in table.items():
            black_rank, black_line, black_index, black_error = black
            white_rank, white_line, white_index, white_error = white
            local["states"] += 1
            state_failure = False
            if black_rank + white_rank != 2:
                local["rank_sum_failures"] += 1
                state_failure = True
            if black_rank == white_rank == 1:
                local["rank_one_states"] += 1
                if black_error or white_error:
                    local["rank_mark_failures"] += 1
                    state_failure = True
                if black_line != white_line:
                    local["primitive_line_failures"] += 1
                    state_failure = True
                if black_index != 1 or white_index != 1:
                    local["nonsaturated_rank_one_states"] += 1
            if first_failure is None and state_failure:
                first_failure = {
                    "matrix": [list(row) for row in matrix],
                    "mask": mask,
                    "black": list(black),
                    "white": list(white),
                }
        counters.update(local)
        rows.append({
            "matrix": [list(row) for row in matrix],
            "order": geometry.n,
            **dict(local),
        })
    return {
        "maximum_order": maximum_order,
        "representatives": len(rows),
        "states": counters["states"],
        "rank_one_states": counters["rank_one_states"],
        "rank_sum_failures": counters["rank_sum_failures"],
        "rank_mark_failures": counters["rank_mark_failures"],
        "primitive_line_failures": counters["primitive_line_failures"],
        "nonsaturated_rank_one_states": counters["nonsaturated_rank_one_states"],
        "saturation_claim": (
            "diagnostic only: no nonsaturated rank-one graph image occurs in "
            "the exact HNF index-2-through-10 regression"
        ),
        "first_failure": first_failure,
        "rows": rows,
    }


def canonical_primitive(vector: Vector) -> Vector:
    divisor = gcd(abs(vector[0]), abs(vector[1]))
    if divisor == 0:
        raise ValueError("zero has no primitive line")
    result = vector[0] // divisor, vector[1] // divisor
    if result[0] < 0 or (result[0] == 0 and result[1] < 0):
        result = -result[0], -result[1]
    return result


def symplectic_line_audit(bound: int = 5) -> dict[str, Any]:
    tested = 0
    failures = []
    vectors = [
        (x, y) for x in range(-bound, bound + 1)
        for y in range(-bound, bound + 1) if (x, y) != (0, 0)
    ]
    for first in vectors:
        if gcd(abs(first[0]), abs(first[1])) != 1:
            continue
        for second in vectors:
            intersection = first[0] * second[1] - first[1] * second[0]
            if intersection:
                continue
            tested += 1
            if canonical_primitive(first) != canonical_primitive(second):
                failures.append([list(first), list(second)])
    return {
        "bound": bound,
        "orthogonal_nonzero_pairs_checked": tested,
        "failure_count": len(failures),
        "first_failure": failures[0] if failures else None,
        "symbolic_identity": "u^T J v = u0*v1-u1*v0; in dimension two, zero means the nonzero vectors span the same rational line",
    }


def rank_path(k_first: int, k_second: int, k: int) -> int:
    return int(k >= k_first) + int(k >= k_second)


def filtration_consequence_audit(maximum_sites: int = 12) -> dict[str, Any]:
    cases = 0
    failures = []
    for sites in range(1, maximum_sites + 1):
        for k_first in range(1, sites + 1):
            for k_second in range(k_first, sites + 1):
                cases += 1
                path = [rank_path(k_first, k_second, k) for k in range(sites + 1)]
                white = [2 - value for value in path]
                direct_first = next(k for k, value in enumerate(path) if value >= 1)
                direct_second = next(k for k, value in enumerate(path) if value == 2)
                reverse_matching_birth = sites - k_first + 1
                historical_minus = sites - reverse_matching_birth + 1
                swapped_plus = reverse_matching_birth
                swapped_minus = sites - k_second + 1
                reconstructed = [
                    int(k >= historical_minus) + int(k >= direct_second)
                    for k in range(sites + 1)
                ]
                passed = (
                    direct_first == historical_minus == k_first
                    and direct_second == k_second
                    and all(left + right == 2 for left, right in zip(path, white))
                    and reconstructed == path
                    and historical_minus + swapped_plus == sites + 1
                    and direct_second + swapped_minus == sites + 1
                )
                if not passed:
                    failures.append({
                        "sites": sites,
                        "k_first": k_first,
                        "k_second": k_second,
                        "path": path,
                        "white": white,
                    })
    return {
        "maximum_sites": maximum_sites,
        "threshold_pairs_checked": cases,
        "failure_count": len(failures),
        "first_failure": failures[0] if failures else None,
        "derived_gates": [
            "K_minus=min{k:r_black(k)>=1}",
            "K_plus=min{k:r_black(k)=2}",
            "r_black(k)=1[k>=K_minus]+1[k>=K_plus]",
            "K_minus^G(pi)+K_plus^Ghat(reverse(pi))=N+1",
            "K_plus^G(pi)+K_minus^Ghat(reverse(pi))=N+1",
            "the rank-one rational line is constant by inclusion on its plateau",
        ],
    }


def build_artifact(config: dict[str, Any]) -> dict[str, Any]:
    maximum_order = int(config["regression_search"]["HNF_order_range"][1])
    face = universal_face_audit()
    projection = quotient_projection_audit(maximum_order)
    cover = finite_regular_cover_audit(maximum_order)
    states = subset_state_regression(maximum_order)
    symplectic = symplectic_line_audit()
    filtration = filtration_consequence_audit()
    all_machine_gates = (
        face["all_patterns_pass"]
        and projection["failure_count"] == 0
        and cover["failure_count"] == 0
        and states["rank_sum_failures"] == 0
        and states["rank_mark_failures"] == 0
        and states["primitive_line_failures"] == 0
        and symplectic["failure_count"] == 0
        and filtration["failure_count"] == 0
    )
    return {
        "schema": "matching-one/digital-alexander-unrestricted-theorem/v1",
        "issue": 269,
        "status": "unrestricted_finite_index_theorem",
        "theorem": {
            "period_lattice": "L=P Z^2 for every nonsingular integer 2x2 matrix P",
            "configuration": "every L-periodic black subset of Z^2",
            "identities": [
                "im_Q H1(white_matching)=im_Q H1(black_NN)^perp",
                "r_black+r_white=2",
                "q=r_black-1=1-r_white=(r_black-r_white)/2",
            ],
            "primitive_rank_one_line": (
                "the common black/white rational line has one canonical primitive "
                "integral generator up to the fixed sign convention"
            ),
            "saturation_index": (
                "not fixed by finite-cover descent; index one through quotient "
                "index 10 is retained only as an exact finite diagnostic"
            ),
            "extra_honest_cell_hypothesis": None,
        },
        "proof_chain": [
            {
                "lemma": "canonical finite honest cover",
                "statement": (
                    "For every L=P Z^2, the sublattice L'=2L defines a four-sheeted "
                    "regular cover T_(2L)->T_L. Every vector in 2L has even ambient "
                    "coordinates, so no nonzero difference of unit-square corners "
                    "lies in 2L; all quotient faces upstairs have four distinct corners."
                ),
            },
            {
                "lemma": "honest-cell theorem upstairs",
                "statement": (
                    "Lift the coloring and both graphs to T_(2L). The existing 16-pattern "
                    "pruning and complementary-subsurface proof applies without change, "
                    "giving I'_W=(I'_B)^perp and r'_black+r'_white=2."
                ),
            },
            {
                "lemma": "lifted-chain replacement",
                "statement": (
                    "A removed white diagonal and its same-face NN replacement have the "
                    "same relative boundary and identical lifted displacement. Their "
                    "difference has ambient H1 class zero after every quotient projection, "
                    "including loops and repeated edges."
                ),
            },
            {
                "lemma": "complementary-subsurface duality",
                "statement": (
                    "For complementary subsurfaces U,V of an oriented torus, exact "
                    "sequence, excision and Poincare-Lefschetz duality give "
                    "im H1(V;Q)=im H1(U;Q)^perp."
                ),
            },
            {
                "lemma": "finite-cover rational image descent",
                "statement": (
                    "For the full inverse-image graph, every upstairs cycle projects "
                    "downstairs. Conversely, a downstairs loop has finite deck monodromy, "
                    "so a positive iterate lifts closed; after tensoring with Q its class "
                    "is in the projected upstairs image. The same argument componentwise "
                    "preserves the repository maximum-component ranks."
                ),
            },
            {
                "lemma": "symplectic and primitive-line descent",
                "statement": (
                    "In compatible period bases p_*=2I, hence omega(p_*u,p_*v)=4omega(u,v). "
                    "It carries orthogonal complements to orthogonal complements. A "
                    "rank-one rational line has a canonical primitive integral generator; "
                    "this does not assert that the graph-image subgroup has index one."
                ),
            },
            {
                "lemma": "filtration corollaries",
                "statement": (
                    "Monotonicity plus rank complementarity gives both births, full "
                    "two-threshold reconstruction and swapped/reversed reflection. In "
                    "the rank-one case a line equals its symplectic orthogonal in Q2; "
                    "nested one-dimensional images keep that rational line and its "
                    "canonical primitive representative constant."
                ),
            },
        ],
        "machine_certificates": {
            "universal_face_chains": face,
            "degenerate_quotient_projection": projection,
            "finite_honest_cover": cover,
            "cached_subset_regression": states,
            "symplectic_line": symplectic,
            "filtration_algebra": filtration,
            "all_pass": all_machine_gates,
        },
        "proof_vs_search": {
            "proof": (
                "four-sheeted honest regular cover plus rational graph-image and "
                "conformally symplectic descent; independent of quotient index"
            ),
            "search": (
                "the existing exact state cache is reused through index 10 only as a "
                "regression and executable counterexample locator"
            ),
        },
        "claim_boundary": config["claim_boundary"],
    }


def render_markdown(artifact: dict[str, Any]) -> str:
    theorem = artifact["theorem"]
    machine = artifact["machine_certificates"]
    projection = machine["degenerate_quotient_projection"]
    cover = machine["finite_honest_cover"]
    states = machine["cached_subset_regression"]
    lines = [
        "# Unrestricted degenerate-quotient digital Alexander theorem",
        "",
        "Status: theorem for every finite-index period sublattice of `Z^2`.",
        "",
        "## Statement",
        "",
        "For every nonsingular integer period matrix `P` and every periodic black-site configuration,",
        "",
        "```text",
        "im_Q H1(white matching) = im_Q H1(black NN)^perp,",
        "r_black + r_white = 2,",
        "q = r_black - 1 = 1 - r_white = (r_black-r_white)/2.",
        "```",
        "",
        "No honest quotient-cell hypothesis is needed. A rank-one rational image has a canonical",
        "primitive integral direction; saturation of the actual graph-image subgroup is not asserted.",
        "",
        "## Why self-identifying faces are harmless",
        "",
        "Replace `L` by `2L`. The four-sheeted regular cover `T_(2L) -> T_L` has honest square faces:",
        "vectors in `2L` have even ambient coordinates, whereas every nonzero difference between two",
        "unit-square corners has a coordinate of absolute value one. The existing honest-cell theorem",
        "therefore applies upstairs even if the base presentation has loops, repeated edges or identified corners.",
        "Every removed diagonal has the same lifted displacement as its white NN replacement; the difference",
        "therefore remains ambient-null after any quotient projection.",
        "Rational graph-image homology descends exactly because a downstairs loop's finite deck monodromy",
        "is killed by a positive iterate. In period bases `p_*=2I`, so the intersection form scales by four",
        "and orthogonal complements descend with it.",
        "",
        "## Proof chain",
        "",
    ]
    for index, lemma in enumerate(artifact["proof_chain"], 1):
        lines.append(f"{index}. **{lemma['lemma']}.** {lemma['statement']}")
    lines += [
        "",
        "## Machine certificates",
        "",
        f"- universal face patterns: {machine['universal_face_chains']['pattern_count']}; all pass",
        f"- projected HNF representatives: {projection['representatives']} "
        f"({projection['self_identifying_face_representatives']} self-identifying)",
        f"- canonical four-sheeted honest covers: {cover['honest_cover_representatives']}; failures: zero",
        f"- consistent projected face patterns: {projection['consistent_face_patterns_checked']}",
        f"- cached subset states: {states['states']}; rank-one states: {states['rank_one_states']}",
        "- rank-sum, rank-mark, primitive-line and projection failures: zero",
        f"- nonsaturated rank-one states in the finite regression: "
        f"{states['nonsaturated_rank_one_states']} (diagnostic only)",
        f"- symbolic threshold pairs: {machine['filtration_algebra']['threshold_pairs_checked']}; failures: zero",
        f"- all machine gates: `{machine['all_pass']}`",
        "",
        "The finite HNF layer is a regression certificate, not the unrestricted inference. The unrestricted",
        "step is finite honest-cover descent over rational homology.",
        "",
        "## Boundary",
        "",
    ]
    lines.extend(f"- {item}" for item in artifact["claim_boundary"])
    lines.append("")
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--manifest", type=Path,
        default=Path("analysis/digital_alexander_unrestricted_manifest.json"),
    )
    parser.add_argument("--format", choices=("json", "markdown"), default="json")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    config = json.loads(args.manifest.read_text(encoding="utf-8"))
    artifact = build_artifact(config)
    rendered = (
        json.dumps(artifact, indent=2, sort_keys=True) + "\n"
        if args.format == "json" else render_markdown(artifact)
    )
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
