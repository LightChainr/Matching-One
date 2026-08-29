#!/usr/bin/env python3
"""Exhaust the short-period HNF frontier for digital Alexander filtrations."""

from __future__ import annotations

import argparse
import itertools
import json
from collections import Counter
from math import factorial
from pathlib import Path
from typing import Any, Optional, Sequence, Tuple

from digital_alexander_filtration_oracle import rank_mark
from integer_period_torus import (
    IntegerTorusGeometry,
    classify_configuration,
    determinant,
    integer_torus_geometry,
)


Matrix = Tuple[Tuple[int, int], Tuple[int, int]]


def hnf_matrices(minimum_order: int = 2, maximum_order: int = 6) -> list[Matrix]:
    if minimum_order < 1 or maximum_order < minimum_order:
        raise ValueError("invalid order range")
    matrices = []
    for order in range(minimum_order, maximum_order + 1):
        for first_diagonal in range(1, order + 1):
            if order % first_diagonal:
                continue
            second_diagonal = order // first_diagonal
            for shear in range(first_diagonal):
                matrices.append(((first_diagonal, shear), (0, second_diagonal)))
    return matrices


def has_four_distinct_face_corners(geometry: IntegerTorusGeometry) -> bool:
    corners = ((0, 0), (1, 0), (0, 1), (1, 1))
    return len({geometry.periods.quotient_key(corner) for corner in corners}) == 4


def _rank_only(
    geometry: IntegerTorusGeometry,
    active: Sequence[bool],
    *,
    matching: bool,
) -> int:
    return classify_configuration(geometry, active, matching=matching)[0].max_rank


def _rank_line_safe(
    geometry: IntegerTorusGeometry,
    active: Sequence[bool],
    *,
    matching: bool,
) -> Tuple[int, Optional[Tuple[int, int]], Optional[int], Optional[str]]:
    rank = _rank_only(geometry, active, matching=matching)
    if rank != 1:
        return rank, None, None, None
    try:
        _, line, index = rank_mark(geometry, active, matching=matching)
        return rank, line, index, None
    except AssertionError as error:
        return rank, None, None, str(error)


def _first_rank_two(
    geometry: IntegerTorusGeometry,
    order: Sequence[int],
    *,
    matching: bool,
) -> int:
    active = [False] * geometry.n
    for rank, vertex in enumerate(order, start=1):
        active[vertex] = True
        if _rank_only(geometry, active, matching=matching) == 2:
            return rank
    raise AssertionError("fully occupied quotient never reached ambient rank two")


def _thresholds(
    geometry: IntegerTorusGeometry,
    permutation: Sequence[int],
    *,
    forward_matching: bool,
    complement_matching: bool,
) -> Tuple[int, int]:
    k_plus = _first_rank_two(geometry, permutation, matching=forward_matching)
    reverse_birth = _first_rank_two(
        geometry,
        tuple(reversed(permutation)),
        matching=complement_matching,
    )
    return geometry.n - reverse_birth + 1, k_plus


def analyze_permutation(
    geometry: IntegerTorusGeometry,
    permutation: Sequence[int],
) -> dict[str, Any]:
    active = [False] * geometry.n
    states = []
    for k in range(geometry.n + 1):
        black = _rank_line_safe(geometry, active, matching=False)
        white = _rank_line_safe(
            geometry,
            tuple(not value for value in active),
            matching=True,
        )
        states.append({"k": k, "black": black, "white": white})
        if k < geometry.n:
            active[permutation[k]] = True

    k_first = next(state["k"] for state in states if state["black"][0] >= 1)
    k_second = next(state["k"] for state in states if state["black"][0] == 2)
    k_minus, k_plus = _thresholds(
        geometry,
        permutation,
        forward_matching=False,
        complement_matching=True,
    )
    swapped_minus, swapped_plus = _thresholds(
        geometry,
        tuple(reversed(permutation)),
        forward_matching=True,
        complement_matching=False,
    )

    rank_sum_steps = []
    reconstruction_steps = []
    line_steps = []
    plateau_lines = []
    index_pairs = []
    for state in states:
        k = state["k"]
        black_rank, black_line, black_index, black_error = state["black"]
        white_rank, white_line, white_index, white_error = state["white"]
        if black_rank + white_rank != 2:
            rank_sum_steps.append(k)
        if black_rank != int(k >= k_minus) + int(k >= k_plus):
            reconstruction_steps.append(k)
        if k_minus <= k < k_plus:
            if (
                black_error
                or white_error
                or black_line is None
                or white_line is None
                or black_line != white_line
            ):
                line_steps.append(k)
            else:
                plateau_lines.append(black_line)
                index_pairs.append((black_index, white_index))
    if len(set(plateau_lines)) > 1:
        line_steps.append(-1)

    return {
        "permutation": list(permutation),
        "birth_residuals": [k_first - k_minus, k_second - k_plus],
        "reflection_residuals": [
            k_minus + swapped_plus - (geometry.n + 1),
            k_plus + swapped_minus - (geometry.n + 1),
        ],
        "rank_sum_steps": rank_sum_steps,
        "reconstruction_steps": reconstruction_steps,
        "line_steps": line_steps,
        "plateau_length": k_plus - k_minus,
        "plateau_line": list(plateau_lines[0]) if plateau_lines else None,
        "index_pairs": [list(pair) for pair in index_pairs],
        "index_evolves": len(set(index_pairs)) > 1,
    }


def summarize_matrix(matrix: Matrix) -> dict[str, Any]:
    geometry = integer_torus_geometry(matrix, name="hnf-frontier")
    counters: Counter[str] = Counter()
    line_counts: Counter[str] = Counter()
    examples = {gate: [] for gate in ("birth", "reflection", "rank_sum", "reconstruction", "line")}
    plateau_steps = 0
    maximum_index = 0
    for permutation in itertools.permutations(range(geometry.n)):
        row = analyze_permutation(geometry, permutation)
        failures = {
            "birth": any(row["birth_residuals"]),
            "reflection": any(row["reflection_residuals"]),
            "rank_sum": bool(row["rank_sum_steps"]),
            "reconstruction": bool(row["reconstruction_steps"]),
            "line": bool(row["line_steps"]),
        }
        for gate, failed in failures.items():
            if failed:
                counters[gate + "_failures"] += 1
                if len(examples[gate]) < 4:
                    examples[gate].append(row)
        plateau_steps += row["plateau_length"]
        if row["plateau_line"] is not None:
            line_counts[",".join(str(value) for value in row["plateau_line"])] += 1
        if row["index_evolves"]:
            counters["index_evolution"] += 1
        for pair in row["index_pairs"]:
            maximum_index = max(maximum_index, pair[0], pair[1])
    return {
        "matrix": [list(part) for part in matrix],
        "order": geometry.n,
        "four_distinct_face_corners": has_four_distinct_face_corners(geometry),
        "permutations": factorial(geometry.n),
        "rank_one_plateau_steps": plateau_steps,
        "plateau_line_counts": dict(sorted(line_counts.items())),
        "maximum_saturation_index": maximum_index,
        "permutations_with_index_evolution": counters["index_evolution"],
        "failure_counts": {
            gate: counters[gate + "_failures"]
            for gate in ("birth", "reflection", "rank_sum", "reconstruction", "line")
        },
        "counterexamples": examples,
    }


def _first_counterexample(rows: Sequence[dict[str, Any]], gate: str) -> Optional[dict[str, Any]]:
    for row in rows:
        examples = row["counterexamples"][gate]
        if examples:
            return {
                "matrix": row["matrix"],
                "order": row["order"],
                "four_distinct_face_corners": row["four_distinct_face_corners"],
                "example": examples[0],
            }
    return None


def build_artifact() -> dict[str, Any]:
    matrices = hnf_matrices()
    rows = [summarize_matrix(matrix) for matrix in matrices]
    total_paths = sum(row["permutations"] for row in rows)
    total_failures = {
        gate: sum(row["failure_counts"][gate] for row in rows)
        for gate in ("birth", "reflection", "rank_sum", "reconstruction", "line")
    }
    first_counterexamples = {
        gate: _first_counterexample(rows, gate)
        for gate in ("birth", "reflection", "rank_sum", "reconstruction", "line")
    }
    honest = [row for row in rows if row["four_distinct_face_corners"]]
    degenerate = [row for row in rows if not row["four_distinct_face_corners"]]
    assert len(matrices) == 32
    assert total_paths == 9558
    assert all(abs(determinant(tuple(tuple(part) for part in row["matrix"]))) == row["order"] for row in rows)
    return {
        "schema": "matching-one/digital-alexander-quotient-frontier/v1",
        "issue": 269,
        "status": (
            "counterexample_found"
            if any(total_failures.values())
            else "no_counterexample_through_index_6"
        ),
        "order_range": [2, 6],
        "HNF_representatives": len(rows),
        "honest_face_representatives": len(honest),
        "self_identifying_face_representatives": len(degenerate),
        "filtration_paths": total_paths,
        "total_failure_counts": total_failures,
        "first_counterexamples": first_counterexamples,
        "geometries": rows,
        "claim_boundary": {
            "proved": "all permutations of every HNF quotient of index 2 through 6",
            "not_proved": "index 7 and above or an unrestricted degenerate-quotient theorem",
        },
    }


def render_markdown(artifact: dict[str, Any]) -> str:
    lines = [
        "# Short-period digital Alexander quotient frontier",
        "",
        "Every permutation of every two-dimensional HNF quotient of index 2 through 6 is exhausted.",
        "",
        "- HNF representatives: `%d`;" % artifact["HNF_representatives"],
        "- honest four-corner representatives: `%d`;" % artifact["honest_face_representatives"],
        "- self-identifying face representatives: `%d`;" % artifact["self_identifying_face_representatives"],
        "- complete filtration paths: `%d`." % artifact["filtration_paths"],
        "",
        "## Failure census",
        "",
        "| gate | failing paths | first counterexample |",
        "|---|---:|---|",
    ]
    for gate, count in artifact["total_failure_counts"].items():
        example = artifact["first_counterexamples"][gate]
        label = "none" if example is None else "order %d, matrix `%s`" % (example["order"], example["matrix"])
        lines.append("| %s | %d | %s |" % (gate, count, label))
    lines.extend(
        [
            "",
            "The face-corner classification is reported explicitly so a failure on a self-identifying",
            "quotient is not confused with the already proved regular square-cell theorem.",
            "",
            "## Boundary",
            "",
            "The result is exhaustive only through index 6. It does not prove an unrestricted theorem",
            "for all degenerate quotients or alter production data semantics.",
            "",
        ]
    )
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--format", choices=("json", "markdown"), default="json")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    artifact = build_artifact()
    rendered = (
        json.dumps(artifact, indent=2, sort_keys=True) + "\n"
        if args.format == "json"
        else render_markdown(artifact)
    )
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
