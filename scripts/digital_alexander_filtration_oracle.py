#!/usr/bin/env python3
"""Exact tiny-permutation oracle for the two essential homology births."""

from __future__ import annotations

import argparse
import itertools
import json
from collections import Counter
from math import factorial, gcd
from pathlib import Path
from typing import Any, Iterable, Optional, Sequence, Tuple

from integer_period_torus import (
    IntegerComponentHomology,
    IntegerTorusGeometry,
    axis_integer_torus,
    classify_configuration,
    gaussian_integer_torus,
)


Vector = Tuple[int, int]


def _canonical_primitive(vector: Vector) -> Vector:
    divisor = gcd(abs(vector[0]), abs(vector[1]))
    if divisor == 0:
        raise ValueError("zero vector has no projective line")
    result = vector[0] // divisor, vector[1] // divisor
    if result[0] < 0 or (result[0] == 0 and result[1] < 0):
        result = -result[0], -result[1]
    return result


def _rank_one_mark(components: Iterable[IntegerComponentHomology]) -> Tuple[Vector, int]:
    generators = [
        vector
        for component in components
        for vector in component.generators
        if vector != (0, 0)
    ]
    lines = {_canonical_primitive(vector) for vector in generators}
    if len(lines) != 1:
        raise AssertionError("rank-one state did not have one rational winding line")
    line = next(iter(lines))
    coefficients = []
    for x, y in generators:
        if line[0] != 0:
            if x % line[0] or y != (x // line[0]) * line[1]:
                raise AssertionError("generator is not integral on its primitive line")
            coefficients.append(x // line[0])
        else:
            if y % line[1] or x != 0:
                raise AssertionError("generator is not integral on its primitive line")
            coefficients.append(y // line[1])
    index = 0
    for coefficient in coefficients:
        index = gcd(index, abs(coefficient))
    if index < 1:
        raise AssertionError("rank-one subgroup has no positive saturation index")
    return line, index


def rank_mark(
    geometry: IntegerTorusGeometry,
    active: Sequence[bool],
    *,
    matching: bool,
) -> Tuple[int, Optional[Vector], Optional[int]]:
    channels, components = classify_configuration(geometry, active, matching=matching)
    rank = channels.max_rank
    if rank != 1:
        return rank, None, None
    line, index = _rank_one_mark(components)
    return rank, line, index


def first_rank_two(
    geometry: IntegerTorusGeometry,
    order: Sequence[int],
    *,
    matching: bool,
) -> int:
    active = [False] * geometry.n
    for rank, vertex in enumerate(order, start=1):
        active[vertex] = True
        if rank_mark(geometry, active, matching=matching)[0] == 2:
            return rank
    raise AssertionError("the fully occupied graph never reached ambient rank two")


def historical_thresholds(
    geometry: IntegerTorusGeometry,
    permutation: Sequence[int],
    *,
    forward_matching: bool,
    complement_matching: bool,
) -> Tuple[int, int]:
    k_plus = first_rank_two(geometry, permutation, matching=forward_matching)
    reverse_birth = first_rank_two(
        geometry,
        tuple(reversed(permutation)),
        matching=complement_matching,
    )
    return geometry.n - reverse_birth + 1, k_plus


def direct_filtration(
    geometry: IntegerTorusGeometry,
    permutation: Sequence[int],
) -> list[dict[str, Any]]:
    active = [False] * geometry.n
    rows = []
    for k in range(geometry.n + 1):
        black_rank, black_line, black_index = rank_mark(
            geometry, active, matching=False
        )
        white_rank, white_line, white_index = rank_mark(
            geometry, tuple(not value for value in active), matching=True
        )
        rows.append(
            {
                "k": k,
                "rank_black": black_rank,
                "rank_white": white_rank,
                "black_line": black_line,
                "white_line": white_line,
                "black_index": black_index,
                "white_index": white_index,
            }
        )
        if k < geometry.n:
            active[permutation[k]] = True
    return rows


def analyze_permutation(
    geometry: IntegerTorusGeometry,
    permutation: Sequence[int],
) -> dict[str, Any]:
    rows = direct_filtration(geometry, permutation)
    k_first = next(row["k"] for row in rows if row["rank_black"] >= 1)
    k_second = next(row["k"] for row in rows if row["rank_black"] == 2)
    k_minus, k_plus = historical_thresholds(
        geometry,
        permutation,
        forward_matching=False,
        complement_matching=True,
    )
    reverse_swapped_minus, reverse_swapped_plus = historical_thresholds(
        geometry,
        tuple(reversed(permutation)),
        forward_matching=True,
        complement_matching=False,
    )

    rank_sum_failures = []
    reconstruction_failures = []
    line_failures = []
    plateau_lines = []
    index_pairs = []
    for row in rows:
        k = row["k"]
        if row["rank_black"] + row["rank_white"] != 2:
            rank_sum_failures.append(k)
        reconstructed = int(k >= k_minus) + int(k >= k_plus)
        if row["rank_black"] != reconstructed:
            reconstruction_failures.append(k)
        if k_minus <= k < k_plus:
            black_line = row["black_line"]
            white_line = row["white_line"]
            if (
                black_line is None
                or white_line is None
                or black_line != white_line
                or black_line[0] * white_line[1] - black_line[1] * white_line[0] != 0
            ):
                line_failures.append(k)
            else:
                plateau_lines.append(black_line)
                index_pairs.append((row["black_index"], row["white_index"]))

    if len(set(plateau_lines)) > 1:
        line_failures.append(-1)
    return {
        "permutation": list(permutation),
        "k_first_direct": k_first,
        "k_second_direct": k_second,
        "k_minus_historical": k_minus,
        "k_plus_historical": k_plus,
        "endpoint_reflection_residuals": [
            k_minus + reverse_swapped_plus - (geometry.n + 1),
            k_plus + reverse_swapped_minus - (geometry.n + 1),
        ],
        "rank_sum_failures": rank_sum_failures,
        "reconstruction_failures": reconstruction_failures,
        "line_failures": line_failures,
        "plateau_length": k_plus - k_minus,
        "plateau_line": list(plateau_lines[0]) if plateau_lines else None,
        "index_pairs": [list(pair) for pair in index_pairs],
        "index_evolves": len(set(index_pairs)) > 1,
    }


def geometry_specs() -> list[Tuple[str, IntegerTorusGeometry]]:
    return [
        ("axis-L2", axis_integer_torus(2)),
        ("gaussian-2-1", gaussian_integer_torus(2, 1)),
    ]


def summarize_geometry(name: str, geometry: IntegerTorusGeometry) -> dict[str, Any]:
    birth_failures = []
    reflection_failures = []
    rank_sum_failures = []
    reconstruction_failures = []
    line_failures = []
    line_counts: Counter[str] = Counter()
    plateau_steps = 0
    positive_plateaus = 0
    index_evolution_count = 0
    maximum_index = 0
    representative = None

    for permutation in itertools.permutations(range(geometry.n)):
        row = analyze_permutation(geometry, permutation)
        payload = {"permutation": row["permutation"]}
        if (
            row["k_first_direct"] != row["k_minus_historical"]
            or row["k_second_direct"] != row["k_plus_historical"]
        ):
            birth_failures.append({**payload, "row": row})
        if any(row["endpoint_reflection_residuals"]):
            reflection_failures.append({**payload, "residuals": row["endpoint_reflection_residuals"]})
        if row["rank_sum_failures"]:
            rank_sum_failures.append({**payload, "steps": row["rank_sum_failures"]})
        if row["reconstruction_failures"]:
            reconstruction_failures.append({**payload, "steps": row["reconstruction_failures"]})
        if row["line_failures"]:
            line_failures.append({**payload, "steps": row["line_failures"]})
        plateau_steps += row["plateau_length"]
        if row["plateau_length"]:
            positive_plateaus += 1
            line_counts[",".join(str(value) for value in row["plateau_line"])] += 1
            if representative is None:
                representative = row
        if row["index_evolves"]:
            index_evolution_count += 1
        for pair in row["index_pairs"]:
            maximum_index = max(maximum_index, pair[0], pair[1])

    permutation_count = factorial(geometry.n)
    for failures in (
        birth_failures,
        reflection_failures,
        rank_sum_failures,
        reconstruction_failures,
        line_failures,
    ):
        del failures[8:]
    return {
        "id": name,
        "N": geometry.n,
        "period_matrix": [list(row) for row in geometry.periods.matrix],
        "permutations": permutation_count,
        "positive_rank_one_plateaus": positive_plateaus,
        "rank_one_plateau_steps": plateau_steps,
        "plateau_line_counts": dict(sorted(line_counts.items())),
        "permutations_with_index_evolution": index_evolution_count,
        "maximum_saturation_index": maximum_index,
        "birth_failure_count_capped": len(birth_failures),
        "reflection_failure_count_capped": len(reflection_failures),
        "rank_sum_failure_count_capped": len(rank_sum_failures),
        "reconstruction_failure_count_capped": len(reconstruction_failures),
        "line_failure_count_capped": len(line_failures),
        "birth_counterexamples": birth_failures,
        "reflection_counterexamples": reflection_failures,
        "rank_sum_counterexamples": rank_sum_failures,
        "reconstruction_counterexamples": reconstruction_failures,
        "line_counterexamples": line_failures,
        "representative_positive_plateau": representative,
    }


def build_artifact() -> dict[str, Any]:
    geometries = [summarize_geometry(name, geometry) for name, geometry in geometry_specs()]
    total = sum(row["permutations"] for row in geometries)
    assert total == 144
    assert all(not row["birth_counterexamples"] for row in geometries)
    assert all(not row["reflection_counterexamples"] for row in geometries)
    assert all(not row["rank_sum_counterexamples"] for row in geometries)
    assert all(not row["reconstruction_counterexamples"] for row in geometries)
    assert all(not row["line_counterexamples"] for row in geometries)
    return {
        "schema": "matching-one/digital-alexander-filtration-oracle/v1",
        "issue": 269,
        "status": "tiny_exact_filtration_semantics",
        "geometries": geometries,
        "totals": {
            "permutations": total,
            "birth_failures": sum(row["birth_failure_count_capped"] for row in geometries),
            "reflection_failures": sum(row["reflection_failure_count_capped"] for row in geometries),
            "rank_reconstruction_failures": sum(
                row["reconstruction_failure_count_capped"] for row in geometries
            ),
            "line_failures": sum(row["line_failure_count_capped"] for row in geometries),
        },
        "exact_identities": [
            "K_minus=min{k:R_k>=1}",
            "K_plus=min{k:R_k=2}",
            "R_k=1[k>=K_minus]+1[k>=K_plus]",
            "K_minus^G(pi)+K_plus^Ghat(reverse(pi))=N+1",
            "K_plus^G(pi)+K_minus^Ghat(reverse(pi))=N+1",
        ],
        "claim_boundary": {
            "proved": "all permutations of the two declared tiny geometries",
            "not_proved": (
                "all regular tori, degenerate short-period quotients, production-stream lineage, "
                "or a continuum interpretation"
            ),
        },
    }


def render_markdown(artifact: dict[str, Any]) -> str:
    lines = [
        "# Essential-H1 filtration oracle",
        "",
        "Every site permutation is exhausted on the two declared tiny tori. Direct rank paths are",
        "recomputed independently of the historical reverse-sweep threshold convention.",
        "",
        "| geometry | N | permutations | positive rank-one plateaus | plateau steps | index evolution | max index |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for row in artifact["geometries"]:
        lines.append(
            "| %s | %d | %d | %d | %d | %d | %d |"
            % (
                row["id"], row["N"], row["permutations"],
                row["positive_rank_one_plateaus"], row["rank_one_plateau_steps"],
                row["permutations_with_index_evolution"], row["maximum_saturation_index"],
            )
        )
    lines.extend(
        [
            "",
            "Across `%d` permutations there are zero birth, rank-reconstruction, endpoint-reflection,"
            % artifact["totals"]["permutations"],
            "or rational winding-line failures. Black and white share one constant saturated primitive",
            "line throughout every nonempty rank-one plateau. Integral saturation indices are recorded",
            "separately and are not inferred from rational-line equality.",
            "",
            "## Boundary",
            "",
            "This is a tiny exact semantics gate. It does not modify production streams or prove the",
            "identities for every regular or degenerate quotient, and it assigns no CFT interpretation.",
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
