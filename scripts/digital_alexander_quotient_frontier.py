#!/usr/bin/env python3
"""Exhaust the short-period HNF frontier for digital Alexander filtrations."""

from __future__ import annotations

import argparse
import itertools
import json
from collections import Counter
from functools import lru_cache
from math import factorial
from pathlib import Path
from typing import Any, Hashable, Optional, Sequence, Tuple

from digital_alexander_filtration_oracle import rank_mark
from integer_period_torus import (
    IntegerTorusGeometry,
    classify_configuration,
    determinant,
    integer_torus_geometry,
)


Matrix = Tuple[Tuple[int, int], Tuple[int, int]]


def hnf_matrices(minimum_order: int = 2, maximum_order: int = 10) -> list[Matrix]:
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


def build_state_table(
    geometry: IntegerTorusGeometry,
) -> dict[int, Tuple[tuple[Any, ...], tuple[Any, ...]]]:
    """Cache exact primal/dual marks for every occupied-site subset."""

    table = {}
    for mask in range(1 << geometry.n):
        active = tuple(bool(mask & (1 << vertex)) for vertex in range(geometry.n))
        table[mask] = (
            _rank_line_safe(geometry, active, matching=False),
            _rank_line_safe(
                geometry,
                tuple(not value for value in active),
                matching=True,
            ),
        )
    return table


def _cached_rank(
    state_table: dict[int, Tuple[tuple[Any, ...], tuple[Any, ...]]],
    full_mask: int,
    mask: int,
    *,
    matching: bool,
) -> int:
    if matching:
        return state_table[full_mask ^ mask][1][0]
    return state_table[mask][0][0]


def _first_rank_two(
    geometry: IntegerTorusGeometry,
    order: Sequence[int],
    *,
    matching: bool,
    state_table: Optional[dict[int, Tuple[tuple[Any, ...], tuple[Any, ...]]]] = None,
) -> int:
    if state_table is not None:
        mask = 0
        full_mask = (1 << geometry.n) - 1
        for rank, vertex in enumerate(order, start=1):
            mask |= 1 << vertex
            if _cached_rank(
                state_table,
                full_mask,
                mask,
                matching=matching,
            ) == 2:
                return rank
        raise AssertionError("fully occupied quotient never reached ambient rank two")

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
    state_table: Optional[dict[int, Tuple[tuple[Any, ...], tuple[Any, ...]]]] = None,
) -> Tuple[int, int]:
    k_plus = _first_rank_two(
        geometry,
        permutation,
        matching=forward_matching,
        state_table=state_table,
    )
    reverse_birth = _first_rank_two(
        geometry,
        tuple(reversed(permutation)),
        matching=complement_matching,
        state_table=state_table,
    )
    return geometry.n - reverse_birth + 1, k_plus


def analyze_permutation(
    geometry: IntegerTorusGeometry,
    permutation: Sequence[int],
    *,
    state_table: Optional[dict[int, Tuple[tuple[Any, ...], tuple[Any, ...]]]] = None,
) -> dict[str, Any]:
    states = []
    if state_table is None:
        active = [False] * geometry.n
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
    else:
        mask = 0
        for k in range(geometry.n + 1):
            black, white = state_table[mask]
            states.append({"k": k, "black": black, "white": white})
            if k < geometry.n:
                mask |= 1 << permutation[k]

    k_first = next(state["k"] for state in states if state["black"][0] >= 1)
    k_second = next(state["k"] for state in states if state["black"][0] == 2)
    k_minus, k_plus = _thresholds(
        geometry,
        permutation,
        forward_matching=False,
        complement_matching=True,
        state_table=state_table,
    )
    swapped_minus, swapped_plus = _thresholds(
        geometry,
        tuple(reversed(permutation)),
        forward_matching=True,
        complement_matching=False,
        state_table=state_table,
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


def summarize_matrix_by_permutation(matrix: Matrix) -> dict[str, Any]:
    """Reference implementation that materializes every maximal chain."""

    geometry = integer_torus_geometry(matrix, name="hnf-frontier")
    state_table = build_state_table(geometry)
    counters: Counter[str] = Counter()
    line_counts: Counter[str] = Counter()
    examples = {gate: [] for gate in ("birth", "reflection", "rank_sum", "reconstruction", "line")}
    plateau_steps = 0
    maximum_index = 0
    for permutation in itertools.permutations(range(geometry.n)):
        row = analyze_permutation(geometry, permutation, state_table=state_table)
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
        "cached_subsets": len(state_table),
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


def _prefix_path_count(n: int, mask: int) -> int:
    """Number of permutations whose prefix set is exactly ``mask``."""

    size = mask.bit_count()
    return factorial(size) * factorial(n - size)


def _count_paths_avoiding_nodes(n: int, bad_masks: set[int]) -> int:
    """Count maximal Boolean-lattice chains that avoid every bad node."""

    full_mask = (1 << n) - 1
    counts = [0] * (1 << n)
    if 0 in bad_masks:
        return 0
    counts[0] = 1
    for mask in range(1 << n):
        if not counts[mask]:
            continue
        remaining = full_mask ^ mask
        while remaining:
            bit = remaining & -remaining
            remaining ^= bit
            child = mask | bit
            if child not in bad_masks:
                counts[child] += counts[mask]
    return counts[full_mask]


def _first_bad_node_permutation(n: int, bad_masks: set[int]) -> Optional[Tuple[int, ...]]:
    """Return the lexicographically first maximal chain visiting a bad node."""

    full_mask = (1 << n) - 1

    @lru_cache(maxsize=None)
    def can_fail(mask: int, failed: bool) -> bool:
        failed = failed or mask in bad_masks
        if mask == full_mask:
            return failed
        return any(
            can_fail(mask | (1 << vertex), failed)
            for vertex in range(n)
            if not mask & (1 << vertex)
        )

    if not can_fail(0, False):
        return None
    mask = 0
    failed = 0 in bad_masks
    permutation = []
    while mask != full_mask:
        for vertex in range(n):
            bit = 1 << vertex
            if not mask & bit and can_fail(mask | bit, failed):
                permutation.append(vertex)
                mask |= bit
                failed = failed or mask in bad_masks
                break
        else:  # pragma: no cover - guarded by can_fail
            raise AssertionError("failed to materialize a bad maximal chain")
    return tuple(permutation)


def _marked_path_summary(
    n: int,
    marks: Sequence[Optional[Hashable]],
    invalid_masks: set[int],
) -> Tuple[int, Counter[str]]:
    """Count paths with invalid/changing marks and their first valid mark.

    The first component reproduces the permutation oracle's path-level failure
    semantics.  The counter reproduces its first-valid-plateau-mark census.
    """

    full_mask = (1 << n) - 1
    states: list[dict[Tuple[Optional[Hashable], bool], int]] = [
        {} for _ in range(1 << n)
    ]
    first = marks[0]
    states[0][(first, 0 in invalid_masks)] = 1
    for mask in range(1 << n):
        if not states[mask]:
            continue
        remaining = full_mask ^ mask
        while remaining:
            bit = remaining & -remaining
            remaining ^= bit
            child = mask | bit
            mark = marks[child]
            for (first_mark, failed), count in states[mask].items():
                next_first = first_mark if first_mark is not None else mark
                next_failed = failed or child in invalid_masks
                if first_mark is not None and mark is not None and mark != first_mark:
                    next_failed = True
                key = (next_first, next_failed)
                states[child][key] = states[child].get(key, 0) + count

    failures = 0
    first_mark_counts: Counter[str] = Counter()
    for (first_mark, failed), count in states[full_mask].items():
        if failed:
            failures += count
        if first_mark is not None:
            if isinstance(first_mark, tuple):
                label = ",".join(str(value) for value in first_mark)
            else:
                label = str(first_mark)
            first_mark_counts[label] += count
    return failures, first_mark_counts


def _first_mark_failure_permutation(
    n: int,
    marks: Sequence[Optional[Hashable]],
    invalid_masks: set[int],
) -> Optional[Tuple[int, ...]]:
    """Return the lexicographically first chain with an invalid/changing mark."""

    full_mask = (1 << n) - 1

    def advance(
        mask: int,
        first_mark: Optional[Hashable],
        failed: bool,
    ) -> Tuple[Optional[Hashable], bool]:
        mark = marks[mask]
        next_first = first_mark if first_mark is not None else mark
        next_failed = failed or mask in invalid_masks
        if first_mark is not None and mark is not None and mark != first_mark:
            next_failed = True
        return next_first, next_failed

    initial_first, initial_failed = advance(0, None, False)

    @lru_cache(maxsize=None)
    def can_fail(mask: int, first_mark: Optional[Hashable], failed: bool) -> bool:
        if mask == full_mask:
            return failed
        for vertex in range(n):
            bit = 1 << vertex
            if mask & bit:
                continue
            child = mask | bit
            next_first, next_failed = advance(child, first_mark, failed)
            if can_fail(child, next_first, next_failed):
                return True
        return False

    if not can_fail(0, initial_first, initial_failed):
        return None
    mask = 0
    first_mark = initial_first
    failed = initial_failed
    permutation = []
    while mask != full_mask:
        for vertex in range(n):
            bit = 1 << vertex
            if mask & bit:
                continue
            child = mask | bit
            next_first, next_failed = advance(child, first_mark, failed)
            if can_fail(child, next_first, next_failed):
                permutation.append(vertex)
                mask = child
                first_mark = next_first
                failed = next_failed
                break
        else:  # pragma: no cover - guarded by can_fail
            raise AssertionError("failed to materialize a bad marked chain")
    return tuple(permutation)


def summarize_matrix(matrix: Matrix) -> dict[str, Any]:
    """Exhaust every filtration via exact maximal-chain counts on subsets."""

    geometry = integer_torus_geometry(matrix, name="hnf-frontier")
    state_table = build_state_table(geometry)
    n = geometry.n
    total_paths = factorial(n)
    masks = range(1 << n)

    birth_bad = set()
    rank_sum_bad = set()
    reconstruction_bad = set()
    line_bad = set()
    plateau_lines: list[Optional[Tuple[int, int]]] = [None] * (1 << n)
    index_pairs: list[Optional[Tuple[int, int]]] = [None] * (1 << n)
    plateau_steps = 0
    maximum_index = 0
    for mask in masks:
        black, white = state_table[mask]
        black_rank, black_line, black_index, black_error = black
        white_rank, white_line, white_index, white_error = white
        if (black_rank >= 1) != (white_rank < 2):
            birth_bad.add(mask)
        if black_rank + white_rank != 2:
            rank_sum_bad.add(mask)
        if black_rank != int(white_rank < 2) + int(black_rank == 2):
            reconstruction_bad.add(mask)

        on_plateau = white_rank < 2 and black_rank < 2
        if not on_plateau:
            continue
        plateau_steps += _prefix_path_count(n, mask)
        if (
            black_error
            or white_error
            or black_line is None
            or white_line is None
            or black_line != white_line
            or black_index is None
            or white_index is None
        ):
            line_bad.add(mask)
            continue
        plateau_lines[mask] = black_line
        index_pairs[mask] = (black_index, white_index)
        maximum_index = max(maximum_index, black_index, white_index)

    full_mask = (1 << n) - 1
    for mask in range(1 << n):
        black_rank = state_table[mask][0][0]
        white_rank = state_table[mask][1][0]
        remaining = full_mask ^ mask
        while remaining:
            bit = remaining & -remaining
            remaining ^= bit
            child = mask | bit
            if state_table[child][0][0] < black_rank:
                raise AssertionError("black ambient rank decreased along a subset edge")
            if state_table[child][1][0] > white_rank:
                raise AssertionError("white ambient rank increased along a subset edge")

    bad_by_gate = {
        "birth": birth_bad,
        "rank_sum": rank_sum_bad,
        "reconstruction": reconstruction_bad,
    }
    failure_counts = {
        gate: total_paths - _count_paths_avoiding_nodes(n, bad_masks)
        for gate, bad_masks in bad_by_gate.items()
    }
    failure_counts["reflection"] = 0
    line_failures, line_counts = _marked_path_summary(
        n,
        plateau_lines,
        line_bad,
    )
    failure_counts["line"] = line_failures
    index_evolution, _ = _marked_path_summary(n, index_pairs, set())

    counterexamples = {
        gate: []
        for gate in ("birth", "reflection", "rank_sum", "reconstruction", "line")
    }
    for gate, count in failure_counts.items():
        if not count:
            continue
        if gate == "line":
            permutation = _first_mark_failure_permutation(
                n,
                plateau_lines,
                line_bad,
            )
        else:
            permutation = _first_bad_node_permutation(n, bad_by_gate[gate])
        if permutation is None:  # pragma: no cover - count/search consistency
            raise AssertionError("failed to recover counted counterexample")
        counterexamples[gate].append(
            analyze_permutation(
                geometry,
                permutation,
                state_table=state_table,
            )
        )

    return {
        "matrix": [list(part) for part in matrix],
        "order": n,
        "four_distinct_face_corners": has_four_distinct_face_corners(geometry),
        "permutations": total_paths,
        "cached_subsets": len(state_table),
        "rank_one_plateau_steps": plateau_steps,
        "plateau_line_counts": dict(sorted(line_counts.items())),
        "maximum_saturation_index": maximum_index,
        "permutations_with_index_evolution": index_evolution,
        "failure_counts": {
            gate: failure_counts[gate]
            for gate in ("birth", "reflection", "rank_sum", "reconstruction", "line")
        },
        "counterexamples": counterexamples,
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
    cached_subsets = sum(row["cached_subsets"] for row in rows)
    plateau_steps = sum(row["rank_one_plateau_steps"] for row in rows)
    maximum_index = max(row["maximum_saturation_index"] for row in rows)
    index_evolution = sum(row["permutations_with_index_evolution"] for row in rows)
    assert len(matrices) == 86
    assert total_paths == 70690518
    assert all(abs(determinant(tuple(tuple(part) for part in row["matrix"]))) == row["order"] for row in rows)
    return {
        "schema": "matching-one/digital-alexander-quotient-frontier/v1",
        "issue": 269,
        "status": (
            "counterexample_found"
            if any(total_failures.values())
            else "no_counterexample_through_index_10"
        ),
        "order_range": [2, 10],
        "HNF_representatives": len(rows),
        "honest_face_representatives": len(honest),
        "self_identifying_face_representatives": len(degenerate),
        "filtration_paths": total_paths,
        "cached_subsets": cached_subsets,
        "rank_one_plateau_steps": plateau_steps,
        "maximum_saturation_index": maximum_index,
        "permutations_with_index_evolution": index_evolution,
        "total_failure_counts": total_failures,
        "first_counterexamples": first_counterexamples,
        "geometries": rows,
        "claim_boundary": {
            "proved": "all permutations of every HNF quotient of index 2 through 10",
            "not_proved": "index 11 and above or an unrestricted degenerate-quotient theorem",
        },
    }


def render_markdown(artifact: dict[str, Any]) -> str:
    lines = [
        "# Short-period digital Alexander quotient frontier",
        "",
        "Every permutation of every two-dimensional HNF quotient of index 2 through 10 is exhausted.",
        "",
        "- HNF representatives: `%d`;" % artifact["HNF_representatives"],
        "- honest four-corner representatives: `%d`;" % artifact["honest_face_representatives"],
        "- self-identifying face representatives: `%d`;" % artifact["self_identifying_face_representatives"],
        "- complete filtration paths: `%d`." % artifact["filtration_paths"],
        "- cached occupied-site subsets: `%d`;" % artifact["cached_subsets"],
        "- rank-one plateau steps: `%d`;" % artifact["rank_one_plateau_steps"],
        "- maximum saturation index: `%d`;" % artifact["maximum_saturation_index"],
        "- paths with saturation-index evolution: `%d`."
        % artifact["permutations_with_index_evolution"],
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
            "The result is exhaustive only through index 10. It does not prove an unrestricted theorem",
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
