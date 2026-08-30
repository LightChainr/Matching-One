#!/usr/bin/env python3
"""Exact Gram/source intersection of the P321 affine sigma-Hom space."""

from __future__ import annotations

import argparse
from fractions import Fraction
from functools import reduce
from math import gcd, lcm
import json
from pathlib import Path
from typing import Any, Iterable, Sequence

from noncrossing_connectivity_codec import noncrossing_states
from p321_graded_closure_extension import intertwiner_constraints
from p321_homology_trace_certificate import (
    Matrix,
    action_matrix,
    exact_rank,
    join_adjacent,
    rotate_state,
)


SCHEMA = "matching-one/p333-gram-source-intertwiner/v1"
ROOT = Path(__file__).resolve().parents[1]
PROTOCOL = ROOT / "analysis/p333_gram_source_intertwiner_protocol.json"
DEFAULT_RESULT = ROOT / "results/p333-gram-source-intertwiner/latest.json"


def join_block_count(left: Sequence[int], right: Sequence[int]) -> int:
    """Number of blocks in the join of two connectivity partitions."""

    width = len(left)
    parent = list(range(width))

    def find(value: int) -> int:
        while parent[value] != value:
            parent[value] = parent[parent[value]]
            value = parent[value]
        return value

    def union(a: int, b: int) -> None:
        a_root, b_root = find(a), find(b)
        if a_root != b_root:
            parent[b_root] = a_root

    for labels in (left, right):
        representatives: dict[int, int] = {}
        for point, label in enumerate(labels):
            if label in representatives:
                union(point, representatives[label])
            else:
                representatives[label] = point
    return len({find(point) for point in range(width)})


def transpose(matrix: Sequence[Sequence[Fraction]]) -> list[list[Fraction]]:
    return [list(row) for row in zip(*matrix)]


def multiply(
    left: Sequence[Sequence[Fraction]], right: Sequence[Sequence[Fraction]]
) -> list[list[Fraction]]:
    columns = list(zip(*right))
    return [
        [sum(a * b for a, b in zip(row, column)) for column in columns]
        for row in left
    ]


def subtract(
    left: Sequence[Sequence[Fraction]], right: Sequence[Sequence[Fraction]]
) -> list[list[Fraction]]:
    return [
        [a - b for a, b in zip(left_row, right_row)]
        for left_row, right_row in zip(left, right)
    ]


def first_jet_radical_gram(states: Sequence[Sequence[int]]) -> list[list[Fraction]]:
    size = len(states)
    g1 = [
        [Fraction(join_block_count(left, right)) for right in states]
        for left in states
    ]
    basis = [
        [Fraction((row == column) - (row == size - 1)) for column in range(size - 1)]
        for row in range(size)
    ]
    return multiply(transpose(basis), multiply(g1, basis))


def rref_solve(
    coefficients: Sequence[Sequence[int | Fraction]],
    rhs: Sequence[int | Fraction],
    variables: int,
) -> dict[str, Any]:
    rows = [
        [Fraction(value) for value in row] + [Fraction(target)]
        for row, target in zip(coefficients, rhs)
    ]
    pivot_row = 0
    pivots: list[int] = []
    for column in range(variables):
        pivot = next(
            (index for index in range(pivot_row, len(rows)) if rows[index][column]),
            None,
        )
        if pivot is None:
            continue
        rows[pivot_row], rows[pivot] = rows[pivot], rows[pivot_row]
        scale = rows[pivot_row][column]
        rows[pivot_row] = [value / scale for value in rows[pivot_row]]
        for index in range(len(rows)):
            if index == pivot_row or not rows[index][column]:
                continue
            scale = rows[index][column]
            rows[index] = [
                value - scale * pivot_value
                for value, pivot_value in zip(rows[index], rows[pivot_row])
            ]
        pivots.append(column)
        pivot_row += 1
    inconsistent = any(
        not any(row[:variables]) and row[variables]
        for row in rows
    )
    if inconsistent:
        return {
            "consistent": False,
            "rank": len(pivots),
            "augmented_rank": len(pivots) + 1,
            "dimension": None,
            "particular": None,
            "nullspace": [],
        }

    free = [column for column in range(variables) if column not in pivots]
    particular = [Fraction(0) for _ in range(variables)]
    for row_index, pivot in enumerate(pivots):
        particular[pivot] = rows[row_index][variables]
    nullspace = []
    for free_column in free:
        vector = [Fraction(0) for _ in range(variables)]
        vector[free_column] = Fraction(1)
        for row_index, pivot in enumerate(pivots):
            vector[pivot] = -rows[row_index][free_column]
        nullspace.append(vector)
    return {
        "consistent": True,
        "rank": len(pivots),
        "augmented_rank": len(pivots),
        "dimension": len(free),
        "particular": particular,
        "nullspace": nullspace,
    }


def endpoint_constraints(size: int) -> tuple[list[list[int]], list[int]]:
    rows, rhs = [], []
    for column in range(size):
        row = [0] * (size * size)
        for output in range(size):
            row[output * size + column] = 1
        rows.append(row)
        rhs.append(1)
    return rows, rhs


def source_constraints(size: int, source: int) -> tuple[list[list[int]], list[int]]:
    rows, rhs = [], []
    for output in range(size):
        row = [0] * (size * size)
        row[output * size + source] = 1
        rows.append(row)
        rhs.append(int(output == source))
    return rows, rhs


def gram_constraints(h: Sequence[Sequence[Fraction]], size: int) -> tuple[list[list[Fraction]], list[int]]:
    radical_size = size - 1
    rows, rhs = [], []
    for left in range(radical_size):
        for right in range(left + 1, radical_size):
            row = [Fraction(0)] * (size * size)
            for pivot in range(radical_size):
                # (H A)_(left,right), A_(pivot,right)=X_(pivot,right)-X_(pivot,last)
                row[pivot * size + right] += h[left][pivot]
                row[pivot * size + (size - 1)] -= h[left][pivot]
                # -(A^T H)_(left,right), A_(pivot,left)=X_(pivot,left)-X_(pivot,last)
                row[pivot * size + left] -= h[pivot][right]
                row[pivot * size + (size - 1)] += h[pivot][right]
            rows.append(row)
            rhs.append(0)
    return rows, rhs


def reshape(vector: Sequence[Fraction], size: int) -> list[list[Fraction]]:
    return [list(vector[row * size : (row + 1) * size]) for row in range(size)]


def radical_action(matrix: Sequence[Sequence[Fraction]]) -> list[list[Fraction]]:
    size = len(matrix)
    return [
        [matrix[row][column] - matrix[row][size - 1] for column in range(size - 1)]
        for row in range(size - 1)
    ]


def primitive_integer_vector(vector: Sequence[Fraction]) -> list[int]:
    denominator = reduce(lcm, (value.denominator for value in vector), 1)
    integers = [value.numerator * (denominator // value.denominator) for value in vector]
    divisor = reduce(gcd, (abs(value) for value in integers if value), 0) or 1
    integers = [value // divisor for value in integers]
    first = next((value for value in integers if value), 1)
    return [-value for value in integers] if first < 0 else integers


def encode_fraction(value: Fraction) -> int | str:
    return value.numerator if value.denominator == 1 else str(value)


def encode_matrix(matrix: Sequence[Sequence[Fraction]]) -> list[list[int | str]]:
    return [[encode_fraction(value) for value in row] for row in matrix]


def matrix_residual_rank(matrix: Sequence[Sequence[Fraction]]) -> int:
    return exact_rank(tuple(tuple(value for value in row) for row in matrix))


def constraint_residual(
    coefficients: Sequence[Sequence[int | Fraction]],
    rhs: Sequence[int | Fraction],
    vector: Sequence[Fraction],
) -> list[Fraction]:
    return [
        sum(Fraction(value) * coordinate for value, coordinate in zip(row, vector))
        - Fraction(target)
        for row, target in zip(coefficients, rhs)
    ]


def residual_summary(residual: Sequence[Fraction]) -> dict[str, Any]:
    return {
        "nonzero_equations": sum(bool(value) for value in residual),
        "values_sha256": __import__("hashlib").sha256(
            json.dumps([str(value) for value in residual], separators=(",", ":")).encode()
        ).hexdigest(),
    }


def encoded_parameterization(solution: dict[str, Any], size: int) -> tuple[dict[str, Any] | None, list[list[Fraction]]]:
    if not solution["consistent"]:
        return None, []
    particular = reshape(solution["particular"], size)
    integer_vectors = [
        [Fraction(value) for value in primitive_integer_vector(vector)]
        for vector in solution["nullspace"]
    ]
    tangent_basis = [reshape(vector, size) for vector in integer_vectors]
    return (
        {
            "formula": "X=X0+sum_a theta_a Y_a",
            "particular_X0": encode_matrix(particular),
            "primitive_integer_tangent_basis": [
                encode_matrix(matrix) for matrix in tangent_basis
            ],
        },
        integer_vectors,
    )


def affine_inconsistency_witness(
    coefficients: Sequence[Sequence[Fraction]], rhs: Sequence[Fraction]
) -> dict[str, Any] | None:
    """Return y with y^T C=0 and y^T rhs=1 when C theta=rhs is inconsistent."""

    row_count = len(coefficients)
    column_count = len(coefficients[0]) if coefficients else 0
    transposed = [
        [coefficients[row][column] for row in range(row_count)]
        for column in range(column_count)
    ]
    left_null = rref_solve(transposed, [0] * column_count, row_count)
    if not left_null["consistent"]:
        raise AssertionError("a homogeneous left-null system cannot be inconsistent")
    for candidate in left_null["nullspace"]:
        pairing = sum(value * target for value, target in zip(candidate, rhs))
        if not pairing:
            continue
        witness = [value / pairing for value in candidate]
        residual = [
            sum(witness[row] * coefficients[row][column] for row in range(row_count))
            for column in range(column_count)
        ]
        if any(residual) or sum(value * target for value, target in zip(witness, rhs)) != 1:
            raise AssertionError("invalid affine inconsistency witness")
        return {
            "left_multiplier": [encode_fraction(value) for value in witness],
            "left_times_parameter_matrix": [encode_fraction(value) for value in residual],
            "left_times_rhs": 1,
            "identity": "y^T C=0 but y^T b=1",
        }
    return None


def width_result(width: int) -> dict[str, Any]:
    states = noncrossing_states(width)
    size = len(states)
    translation = action_matrix(width, lambda state: rotate_state(state, 1))
    joins = tuple(
        action_matrix(width, lambda state, site=site: join_adjacent(state, site))
        for site in range(width)
    )
    affine_rows = [
        list(row)
        for row in (
            intertwiner_constraints(joins, joins[1:] + joins[:1])
            + intertwiner_constraints((translation,), (translation,))
        )
    ]
    affine_rhs = [0] * len(affine_rows)
    endpoint_rows, endpoint_rhs = endpoint_constraints(size)
    h = first_jet_radical_gram(states)
    gram_rows, gram_rhs = gram_constraints(h, size)
    source_state = tuple(range(width))
    source_index = states.index(source_state)
    source_rows, source_rhs = source_constraints(size, source_index)
    stages = [
        ("affine_sigma", affine_rows, affine_rhs),
        (
            "endpoint_normalized",
            affine_rows + endpoint_rows,
            affine_rhs + endpoint_rhs,
        ),
        (
            "gram_self_adjoint",
            affine_rows + endpoint_rows + gram_rows,
            affine_rhs + endpoint_rhs + gram_rhs,
        ),
        (
            "source_normalized",
            affine_rows + endpoint_rows + gram_rows + source_rows,
            affine_rhs + endpoint_rhs + gram_rhs + source_rhs,
        ),
    ]

    translation_vector = [Fraction(value) for row in translation for value in row]
    stage_results = {}
    solved = None
    solutions = {}
    for name, rows, rhs in stages:
        solution = rref_solve(rows, rhs, size * size)
        solutions[name] = solution
        stage_results[name] = {
            "equations": len(rows),
            "consistent": solution["consistent"],
            "rank": solution["rank"],
            "affine_tangent_dimension": solution["dimension"],
            "canonical_X_equals_T": residual_summary(
                constraint_residual(rows, rhs, translation_vector)
            ),
        }
        if name == "source_normalized":
            solved = solution

    assert solved is not None
    canonical_radical = radical_action(
        [[Fraction(value) for value in row] for row in translation]
    )
    canonical_gram_residual = subtract(
        multiply(h, canonical_radical),
        multiply(transpose(canonical_radical), h),
    )
    canonical_in_final = (
        stage_results["source_normalized"]["canonical_X_equals_T"]["nonzero_equations"]
        == 0
    )

    parameterization, _final_integer_vectors = encoded_parameterization(solved, size)
    gram_parameterization, gram_integer_vectors = encoded_parameterization(
        solutions["gram_self_adjoint"], size
    )
    gram_particular = solutions["gram_self_adjoint"]["particular"]
    source_parameter_rows = []
    source_parameter_rhs = []
    if gram_particular is not None:
        for row, target in zip(source_rows, source_rhs):
            source_parameter_rows.append(
                [
                    sum(Fraction(value) * coordinate for value, coordinate in zip(row, vector))
                    for vector in gram_integer_vectors
                ]
            )
            source_parameter_rhs.append(
                Fraction(target)
                - sum(
                    Fraction(value) * coordinate
                    for value, coordinate in zip(row, gram_particular)
                )
            )
        source_parameter_solution = rref_solve(
            source_parameter_rows,
            source_parameter_rhs,
            len(gram_integer_vectors),
        )
        source_image_x0 = [
            gram_particular[output * size + source_index]
            for output in range(size)
        ]
        source_image_tangents = [
            [vector[output * size + source_index] for output in range(size)]
            for vector in gram_integer_vectors
        ]
        witness = (
            affine_inconsistency_witness(
                source_parameter_rows, source_parameter_rhs
            )
            if not source_parameter_solution["consistent"]
            else None
        )
        if witness is not None:
            witness["nonzero_state_coefficients"] = [
                {
                    "state": list(states[index]),
                    "coefficient": coefficient,
                }
                for index, coefficient in enumerate(witness["left_multiplier"])
                if coefficient != 0
            ]
        source_restriction = {
            "parameter_count_before_source": len(gram_integer_vectors),
            "coefficient_rank": source_parameter_solution["rank"],
            "augmented_rank": source_parameter_solution["augmented_rank"],
            "consistent": source_parameter_solution["consistent"],
            "remaining_dimension": source_parameter_solution["dimension"],
            "source_target": [int(output == source_index) for output in range(size)],
            "source_image_X0": [encode_fraction(value) for value in source_image_x0],
            "source_images_tangent_basis": [
                [encode_fraction(value) for value in vector]
                for vector in source_image_tangents
            ],
            "inconsistency_witness": witness,
        }
    else:
        source_restriction = None

    if not solved["consistent"]:
        decision = "empty_intersection"
    elif canonical_in_final and solved["dimension"] == 0:
        decision = "canonical_unique"
    elif canonical_in_final:
        decision = "canonical_not_unique"
    else:
        decision = "canonical_excluded"

    return {
        "width": width,
        "module_dimension": size,
        "radical_dimension": size - 1,
        "source_state": list(source_state),
        "source_index": source_index,
        "first_jet_radical_gram": {
            "rank": matrix_residual_rank(h),
            "matrix": encode_matrix(h),
        },
        "stages": stage_results,
        "canonical_translation": {
            "matrix": [list(row) for row in translation],
            "is_identity": all(
                translation[row][column] == int(row == column)
                for row in range(size)
                for column in range(size)
            ),
            "gram_self_adjoint_residual_rank": matrix_residual_rank(
                canonical_gram_residual
            ),
            "gram_self_adjoint_residual": encode_matrix(canonical_gram_residual),
            "belongs_to_final_intersection": canonical_in_final,
        },
        "gram_stage_parameterization": gram_parameterization,
        "source_restriction_on_gram_moduli": source_restriction,
        "final_parameterization": parameterization,
        "decision": decision,
    }


def build_result() -> dict[str, Any]:
    protocol = json.loads(PROTOCOL.read_text(encoding="utf-8"))
    if protocol["status"] != "frozen_before_exact_intersection":
        raise AssertionError("protocol is not frozen")
    widths = [width_result(width) for width in protocol["spaces"]["widths"]]
    return {
        "schema": SCHEMA,
        "status": "exact_rational_intersection_certificate",
        "issues": [333, 370],
        "protocol": {
            "path": str(PROTOCOL.relative_to(ROOT)),
            "sha256": __import__("hashlib").sha256(PROTOCOL.read_bytes()).hexdigest(),
            "commit": "023d4de",
        },
        "widths": widths,
        "decision_summary": {
            str(row["width"]): row["decision"] for row in widths
        },
        "global_decision": (
            "Width two selects X=T only because translation is identity in the two-state quotient. "
            "At the first nondegenerate widths three and four, T is not first-jet-Gram self-adjoint, "
            "and the entire Gram-compatible affine sigma-Hom slice is incompatible with fixing the "
            "all-singleton source. The frozen physical intersection is empty there."
        ),
        "deduplication": protocol["deduplication"],
        "claim_boundary": protocol["claim_boundary"],
    }


def render_markdown(result: dict[str, Any]) -> str:
    lines = [
        "# P333/P370 Gram-source affine intertwiner intersection",
        "",
        "## Exact dimensions",
        "",
        "| width | dim V | affine Hom | endpoint-normalized | + Gram | + source | X=T in final | decision |",
        "|---:|---:|---:|---:|---:|---:|:---:|---|",
    ]
    for row in result["widths"]:
        stages = row["stages"]
        source_dimension = stages["source_normalized"]["affine_tangent_dimension"]
        lines.append(
            f"| {row['width']} | {row['module_dimension']} | "
            f"{stages['affine_sigma']['affine_tangent_dimension']} | "
            f"{stages['endpoint_normalized']['affine_tangent_dimension']} | "
            f"{stages['gram_self_adjoint']['affine_tangent_dimension']} | "
            f"{source_dimension if source_dimension is not None else '—'} | "
            f"{'yes' if row['canonical_translation']['belongs_to_final_intersection'] else 'no'} | "
            f"`{row['decision']}` |"
        )
    lines.extend(
        [
            "",
            "The reported normalized dimensions are affine tangent dimensions. Exact rational particular solutions and primitive-integer tangent bases are stored in `latest.json`.",
            "",
            "## Interpretation",
            "",
            result["global_decision"],
            "",
        ]
    )
    for row in result["widths"]:
        lines.append(
            f"- Width {row['width']}: canonical `T` has Gram residual rank "
            f"{row['canonical_translation']['gram_self_adjoint_residual_rank']}; "
            f"final decision `{row['decision']}`."
        )
        witness = row["source_restriction_on_gram_moduli"]["inconsistency_witness"]
        if witness is not None:
            lines.append(
                f"  Exact inconsistency witness: `{witness['identity']}` with nonzero state "
                f"coefficients `{witness['nonzero_state_coefficients']}`."
            )
    lines.extend(
        [
            "",
            "This intersection is distinct from the signed-history nilpotent `K=(D-J)^2`: K lives in the three-mark endpoint radical, whereas this calculation selects off-diagonal maps between crossed/trivial affine closure modules.",
            "It is also distinct from PR #385: an exact selected line does not by itself make Jordan behavior statistically separated from near-colliding ordinary models.",
            "",
            "## Boundary",
            "",
        ]
    )
    lines.extend(f"- {item}" for item in result["claim_boundary"])
    lines.append("")
    return "\n".join(lines)


def render_card(result: dict[str, Any]) -> str:
    decisions = ", ".join(
        f"w{row['width']}={row['decision']}"
        for row in result["widths"]
    )
    return "\n".join(
        [
            "# Scientific card: Gram/source affine intertwiner gate",
            "",
            f"- **Mechanism space changed:** `{decisions}` under the frozen exact intersection.",
            f"- **Result:** {result['global_decision']}",
            "- **Observer/sector/source:** crossed-to-trivial affine closure intertwiner | first-jet connectivity radical | all-singleton connectivity source.",
            "- **Not proved:** no physical thermal insertion, generic-Q detach algebra, continuum Jordan field or noise-level identifiability.",
            "- **Deduplication:** d0fca79 supplies unconstrained affine Hom; 6c60b0e supplies a separate signed-history Gram control; PR #385 supplies the finite-noise closure theorem.",
            "- **Next upweight:** add one declared generic-Q detach/loop-weight generator and rerun the same exact intersection without changing Gram or source conventions.",
            "",
        ]
    )


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", type=Path, default=DEFAULT_RESULT)
    parser.add_argument("--markdown", type=Path, required=True)
    parser.add_argument("--card", type=Path, required=True)
    args = parser.parse_args(argv)
    result = build_result()
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.markdown.parent.mkdir(parents=True, exist_ok=True)
    args.card.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    args.markdown.write_text(render_markdown(result), encoding="utf-8")
    args.card.write_text(render_card(result), encoding="utf-8")
    print(json.dumps(result["decision_summary"], indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
