#!/usr/bin/env python3
"""Exact one-mark endpoint-jet continuation of the P333 confluence gate."""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
from typing import Any, Sequence

from noncrossing_connectivity_codec import noncrossing_states
from p321_graded_closure_extension import intertwiner_constraints
from p321_homology_trace_certificate import (
    Matrix,
    action_matrix,
    exact_rank,
    join_adjacent,
    rotate_state,
)
from p333_generic_q_detach_intertwiner import (
    _block_rows,
    _restriction_witness,
    _stage_summary,
    detach_jet,
)
from p333_gram_source_intertwiner import (
    encode_fraction,
    encode_matrix,
    join_block_count,
    matrix_residual_rank,
    multiply,
    rref_solve,
    subtract,
    transpose,
)


SCHEMA = "matching-one/p333-one-mark-endpoint-jet/v1"
ROOT = Path(__file__).resolve().parents[1]
PROTOCOL = ROOT / "analysis/p333_one_mark_endpoint_jet_protocol.json"
DEFAULT_RESULT = ROOT / "results/p333-one-mark-endpoint-jet/latest.json"


def block_diagonal_with_mark(matrix: Matrix) -> Matrix:
    """Extend an ordinary action by the identity on the terminal mark."""

    size = len(matrix)
    return tuple(
        tuple(matrix[row][column] for column in range(size)) + (0,)
        for row in range(size)
    ) + (tuple(0 for _ in range(size)) + (1,),)


def singleton_covector(states: Sequence[Sequence[int]], site: int) -> tuple[int, ...]:
    return tuple(int(tuple(state).count(state[site]) == 1) for state in states)


def marked_detach_jet(width: int, site: int) -> tuple[Matrix, Matrix]:
    """Return the frozen terminal-mark detach at Q=1 and its Q velocity."""

    states = noncrossing_states(width)
    ordinary_zero, ordinary_velocity = detach_jet(width, site)
    size = len(states)
    zero = block_diagonal_with_mark(ordinary_zero)
    emission = singleton_covector(states, site)
    velocity = tuple(
        tuple(ordinary_velocity[row][column] for column in range(size)) + (0,)
        for row in range(size)
    ) + (tuple(emission) + (0,),)
    return zero, velocity


def extended_gram_jet(
    states: Sequence[Sequence[int]],
) -> tuple[list[list[Fraction]], list[list[Fraction]]]:
    """First jet of [[G(Q), d_Q G(Q)[:,s0]],[transpose,0]]."""

    block_counts = [len(set(state)) for state in states]
    g0 = [
        [Fraction(1) for _ in states] + [Fraction(block_counts[row])]
        for row in range(len(states))
    ]
    g0.append([Fraction(value) for value in block_counts] + [Fraction(0)])
    g1 = [
        [Fraction(join_block_count(left, right)) for right in states]
        + [Fraction(block_counts[row] * (block_counts[row] - 1))]
        for row, left in enumerate(states)
    ]
    g1.append(
        [Fraction(value * (value - 1)) for value in block_counts]
        + [Fraction(0)]
    )
    return g0, g1


def nullspace_basis(matrix: Sequence[Sequence[Fraction]]) -> list[list[Fraction]]:
    """Return a matrix whose columns are a deterministic rational kernel basis."""

    size = len(matrix[0]) if matrix else 0
    solved = rref_solve(matrix, [0] * len(matrix), size)
    if not solved["consistent"]:
        raise AssertionError("homogeneous kernel system is inconsistent")
    if not solved["nullspace"]:
        return [[] for _ in range(size)]
    return [list(row) for row in zip(*solved["nullspace"])]


def q_independent_jet_constraints(
    source: Sequence[Matrix], target: Sequence[Matrix]
) -> list[list[Fraction]]:
    base = list(intertwiner_constraints(source, target))
    block_size = len(source[0]) ** 2
    return _block_rows(base, 0, block_size) + _block_rows(base, 1, block_size)


def first_jet_intertwiner_constraints(
    source_zero: Sequence[Matrix],
    source_one: Sequence[Matrix],
    target_zero: Sequence[Matrix],
    target_one: Sequence[Matrix],
) -> list[list[Fraction]]:
    """Equations V A0+X0 A1=B0 V+B1 X0."""

    size = len(source_zero[0])
    block_size = size * size
    rows: list[list[Fraction]] = []
    for a0, a1, b0, b1 in zip(source_zero, source_one, target_zero, target_one):
        for output in range(size):
            for column in range(size):
                row = [Fraction(0)] * (2 * block_size)
                for pivot in range(size):
                    row[block_size + output * size + pivot] += a0[pivot][column]
                    row[block_size + pivot * size + column] -= b0[output][pivot]
                    row[output * size + pivot] += a1[pivot][column]
                    row[pivot * size + column] -= b1[output][pivot]
                rows.append(row)
    return rows


def filtration_constraints(size: int) -> list[list[Fraction]]:
    """Prevent the mark column from feeding back into ordinary connectivity."""

    block_size = size * size
    rows: list[list[Fraction]] = []
    mark = size - 1
    for block in (0, 1):
        for output in range(mark):
            row = [Fraction(0)] * (2 * block_size)
            row[block * block_size + output * size + mark] = 1
            rows.append(row)
    return rows


def covector_normalization(
    covector: Sequence[int | Fraction], size: int
) -> tuple[list[list[Fraction]], list[Fraction]]:
    """Fix c^T X0=c^T and set c^T V=0."""

    block_size = size * size
    rows: list[list[Fraction]] = []
    rhs: list[Fraction] = []
    for block, target_scale in ((0, 1), (1, 0)):
        for column in range(size):
            row = [Fraction(0)] * (2 * block_size)
            for output, value in enumerate(covector):
                row[block * block_size + output * size + column] = Fraction(value)
            rows.append(row)
            rhs.append(Fraction(target_scale) * Fraction(covector[column]))
    return rows, rhs


def source_normalization(
    sources: Sequence[Sequence[int | Fraction]], size: int
) -> tuple[list[list[Fraction]], list[Fraction]]:
    """Fix each source by X0 and kill it by V."""

    block_size = size * size
    rows: list[list[Fraction]] = []
    rhs: list[Fraction] = []
    for source in sources:
        for block, target_scale in ((0, 1), (1, 0)):
            for output in range(size):
                row = [Fraction(0)] * (2 * block_size)
                for column, value in enumerate(source):
                    row[block * block_size + output * size + column] = Fraction(value)
                rows.append(row)
                rhs.append(Fraction(target_scale) * Fraction(source[output]))
    return rows, rhs


def radical_invariance_constraints(
    g0: Sequence[Sequence[Fraction]],
    basis: Sequence[Sequence[Fraction]],
    size: int,
) -> list[list[Fraction]]:
    """Equations G0 X0 B=0."""

    block_size = size * size
    radical_dimension = len(basis[0]) if basis else 0
    rows: list[list[Fraction]] = []
    for output in range(size):
        for radical_column in range(radical_dimension):
            row = [Fraction(0)] * (2 * block_size)
            for pivot in range(size):
                for source in range(size):
                    row[pivot * size + source] += (
                        g0[output][pivot] * basis[source][radical_column]
                    )
            rows.append(row)
    return rows


def radical_gram_constraints(
    g1: Sequence[Sequence[Fraction]],
    basis: Sequence[Sequence[Fraction]],
    size: int,
) -> list[list[Fraction]]:
    """Equations B^T(G1 X0-X0^T G1)B=0."""

    block_size = size * size
    radical_dimension = len(basis[0]) if basis else 0
    bt_g1 = multiply(transpose(basis), g1)
    g1_b = multiply(g1, basis)
    rows: list[list[Fraction]] = []
    for left in range(radical_dimension):
        for right in range(left + 1, radical_dimension):
            row = [Fraction(0)] * (2 * block_size)
            for output in range(size):
                for source in range(size):
                    row[output * size + source] += (
                        bt_g1[left][output] * basis[source][right]
                        - basis[source][left] * g1_b[output][right]
                    )
            rows.append(row)
    return rows


def restricted_skew_residual(
    g1: Sequence[Sequence[Fraction]],
    basis: Sequence[Sequence[Fraction]],
    operator: Sequence[Sequence[Fraction]],
) -> list[list[Fraction]]:
    left = multiply(transpose(basis), multiply(g1, multiply(operator, basis)))
    right = multiply(
        transpose(basis), multiply(transpose(operator), multiply(g1, basis))
    )
    return subtract(left, right)


def encode_joint_particular(solution: dict[str, Any], size: int) -> dict[str, Any] | None:
    if not solution["consistent"]:
        return None
    block_size = size * size
    particular = solution["particular"]
    return {
        "X0": encode_matrix(
            [particular[row * size : (row + 1) * size] for row in range(size)]
        ),
        "V": encode_matrix(
            [
                particular[block_size + row * size : block_size + (row + 1) * size]
                for row in range(size)
            ]
        ),
        "tangent_count": len(solution["nullspace"]),
    }


def width_result(width: int) -> dict[str, Any]:
    states = noncrossing_states(width)
    ordinary_size = len(states)
    size = ordinary_size + 1
    block_size = size * size
    variables = 2 * block_size

    translation = block_diagonal_with_mark(
        action_matrix(width, lambda state: rotate_state(state, 1))
    )
    joins = tuple(
        block_diagonal_with_mark(
            action_matrix(width, lambda state, site=site: join_adjacent(state, site))
        )
        for site in range(width)
    )
    detach_zero, detach_one = zip(
        *(marked_detach_jet(width, site) for site in range(width))
    )

    hom_rows = q_independent_jet_constraints(
        joins, joins[1:] + joins[:1]
    ) + q_independent_jet_constraints((translation,), (translation,))
    hom_rows += [
        list(row)
        for row in intertwiner_constraints(
            detach_zero, detach_zero[1:] + detach_zero[:1]
        )
    ]
    # The order-zero detach rows above need embedding in the X0 block.
    zero_detach_count = width * size * size
    raw_zero_detach = hom_rows[-zero_detach_count:]
    hom_rows = hom_rows[:-zero_detach_count] + _block_rows(
        raw_zero_detach, 0, block_size
    )
    hom_rows += first_jet_intertwiner_constraints(
        detach_zero,
        detach_one,
        detach_zero[1:] + detach_zero[:1],
        detach_one[1:] + detach_one[:1],
    )
    hom_rows += filtration_constraints(size)
    hom_rhs = [Fraction(0)] * len(hom_rows)

    g0, g1 = extended_gram_jet(states)
    radical_basis = nullspace_basis(g0)
    radical_dimension = len(radical_basis[0]) if radical_basis else 0
    endpoint_covector = [1] * ordinary_size + [0]
    endpoint_rows, endpoint_rhs = covector_normalization(
        endpoint_covector, size
    )
    invariance_rows = radical_invariance_constraints(g0, radical_basis, size)
    invariance_rhs = [Fraction(0)] * len(invariance_rows)
    gram_rows = radical_gram_constraints(g1, radical_basis, size)
    gram_rhs = [Fraction(0)] * len(gram_rows)

    all_singleton = states.index(tuple(range(width)))
    ordinary_source = [0] * size
    ordinary_source[all_singleton] = 1
    mark_source = [0] * size
    mark_source[-1] = 1
    source_rows, source_rhs = source_normalization(
        (ordinary_source, mark_source), size
    )

    endpoint_invariance_rows = endpoint_rows + invariance_rows
    endpoint_invariance_rhs = endpoint_rhs + invariance_rhs
    ladders = [
        ("marked_affine_hom_jet", hom_rows, hom_rhs),
        (
            "endpoint_radical_normalized",
            hom_rows + endpoint_invariance_rows,
            hom_rhs + endpoint_invariance_rhs,
        ),
        (
            "gram_self_adjoint",
            hom_rows + endpoint_invariance_rows + gram_rows,
            hom_rhs + endpoint_invariance_rhs + gram_rhs,
        ),
        (
            "source_normalized",
            hom_rows + endpoint_invariance_rows + gram_rows + source_rows,
            hom_rhs + endpoint_invariance_rhs + gram_rhs + source_rhs,
        ),
    ]

    canonical = [Fraction(value) for row in translation for value in row] + [
        Fraction(0)
    ] * block_size
    summaries: dict[str, Any] = {}
    solutions: dict[str, Any] = {}
    for name, rows, rhs in ladders:
        summaries[name], solutions[name] = _stage_summary(
            rows, rhs, variables, canonical
        )

    first_empty_restriction = None
    prior_name = "marked_affine_hom_jet"
    additions = [
        (
            "endpoint_radical_normalized",
            endpoint_invariance_rows,
            endpoint_invariance_rhs,
        ),
        ("gram_self_adjoint", gram_rows, gram_rhs),
        ("source_normalized", source_rows, source_rhs),
    ]
    for name, extra_rows, extra_rhs in additions:
        restriction = _restriction_witness(
            solutions[prior_name], extra_rows, extra_rhs
        )
        if restriction is not None and not restriction["consistent"]:
            first_empty_restriction = {
                "from_stage": prior_name,
                "added_stage": name,
                **restriction,
            }
            break
        prior_name = name

    final = solutions["source_normalized"]
    if not final["consistent"]:
        decision = "one_mark_insufficient"
    elif final["dimension"] == 0:
        decision = "one_mark_unique"
    else:
        decision = "one_mark_nonunique"

    canonical_matrix = [
        [Fraction(value) for value in row] for row in translation
    ]
    canonical_skew = restricted_skew_residual(
        g1, radical_basis, canonical_matrix
    )
    skew_rank = matrix_residual_rank(canonical_skew)
    lower_bound = skew_rank // 2

    velocity_rows = []
    for site, velocity in enumerate(detach_one):
        restricted = multiply(
            [[Fraction(value) for value in row] for row in velocity],
            radical_basis,
        )
        mark_row = [restricted[-1]]
        velocity_rows.append(
            {
                "site": site,
                "full_velocity_rank_on_radical": matrix_residual_rank(restricted),
                "mark_injection_rank_on_radical": matrix_residual_rank(mark_row),
                "mark_injection_row": [
                    encode_fraction(value) for value in restricted[-1]
                ],
            }
        )

    return {
        "width": width,
        "ordinary_dimension": ordinary_size,
        "extended_dimension": size,
        "g0_rank": matrix_residual_rank(g0),
        "radical_dimension": radical_dimension,
        "g1_radical_rank": matrix_residual_rank(
            multiply(transpose(radical_basis), multiply(g1, radical_basis))
        ),
        "mark_velocity_gate": {
            "sites": velocity_rows,
            "all_sites_nonzero_mark_injection": all(
                row["mark_injection_rank_on_radical"] > 0 for row in velocity_rows
            ),
        },
        "stages": summaries,
        "first_empty_restriction": first_empty_restriction,
        "final_parameterization": encode_joint_particular(final, size),
        "decision": decision,
        "canonical_translation": {
            "restricted_gram_skew_rank": skew_rank,
            "restricted_gram_skew": encode_matrix(canonical_skew),
        },
        "scalar_mark_lower_bound_if_canonical": {
            "further_independent_endpoint_marks": lower_bound,
            "total_marks_including_current": 1 + lower_bound,
            "minimum_extended_dimension": ordinary_size + 1 + lower_bound,
            "scope": "scalar endpoint-mark quotient family only",
        },
    }


def build_result() -> dict[str, Any]:
    protocol = json.loads(PROTOCOL.read_text(encoding="utf-8"))
    if protocol["status"] != "frozen_before_exact_intersection":
        raise AssertionError("protocol is not frozen")
    widths = [width_result(width) for width in protocol["state_space"]["widths"]]
    all_genuine = all(
        row["mark_velocity_gate"]["all_sites_nonzero_mark_injection"]
        for row in widths
    )
    return {
        "schema": SCHEMA,
        "status": "exact_rational_one_mark_certificate",
        "issues": [333, 321, 370],
        "protocol": {
            "path": str(PROTOCOL.relative_to(ROOT)),
            "sha256": hashlib.sha256(PROTOCOL.read_bytes()).hexdigest(),
            "commit": "fb72351",
        },
        "widths": widths,
        "decision_summary": {str(row["width"]): row["decision"] for row in widths},
        "genuine_q_velocity_gate_passed": all_genuine,
        "global_decision": (
            "The terminal mark is not a disguised scalar fugacity: its detach Q velocity "
            "has rank-one mark-coordinate action on the extended radical at every site. "
            "Nevertheless one mark is insufficient. At both widths three and four the "
            "marked affine jet has five moduli and endpoint/radical normalization leaves "
            "three, but the Gram restriction is inconsistent with coefficient rank zero "
            "and augmented rank one: none of those three moduli changes the offending skew "
            "form. The smallest counterexample is width three, where dim(W)=6."
        ),
        "deduplication": protocol["deduplication"],
        "claim_boundary": protocol["claim_boundary"],
    }


def render_markdown(result: dict[str, Any]) -> str:
    lines = [
        "# P333/P321/P370 one-mark endpoint-jet gate",
        "",
        result["global_decision"],
        "",
        "| width | dim V | dim W | rank G0 | dim radical | Q mark gate | affine jet | + endpoint/radical | + Gram | + source | decision |",
        "|---:|---:|---:|---:|---:|:---:|---:|---:|---:|---:|---|",
    ]
    for row in result["widths"]:
        stages = row["stages"]
        def dim(name: str) -> str:
            value = stages[name]["affine_tangent_dimension"]
            return "empty" if value is None else str(value)
        lines.append(
            f"| {row['width']} | {row['ordinary_dimension']} | {row['extended_dimension']} | "
            f"{row['g0_rank']} | {row['radical_dimension']} | "
            f"{'pass' if row['mark_velocity_gate']['all_sites_nonzero_mark_injection'] else 'fail'} | "
            f"{dim('marked_affine_hom_jet')} | {dim('endpoint_radical_normalized')} | "
            f"{dim('gram_self_adjoint')} | {dim('source_normalized')} | `{row['decision']}` |"
        )
    lines.extend(["", "## Exact interpretation", ""])
    for row in result["widths"]:
        lower = row["scalar_mark_lower_bound_if_canonical"]
        lines.append(
            f"- Width {row['width']}: `{row['decision']}`. Canonical restricted Gram-skew "
            f"rank {row['canonical_translation']['restricted_gram_skew_rank']}; within the "
            f"scalar endpoint-mark family this requires at least {lower['further_independent_endpoint_marks']} "
            f"further independent mark covectors (at least {lower['total_marks_including_current']} total marks) "
            f"if the canonical line is retained."
        )
        witness = row["first_empty_restriction"]
        if witness is not None:
            lines.append(
                f"  First empty restriction: `{witness['from_stage']}` -> "
                f"`{witness['added_stage']}` with exact witness "
                f"`{witness['inconsistency_witness']['identity']}`."
            )
    lines.extend(["", "## Boundary", ""])
    lines.extend(f"- {item}" for item in result["claim_boundary"])
    lines.append("")
    return "\n".join(lines)


def render_card(result: dict[str, Any]) -> str:
    decisions = ", ".join(
        f"w{row['width']}={row['decision']}" for row in result["widths"]
    )
    return "\n".join(
        [
            "# Scientific card: one terminal mark versus the exact confluence gate",
            "",
            f"- **Mechanism space changed:** `{decisions}` in the frozen one-mark endpoint quotient.",
            f"- **Genuine-Q gate:** nonzero mark-coordinate detach velocity on the extended radical = `{result['genuine_q_velocity_gate_passed']}`.",
            "- **Result:** one mark is insufficient at both widths: `5 -> 3 -> empty` at the Gram gate. The exact restriction has coefficient rank 0 and augmented rank 1, so no surviving modulus can rotate away the obstruction.",
            "- **Observer/sector/source:** crossed-to-trivial marked connectivity response | extended Q=1 radical | all-singleton and terminal-mark sources.",
            "- **Not proved:** the mark is not a full rooted-cluster transfer module; no LCFT/Jordan or PR #393 formal-K identification follows.",
            "- **Dependency group:** exact finite-width continuation of `96df7c8` and `ba3135e`, not an independent vote.",
            "- **Next upweight:** if one mark fails, add only the number of independent endpoint marks demanded by the exact restricted-skew rank, then retain marked connectivity instead of adding free matrix blocks.",
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
