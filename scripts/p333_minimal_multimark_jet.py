#!/usr/bin/env python3
"""Exact lower-bound-attainment test for scalar endpoint marks."""

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
    subtract,
    transpose,
)
from p333_one_mark_endpoint_jet import (
    covector_normalization,
    first_jet_intertwiner_constraints,
    nullspace_basis,
    q_independent_jet_constraints,
    radical_gram_constraints,
    radical_invariance_constraints,
    restricted_skew_residual,
    source_normalization,
)


SCHEMA = "matching-one/p333-minimal-multimark-jet/v1"
ROOT = Path(__file__).resolve().parents[1]
PROTOCOL = ROOT / "analysis/p333_minimal_multimark_jet_protocol.json"
DEFAULT_RESULT = ROOT / "results/p333-minimal-multimark-jet/latest.json"


def falling_factorial(value: int, order: int) -> int:
    output = 1
    for step in range(order):
        output *= value - step
    return output


def block_diagonal_with_marks(matrix: Matrix, marks: int) -> Matrix:
    ordinary = len(matrix)
    return tuple(
        tuple(matrix[row][column] for column in range(ordinary))
        + tuple(0 for _ in range(marks))
        for row in range(ordinary)
    ) + tuple(
        tuple(0 for _ in range(ordinary))
        + tuple(int(row == column) for column in range(marks))
        for row in range(marks)
    )


def multimark_detach_jet(
    width: int, site: int, marks: int
) -> tuple[Matrix, Matrix]:
    states = noncrossing_states(width)
    ordinary_zero, ordinary_velocity = detach_jet(width, site)
    ordinary = len(states)
    zero = block_diagonal_with_marks(ordinary_zero, marks)
    velocity_top = tuple(
        tuple(ordinary_velocity[row][column] for column in range(ordinary))
        + tuple(0 for _ in range(marks))
        for row in range(ordinary)
    )
    velocity_marks = []
    for order in range(1, marks + 1):
        row = []
        for state in states:
            singleton = tuple(state).count(state[site]) == 1
            blocks_after_closure = len(set(state)) - 1
            row.append(
                int(singleton)
                * falling_factorial(blocks_after_closure, order - 1)
            )
        velocity_marks.append(
            tuple(row) + tuple(0 for _ in range(marks))
        )
    return zero, velocity_top + tuple(velocity_marks)


def extended_multimark_gram_jet(
    states: Sequence[Sequence[int]], marks: int
) -> tuple[list[list[Fraction]], list[list[Fraction]]]:
    block_counts = [len(set(state)) for state in states]
    ordinary = len(states)
    g0 = []
    g1 = []
    for left_index, left in enumerate(states):
        g0.append(
            [Fraction(1) for _ in states]
            + [
                Fraction(falling_factorial(block_counts[left_index], order))
                for order in range(1, marks + 1)
            ]
        )
        g1.append(
            [Fraction(join_block_count(left, right)) for right in states]
            + [
                Fraction(falling_factorial(block_counts[left_index], order + 1))
                for order in range(1, marks + 1)
            ]
        )
    for order in range(1, marks + 1):
        g0.append(
            [
                Fraction(falling_factorial(value, order))
                for value in block_counts
            ]
            + [Fraction(0) for _ in range(marks)]
        )
        g1.append(
            [
                Fraction(falling_factorial(value, order + 1))
                for value in block_counts
            ]
            + [Fraction(0) for _ in range(marks)]
        )
    if len(g0) != ordinary + marks:
        raise AssertionError("extended Gram size mismatch")
    return g0, g1


def multimark_filtration_constraints(
    ordinary: int, marks: int
) -> list[list[Fraction]]:
    size = ordinary + marks
    block_size = size * size
    rows: list[list[Fraction]] = []
    for block in (0, 1):
        for output in range(ordinary):
            for mark in range(marks):
                row = [Fraction(0)] * (2 * block_size)
                row[
                    block * block_size + output * size + ordinary + mark
                ] = 1
                rows.append(row)
    return rows


def encode_final_solution(
    solution: dict[str, Any], size: int, canonical: Sequence[Fraction]
) -> dict[str, Any] | None:
    if not solution["consistent"]:
        return None
    block_size = size * size
    particular = solution["particular"]
    velocity = particular[block_size:]
    canonical_in = all(
        coordinate == target
        for coordinate, target in zip(particular, canonical)
    )
    # Solve directly whether canonical-particular lies in the tangent span.
    if solution["nullspace"]:
        from p333_gram_source_intertwiner import rref_solve

        basis_matrix = [
            [tangent[index] for tangent in solution["nullspace"]]
            for index in range(2 * block_size)
        ]
        target = [canonical[index] - particular[index] for index in range(2 * block_size)]
        canonical_in = rref_solve(
            basis_matrix, target, len(solution["nullspace"])
        )["consistent"]
    return {
        "particular_X0": encode_matrix(
            [particular[row * size : (row + 1) * size] for row in range(size)]
        ),
        "particular_V": encode_matrix(
            [
                particular[block_size + row * size : block_size + (row + 1) * size]
                for row in range(size)
            ]
        ),
        "affine_tangent_dimension": solution["dimension"],
        "contains_canonical_X0_Tbar_V0": canonical_in,
        "particular_velocity_is_zero": all(value == 0 for value in velocity),
        "all_tangent_velocities_are_zero": all(
            all(value == 0 for value in tangent[block_size:])
            for tangent in solution["nullspace"]
        ),
        "velocity_tangent_dimension": matrix_residual_rank(
            [list(tangent[block_size:]) for tangent in solution["nullspace"]]
        )
        if solution["nullspace"]
        else 0,
    }


def case_result(width: int, marks: int) -> dict[str, Any]:
    states = noncrossing_states(width)
    ordinary = len(states)
    size = ordinary + marks
    block_size = size * size
    variables = 2 * block_size

    translation = block_diagonal_with_marks(
        action_matrix(width, lambda state: rotate_state(state, 1)), marks
    )
    joins = tuple(
        block_diagonal_with_marks(
            action_matrix(width, lambda state, site=site: join_adjacent(state, site)),
            marks,
        )
        for site in range(width)
    )
    detach_zero, detach_one = zip(
        *(multimark_detach_jet(width, site, marks) for site in range(width))
    )

    hom_rows = q_independent_jet_constraints(
        joins, joins[1:] + joins[:1]
    ) + q_independent_jet_constraints((translation,), (translation,))
    raw_zero_detach = list(
        intertwiner_constraints(
            detach_zero, detach_zero[1:] + detach_zero[:1]
        )
    )
    hom_rows += _block_rows(raw_zero_detach, 0, block_size)
    hom_rows += first_jet_intertwiner_constraints(
        detach_zero,
        detach_one,
        detach_zero[1:] + detach_zero[:1],
        detach_one[1:] + detach_one[:1],
    )
    hom_rows += multimark_filtration_constraints(ordinary, marks)
    hom_rhs = [Fraction(0)] * len(hom_rows)

    g0, g1 = extended_multimark_gram_jet(states, marks)
    basis = nullspace_basis(g0)
    radical_dimension = len(basis[0]) if basis else 0
    endpoint = [1] * ordinary + [0] * marks
    endpoint_rows, endpoint_rhs = covector_normalization(endpoint, size)
    invariance_rows = radical_invariance_constraints(g0, basis, size)
    invariance_rhs = [Fraction(0)] * len(invariance_rows)
    endpoint_invariance_rows = endpoint_rows + invariance_rows
    endpoint_invariance_rhs = endpoint_rhs + invariance_rhs
    gram_rows = radical_gram_constraints(g1, basis, size)
    gram_rhs = [Fraction(0)] * len(gram_rows)

    ordinary_source = [0] * size
    ordinary_source[states.index(tuple(range(width)))] = 1
    sources = [ordinary_source]
    for mark in range(marks):
        source = [0] * size
        source[ordinary + mark] = 1
        sources.append(source)
    source_rows, source_rhs = source_normalization(sources, size)

    ladders = [
        ("affine_q_jet", hom_rows, hom_rhs),
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

    first_empty = None
    previous = "affine_q_jet"
    for name, rows, rhs in (
        (
            "endpoint_radical_normalized",
            endpoint_invariance_rows,
            endpoint_invariance_rhs,
        ),
        ("gram_self_adjoint", gram_rows, gram_rhs),
        ("source_normalized", source_rows, source_rhs),
    ):
        restriction = _restriction_witness(solutions[previous], rows, rhs)
        if restriction is not None and not restriction["consistent"]:
            first_empty = {
                "from_stage": previous,
                "added_stage": name,
                **restriction,
            }
            break
        previous = name

    final = solutions["source_normalized"]
    if not final["consistent"]:
        decision = "lower_bound_not_attained"
    elif final["dimension"] == 0:
        decision = "lower_bound_attained_unique"
    else:
        decision = "lower_bound_attained_nonunique"

    canonical_matrix = [
        [Fraction(value) for value in row] for row in translation
    ]
    canonical_skew = restricted_skew_residual(g1, basis, canonical_matrix)
    skew_rank = matrix_residual_rank(canonical_skew)
    further_marks = skew_rank // 2
    raised_total = marks + further_marks
    maximum_nonzero_marks = width

    velocity_sites = []
    for site, velocity in enumerate(detach_one):
        restricted = multiply(
            [[Fraction(value) for value in row] for row in velocity], basis
        )
        mark_block = restricted[ordinary:]
        velocity_sites.append(
            {
                "site": site,
                "full_velocity_rank_on_radical": matrix_residual_rank(restricted),
                "mark_emission_rank_on_radical": matrix_residual_rank(mark_block),
                "mark_emission_block": encode_matrix(mark_block),
            }
        )

    return {
        "width": width,
        "marks": marks,
        "ordinary_dimension": ordinary,
        "extended_dimension": size,
        "g0_rank": matrix_residual_rank(g0),
        "radical_dimension": radical_dimension,
        "g1_radical_rank": matrix_residual_rank(
            multiply(transpose(basis), multiply(g1, basis))
        ),
        "mark_velocity_gate": {
            "sites": velocity_sites,
            "minimum_mark_emission_rank": min(
                row["mark_emission_rank_on_radical"] for row in velocity_sites
            ),
        },
        "stages": summaries,
        "first_empty_restriction": first_empty,
        "decision": decision,
        "final_solution": encode_final_solution(final, size, canonical),
        "canonical_translation": {
            "restricted_gram_skew_rank": skew_rank,
            "restricted_gram_skew": encode_matrix(canonical_skew),
        },
        "raised_bound_if_empty": {
            "further_independent_scalar_marks": further_marks,
            "total_scalar_marks": raised_total,
            "minimum_extended_dimension": ordinary + raised_total,
            "maximum_nonzero_falling_factorial_marks": maximum_nonzero_marks,
            "exceeds_available_independent_marks": raised_total > maximum_nonzero_marks,
            "scope": "falling-factorial scalar endpoint-mark quotient family",
        }
        if not final["consistent"]
        else None,
    }


def build_result() -> dict[str, Any]:
    protocol = json.loads(PROTOCOL.read_text(encoding="utf-8"))
    if protocol["status"] != "frozen_before_exact_intersection":
        raise AssertionError("protocol is not frozen")
    cases = [
        case_result(row["width"], row["marks"])
        for row in protocol["state_space"]["frozen_cases"]
    ]
    return {
        "schema": SCHEMA,
        "status": "exact_rational_minimal_multimark_certificate",
        "issues": [333, 321, 370],
        "protocol": {
            "path": str(PROTOCOL.relative_to(ROOT)),
            "sha256": hashlib.sha256(PROTOCOL.read_bytes()).hexdigest(),
            "commit": "b4cf829",
        },
        "cases": cases,
        "decision_summary": {
            f"w{row['width']}_k{row['marks']}": row["decision"] for row in cases
        },
        "global_decision": (
            "Neither lower bound from 081a5ed is attained by the frozen falling-factorial "
            "marks. Width three with two marks has ladder 10 -> 8 -> empty; width four "
            "with three marks has 17 -> 15 -> empty. In both cases the Gram restriction "
            "has coefficient rank zero and augmented rank one, so every surviving affine "
            "modulus is invisible to the obstruction. The family-specific bounds rise to "
            "three total marks at width three and five at width four. The width-four bound "
            "exceeds its four nonzero falling-factorial responses, exhausting that scalar "
            "family without testing another case."
        ),
        "deduplication": protocol["deduplication"],
        "claim_boundary": protocol["claim_boundary"],
    }


def render_markdown(result: dict[str, Any]) -> str:
    lines = [
        "# P333/P321/P370 minimal multimark lower-bound gate",
        "",
        result["global_decision"],
        "",
        "| width | marks | dim W | rank G0 | dim radical | mark velocity rank | affine jet | + endpoint/radical | + Gram | + source | decision |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|",
    ]
    for row in result["cases"]:
        stages = row["stages"]
        def dim(name: str) -> str:
            value = stages[name]["affine_tangent_dimension"]
            return "empty" if value is None else str(value)
        lines.append(
            f"| {row['width']} | {row['marks']} | {row['extended_dimension']} | "
            f"{row['g0_rank']} | {row['radical_dimension']} | "
            f"{row['mark_velocity_gate']['minimum_mark_emission_rank']} | "
            f"{dim('affine_q_jet')} | {dim('endpoint_radical_normalized')} | "
            f"{dim('gram_self_adjoint')} | {dim('source_normalized')} | `{row['decision']}` |"
        )
    lines.extend(["", "## Exact interpretation", ""])
    for row in result["cases"]:
        lines.append(
            f"- Width {row['width']} with exactly {row['marks']} marks: "
            f"`{row['decision']}`; canonical Gram-skew rank "
            f"{row['canonical_translation']['restricted_gram_skew_rank']}."
        )
        if row["final_solution"] is not None:
            solution = row["final_solution"]
            lines.append(
                f"  Final tangent dimension {solution['affine_tangent_dimension']}; "
                f"contains canonical `Tbar,0`={solution['contains_canonical_X0_Tbar_V0']}; "
                f"all tangent velocities zero={solution['all_tangent_velocities_are_zero']}."
            )
        else:
            witness = row["first_empty_restriction"]
            raised = row["raised_bound_if_empty"]
            lines.append(
                f"  First empty restriction `{witness['from_stage']} -> "
                f"{witness['added_stage']}` with `{witness['inconsistency_witness']['identity']}`. "
                f"The family-specific lower bound rises to {raised['total_scalar_marks']} total marks; "
                f"available nonzero falling-factorial marks={raised['maximum_nonzero_falling_factorial_marks']}, "
                f"family exhausted={raised['exceeds_available_independent_marks']}. The larger case was not tested."
            )
    lines.extend(["", "## Boundary", ""])
    lines.extend(f"- {item}" for item in result["claim_boundary"])
    lines.append("")
    return "\n".join(lines)


def render_card(result: dict[str, Any]) -> str:
    decisions = ", ".join(
        f"w{row['width']}/k{row['marks']}={row['decision']}"
        for row in result["cases"]
    )
    return "\n".join(
        [
            "# Scientific card: exact attainment of the scalar-mark lower bound",
            "",
            f"- **Mechanism space changed:** `{decisions}` for the only two frozen cases.",
            "- **Construction:** marks are the first k falling-factorial endpoint responses `(b)_r`; no mark direction was selected from the outcomes.",
            "- **Result:** neither frozen lower bound is attained. The exact ladders are `10 -> 8 -> empty` at width 3/k=2 and `17 -> 15 -> empty` at width 4/k=3; both Gram restrictions have coefficient rank 0 and augmented rank 1.",
            "- **Raised boundary:** width 3 rises to at least 3 total marks. Width 4 rises to at least 5, exceeding the 4 nonzero falling-factorial block-count responses and thereby exhausting this scalar family.",
            "- **Observer/sector/source:** crossed-to-trivial scalar marked endpoint quotient | extended Q=1 radical | all-singleton plus every frozen mark source.",
            "- **Not proved:** no rooted-cluster transfer module, LCFT/Jordan identity or formal Gaussian-K identification.",
            "- **Dependency group:** exact continuation of `081a5ed` and `ba3135e`, not an independent vote.",
            "- **Stop rule:** no mark count beyond the frozen two cases was tested; an empty result only raises the family-specific lower bound.",
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
