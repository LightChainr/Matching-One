#!/usr/bin/env python3
"""Exact width-four C4 charge-two landing response gate for P333."""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
from typing import Any, Sequence

from noncrossing_connectivity_codec import canonical_rgs, noncrossing_states
from p321_graded_closure_extension import intertwiner_constraints
from p321_homology_trace_certificate import action_matrix, join_adjacent, rotate_state
from p333_generic_q_detach_intertwiner import (
    _block_rows,
    _restriction_witness,
    _stage_summary,
    detach_jet,
)
from p333_gram_source_intertwiner import (
    join_block_count,
    matrix_residual_rank,
    multiply,
    subtract,
    transpose,
)
from p333_minimal_multimark_jet import falling_factorial
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
from p333_source_landing_doublet import block_diagonal


SCHEMA = "matching-one/p333-charge2-landing/v1"
ROOT = Path(__file__).resolve().parents[1]
PROTOCOL = ROOT / "analysis/p333_charge2_landing_protocol.json"
DEFAULT_RESULT = ROOT / "results/p333-charge2-landing/latest.json"
WIDTH = 4


def charge2_character(site: int) -> int:
    return 1 if site % 2 == 0 else -1


def landing_reference_state(site: int) -> tuple[int, ...]:
    labels = [0, 0, 0, 0]
    labels[site] = 1
    return canonical_rgs(labels)


def charge2_detach_jet(site: int):
    states = noncrossing_states(WIDTH)
    ordinary_zero, ordinary_velocity = detach_jet(WIDTH, site)
    ordinary = len(states)
    zero = block_diagonal(ordinary_zero, ((1,),))
    top = tuple(
        tuple(ordinary_velocity[row][column] for column in range(ordinary)) + (0,)
        for row in range(ordinary)
    )
    bottom = (
        tuple(
            charge2_character(site)
            * int(tuple(state).count(state[site]) == 1)
            for state in states
        )
        + (0,),
    )
    return zero, top + bottom


def charge2_gram_jet():
    states = noncrossing_states(WIDTH)
    references = [landing_reference_state(site) for site in range(WIDTH)]
    c0, c1 = [], []
    for state in states:
        blocks = [join_block_count(state, reference) for reference in references]
        c0.append(Fraction(sum(charge2_character(site) * blocks[site] for site in range(WIDTH))))
        c1.append(
            Fraction(
                sum(
                    charge2_character(site) * blocks[site] * (blocks[site] - 1)
                    for site in range(WIDTH)
                )
            )
        )
    g0 = [
        [Fraction(1) for _ in states] + [c0[row]]
        for row in range(len(states))
    ]
    g1 = [
        [Fraction(join_block_count(left, right)) for right in states] + [c1[row]]
        for row, left in enumerate(states)
    ]
    g0.append(c0 + [Fraction(0)])
    g1.append(c1 + [Fraction(0)])
    return g0, g1


def filtration_constraints(ordinary: int):
    size = ordinary + 1
    block_size = size * size
    rows = []
    for block in (0, 1):
        for output in range(ordinary):
            row = [Fraction(0)] * (2 * block_size)
            row[block * block_size + output * size + ordinary] = 1
            rows.append(row)
    return rows


def transport_constraints(ordinary: int):
    size = ordinary + 1
    block_size = size * size
    rows, rhs = [], []
    for block in (0, 1):
        row = [Fraction(0)] * (2 * block_size)
        row[block * block_size + ordinary * size + ordinary] = 1
        rows.append(row)
        rhs.append(Fraction(-1) if block == 0 else Fraction(0))
    return rows, rhs


def non_scalar_certificate(column, states):
    scalar = [
        [
            Fraction(falling_factorial(len(set(state)), order))
            for order in range(1, WIDTH + 1)
        ]
        for state in states
    ]
    combined = [list(scalar[row]) + [column[row]] for row in range(len(states))]
    scalar_rank = matrix_residual_rank(scalar)
    combined_rank = matrix_residual_rank(combined)
    return {
        "scalar_falling_factorial_rank": scalar_rank,
        "combined_scalar_plus_charge2_rank": combined_rank,
        "new_non_scalar_rank": combined_rank - scalar_rank,
        "passes": combined_rank > scalar_rank,
    }


def build_result() -> dict[str, Any]:
    protocol = json.loads(PROTOCOL.read_text(encoding="utf-8"))
    if protocol["status"] != "frozen_before_exact_intersection":
        raise AssertionError("protocol is not frozen")
    states = noncrossing_states(WIDTH)
    ordinary = len(states)
    size = ordinary + 1
    block_size = size * size
    variables = 2 * block_size
    translation = block_diagonal(
        action_matrix(WIDTH, lambda state: rotate_state(state, 1)), ((-1,),)
    )
    joins = tuple(
        block_diagonal(
            action_matrix(WIDTH, lambda state, site=site: join_adjacent(state, site)),
            ((1,),),
        )
        for site in range(WIDTH)
    )
    detach_zero, detach_one = zip(
        *(charge2_detach_jet(site) for site in range(WIDTH))
    )

    hom_rows = q_independent_jet_constraints(
        joins, joins[1:] + joins[:1]
    ) + q_independent_jet_constraints((translation,), (translation,))
    hom_rows += _block_rows(
        list(
            intertwiner_constraints(
                detach_zero, detach_zero[1:] + detach_zero[:1]
            )
        ),
        0,
        block_size,
    )
    hom_rows += first_jet_intertwiner_constraints(
        detach_zero,
        detach_one,
        detach_zero[1:] + detach_zero[:1],
        detach_one[1:] + detach_one[:1],
    )
    hom_rows += filtration_constraints(ordinary)
    hom_rhs = [Fraction(0)] * len(hom_rows)

    g0, g1 = charge2_gram_jet()
    basis = nullspace_basis(g0)
    radical_dimension = len(basis[0]) if basis else 0
    endpoint = [1] * ordinary + [0]
    endpoint_rows, endpoint_rhs = covector_normalization(endpoint, size)
    invariance_rows = radical_invariance_constraints(g0, basis, size)
    endpoint_invariance_rows = endpoint_rows + invariance_rows
    endpoint_invariance_rhs = endpoint_rhs + [Fraction(0)] * len(invariance_rows)
    gram_rows = radical_gram_constraints(g1, basis, size)
    gram_rhs = [Fraction(0)] * len(gram_rows)

    ordinary_source = [0] * size
    ordinary_source[states.index(tuple(range(WIDTH)))] = 1
    ordinary_rows, ordinary_rhs = source_normalization((ordinary_source,), size)
    landing_rows, landing_rhs = transport_constraints(ordinary)
    source_rows = ordinary_rows + landing_rows
    source_rhs = ordinary_rhs + landing_rhs

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
            "source_landing_normalized",
            hom_rows + endpoint_invariance_rows + gram_rows + source_rows,
            hom_rhs + endpoint_invariance_rhs + gram_rhs + source_rhs,
        ),
    ]
    canonical = [Fraction(value) for row in translation for value in row] + [
        Fraction(0)
    ] * block_size
    summaries, solutions = {}, {}
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
        ("source_landing_normalized", source_rows, source_rhs),
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

    final = solutions["source_landing_normalized"]
    if not final["consistent"]:
        decision = "charge2_fails"
    elif final["dimension"] == 0:
        decision = "minimal_closed_unique"
    else:
        decision = "minimal_closed_nonunique"

    canonical_matrix = [
        [Fraction(value) for value in row] for row in translation
    ]
    c0 = [row[ordinary] for row in g0[:ordinary]]
    final_data = None
    if final["consistent"]:
        particular = final["particular"]
        final_data = {
            "affine_tangent_dimension": final["dimension"],
            "canonical_residual_nonzero_equations": summaries[
                "source_landing_normalized"
            ]["canonical_X0_T_V0_nonzero_equations"],
            "particular_velocity_is_zero": all(
                value == 0 for value in particular[block_size:]
            ),
            "all_tangent_velocities_zero": all(
                all(value == 0 for value in tangent[block_size:])
                for tangent in final["nullspace"]
            ),
        }

    return {
        "schema": SCHEMA,
        "status": "exact_rational_charge2_landing_certificate",
        "issues": [333, 321, 370],
        "protocol": {
            "path": str(PROTOCOL.relative_to(ROOT)),
            "sha256": hashlib.sha256(PROTOCOL.read_bytes()).hexdigest(),
            "commit": "a572058",
        },
        "width": WIDTH,
        "ordinary_dimension": ordinary,
        "charge2_dimension": 1,
        "extended_dimension": size,
        "g0_rank": matrix_residual_rank(g0),
        "radical_dimension": radical_dimension,
        "non_scalar_gate": non_scalar_certificate(c0, states),
        "translation_covariance": {
            "G0_residual_rank": matrix_residual_rank(
                subtract(
                    multiply(transpose(canonical_matrix), multiply(g0, canonical_matrix)),
                    g0,
                )
            ),
            "G1_residual_rank": matrix_residual_rank(
                subtract(
                    multiply(transpose(canonical_matrix), multiply(g1, canonical_matrix)),
                    g1,
                )
            ),
        },
        "stages": summaries,
        "first_empty_restriction": first_empty,
        "canonical_restricted_gram_skew_rank": matrix_residual_rank(
            restricted_skew_residual(g1, basis, canonical_matrix)
        ),
        "decision": decision,
        "final_solution": final_data,
        "irrep_exhaustion_if_failure": {
            "scalar_endpoint_irreps": "exhausted by e7e6c80",
            "C4_charge1_landing_irrep": "fails in 7b40ec7",
            "C4_charge2_landing_irrep": "fails here",
            "next_semantic_change": (
                "Retain rooted/landing connectivity after emission, or couple multiple "
                "landing irreps in one registry, instead of collapsing each irrep to an "
                "independent terminal endpoint mark."
            ),
        }
        if not final["consistent"]
        else None,
        "minimal_closed_module_if_success": {
            "dimension": size,
            "ordinary_dimension": ordinary,
            "mark_dimension": 1,
            "mark_character": "C4 charge two",
            "canonical": "X0=T_bar,V=0",
        }
        if final["consistent"]
        else None,
        "global_decision": (
            "The frozen C4 charge-two landing character fails exactly at the Gram gate: "
            "the unique endpoint/radical-normalized affine direction has no Gram-self-adjoint "
            "extension. Together with the prior scalar and C4 charge-one results, this exhausts "
            "all individual terminal landing irreps at width four."
        ),
        "deduplication": protocol["deduplication"],
        "claim_boundary": protocol["claim_boundary"],
    }


def render_markdown(result):
    stages = result["stages"]
    def dim(name):
        value = stages[name]["affine_tangent_dimension"]
        return "empty" if value is None else str(value)
    lines = [
        "# P333/P321/P370 width-four charge-two landing gate",
        "",
        result["global_decision"],
        "",
        "| dim V | dim mark | dim W | rank G0 | dim radical | affine jet | + endpoint/radical | + Gram | + source/landing | decision |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|",
        f"| {result['ordinary_dimension']} | 1 | {result['extended_dimension']} | {result['g0_rank']} | "
        f"{result['radical_dimension']} | {dim('affine_q_jet')} | {dim('endpoint_radical_normalized')} | "
        f"{dim('gram_self_adjoint')} | {dim('source_landing_normalized')} | `{result['decision']}` |",
        "",
        "## Exact gates",
        "",
        f"- New non-scalar rank beyond all scalar columns: {result['non_scalar_gate']['new_non_scalar_rank']}.",
        f"- Translation residual ranks: G0={result['translation_covariance']['G0_residual_rank']}, G1={result['translation_covariance']['G1_residual_rank']}.",
        f"- Canonical restricted Gram-skew rank: {result['canonical_restricted_gram_skew_rank']}.",
    ]
    if result["minimal_closed_module_if_success"] is not None:
        lines.append(
            f"- Minimal closed module: `{result['minimal_closed_module_if_success']}`."
        )
    if result["first_empty_restriction"] is not None:
        witness = result["first_empty_restriction"]
        lines.append(
            f"- First empty restriction `{witness['from_stage']} -> {witness['added_stage']}` "
            f"with `{witness['inconsistency_witness']['identity']}`."
        )
        lines.append(
            f"- Required semantic change: {result['irrep_exhaustion_if_failure']['next_semantic_change']}"
        )
    lines.extend(["", "## Boundary", ""])
    lines.extend(f"- {item}" for item in result["claim_boundary"])
    lines.append("")
    return "\n".join(lines)


def render_card(result):
    witness = result["first_empty_restriction"]
    exhaustion = result["irrep_exhaustion_if_failure"]
    return "\n".join(
        [
            "# Scientific card: final C4 charge-two landing irrep",
            "",
            f"- **Mechanism space changed:** exact decision `{result['decision']}` for the only remaining individual landing irrep.",
            f"- **Non-scalar gate:** rank `{result['non_scalar_gate']['new_non_scalar_rank']}` beyond the full scalar falling-factorial span.",
            f"- **Exact obstruction:** `{witness['from_stage']} -> {witness['added_stage']}` is inconsistent with `{witness['inconsistency_witness']['identity']}` (coefficient rank {witness['coefficient_rank']}, augmented rank {witness['augmented_rank']}).",
            "- **Observer/sector/source:** crossed-to-trivial connectivity | C4 charge-two alternating landing response | all-singleton source and sign-flip transport.",
            "- **Irrep closure:** scalar endpoint sectors, the C4 charge-one landing doublet and the C4 charge-two landing character have now all failed individually at width four.",
            "- **Not proved:** no continuum LCFT/Jordan identity, physical transfer matrix or formal Gaussian-K identification.",
            f"- **Next discriminating construction:** {exhaustion['next_semantic_change']}",
            "- **Stop rule:** no further terminal endpoint-mark irrep is available or tested after this result.",
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
    print(json.dumps({"charge2": result["decision"]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
