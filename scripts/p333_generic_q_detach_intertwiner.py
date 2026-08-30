#!/usr/bin/env python3
"""Exact generic-Q detach lift of the P333 affine-intertwiner gate."""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
from typing import Any, Sequence

from noncrossing_connectivity_codec import canonical_rgs, noncrossing_states
from p321_graded_closure_extension import intertwiner_constraints
from p321_homology_trace_certificate import (
    Matrix,
    action_matrix,
    exact_rank,
    join_adjacent,
    rotate_state,
)
from p333_gram_source_intertwiner import (
    affine_inconsistency_witness,
    encode_fraction,
    encode_matrix,
    endpoint_constraints,
    first_jet_radical_gram,
    gram_constraints,
    join_block_count,
    matrix_residual_rank,
    multiply,
    radical_action,
    rref_solve,
    source_constraints,
    subtract,
    transpose,
)


SCHEMA = "matching-one/p333-generic-q-detach-intertwiner/v1"
ROOT = Path(__file__).resolve().parents[1]
PROTOCOL = ROOT / "analysis/p333_generic_q_detach_intertwiner_protocol.json"
DEFAULT_RESULT = ROOT / "results/p333-generic-q-detach-intertwiner/latest.json"
OLD_RESULT = ROOT / "results/p333-gram-source-intertwiner/latest.json"


def detach_state(state: Sequence[int], site: int) -> tuple[int, ...]:
    """Split ``site`` into a singleton; an existing singleton is unchanged."""

    state = tuple(state)
    label = state[site]
    if state.count(label) == 1:
        return state
    fresh = max(state) + 1
    output = list(state)
    output[site] = fresh
    return canonical_rgs(output)


def singleton_projector(states: Sequence[Sequence[int]], site: int) -> Matrix:
    """Coefficient of epsilon in D_i(1+epsilon)."""

    size = len(states)
    return tuple(
        tuple(
            int(row == column and tuple(states[column]).count(states[column][site]) == 1)
            for column in range(size)
        )
        for row in range(size)
    )


def detach_jet(width: int, site: int) -> tuple[Matrix, Matrix]:
    """Return D_i(1) and dD_i/dQ at Q=1 in the frozen normalization."""

    states = noncrossing_states(width)
    constant = action_matrix(width, lambda state: detach_state(state, site))
    velocity = singleton_projector(states, site)
    return constant, velocity


def _block_rows(
    rows: Sequence[Sequence[int | Fraction]],
    block: int,
    variables_per_block: int,
) -> list[list[Fraction]]:
    output = []
    for row in rows:
        embedded = [Fraction(0)] * (2 * variables_per_block)
        offset = block * variables_per_block
        for index, value in enumerate(row):
            embedded[offset + index] = Fraction(value)
        output.append(embedded)
    return output


def detach_velocity_constraints(
    constants: Sequence[Matrix], velocities: Sequence[Matrix]
) -> list[list[Fraction]]:
    """Coefficient of epsilon in X D_i = D_(i+1) X."""

    size = len(constants[0])
    block_size = size * size
    rows: list[list[Fraction]] = []
    for site, (left, left_velocity) in enumerate(zip(constants, velocities)):
        right = constants[(site + 1) % len(constants)]
        right_velocity = velocities[(site + 1) % len(velocities)]
        for output in range(size):
            for column in range(size):
                row = [Fraction(0)] * (2 * block_size)
                for pivot in range(size):
                    # V D_i(1) - D_(i+1)(1) V
                    row[block_size + output * size + pivot] += left[pivot][column]
                    row[block_size + pivot * size + column] -= right[output][pivot]
                    # X0 P_i - P_(i+1) X0
                    row[output * size + pivot] += left_velocity[pivot][column]
                    row[pivot * size + column] -= right_velocity[output][pivot]
                rows.append(row)
    return rows


def full_gram_jet_constraints(
    g1: Sequence[Sequence[Fraction]], size: int
) -> list[list[Fraction]]:
    """Unprojected epsilon coefficient of G X = X^T G."""

    block_size = size * size
    rows: list[list[Fraction]] = []
    # G0 is the all-ones matrix at Q=1.  Off-diagonal equations suffice.
    for left in range(size):
        for right in range(left + 1, size):
            row = [Fraction(0)] * (2 * block_size)
            for pivot in range(size):
                # (G1 X0)_(left,right) - (X0^T G1)_(left,right)
                row[pivot * size + right] += g1[left][pivot]
                row[pivot * size + left] -= g1[pivot][right]
                # (G0 V)_(left,right) - (V^T G0)_(left,right)
                row[block_size + pivot * size + right] += 1
                row[block_size + pivot * size + left] -= 1
            rows.append(row)
    return rows


def _stage_summary(
    rows: Sequence[Sequence[Fraction]],
    rhs: Sequence[int | Fraction],
    variables: int,
    canonical: Sequence[Fraction],
) -> tuple[dict[str, Any], dict[str, Any]]:
    solved = rref_solve(rows, rhs, variables)
    residual = [
        sum(Fraction(value) * coordinate for value, coordinate in zip(row, canonical))
        - Fraction(target)
        for row, target in zip(rows, rhs)
    ]
    return (
        {
            "equations": len(rows),
            "consistent": solved["consistent"],
            "rank": solved["rank"],
            "affine_tangent_dimension": solved["dimension"],
            "canonical_X0_T_V0_nonzero_equations": sum(bool(value) for value in residual),
        },
        solved,
    )


def _joint_parameterization(solution: dict[str, Any], size: int) -> dict[str, Any] | None:
    if not solution["consistent"]:
        return None
    block_size = size * size
    particular = solution["particular"]
    return {
        "formula": "(X0,V)=(X0*,V*)+sum_a theta_a (Y0_a,Y1_a)",
        "particular_X0": encode_matrix(
            [particular[row * size : (row + 1) * size] for row in range(size)]
        ),
        "particular_V": encode_matrix(
            [
                particular[block_size + row * size : block_size + (row + 1) * size]
                for row in range(size)
            ]
        ),
        "tangent_count": len(solution["nullspace"]),
    }


def _restriction_witness(
    base_solution: dict[str, Any],
    extra_rows: Sequence[Sequence[Fraction]],
    extra_rhs: Sequence[int | Fraction],
) -> dict[str, Any] | None:
    if not base_solution["consistent"]:
        return None
    particular = base_solution["particular"]
    tangents = base_solution["nullspace"]
    coefficients = [
        [sum(value * coordinate for value, coordinate in zip(row, tangent)) for tangent in tangents]
        for row in extra_rows
    ]
    rhs = [
        Fraction(target)
        - sum(value * coordinate for value, coordinate in zip(row, particular))
        for row, target in zip(extra_rows, extra_rhs)
    ]
    restricted = rref_solve(coefficients, rhs, len(tangents))
    witness = (
        affine_inconsistency_witness(coefficients, rhs)
        if not restricted["consistent"]
        else None
    )
    return {
        "parameter_count_before_restriction": len(tangents),
        "coefficient_rank": restricted["rank"],
        "augmented_rank": restricted["augmented_rank"],
        "consistent": restricted["consistent"],
        "remaining_dimension": restricted["dimension"],
        "inconsistency_witness": witness,
    }


def width_result(width: int, old_row: dict[str, Any]) -> dict[str, Any]:
    states = noncrossing_states(width)
    size = len(states)
    block_size = size * size
    variables = 2 * block_size
    translation = action_matrix(width, lambda state: rotate_state(state, 1))
    joins = tuple(
        action_matrix(width, lambda state, site=site: join_adjacent(state, site))
        for site in range(width)
    )
    detach_constants, detach_velocities = zip(
        *(detach_jet(width, site) for site in range(width))
    )

    join_translation = list(
        intertwiner_constraints(joins, joins[1:] + joins[:1])
        + intertwiner_constraints((translation,), (translation,))
    )
    detach_zero = list(
        intertwiner_constraints(
            detach_constants,
            detach_constants[1:] + detach_constants[:1],
        )
    )
    hom_rows = (
        _block_rows(join_translation + detach_zero, 0, block_size)
        + _block_rows(join_translation, 1, block_size)
        + detach_velocity_constraints(detach_constants, detach_velocities)
    )
    hom_rhs = [0] * len(hom_rows)

    endpoint_zero_rows, endpoint_zero_rhs = endpoint_constraints(size)
    endpoint_velocity_rows, _ = endpoint_constraints(size)
    endpoint_rows = _block_rows(endpoint_zero_rows, 0, block_size) + _block_rows(
        endpoint_velocity_rows, 1, block_size
    )
    endpoint_rhs = endpoint_zero_rhs + [0] * len(endpoint_velocity_rows)

    h = first_jet_radical_gram(states)
    gram_zero_rows, gram_rhs = gram_constraints(h, size)
    gram_rows = _block_rows(gram_zero_rows, 0, block_size)

    source_state = tuple(range(width))
    source_index = states.index(source_state)
    source_zero_rows, source_zero_rhs = source_constraints(size, source_index)
    source_velocity_rows, _ = source_constraints(size, source_index)
    source_rows = _block_rows(source_zero_rows, 0, block_size) + _block_rows(
        source_velocity_rows, 1, block_size
    )
    source_rhs = source_zero_rhs + [0] * len(source_velocity_rows)

    g1 = [
        [Fraction(join_block_count(left, right)) for right in states]
        for left in states
    ]
    full_gram_rows = full_gram_jet_constraints(g1, size)
    full_gram_rhs = [0] * len(full_gram_rows)

    canonical = [Fraction(value) for row in translation for value in row] + [
        Fraction(0)
    ] * block_size
    ladders = [
        ("generic_q_affine_hom_jet", hom_rows, hom_rhs),
        (
            "endpoint_normalized",
            hom_rows + endpoint_rows,
            hom_rhs + endpoint_rhs,
        ),
        (
            "gram_radical_self_adjoint",
            hom_rows + endpoint_rows + gram_rows,
            hom_rhs + endpoint_rhs + gram_rhs,
        ),
        (
            "source_normalized",
            hom_rows + endpoint_rows + gram_rows + source_rows,
            hom_rhs + endpoint_rhs + gram_rhs + source_rhs,
        ),
    ]
    summaries: dict[str, Any] = {}
    solutions: dict[str, Any] = {}
    for name, rows, rhs in ladders:
        summaries[name], solutions[name] = _stage_summary(
            rows, rhs, variables, canonical
        )

    full_rows = hom_rows + endpoint_rows + gram_rows + full_gram_rows + source_rows
    full_rhs = hom_rhs + endpoint_rhs + gram_rhs + full_gram_rhs + source_rhs
    full_summary, full_solution = _stage_summary(full_rows, full_rhs, variables, canonical)

    gram_solution = solutions["gram_radical_self_adjoint"]
    endpoint_solution = solutions["endpoint_normalized"]
    gram_restriction = _restriction_witness(
        endpoint_solution, gram_rows, gram_rhs
    )
    source_restriction = _restriction_witness(gram_solution, source_rows, source_rhs)
    final_solution = solutions["source_normalized"]
    if not final_solution["consistent"]:
        decision = "remains_empty"
    elif final_solution["dimension"] == 0:
        decision = "reopened_unique"
    else:
        decision = "reopened_moduli"

    velocity_overlap = None
    if final_solution["consistent"]:
        velocity = [
            final_solution["particular"][
                block_size + row * size : block_size + (row + 1) * size
            ]
            for row in range(size)
        ]
        induced = radical_action(velocity)
        velocity_overlap = {
            "particular_induced_radical_velocity": encode_matrix(induced),
            "particular_rank": matrix_residual_rank(induced),
            "all_final_tangents_have_zero_radical_velocity": all(
                matrix_residual_rank(
                    radical_action(
                        [
                            tangent[
                                block_size + row * size : block_size + (row + 1) * size
                            ]
                            for row in range(size)
                        ]
                    )
                )
                == 0
                for tangent in final_solution["nullspace"]
            ),
        }

    detach_certificate = []
    for site, (constant, velocity) in enumerate(
        zip(detach_constants, detach_velocities)
    ):
        detach_certificate.append(
            {
                "site": site,
                "D_at_Q1": [list(row) for row in constant],
                "dD_dQ_at_Q1": [list(row) for row in velocity],
                "D_at_Q1_rank": exact_rank(constant),
                "dD_dQ_rank": exact_rank(velocity),
            }
        )

    inherited_no_go = (
        old_row["decision"] == "empty_intersection"
        and not final_solution["consistent"]
    )
    endpoint_uniquely_canonical = (
        solutions["endpoint_normalized"]["consistent"]
        and solutions["endpoint_normalized"]["dimension"] == 0
        and summaries["endpoint_normalized"][
            "canonical_X0_T_V0_nonzero_equations"
        ]
        == 0
    )
    return {
        "width": width,
        "module_dimension": size,
        "joint_variable_count": variables,
        "source_state": list(source_state),
        "source_index": source_index,
        "detach_polynomial_certificate": detach_certificate,
        "stages": summaries,
        "secondary_full_gram_source": full_summary,
        "gram_restriction_on_endpoint_moduli": gram_restriction,
        "source_restriction_on_gram_moduli": source_restriction,
        "final_parameterization": _joint_parameterization(final_solution, size),
        "velocity_radical_overlap": velocity_overlap,
        "decision": decision,
        "endpoint_uniquely_selects_X0_T_V0": endpoint_uniquely_canonical,
        "zeroth_order_projection": {
            "join_only_96df7c8_decision": old_row["decision"],
            "inherited_no_go": inherited_no_go,
            "logic": (
                "Every full-Q jet solution projects to an X0 satisfying the 96df7c8 "
                "join/endpoint/radical-Gram/source equations; therefore an empty old "
                "intersection cannot be repaired by V."
            ),
        },
    }


def build_result() -> dict[str, Any]:
    protocol = json.loads(PROTOCOL.read_text(encoding="utf-8"))
    if protocol["status"] != "frozen_before_exact_intersection":
        raise AssertionError("protocol is not frozen")
    old = json.loads(OLD_RESULT.read_text(encoding="utf-8"))
    old_by_width = {row["width"]: row for row in old["widths"]}
    widths = [
        width_result(width, old_by_width[width])
        for width in protocol["representation"]["widths"]
    ]
    nondegenerate = [row for row in widths if row["width"] >= 3]
    return {
        "schema": SCHEMA,
        "status": "exact_rational_first_jet_certificate",
        "issues": [333, 321, 370],
        "protocol": {
            "path": str(PROTOCOL.relative_to(ROOT)),
            "sha256": hashlib.sha256(PROTOCOL.read_bytes()).hexdigest(),
            "commit": "6592d11",
        },
        "widths": widths,
        "decision_summary": {str(row["width"]): row["decision"] for row in widths},
        "global_decision": (
            "The standard generic-Q detach first makes the algebra sharper: its affine-Hom "
            "jet has two moduli at every tested width, and endpoint normalization uniquely "
            "selects X0=T,V=0. That selected line survives at translation-degenerate width "
            "two but fails the radical Gram equation at widths three and four, so the full "
            "physical intersection does not reopen. This is also a zeroth-order obstruction: "
            "every full-Q jet projects to the already inconsistent join/endpoint/radical-"
            "Gram/source system, hence no detach velocity V can repair it. A physical "
            "confluence requires a larger marked or direct-sum module, not merely the "
            "missing scalar loop weight."
        ),
        "nondegenerate_widths_all_remain_empty": all(
            row["decision"] == "remains_empty" for row in nondegenerate
        ),
        "endpoint_uniquely_selects_canonical_all_widths": all(
            row["endpoint_uniquely_selects_X0_T_V0"] for row in widths
        ),
        "deduplication": protocol["deduplication"],
        "claim_boundary": protocol["claim_boundary"],
    }


def render_markdown(result: dict[str, Any]) -> str:
    lines = [
        "# P333/P321/P370 generic-Q detach intertwiner gate",
        "",
        "The frozen generator is the standard FK detach: split a non-singleton site with coefficient 1; detaching an existing singleton multiplies the state by `Q`. The exact calculation uses `Q=1+epsilon` and `X=X0+epsilon V` modulo `epsilon^2`.",
        "",
        "| width | dim V_w | affine Hom jet | + endpoint | + radical Gram | + source | full Gram + source | decision |",
        "|---:|---:|---:|---:|---:|---:|---:|---|",
    ]
    for row in result["widths"]:
        stages = row["stages"]
        def dim(name: str) -> str:
            value = stages[name]["affine_tangent_dimension"]
            return "empty" if value is None else str(value)
        full = row["secondary_full_gram_source"]["affine_tangent_dimension"]
        lines.append(
            f"| {row['width']} | {row['module_dimension']} | "
            f"{dim('generic_q_affine_hom_jet')} | {dim('endpoint_normalized')} | "
            f"{dim('gram_radical_self_adjoint')} | {dim('source_normalized')} | "
            f"{'empty' if full is None else full} | `{row['decision']}` |"
        )
    lines.extend(["", "## Exact decision", "", result["global_decision"], ""])
    for row in result["widths"]:
        overlap = row["velocity_radical_overlap"]
        if overlap is None:
            detail = "no final velocity exists"
        else:
            detail = (
                f"particular radical-velocity rank {overlap['particular_rank']}; "
                f"all tangent radical velocities zero={overlap['all_final_tangents_have_zero_radical_velocity']}"
            )
        lines.append(
            f"- Width {row['width']}: `{row['decision']}`; {detail}; "
            f"inherited zeroth-order no-go={row['zeroth_order_projection']['inherited_no_go']}."
        )
    lines.extend(
        [
            "",
            "The stronger unprojected first-jet Gram equation is reported only as a secondary check. It cannot rescue a failure of the frozen radical condition.",
            "",
            "## Boundary",
            "",
        ]
    )
    lines.extend(f"- {item}" for item in result["claim_boundary"])
    lines.append("")
    return "\n".join(lines)


def render_card(result: dict[str, Any]) -> str:
    dims = ", ".join(
        f"w{row['width']}={row['decision']}" for row in result["widths"]
    )
    return "\n".join(
        [
            "# Scientific card: minimal generic-Q detach does not repair the confluence",
            "",
            f"- **Mechanism space changed:** `{dims}` after adding the standard loop-weight detach and its exact Q velocity. Before Gram compatibility, endpoint normalization uniquely selects `X0=T,V=0` at all three widths.",
            f"- **Result:** {result['global_decision']}",
            "- **Observer/sector/source:** crossed-to-trivial affine connectivity map | Q=1 first-jet radical | all-singleton source.",
            "- **What this rules out:** scalar loop fugacity as the missing datum that could repair the width-3/4 join-only Gram/source obstruction; it selects the canonical translation but cannot make that line Gram-compatible.",
            "- **What it does not prove:** no LCFT/Jordan identification, no exclusion of marked-cluster or direct-sum extensions, and no relation to the formal K in PR #393.",
            "- **Dependency group:** exact finite-width connectivity representation; this is the declared continuation of `96df7c8`, not a new independent vote.",
            "- **Next discriminant:** add one physically named mark/closure block whose Q derivative has a nonzero radical projection; another scalar Q weight cannot change the zeroth-order no-go.",
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
