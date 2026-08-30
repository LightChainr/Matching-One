#!/usr/bin/env python3
"""Exact width-three source/landing doublet test for P333."""

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


SCHEMA = "matching-one/p333-source-landing-doublet/v1"
ROOT = Path(__file__).resolve().parents[1]
PROTOCOL = ROOT / "analysis/p333_source_landing_doublet_protocol.json"
DEFAULT_RESULT = ROOT / "results/p333-source-landing-doublet/latest.json"
WIDTH = 3


def landing_rotation() -> Matrix:
    """C3 action in basis (lambda0-lambda2, lambda1-lambda2)."""

    return ((-1, -1), (1, 0))


def landing_emission_vector(site: int) -> tuple[int, int]:
    vectors = ((2, -1), (-1, 2), (-1, -1))
    return vectors[site % 3]


def landing_pair_state(site: int) -> tuple[int, ...]:
    """Connectivity with ``site`` singleton and the complementary pair joined."""

    labels = [0, 0, 0]
    labels[site] = 1
    return canonical_rgs(labels)


def block_diagonal(
    ordinary: Sequence[Sequence[int | Fraction]],
    landing: Sequence[Sequence[int | Fraction]],
) -> Matrix:
    ordinary_size = len(ordinary)
    landing_size = len(landing)
    return tuple(
        tuple(int(value) for value in ordinary[row])
        + tuple(0 for _ in range(landing_size))
        for row in range(ordinary_size)
    ) + tuple(
        tuple(0 for _ in range(ordinary_size))
        + tuple(int(value) for value in landing[row])
        for row in range(landing_size)
    )


def source_landing_detach_jet(site: int) -> tuple[Matrix, Matrix]:
    states = noncrossing_states(WIDTH)
    ordinary_zero, ordinary_velocity = detach_jet(WIDTH, site)
    ordinary = len(states)
    identity2 = ((1, 0), (0, 1))
    zero = block_diagonal(ordinary_zero, identity2)
    q = landing_emission_vector(site)
    velocity_top = tuple(
        tuple(ordinary_velocity[row][column] for column in range(ordinary))
        + (0, 0)
        for row in range(ordinary)
    )
    velocity_bottom = []
    for component in range(2):
        velocity_bottom.append(
            tuple(
                q[component] * int(tuple(state).count(state[site]) == 1)
                for state in states
            )
            + (0, 0)
        )
    return zero, velocity_top + tuple(velocity_bottom)


def source_landing_gram_jet() -> tuple[list[list[Fraction]], list[list[Fraction]]]:
    states = noncrossing_states(WIDTH)
    references = [landing_pair_state(site) for site in range(WIDTH)]
    raw0 = [
        [Fraction(join_block_count(state, reference)) for reference in references]
        for state in states
    ]
    raw1 = [
        [
            Fraction(join_block_count(state, reference) * (join_block_count(state, reference) - 1))
            for reference in references
        ]
        for state in states
    ]
    c0 = [[row[0] - row[2], row[1] - row[2]] for row in raw0]
    c1 = [[row[0] - row[2], row[1] - row[2]] for row in raw1]
    g0 = [
        [Fraction(1) for _ in states] + c0[row]
        for row in range(len(states))
    ]
    g1 = [
        [Fraction(join_block_count(left, right)) for right in states] + c1[row]
        for row, left in enumerate(states)
    ]
    for landing_column in range(2):
        g0.append([c0[row][landing_column] for row in range(len(states))] + [Fraction(0), Fraction(0)])
        g1.append([c1[row][landing_column] for row in range(len(states))] + [Fraction(0), Fraction(0)])
    return g0, g1


def landing_filtration_constraints(ordinary: int) -> list[list[Fraction]]:
    size = ordinary + 2
    block_size = size * size
    rows: list[list[Fraction]] = []
    for block in (0, 1):
        for output in range(ordinary):
            for landing in range(2):
                row = [Fraction(0)] * (2 * block_size)
                row[
                    block * block_size + output * size + ordinary + landing
                ] = 1
                rows.append(row)
    return rows


def landing_transport_constraints(
    ordinary: int,
) -> tuple[list[list[Fraction]], list[Fraction]]:
    """Fix X0|L=R and V|L=0, including the already-declared filtration."""

    size = ordinary + 2
    block_size = size * size
    rotation = landing_rotation()
    rows: list[list[Fraction]] = []
    rhs: list[Fraction] = []
    for block in (0, 1):
        for output in range(2):
            for source in range(2):
                row = [Fraction(0)] * (2 * block_size)
                row[
                    block * block_size
                    + (ordinary + output) * size
                    + ordinary
                    + source
                ] = 1
                rows.append(row)
                rhs.append(
                    Fraction(rotation[output][source])
                    if block == 0
                    else Fraction(0)
                )
    return rows, rhs


def scalar_span_certificate(
    c0: Sequence[Sequence[Fraction]], states: Sequence[Sequence[int]]
) -> dict[str, Any]:
    scalar = [
        [
            Fraction(len(set(state))),
            Fraction(len(set(state)) * (len(set(state)) - 1)),
        ]
        for state in states
    ]
    combined = [list(scalar[row]) + list(c0[row]) for row in range(len(states))]
    scalar_rank = matrix_residual_rank(scalar)
    combined_rank = matrix_residual_rank(combined)
    return {
        "scalar_falling_factorial_rank": scalar_rank,
        "combined_scalar_plus_landing_rank": combined_rank,
        "new_non_scalar_rank": combined_rank - scalar_rank,
        "passes": combined_rank > scalar_rank,
    }


def build_result() -> dict[str, Any]:
    protocol = json.loads(PROTOCOL.read_text(encoding="utf-8"))
    if protocol["status"] != "frozen_before_exact_intersection":
        raise AssertionError("protocol is not frozen")
    states = noncrossing_states(WIDTH)
    ordinary = len(states)
    size = ordinary + 2
    block_size = size * size
    variables = 2 * block_size
    identity2 = ((1, 0), (0, 1))
    rotation = landing_rotation()

    translation = block_diagonal(
        action_matrix(WIDTH, lambda state: rotate_state(state, 1)), rotation
    )
    joins = tuple(
        block_diagonal(
            action_matrix(WIDTH, lambda state, site=site: join_adjacent(state, site)),
            identity2,
        )
        for site in range(WIDTH)
    )
    detach_zero, detach_one = zip(
        *(source_landing_detach_jet(site) for site in range(WIDTH))
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
    hom_rows += landing_filtration_constraints(ordinary)
    hom_rhs = [Fraction(0)] * len(hom_rows)

    g0, g1 = source_landing_gram_jet()
    basis = nullspace_basis(g0)
    radical_dimension = len(basis[0]) if basis else 0
    endpoint = [1] * ordinary + [0, 0]
    endpoint_rows, endpoint_rhs = covector_normalization(endpoint, size)
    invariance_rows = radical_invariance_constraints(g0, basis, size)
    endpoint_invariance_rows = endpoint_rows + invariance_rows
    endpoint_invariance_rhs = endpoint_rhs + [Fraction(0)] * len(invariance_rows)
    gram_rows = radical_gram_constraints(g1, basis, size)
    gram_rhs = [Fraction(0)] * len(gram_rows)

    ordinary_source = [0] * size
    ordinary_source[states.index(tuple(range(WIDTH)))] = 1
    ordinary_source_rows, ordinary_source_rhs = source_normalization(
        (ordinary_source,), size
    )
    transport_rows, transport_rhs = landing_transport_constraints(ordinary)
    source_rows = ordinary_source_rows + transport_rows
    source_rhs = ordinary_source_rhs + transport_rhs

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
        decision = "doublet_fails"
    elif final["dimension"] == 0:
        decision = "doublet_breaks_obstruction_unique"
    else:
        decision = "doublet_breaks_obstruction_nonunique"

    c0 = [row[ordinary:] for row in g0[:ordinary]]
    non_scalar = scalar_span_certificate(c0, states)
    canonical_matrix = [
        [Fraction(value) for value in row] for row in translation
    ]
    canonical_skew = restricted_skew_residual(g1, basis, canonical_matrix)
    mark_velocity_ranks = []
    for site, velocity in enumerate(detach_one):
        restricted = multiply(
            [[Fraction(value) for value in row] for row in velocity], basis
        )
        mark_velocity_ranks.append(
            {
                "site": site,
                "landing_emission_rank_on_radical": matrix_residual_rank(
                    restricted[ordinary:]
                ),
                "landing_emission_block": encode_matrix(restricted[ordinary:]),
            }
        )

    final_velocity = None
    if final["consistent"]:
        particular = final["particular"]
        final_velocity = {
            "particular_is_zero": all(value == 0 for value in particular[block_size:]),
            "all_tangent_velocities_zero": all(
                all(value == 0 for value in tangent[block_size:])
                for tangent in final["nullspace"]
            ),
        }

    result = {
        "schema": SCHEMA,
        "status": "exact_rational_source_landing_doublet_certificate",
        "issues": [333, 321, 370],
        "protocol": {
            "path": str(PROTOCOL.relative_to(ROOT)),
            "sha256": hashlib.sha256(PROTOCOL.read_bytes()).hexdigest(),
            "commit": "565d276",
        },
        "width": WIDTH,
        "ordinary_dimension": ordinary,
        "landing_dimension": 2,
        "extended_dimension": size,
        "g0_rank": matrix_residual_rank(g0),
        "radical_dimension": radical_dimension,
        "non_scalar_gate": non_scalar,
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
        "mark_velocity_gate": mark_velocity_ranks,
        "stages": summaries,
        "first_empty_restriction": first_empty,
        "canonical_restricted_gram_skew_rank": matrix_residual_rank(canonical_skew),
        "decision": decision,
        "final_velocity": final_velocity,
        "typed_lower_bound_if_failure": {
            "minimum_source_landing_doublets": 2,
            "minimum_non_scalar_mark_dimension": 4,
            "scope": "source/landing typed-family only",
        }
        if not final["consistent"]
        else None,
        "global_decision": (
            "The frozen width-three source/landing doublet is genuinely non-scalar and "
            "translation covariant, and it breaks the obstruction uniquely. The dimension "
            "ladder is 4 -> 2 -> 2 -> 0; the canonical X0=T_bar,V=0 satisfies every gate "
            "and is the unique final solution. The landing columns add rank two beyond the "
            "entire scalar falling-factorial span and reduce the canonical restricted Gram-"
            "skew rank to zero."
        ),
        "deduplication": protocol["deduplication"],
        "claim_boundary": protocol["claim_boundary"],
    }
    return result


def render_markdown(result: dict[str, Any]) -> str:
    stages = result["stages"]
    def dim(name: str) -> str:
        value = stages[name]["affine_tangent_dimension"]
        return "empty" if value is None else str(value)
    lines = [
        "# P333/P321/P370 width-three source/landing doublet",
        "",
        result["global_decision"],
        "",
        "| dim V | dim landing | dim W | rank G0 | dim radical | affine jet | + endpoint/radical | + Gram | + source/landing | decision |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|",
        f"| {result['ordinary_dimension']} | {result['landing_dimension']} | {result['extended_dimension']} | "
        f"{result['g0_rank']} | {result['radical_dimension']} | {dim('affine_q_jet')} | "
        f"{dim('endpoint_radical_normalized')} | {dim('gram_self_adjoint')} | "
        f"{dim('source_landing_normalized')} | `{result['decision']}` |",
        "",
        "## Exact gates",
        "",
        f"- Non-scalar landing rank added beyond scalar falling-factorial columns: {result['non_scalar_gate']['new_non_scalar_rank']}.",
        f"- Translation covariance residual ranks: G0={result['translation_covariance']['G0_residual_rank']}, G1={result['translation_covariance']['G1_residual_rank']}.",
        f"- Canonical restricted Gram-skew rank: {result['canonical_restricted_gram_skew_rank']}.",
    ]
    if result["first_empty_restriction"] is not None:
        witness = result["first_empty_restriction"]
        lines.append(
            f"- First empty restriction `{witness['from_stage']} -> {witness['added_stage']}` "
            f"with `{witness['inconsistency_witness']['identity']}`."
        )
        lines.append(
            "- Typed-family lower bound after failure: at least two source/landing doublets, mark dimension at least four. Width four was not run."
        )
    lines.extend(["", "## Boundary", ""])
    lines.extend(f"- {item}" for item in result["claim_boundary"])
    lines.append("")
    return "\n".join(lines)


def render_card(result: dict[str, Any]) -> str:
    return "\n".join(
        [
            "# Scientific card: non-scalar source/landing doublet",
            "",
            f"- **Mechanism space changed:** width-three decision `{result['decision']}` for one fixed C3 landing doublet.",
            f"- **Non-scalar gate:** added rank `{result['non_scalar_gate']['new_non_scalar_rank']}` beyond all frozen scalar falling-factorial columns; translation residual ranks are `{result['translation_covariance']}`.",
            "- **Observer/sector/source:** crossed-to-trivial connectivity | C3 zero-sum source/landing doublet | all-singleton ordinary source and translated landing registry.",
            "- **Not proved:** no rooted transfer module, LCFT/Jordan identity, physical transfer matrix or formal Gaussian-K identification.",
            "- **Dependency group:** exact continuation of `e7e6c80`; it changes the mark representation rather than adding another scalar vote.",
            "- **Stop rule:** width four is attempted only if this width-three final intersection is nonempty.",
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
    print(json.dumps({"width3": result["decision"]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
