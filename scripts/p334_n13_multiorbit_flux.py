#!/usr/bin/env python3
"""Exact N=13 multi-orbit projective birth/exit flux for Issue #334."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from decimal import Decimal, localcontext
from fractions import Fraction
import json
from pathlib import Path
from typing import Iterable, Optional, Sequence

from integer_period_torus import gaussian_integer_torus
from projective_essential_birth_oracle import chi4, subset_marks


ComplexQ = tuple[Fraction, Fraction]
Vector = tuple[int, int]


P_REF = Fraction(59274605079, 100000000000)


def _active(mask: int, n: int) -> tuple[bool, ...]:
    return tuple(bool(mask & (1 << vertex)) for vertex in range(n))


def _qadd(left: ComplexQ, right: ComplexQ) -> ComplexQ:
    return left[0] + right[0], left[1] + right[1]


def _qsub(left: ComplexQ, right: ComplexQ) -> ComplexQ:
    return left[0] - right[0], left[1] - right[1]


def _qscale(value: ComplexQ, scalar: Fraction | int) -> ComplexQ:
    return value[0] * scalar, value[1] * scalar


def _qpayload(value: ComplexQ) -> dict[str, str]:
    return {"real": str(value[0]), "imag": str(value[1])}


def _decimal(value: Fraction, digits: int = 18) -> str:
    with localcontext() as context:
        context.prec = digits + 8
        return format(Decimal(value.numerator) / Decimal(value.denominator), f".{digits}g")


def _complex_decimal(value: ComplexQ) -> dict[str, str]:
    return {"real": _decimal(value[0]), "imag": _decimal(value[1])}


def orbit_label(line: Optional[Vector]) -> str:
    if line in {(1, 0), (0, 1)}:
        return "axis_orbit"
    if line in {(1, 1), (1, -1)}:
        return "diagonal_orbit"
    raise AssertionError(f"unexpected rank-one line {line}")


def exact_census(
    a: int = 3,
    b: int = 2,
    *,
    include_direct_rank2: bool = False,
) -> dict[str, object]:
    """Resolve the two declared line orbits on a modest Gaussian quotient.

    The defaults preserve the checked-in N=13 certificate.  The optional
    parameters let a stacked exact comparison reuse the identical orbit basis
    and Bernstein source/sink convention without copying the census engine.
    """

    geometry = gaussian_integer_torus(a, b)
    marks = subset_marks(geometry, matching=False)
    n = geometry.n
    state_counts: Counter[tuple[int, str]] = Counter()
    state_characters: dict[str, ComplexQ] = {}
    state_lines: dict[str, set[Vector]] = defaultdict(set)
    birth: Counter[tuple[int, str]] = Counter()
    exit_flux: Counter[tuple[int, str]] = Counter()
    total_edges = 0
    direct_rank2 = 0
    rank_decrease = 0
    rank_one_line_change = 0

    for mask, (old_rank, old_line, _) in enumerate(marks):
        lower_k = mask.bit_count()
        if old_rank == 1:
            assert old_line is not None
            label = orbit_label(old_line)
            value = chi4(geometry.periods.matrix, old_line)
            state_counts[(lower_k, label)] += 1
            state_lines[label].add(old_line)
            if label in state_characters and state_characters[label] != value:
                raise AssertionError("one D4 orbit acquired two chi4 values")
            state_characters[label] = value
        for vertex in range(n):
            bit = 1 << vertex
            if mask & bit:
                continue
            total_edges += 1
            new_rank, new_line, _ = marks[mask | bit]
            if new_rank < old_rank:
                rank_decrease += 1
            if old_rank == new_rank == 1 and old_line != new_line:
                rank_one_line_change += 1
            if old_rank == 0 and new_rank == 1:
                assert new_line is not None
                birth[(lower_k, orbit_label(new_line))] += 1
            if old_rank == 0 and new_rank == 2:
                direct_rank2 += 1
            if old_rank == 1 and new_rank == 2:
                assert old_line is not None
                exit_flux[(lower_k, orbit_label(old_line))] += 1

    labels = sorted(state_characters)
    coefficient_rows = []
    coefficient_failures = 0
    for lower_k in range(n):
        row: dict[str, object] = {"lower_subset_size": lower_k}
        for label in labels:
            character = state_characters[label]
            state_k = state_counts[(lower_k, label)]
            state_k1 = state_counts[(lower_k + 1, label)]
            bernstein_derivative = _qscale(
                character,
                (lower_k + 1) * state_k1 - (n - lower_k) * state_k,
            )
            transition = _qscale(
                character,
                birth[(lower_k, label)] - exit_flux[(lower_k, label)],
            )
            if transition != bernstein_derivative:
                coefficient_failures += 1
            row[label] = {
                "rank_one_states_at_k": state_k,
                "rank_one_states_at_k_plus_1": state_k1,
                "birth_edges": birth[(lower_k, label)],
                "exit_edges": exit_flux[(lower_k, label)],
                "birth_minus_exit_edges": (
                    birth[(lower_k, label)] - exit_flux[(lower_k, label)]
                ),
                "complex_derivative_coefficient": _qpayload(transition),
                "coefficient_identity_pass": transition == bernstein_derivative,
            }
        coefficient_rows.append(row)

    characters_are_opposite = (
        len(labels) == 2
        and _qadd(state_characters[labels[0]], state_characters[labels[1]])
        == (Fraction(0), Fraction(0))
    )
    result = {
        "geometry": {
            "id": f"gaussian-{a}-plus-{b}i",
            "N": n,
            "period_matrix": [list(row) for row in geometry.periods.matrix],
            "subset_states": 1 << n,
            "directed_addition_edges": total_edges,
        },
        "orbits": {
            label: {
                "primitive_lines": [list(line) for line in sorted(state_lines[label])],
                "chi4": _qpayload(state_characters[label]),
                "rank_one_state_count": sum(
                    count for (k, current), count in state_counts.items() if current == label
                ),
                "birth_edge_count": sum(
                    count for (k, current), count in birth.items() if current == label
                ),
                "exit_edge_count": sum(
                    count for (k, current), count in exit_flux.items() if current == label
                ),
            }
            for label in labels
        },
        "coefficient_rows": coefficient_rows,
        "gates": {
            "orbit_count": len(labels),
            "characters_are_exact_opposites": characters_are_opposite,
            "rank_never_decreases": rank_decrease == 0,
            "rank_one_line_never_changes": rank_one_line_change == 0,
            "coefficientwise_dA4_equals_birth_minus_exit": coefficient_failures == 0,
            "coefficient_failures": coefficient_failures,
            "all_pass": (
                len(labels) == 2
                and characters_are_opposite
                and rank_decrease == 0
                and rank_one_line_change == 0
                and coefficient_failures == 0
                and total_edges == n * (1 << (n - 1))
            ),
        },
    }
    if include_direct_rank2:
        result["direct_rank2_edge_count"] = direct_rank2
    return result


def _evaluate_incidence(
    census: dict[str, object], label: str, source: str, p: Fraction
) -> Fraction:
    n = census["geometry"]["N"]
    total = Fraction(0)
    key = "birth_edges" if source == "birth" else "exit_edges"
    for row in census["coefficient_rows"]:
        lower_k = row["lower_subset_size"]
        count = row[label][key]
        total += count * p**lower_k * (1 - p) ** (n - lower_k - 1)
    return total


def evaluate_reference(census: dict[str, object], p: Fraction = P_REF) -> dict[str, object]:
    labels = sorted(census["orbits"])
    orbit_rows = []
    total_birth: ComplexQ = (Fraction(0), Fraction(0))
    total_exit: ComplexQ = (Fraction(0), Fraction(0))
    derivatives: dict[str, ComplexQ] = {}
    for label in labels:
        character_payload = census["orbits"][label]["chi4"]
        character = (
            Fraction(character_payload["real"]),
            Fraction(character_payload["imag"]),
        )
        birth_incidence = _evaluate_incidence(census, label, "birth", p)
        exit_incidence = _evaluate_incidence(census, label, "exit", p)
        birth_value = _qscale(character, birth_incidence)
        exit_value = _qscale(character, exit_incidence)
        derivative = _qsub(birth_value, exit_value)
        derivatives[label] = derivative
        total_birth = _qadd(total_birth, birth_value)
        total_exit = _qadd(total_exit, exit_value)
        orbit_rows.append(
            {
                "orbit": label,
                "birth_incidence": _decimal(birth_incidence),
                "exit_incidence": _decimal(exit_incidence),
                "birth_minus_exit_incidence": _decimal(birth_incidence - exit_incidence),
                "complex_birth_flux": _complex_decimal(birth_value),
                "complex_exit_flux": _complex_decimal(exit_value),
                "complex_dA4_contribution": _complex_decimal(derivative),
            }
        )
    total_derivative = _qsub(total_birth, total_exit)
    # Every contribution is collinear here.  Use the real component to form an
    # exact signed share; both shares are positive because the second orbit has
    # both opposite chi4 and opposite birth-minus-exit incidence.
    shares = {
        label: derivatives[label][0] / total_derivative[0] for label in labels
    }
    return {
        "p_ref": str(p),
        "p_ref_decimal": _decimal(p),
        "orbit_flux": orbit_rows,
        "total_complex_birth_flux": _complex_decimal(total_birth),
        "total_complex_exit_flux": _complex_decimal(total_exit),
        "total_complex_dA4_dp": _complex_decimal(total_derivative),
        "signed_collinear_share": {
            label: _decimal(shares[label]) for label in labels
        },
        "shares_sum_to_one": sum(shares.values()) == 1,
        "both_orbits_reinforce_total": all(value > 0 for value in shares.values()),
    }


def build_certificate() -> dict[str, object]:
    census = exact_census()
    reference = evaluate_reference(census)
    return {
        "schema": "matching-one/p334-n13-multiorbit-flux/v1",
        "issue": 334,
        "status": "exact_multiorbit_birth_exit_flux",
        "result": {
            "statement": (
                "Gaussian quotient 3+2i (N=13) has two inequivalent primitive-line "
                "orbits with exactly opposite chi4. Both have nonzero birth and exit "
                "flux, and their orbit-resolved birth-minus-exit imbalances reinforce "
                "the same total dA4/dp direction at p_ref."
            ),
            "path_enumeration_avoided": (
                "The exact certificate visits 2^13 states and 13*2^12 directed "
                "subset-boundary edges; it does not enumerate 13! paths."
            ),
            "identity": (
                "For every orbit and lower size k, J_birth(k)-J_exit(k)="
                "(k+1)A_(k+1)-(N-k)A_k, coefficientwise in the Bernstein basis."
            ),
        },
        "census": census,
        "reference_evaluation": reference,
        "claim_boundary": [
            "This localizes the exact finite-volume H4 derivative into line orbits; it does not identify a continuum field.",
            "The p_ref shares are an exact evaluation coordinate, not a fitted asymptotic amplitude.",
            "No 13! ordering enumeration, Monte Carlo sample, or new Huawei production is used.",
        ],
    }


def render_markdown(payload: dict[str, object]) -> str:
    census = payload["census"]
    reference = payload["reference_evaluation"]
    lines = [
        "# Exact N=13 multi-orbit projective birth/exit flux",
        "",
        f"Status: `{payload['status']}`.",
        "",
        payload["result"]["statement"],
        "",
        "## Exact frontier",
        "",
        f"- subset states: {census['geometry']['subset_states']:,}",
        f"- directed addition edges: {census['geometry']['directed_addition_edges']:,}",
        f"- orbit count: {census['gates']['orbit_count']}",
        f"- coefficientwise source/sink identity: `{census['gates']['coefficientwise_dA4_equals_birth_minus_exit']}`",
        "",
        "| orbit | primitive lines | chi4 | birth edges | exit edges |",
        "|---|---|---|---:|---:|",
    ]
    for label, row in census["orbits"].items():
        lines.append(
            f"| {label} | {row['primitive_lines']} | ({row['chi4']['real']}, {row['chi4']['imag']}) | "
            f"{row['birth_edge_count']} | {row['exit_edge_count']} |"
        )
    lines += [
        "",
        f"## At p_ref = {reference['p_ref_decimal']}",
        "",
        f"Total dA4/dp = `{reference['total_complex_dA4_dp']}`.",
        "",
    ]
    for row in reference["orbit_flux"]:
        share = reference["signed_collinear_share"][row["orbit"]]
        lines.append(
            f"- {row['orbit']}: birth-exit incidence {row['birth_minus_exit_incidence']}; "
            f"signed reinforcing share {share}."
        )
    lines += [
        "",
        "The diagonal orbit has the opposite chi4 and the opposite incidence imbalance, so its contribution reinforces rather than cancels the axis-orbit H4 derivative.",
        "",
        "## Boundary",
        "",
    ]
    lines.extend(f"- {item}" for item in payload["claim_boundary"])
    return "\n".join(lines)


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--format", choices=("json", "markdown"), default="json")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args(argv)
    payload = build_certificate()
    rendered = (
        json.dumps(payload, indent=2, sort_keys=True) + "\n"
        if args.format == "json"
        else render_markdown(payload) + "\n"
    )
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
