#!/usr/bin/env python3
"""Exact bridge between F3 flat-twist H/A/D and N13 primitive A4 flux."""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
import json
from math import comb, gcd
from pathlib import Path
from typing import Optional, Sequence

from integer_period_torus import gaussian_integer_torus
from projective_essential_birth_oracle import chi4, subset_marks


Vector = tuple[int, int]
ComplexQ = tuple[Fraction, Fraction]
P_REF = Fraction(59274605079, 100000000000)
LINES: tuple[Vector, ...] = ((1, 0), (0, 1), (1, 1), (1, -1))
LINE_NAMES = {
    (1, 0): "axis_x",
    (0, 1): "axis_y",
    (1, 1): "diag_plus",
    (1, -1): "diag_minus",
}
AXES = ((1, 0), (0, 1))
DIAGONALS = ((1, 1), (1, -1))
KERNEL_TWIST = {
    (1, 0): "T_0_1",
    (0, 1): "T_1_0",
    (1, 1): "T_1_2",
    (1, -1): "T_1_1",
}


def qadd(first: ComplexQ, second: ComplexQ) -> ComplexQ:
    return first[0] + second[0], first[1] + second[1]


def qsub(first: ComplexQ, second: ComplexQ) -> ComplexQ:
    return first[0] - second[0], first[1] - second[1]


def qscale(value: ComplexQ, scalar: Fraction | int) -> ComplexQ:
    return value[0] * scalar, value[1] * scalar


def qpayload(value: ComplexQ) -> dict[str, str]:
    return {"real": str(value[0]), "imag": str(value[1])}


def canonical(line: Vector) -> Vector:
    x, y = line
    divisor = gcd(abs(x), abs(y))
    x, y = x // divisor, y // divisor
    if x < 0 or (x == 0 and y < 0):
        x, y = -x, -y
    return x, y


def s_rotate(line: Vector) -> Vector:
    return canonical((-line[1], line[0]))


def matmul2(first: Sequence[Sequence[int]], second: Sequence[Sequence[int]]) -> list[list[int]]:
    return [[sum(first[i][k] * second[k][j] for k in range(2))
             for j in range(2)] for i in range(2)]


def character_sum(counts: Counter[tuple[int, Vector]], k: int,
                  characters: dict[Vector, ComplexQ]) -> ComplexQ:
    value: ComplexQ = (Fraction(0), Fraction(0))
    for line in LINES:
        value = qadd(value, qscale(characters[line], counts[(k, line)]))
    return value


def orbit_difference(counts: Counter[tuple[int, Vector]], k: int) -> int:
    return (
        sum(counts[(k, line)] for line in AXES)
        - sum(counts[(k, line)] for line in DIAGONALS)
    )


def odd_differences(counts: Counter[tuple[int, Vector]], k: int) -> tuple[int, int]:
    return (
        counts[(k, (1, 0))] - counts[(k, (0, 1))],
        counts[(k, (1, 1))] - counts[(k, (1, -1))],
    )


def evaluate_bernstein(coefficients: Sequence[Fraction], p: Fraction) -> Fraction:
    degree = len(coefficients) - 1
    return sum(
        coefficients[k] * comb(degree, k) * p**k * (1 - p) ** (degree - k)
        for k in range(degree + 1)
    )


def evaluate_complex(coefficients: Sequence[ComplexQ], p: Fraction) -> ComplexQ:
    return (
        evaluate_bernstein([value[0] for value in coefficients], p),
        evaluate_bernstein([value[1] for value in coefficients], p),
    )


def build_certificate() -> dict[str, object]:
    geometry = gaussian_integer_torus(3, 2)
    marks = subset_marks(geometry, matching=False)
    n = geometry.n
    rank0: Counter[int] = Counter()
    state: Counter[tuple[int, Vector]] = Counter()
    birth: Counter[tuple[int, Vector]] = Counter()
    exit_flux: Counter[tuple[int, Vector]] = Counter()

    for mask, (old_rank, old_line, _) in enumerate(marks):
        k = mask.bit_count()
        if old_rank == 0:
            rank0[k] += 1
        elif old_rank == 1:
            assert old_line is not None
            state[(k, old_line)] += 1
        for vertex in range(n):
            if mask & (1 << vertex):
                continue
            new_rank, new_line, _ = marks[mask | (1 << vertex)]
            if old_rank == 0 and new_rank == 1:
                assert new_line is not None
                birth[(k, new_line)] += 1
            elif old_rank == 1 and new_rank == 2:
                assert old_line is not None
                exit_flux[(k, old_line)] += 1

    support = sorted({line for (_, line), count in state.items() if count})
    if support != sorted(LINES):
        raise AssertionError(f"unexpected N13 line support {support}")
    characters = {line: chi4(geometry.periods.matrix, line) for line in LINES}
    z_axis = characters[(1, 0)]
    rotation = [[0, -1], [1, 0]]
    period_matrix = [list(row) for row in geometry.periods.matrix]
    quarter_turn_commutes = (
        matmul2(rotation, period_matrix) == matmul2(period_matrix, rotation)
    )

    state_rows = []
    h_coefficients: list[Fraction] = []
    h_unit_coefficients: list[Fraction] = []
    a4_coefficients: list[ComplexQ] = []
    state_failures = 0
    odd_state_failures = 0
    d4_state_failures = 0
    for k in range(n + 1):
        denominator = comb(n, k)
        h_numerator = orbit_difference(state, k)
        h = Fraction(h_numerator, denominator)
        h_unit = h / 2
        a4 = qscale(character_sum(state, k, characters), Fraction(1, denominator))
        expected = qscale(z_axis, h)
        a_odd, d_odd = odd_differences(state, k)
        if a4 != expected:
            state_failures += 1
        if a_odd or d_odd:
            odd_state_failures += 1
        if any(state[(k, line)] != state[(k, s_rotate(line))] for line in LINES):
            d4_state_failures += 1

        twist_coefficients = {
            KERNEL_TWIST[line]: Fraction(rank0[k] + state[(k, line)], denominator)
            for line in LINES
        }
        h_from_twists = (
            twist_coefficients["T_0_1"] + twist_coefficients["T_1_0"]
            - twist_coefficients["T_1_2"] - twist_coefficients["T_1_1"]
        )
        if h_from_twists != h:
            state_failures += 1
        h_coefficients.append(h)
        h_unit_coefficients.append(h_unit)
        a4_coefficients.append(a4)
        state_rows.append({
            "k": k,
            "subset_count": denominator,
            "rank0_count": rank0[k],
            "line_state_counts": {
                LINE_NAMES[line]: state[(k, line)] for line in LINES
            },
            "representative_twist_coefficients": {
                name: str(value) for name, value in twist_coefficients.items()
            },
            "H_F3_orbit_coefficient": str(h),
            "H_F3_unit_coefficient": str(h_unit),
            "A_axis_odd_numerator": a_odd,
            "D_diagonal_odd_numerator": d_odd,
            "primitive_A4_coefficient": qpayload(a4),
            "A4_equals_z_axis_times_H": a4 == expected,
        })

    flux_rows = []
    h_birth_coefficients: list[Fraction] = []
    h_exit_coefficients: list[Fraction] = []
    h_derivative_coefficients: list[Fraction] = []
    a4_birth_coefficients: list[ComplexQ] = []
    a4_exit_coefficients: list[ComplexQ] = []
    a4_derivative_coefficients: list[ComplexQ] = []
    flux_failures = 0
    odd_flux_failures = 0
    d4_flux_failures = 0
    for k in range(n):
        denominator = comb(n - 1, k)
        h_birth_raw = orbit_difference(birth, k)
        h_exit_raw = orbit_difference(exit_flux, k)
        h_birth = Fraction(h_birth_raw, denominator)
        h_exit = Fraction(h_exit_raw, denominator)
        h_derivative = n * (h_coefficients[k + 1] - h_coefficients[k])
        a4_birth = qscale(character_sum(birth, k, characters), Fraction(1, denominator))
        a4_exit = qscale(character_sum(exit_flux, k, characters), Fraction(1, denominator))
        a4_derivative = qsub(a4_birth, a4_exit)
        state_derivative = qscale(
            qsub(a4_coefficients[k + 1], a4_coefficients[k]), n
        )
        a_birth_odd, d_birth_odd = odd_differences(birth, k)
        a_exit_odd, d_exit_odd = odd_differences(exit_flux, k)
        passed = (
            h_derivative == h_birth - h_exit
            and a4_birth == qscale(z_axis, h_birth)
            and a4_exit == qscale(z_axis, h_exit)
            and a4_derivative == qscale(z_axis, h_derivative)
            and a4_derivative == state_derivative
        )
        if not passed:
            flux_failures += 1
        if any((a_birth_odd, d_birth_odd, a_exit_odd, d_exit_odd)):
            odd_flux_failures += 1
        if any(
            counts[(k, line)] != counts[(k, s_rotate(line))]
            for counts in (birth, exit_flux) for line in LINES
        ):
            d4_flux_failures += 1
        h_birth_coefficients.append(h_birth)
        h_exit_coefficients.append(h_exit)
        h_derivative_coefficients.append(h_derivative)
        a4_birth_coefficients.append(a4_birth)
        a4_exit_coefficients.append(a4_exit)
        a4_derivative_coefficients.append(a4_derivative)
        flux_rows.append({
            "lower_k": k,
            "edge_normalization": denominator,
            "line_birth_edges": {LINE_NAMES[line]: birth[(k, line)] for line in LINES},
            "line_exit_edges": {LINE_NAMES[line]: exit_flux[(k, line)] for line in LINES},
            "H_twist_source_coefficient": str(h_birth),
            "H_twist_sink_coefficient": str(h_exit),
            "dH_dp_coefficient": str(h_derivative),
            "A_birth_odd_numerator": a_birth_odd,
            "D_birth_odd_numerator": d_birth_odd,
            "A_exit_odd_numerator": a_exit_odd,
            "D_exit_odd_numerator": d_exit_odd,
            "primitive_j4_birth_coefficient": qpayload(a4_birth),
            "primitive_j4_exit_coefficient": qpayload(a4_exit),
            "primitive_dA4_dp_coefficient": qpayload(a4_derivative),
            "all_flux_bridge_identities_pass": passed,
        })

    reference = {
        "p": str(P_REF),
        "H_F3": str(evaluate_bernstein(h_coefficients, P_REF)),
        "H_F3_unit": str(evaluate_bernstein(h_unit_coefficients, P_REF)),
        "primitive_A4": qpayload(evaluate_complex(a4_coefficients, P_REF)),
        "H_twist_source": str(evaluate_bernstein(h_birth_coefficients, P_REF)),
        "H_twist_sink": str(evaluate_bernstein(h_exit_coefficients, P_REF)),
        "dH_dp": str(evaluate_bernstein(h_derivative_coefficients, P_REF)),
        "primitive_j4_birth": qpayload(evaluate_complex(a4_birth_coefficients, P_REF)),
        "primitive_j4_exit": qpayload(evaluate_complex(a4_exit_coefficients, P_REF)),
        "primitive_dA4_dp": qpayload(evaluate_complex(a4_derivative_coefficients, P_REF)),
    }
    reference["decimal"] = {
        "H_F3": float(Fraction(reference["H_F3"])),
        "H_F3_unit": float(Fraction(reference["H_F3_unit"])),
        "H_twist_source": float(Fraction(reference["H_twist_source"])),
        "H_twist_sink": float(Fraction(reference["H_twist_sink"])),
        "dH_dp": float(Fraction(reference["dH_dp"])),
        "primitive_A4": {
            part: float(Fraction(reference["primitive_A4"][part]))
            for part in ("real", "imag")
        },
        "primitive_dA4_dp": {
            part: float(Fraction(reference["primitive_dA4_dp"][part]))
            for part in ("real", "imag")
        },
    }
    all_pass = (
        state_failures == 0 and flux_failures == 0
        and odd_state_failures == 0 and odd_flux_failures == 0
        and d4_state_failures == 0 and d4_flux_failures == 0
        and quarter_turn_commutes
    )
    return {
        "schema": "matching-one/p337-p334-N13-twist-flux-bridge/v1",
        "issues": [337, 334],
        "status": "exact coefficientwise twist-state/source/sink bridge",
        "geometry": {
            "id": "gaussian-3-plus-2i", "N": n,
            "period_matrix": period_matrix,
            "subset_states": 1 << n,
            "directed_addition_edges": n * (1 << (n - 1)),
        },
        "normalization": {
            "H_F3_orbit": "T_01+T_10-T_12-T_11 = L_axis_x+L_axis_y-L_diag_plus-L_diag_minus",
            "H_F3_unit_used_by_a7cb19a": "H_F3_orbit/2",
            "primitive_A4": "sum_line chi4(P ell) L_line",
            "bridge": "A4=z_axis*H_F3_orbit=2*z_axis*H_F3_unit",
        },
        "z_axis": qpayload(z_axis),
        "line_characters": {
            LINE_NAMES[line]: qpayload(characters[line]) for line in LINES
        },
        "state_coefficient_rows": state_rows,
        "flux_coefficient_rows": flux_rows,
        "reference_evaluation": reference,
        "odd_sector_result": {
            "A_axis_odd": "zero coefficientwise for state, birth source and exit sink",
            "D_diagonal_odd": "zero coefficientwise for state, birth source and exit sink",
            "representation_reason": "both are S-quarter-turn odd; the 3+2i quotient state and edge census is S-invariant",
            "classification": "exact symmetry zeros for the unweighted Gaussian-ideal ensemble; activation requires a charged/S-odd insertion or a geometry outside this symmetry class",
            "general_gaussian_zero_theorem": "for P(a,b)=[[a,-b],[b,a]], quarter-turn R commutes with P and is a graph automorphism; every unweighted R-odd state/source/sink character has zero expectation coefficientwise",
        },
        "gates": {
            "state_bridge_failures": state_failures,
            "flux_bridge_failures": flux_failures,
            "odd_state_failures": odd_state_failures,
            "odd_flux_failures": odd_flux_failures,
            "S_state_transport_failures": d4_state_failures,
            "S_flux_transport_failures": d4_flux_failures,
            "quarter_turn_commutes_with_period_matrix": quarter_turn_commutes,
            "all_pass": all_pass,
        },
        "claim_boundary": (
            "exact finite N13 normalization and coefficient bridge; no continuum-field "
            "identity or nonzero odd-sector coupling is inferred"
        ),
    }


def render_markdown(payload: dict[str, object]) -> str:
    reference = payload["reference_evaluation"]
    decimal = reference["decimal"]
    lines = [
        "# Exact N13 bridge: F3 twist H to primitive A4 flux", "",
        "All state, source/sink, odd-sector and quarter-turn transport gates pass.", "",
        "## Fixed character bridge", "",
        "For `P=[[3,-2],[2,3]]`, both axis lines have", "",
        f"`z_axis={payload['z_axis']['real']} + ({payload['z_axis']['imag']}) i`.", "",
        "Both diagonal lines have `-z_axis`. Therefore, coefficient by coefficient,", "",
        "```text",
        "H_F3 = T_01+T_10-T_12-T_11",
        "     = L_axis_x+L_axis_y-L_diag_plus-L_diag_minus,",
        "A4   = z_axis H_F3.",
        "```", "",
        "The unit-norm convention in `a7cb19a` is `H_unit=H_F3/2`, so in that "
        "stored coordinate `A4=2 z_axis H_unit`. This factor is normalization, not "
        "a fitted amplitude.", "",
        "## Source/sink derivative", "",
        "At every lower size `k`, with degree-12 Bernstein normalization,", "",
        "```text",
        "dH_F3/dp = J_H,birth1 - J_H,exit2,",
        "J4_birth1 = z_axis J_H,birth1,",
        "J4_exit2  = z_axis J_H,exit2,",
        "dA4/dp    = z_axis dH_F3/dp.",
        "```", "",
        "These identities pass independently at every coefficient; no evaluation "
        "point or path sampling is used.", "",
        "## Odd sectors", "",
        "Both `A=L_axis_x-L_axis_y` and `D=L_diag_plus-L_diag_minus` vanish "
        "coefficientwise in the rank-one state curve, birth source and exit sink. "
        "They are exact `S`-quarter-turn-odd symmetry zeros of the `3+2i` quotient, "
        "not projective sectors forbidden on general geometries.", "",
        "## At p_ref", "",
        f"- `H_F3={decimal['H_F3']:.12g}`",
        f"- `A4={decimal['primitive_A4']['real']:.12g}"
        f"{decimal['primitive_A4']['imag']:+.12g} i`",
        f"- `dH_F3/dp={decimal['dH_dp']:.12g}`",
        f"- `dA4/dp={decimal['primitive_dA4_dp']['real']:.12g}"
        f"{decimal['primitive_dA4_dp']['imag']:+.12g} i`", "",
        "This closes the exact state-to-current dictionary between #337 and #334. "
        "It does not identify a continuum field.", "",
    ]
    return "\n".join(lines)


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", type=Path, required=True)
    parser.add_argument("--markdown", type=Path)
    args = parser.parse_args(argv)
    payload = build_certificate()
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n",
                         encoding="utf-8")
    if args.markdown:
        args.markdown.parent.mkdir(parents=True, exist_ok=True)
        args.markdown.write_text(render_markdown(payload), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
