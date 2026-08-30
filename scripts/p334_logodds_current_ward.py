#!/usr/bin/env python3
"""Exact natural-coordinate Ward identity for projective birth/exit current."""

from __future__ import annotations

import argparse
from decimal import Decimal, localcontext
from fractions import Fraction
import json
import math
from pathlib import Path
from typing import Iterable, Mapping, Sequence

from logodds_derivative_chain import p_to_eta_derivatives, polynomial_derivatives


ROOT = Path(__file__).resolve().parents[1]
N13_RESULT = ROOT / "results" / "p334-n13-multiorbit-flux" / "latest.json"
N17_RESULT = ROOT / "results" / "p334-n17-multiorbit-flux" / "latest.json"
P_REF = Fraction(59274605079, 100000000000)
MAXIMUM_ORDER = 6

ComplexQ = tuple[Fraction, Fraction]


def qadd(left: ComplexQ, right: ComplexQ) -> ComplexQ:
    return left[0] + right[0], left[1] + right[1]


def qsub(left: ComplexQ, right: ComplexQ) -> ComplexQ:
    return left[0] - right[0], left[1] - right[1]


def qscale(value: ComplexQ, scalar: Fraction | int) -> ComplexQ:
    return value[0] * scalar, value[1] * scalar


def qpayload(value: ComplexQ) -> dict[str, str]:
    return {"real": str(value[0]), "imag": str(value[1])}


def decimal(value: Fraction, digits: int = 18) -> str:
    with localcontext() as context:
        context.prec = digits + 8
        return format(Decimal(value.numerator) / Decimal(value.denominator), f".{digits}g")


def parse_q(value: Mapping[str, str]) -> ComplexQ:
    return Fraction(value["real"]), Fraction(value["imag"])


def poly_add(left: Sequence[ComplexQ], right: Sequence[ComplexQ]) -> list[ComplexQ]:
    degree = max(len(left), len(right))
    output = [(Fraction(0), Fraction(0)) for _ in range(degree)]
    for index in range(degree):
        if index < len(left):
            output[index] = qadd(output[index], left[index])
        if index < len(right):
            output[index] = qadd(output[index], right[index])
    return output


def poly_sub(left: Sequence[ComplexQ], right: Sequence[ComplexQ]) -> list[ComplexQ]:
    return poly_add(left, [qscale(value, -1) for value in right])


def raw_bernstein_power(values: Sequence[ComplexQ]) -> list[ComplexQ]:
    """Expand sum_k values[k] p^k (1-p)^(degree-k) in ascending powers."""

    degree = len(values) - 1
    output = [(Fraction(0), Fraction(0)) for _ in range(degree + 1)]
    for k, value in enumerate(values):
        for j in range(degree - k + 1):
            output[k + j] = qadd(
                output[k + j], qscale(value, (-1) ** j * math.comb(degree - k, j))
            )
    return output


def poly_derivative(values: Sequence[ComplexQ]) -> list[ComplexQ]:
    return [qscale(values[degree], degree) for degree in range(1, len(values))]


def multiply_p_one_minus_p(values: Sequence[ComplexQ]) -> list[ComplexQ]:
    output = [(Fraction(0), Fraction(0)) for _ in range(len(values) + 2)]
    for degree, value in enumerate(values):
        output[degree + 1] = qadd(output[degree + 1], value)
        output[degree + 2] = qsub(output[degree + 2], value)
    while output and output[-1] == (Fraction(0), Fraction(0)):
        output.pop()
    return output


def evaluate(values: Sequence[ComplexQ], p: Fraction) -> ComplexQ:
    result = (Fraction(0), Fraction(0))
    for value in reversed(values):
        result = qadd(qscale(result, p), value)
    return result


def scalar_jet(values: Sequence[ComplexQ], component: int, p: Fraction, order: int) -> list[Fraction]:
    coefficients = [value[component] for value in values]
    return polynomial_derivatives(coefficients, p, order)


def eta_jet(values: Sequence[ComplexQ], p: Fraction, order: int) -> list[ComplexQ]:
    real = p_to_eta_derivatives(scalar_jet(values, 0, p, order), p, order)
    imag = p_to_eta_derivatives(scalar_jet(values, 1, p, order), p, order)
    return list(zip(real, imag))


def orbit_coefficients(census: Mapping[str, object], label: str) -> dict[str, list[ComplexQ]]:
    n = int(census["geometry"]["N"])
    character = parse_q(census["orbits"][label]["chi4"])
    states: list[ComplexQ | None] = [None for _ in range(n + 1)]
    births: list[ComplexQ] = []
    exits: list[ComplexQ] = []
    for row in census["coefficient_rows"]:
        k = int(row["lower_subset_size"])
        current = row[label]
        state_k = qscale(character, int(current["rank_one_states_at_k"]))
        state_k1 = qscale(character, int(current["rank_one_states_at_k_plus_1"]))
        if states[k] is not None and states[k] != state_k:
            raise AssertionError("state coefficient drifted between adjacent rows")
        if states[k + 1] is not None and states[k + 1] != state_k1:
            raise AssertionError("state coefficient drifted between adjacent rows")
        states[k] = state_k
        states[k + 1] = state_k1
        births.append(qscale(character, int(current["birth_edges"])))
        exits.append(qscale(character, int(current["exit_edges"])))
    if any(value is None for value in states):
        raise AssertionError("state coefficient table is incomplete")
    return {
        "state": raw_bernstein_power([value for value in states if value is not None]),
        "birth": raw_bernstein_power(births),
        "exit": raw_bernstein_power(exits),
    }


def exact_share(value: ComplexQ, total: ComplexQ) -> Fraction:
    candidates = []
    for component in range(2):
        if total[component] != 0:
            candidates.append(value[component] / total[component])
        elif value[component] != 0:
            raise AssertionError("orbit current is not collinear with total")
    if not candidates or any(candidate != candidates[0] for candidate in candidates[1:]):
        raise AssertionError("orbit current is not collinear with total")
    return candidates[0]


def quotient_certificate(identifier: str, census: Mapping[str, object]) -> dict[str, object]:
    labels = sorted(census["orbits"])
    orbit_polynomials = {label: orbit_coefficients(census, label) for label in labels}
    totals = {
        key: [(Fraction(0), Fraction(0))]
        for key in ("state", "birth", "exit")
    }
    rows = []
    for label in labels:
        current = orbit_polynomials[label]
        for key in totals:
            totals[key] = poly_add(totals[key], current[key])
        net = poly_sub(current["birth"], current["exit"])
        weighted_net = multiply_p_one_minus_p(net)
        state_eta = eta_jet(current["state"], P_REF, MAXIMUM_ORDER)
        current_eta = eta_jet(weighted_net, P_REF, MAXIMUM_ORDER - 1)
        shifted_jet_pass = all(
            state_eta[order] == current_eta[order - 1]
            for order in range(1, MAXIMUM_ORDER + 1)
        )
        polynomial_pass = multiply_p_one_minus_p(poly_derivative(current["state"])) == weighted_net
        if not polynomial_pass or not shifted_jet_pass:
            raise AssertionError("projective-current Ward identity failed")
        rows.append(
            {
                "orbit": label,
                "ward_polynomial_exact": polynomial_pass,
                "eta_jet_shift_exact_through_order": MAXIMUM_ORDER,
                "state_at_p_ref": qpayload(evaluate(current["state"], P_REF)),
                "birth_current_at_p_ref": qpayload(evaluate(current["birth"], P_REF)),
                "exit_current_at_p_ref": qpayload(evaluate(current["exit"], P_REF)),
                "dA_dp_at_p_ref": qpayload(evaluate(net, P_REF)),
                "dA_deta_at_p_ref": qpayload(evaluate(weighted_net, P_REF)),
                "eta_derivatives": [qpayload(value) for value in state_eta],
            }
        )

    total_net = poly_sub(totals["birth"], totals["exit"])
    total_weighted = multiply_p_one_minus_p(total_net)
    total_dp = evaluate(total_net, P_REF)
    total_deta = evaluate(total_weighted, P_REF)
    share_rows = []
    for label in labels:
        orbit_net = poly_sub(orbit_polynomials[label]["birth"], orbit_polynomials[label]["exit"])
        orbit_dp = evaluate(orbit_net, P_REF)
        orbit_deta = evaluate(multiply_p_one_minus_p(orbit_net), P_REF)
        share_p = exact_share(orbit_dp, total_dp)
        share_eta = exact_share(orbit_deta, total_deta)
        if share_p != share_eta:
            raise AssertionError("orbit share depends on the thermal coordinate")
        share_rows.append(
            {
                "orbit": label,
                "exact_share": str(share_p),
                "decimal_share": decimal(share_p),
                "same_in_p_and_eta": True,
            }
        )

    boundary_zero = evaluate(totals["state"], Fraction(0)) == (0, 0) and evaluate(
        totals["state"], Fraction(1)
    ) == (0, 0)
    # Since dA/dp=J_birth-J_exit and both endpoints vanish, the signed total
    # source and sink currents have an exact zero integral over [0,1].
    integral_net = sum(
        value[0] / (degree + 1) for degree, value in enumerate(total_net)
    ), sum(value[1] / (degree + 1) for degree, value in enumerate(total_net))
    return {
        "geometry": identifier,
        "N": census["geometry"]["N"],
        "identity": "dA4/deta = p(1-p) (J_birth-J_exit)",
        "orbit_rows": rows,
        "total_dA_dp_at_p_ref": qpayload(total_dp),
        "total_dA_deta_at_p_ref": qpayload(total_deta),
        "coordinate_free_orbit_shares": share_rows,
        "empty_and_full_state_zero": boundary_zero,
        "integrated_net_current_zero": integral_net == (Fraction(0), Fraction(0)),
        "all_exact_gates_pass": boundary_zero
        and integral_net == (Fraction(0), Fraction(0))
        and all(row["ward_polynomial_exact"] for row in rows),
    }


def load_censuses() -> list[tuple[str, Mapping[str, object]]]:
    n13 = json.loads(N13_RESULT.read_text(encoding="utf-8"))["census"]
    n17 = json.loads(N17_RESULT.read_text(encoding="utf-8"))["n17_census"]
    return [(n13["geometry"]["id"], n13), (n17["geometry"]["id"], n17)]


def build_certificate() -> dict[str, object]:
    quotients = [quotient_certificate(identifier, census) for identifier, census in load_censuses()]
    return {
        "schema": "matching-one/p334-logodds-current-ward/v1",
        "issues": [182, 334],
        "status": "exact_finite_natural_coordinate_ward_identity",
        "p_ref": str(P_REF),
        "maximum_eta_derivative_order": MAXIMUM_ORDER,
        "theorem": {
            "raw_coordinate": "dA_O/dp = J_O,birth-J_O,exit",
            "natural_coordinate": "dA_O/deta = p(1-p)(J_O,birth-J_O,exit)",
            "orbit_share": "(dA_orbit/du)/(dA_total/du) is unchanged by every regular scalar thermal coordinate u",
            "stationary_point": "a net-current zero is a stationary point of the corresponding finite-volume character amplitude",
        },
        "quotients": quotients,
        "all_exact_gates_pass": all(row["all_exact_gates_pass"] for row in quotients),
        "claim_boundary": [
            "This is an exact finite-volume Ward/continuity identity, not a continuum Ward identity.",
            "It upgrades orbit shares and net-current zeros to coordinate-free finite observables; it does not assert their asymptotic limits.",
            "The N13 and N17 coefficient tables are reused without Monte Carlo or path enumeration.",
        ],
    }


def render_markdown(payload: Mapping[str, object]) -> str:
    lines = [
        "# Log-odds projective-current Ward identity",
        "",
        "The exact finite-volume identity is",
        "",
        "`dA_O/deta = p(1-p) (J_O,birth - J_O,exit)`.",
        "",
        "It holds as a polynomial for every recorded line orbit and through the complete",
        f"eta jet of order {payload['maximum_eta_derivative_order']} at the frozen p_ref.",
        "",
    ]
    for quotient in payload["quotients"]:
        lines += [
            f"## {quotient['geometry']} (N={quotient['N']})",
            "",
            f"- all exact gates: `{quotient['all_exact_gates_pass']}`",
            f"- empty/full rank-one amplitude vanishes: `{quotient['empty_and_full_state_zero']}`",
            f"- integrated net-current sum rule: `{quotient['integrated_net_current_zero']}`",
            "- coordinate-free shares at p_ref: "
            + ", ".join(
                f"{row['orbit']}={row['decimal_share']}" for row in quotient["coordinate_free_orbit_shares"]
            ),
            "",
        ]
    lines += [
        "## Consequence",
        "",
        "A zero of an orbit-resolved birth-minus-exit current is a stationary point of that orbit's finite H4 character amplitude. A zero of the total current is the corresponding total-amplitude stationary point and a pole of signed orbit shares. These locations survive every regular scalar reparameterization of the thermal coordinate.",
        "",
        "## Boundary",
        "",
    ]
    lines.extend(f"- {item}" for item in payload["claim_boundary"])
    return "\n".join(lines) + "\n"


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--format", choices=("json", "markdown"), default="json")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args(argv)
    payload = build_certificate()
    rendered = (
        json.dumps(payload, indent=2, sort_keys=True) + "\n"
        if args.format == "json"
        else render_markdown(payload)
    )
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
