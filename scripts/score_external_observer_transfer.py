#!/usr/bin/env python3
"""Frozen two-size score for the external Euler/rank-birth experiment."""

from __future__ import annotations

import argparse
import cmath
import json
import math
from pathlib import Path
from typing import Callable, Sequence


CHANNELS = {
    "external_D": "P4_connected_O_ext_J_D",
    "far_D": "P4_connected_O_far_J_D",
    "near_D": "P4_connected_O_near_J_D",
    "external_S_control": "P4_connected_O_ext_J_S",
    "far_S_control": "P4_connected_O_far_J_S",
}


def load(path: Path) -> dict:
    result = json.loads(path.read_text(encoding="utf-8"))
    if result.get("schema") != "matching-one/marked-birth-path-score/v2":
        raise ValueError(f"unexpected score schema: {path}")
    failures = result["complement_validation"]["totals"]
    if any(int(value) for value in failures.values()):
        raise ValueError(f"complement gate failed: {path}: {failures}")
    return result


def point(report: dict, prefix: str) -> complex:
    return complex(
        float(report["P4_point"][prefix + "_re"]),
        float(report["P4_point"][prefix + "_im"]),
    )


def covariance(report: dict, names: Sequence[str]) -> list[list[float]]:
    order = report["covariance_metric_order"]
    source = report["delete_one_covariance"]
    indices = [order.index(name) for name in names]
    return [[float(source[i][j]) for j in indices] for i in indices]


def block_diagonal(first: list[list[float]], second: list[list[float]]) -> list[list[float]]:
    n, m = len(first), len(second)
    result = [[0.0] * (n + m) for _ in range(n + m)]
    for i in range(n):
        for j in range(n):
            result[i][j] = first[i][j]
    for i in range(m):
        for j in range(m):
            result[n + i][n + j] = second[i][j]
    return result


def quadratic(cov: Sequence[Sequence[float]], left: Sequence[float], right: Sequence[float]) -> float:
    return sum(left[i] * cov[i][j] * right[j]
               for i in range(len(left)) for j in range(len(right)))


def mahalanobis(value: complex, cov: Sequence[Sequence[float]]) -> float:
    a, b, c = cov[0][0], cov[0][1], cov[1][1]
    det = a * c - b * b
    if det <= 0:
        raise ValueError("non-positive complex covariance determinant")
    return (c * value.real**2 - 2 * b * value.real * value.imag +
            a * value.imag**2) / det


def ratio_summary(first: complex, second: complex,
                  cov: Sequence[Sequence[float]]) -> dict:
    value = second / first
    # d(z2/z1) = dz2/z1 - (z2/z1) dz1/z1.
    derivatives = (
        -second / (first * first),
        -1j * second / (first * first),
        1 / first,
        1j / first,
    )
    jacobian = ([entry.real for entry in derivatives],
                [entry.imag for entry in derivatives])
    ratio_cov = [[quadratic(cov, jacobian[i], jacobian[j]) for j in range(2)]
                 for i in range(2)]
    amplitude = abs(value)
    amp_grad = [value.real / amplitude, value.imag / amplitude]
    phase_grad = [-value.imag / amplitude**2, value.real / amplitude**2]
    return {
        "complex": [value.real, value.imag],
        "amplitude": amplitude,
        "phase_radians": cmath.phase(value),
        "amplitude_se": math.sqrt(max(0.0, quadratic(ratio_cov, amp_grad, amp_grad))),
        "phase_se": math.sqrt(max(0.0, quadratic(ratio_cov, phase_grad, phase_grad))),
        "covariance_re_im": ratio_cov,
    }


def numerical_gradient(function: Callable[[list[float]], float], values: list[float]) -> list[float]:
    gradient = []
    for index, value in enumerate(values):
        step = 1e-6 * max(1.0, abs(value))
        lower, upper = list(values), list(values)
        lower[index] -= step
        upper[index] += step
        gradient.append((function(upper) - function(lower)) / (2 * step))
    return gradient


def contact_phase_lock(first: dict, second: dict) -> dict:
    prefixes = (CHANNELS["external_D"], CHANNELS["external_S_control"])
    names = [prefix + part for prefix in prefixes for part in ("_re", "_im")]
    cov = block_diagonal(covariance(first, names), covariance(second, names))
    values = []
    for report in (first, second):
        for prefix in prefixes:
            z = point(report, prefix)
            values.extend((z.real, z.imag))

    def phase(row: list[float]) -> float:
        d1, s1 = complex(*row[0:2]), complex(*row[2:4])
        d2, s2 = complex(*row[4:6]), complex(*row[6:8])
        return cmath.phase((d2 / d1) / (s2 / s1))

    value = phase(values)
    gradient = numerical_gradient(phase, values)
    standard_error = math.sqrt(max(0.0, quadratic(cov, gradient, gradient)))
    return {
        "null": "arg((C_D2/C_D1)/(C_S2/C_S1)) = 0 modulo pi; amplitude free",
        "phase_radians": value,
        "phase_se": standard_error,
        "z_to_zero_local_branch": value / standard_error if standard_error else None,
    }


def build(first: dict, second: dict) -> dict:
    sizes = [int(first["N"]), int(second["N"])]
    if sizes[0] >= sizes[1]:
        raise ValueError("first score must be the smaller size")
    channels = {}
    for label, prefix in CHANNELS.items():
        names = [prefix + "_re", prefix + "_im"]
        z1, z2 = point(first, prefix), point(second, prefix)
        c1, c2 = covariance(first, names), covariance(second, names)
        channels[label] = {
            str(sizes[0]): {
                "complex": [z1.real, z1.imag],
                "covariance_re_im": c1,
                "nonzero_mahalanobis_chi2_2d": mahalanobis(z1, c1),
            },
            str(sizes[1]): {
                "complex": [z2.real, z2.imag],
                "covariance_re_im": c2,
                "nonzero_mahalanobis_chi2_2d": mahalanobis(z2, c2),
            },
            "transfer_second_over_first": ratio_summary(
                z1, z2, block_diagonal(c1, c2)
            ),
        }
    selected_names = [prefix + part for prefix in CHANNELS.values()
                      for part in ("_re", "_im")]
    return {
        "schema": "matching-one/external-observer-transfer-score/v1",
        "sizes": sizes,
        "primary": "connected P4 O_ext-J_D complex vector at both sizes",
        "mechanism_gate": "connected P4 O_far-J_D must remain nonzero; O_near is a nuisance, not a promoted field",
        "contact_null": "phase transfer of external JD is locked to external JS",
        "channels": channels,
        "contact_phase_lock": contact_phase_lock(first, second),
        "selected_metric_order_per_size": selected_names,
        "selected_block_covariance": block_diagonal(
            covariance(first, selected_names), covariance(second, selected_names)
        ),
        "dependency": "sizes are independent counter/seed groups; within-size metrics use the full delete-one covariance",
        "claim_boundary": "two sizes test nonzero complex coupling and a frozen transfer/null model; they do not identify an exponent",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--first-score", type=Path, required=True)
    parser.add_argument("--second-score", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = build(load(args.first_score), load(args.second_score))
    text = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")
    else:
        print(text, end="")


if __name__ == "__main__":
    main()
