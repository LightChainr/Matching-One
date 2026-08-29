#!/usr/bin/env python3
"""Exact two-activation/Beta-mixture oracle for Issue #28."""

from __future__ import annotations

import argparse
import itertools
import json
import math
import sys
from fractions import Fraction
from pathlib import Path
from typing import Iterable

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from integer_period_torus import axis_integer_torus, classify_configuration  # noqa: E402
from threshold_rank_nz import enumerate_exact, threshold_ranks  # noqa: E402


def rank_of_active(geometry, active: list[bool]) -> int:
    channels, _ = classify_configuration(geometry, active)
    return int(channels.max_rank)


def rank_trace(geometry, permutation: Iterable[int]) -> list[int]:
    active = [False] * geometry.n
    trace = [0]
    for vertex in permutation:
        active[vertex] = True
        trace.append(rank_of_active(geometry, active))
    return trace


def activations_from_trace(trace: list[int]) -> tuple[int, int]:
    if trace[0] != 0 or trace[-1] != 2:
        raise ValueError("rank trace must run from zero to two")
    if any(left > right for left, right in zip(trace, trace[1:])):
        raise ValueError("rank trace is not monotone")
    k1 = next(index for index, rank in enumerate(trace) if rank >= 1)
    k2 = next(index for index, rank in enumerate(trace) if rank == 2)
    for n, rank in enumerate(trace):
        if rank != int(n >= k1) + int(n >= k2):
            raise AssertionError("two-activation decomposition failed")
    return k1, k2


def binomial_tail(N: int, k: int, p: Fraction) -> Fraction:
    return sum(
        Fraction(math.comb(N, n)) * p**n * (1 - p) ** (N - n)
        for n in range(k, N + 1)
    )


def beta_density(N: int, k: int, p: Fraction) -> Fraction:
    return (
        Fraction(N * math.comb(N - 1, k - 1))
        * p ** (k - 1)
        * (1 - p) ** (N - k)
    )


def power_from_bernstein(values: list[Fraction]) -> list[Fraction]:
    N = len(values) - 1
    coefficients = [Fraction(0) for _ in range(N + 1)]
    for n, value in enumerate(values):
        for j in range(N - n + 1):
            coefficients[n + j] += (
                value
                * math.comb(N, n)
                * math.comb(N - n, j)
                * (-1) ** j
            )
    return coefficients


def evaluate_power(coefficients: list[Fraction], p: Fraction) -> Fraction:
    answer = Fraction(0)
    for coefficient in reversed(coefficients):
        answer = answer * p + coefficient
    return answer


def direct_microcanonical_profile(geometry) -> list[Fraction]:
    sums = [0] * (geometry.n + 1)
    counts = [0] * (geometry.n + 1)
    for mask in range(1 << geometry.n):
        active = [bool(mask & (1 << vertex)) for vertex in range(geometry.n)]
        n = sum(active)
        sums[n] += rank_of_active(geometry, active) - 1
        counts[n] += 1
    return [Fraction(total, count) for total, count in zip(sums, counts)]


def fraction_text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def exact_axis_L2_oracle() -> dict[str, object]:
    geometry = axis_integer_torus(2)
    counts = enumerate_exact(geometry)
    trace_joint: dict[tuple[int, int], int] = {}
    for permutation in itertools.permutations(range(geometry.n)):
        trace = rank_trace(geometry, permutation)
        pair = activations_from_trace(trace)
        if pair != threshold_ranks(geometry, permutation):
            raise AssertionError("rank activations disagree with existing threshold histogram")
        trace_joint[pair] = trace_joint.get(pair, 0) + 1
    if trace_joint != counts.joint:
        raise AssertionError("trace joint counts disagree with enumerate_exact")

    samples = counts.sample_count
    micro_from_thresholds = [
        Fraction(
            sum(count for (k1, _), count in counts.joint.items() if k1 <= n)
            + sum(count for (_, k2), count in counts.joint.items() if k2 <= n),
            samples,
        )
        - 1
        for n in range(geometry.n + 1)
    ]
    micro_direct = direct_microcanonical_profile(geometry)
    if micro_from_thresholds != micro_direct:
        raise AssertionError("threshold and direct microcanonical profiles differ")

    matching_power = power_from_bernstein(micro_from_thresholds)
    density_power = [
        Fraction(degree) * matching_power[degree] / 2
        for degree in range(1, len(matching_power))
    ]
    probes = (Fraction(1, 5), Fraction(1, 2), Fraction(4, 5))
    probe_rows = []
    for p in probes:
        onset = sum(
            Fraction(count, samples) * binomial_tail(geometry.n, k1, p)
            for (k1, _), count in counts.joint.items()
        )
        completion = sum(
            Fraction(count, samples) * binomial_tail(geometry.n, k2, p)
            for (_, k2), count in counts.joint.items()
        )
        matching = -1 + onset + completion
        density = Fraction(1, 2) * sum(
            Fraction(count, samples)
            * (beta_density(geometry.n, k1, p) + beta_density(geometry.n, k2, p))
            for (k1, k2), count in counts.joint.items()
        )
        if matching != evaluate_power(matching_power, p):
            raise AssertionError("Bernoulli convolution polynomial mismatch")
        if density != evaluate_power(density_power, p):
            raise AssertionError("Beta mixture derivative mismatch")
        probe_rows.append(
            {
                "p": fraction_text(p),
                "onset_CDF": fraction_text(onset),
                "completion_CDF": fraction_text(completion),
                "matching_M": fraction_text(matching),
                "density_rho": fraction_text(density),
            }
        )

    midpoint_gap: dict[tuple[Fraction, int], int] = {}
    for (k1, k2), count in counts.joint.items():
        key = (Fraction(k1 + k2, 2), k2 - k1)
        midpoint_gap[key] = midpoint_gap.get(key, 0) + count
    mean_gap = Fraction(
        sum((k2 - k1) * count for (k1, k2), count in counts.joint.items()),
        samples,
    )
    mean_midpoint = Fraction(
        sum(Fraction(k1 + k2, 2) * count for (k1, k2), count in counts.joint.items()),
        samples,
    )

    def balance(p: float) -> float:
        onset = sum(
            count / samples
            * sum(
                math.comb(geometry.n, n) * p**n * (1.0 - p) ** (geometry.n - n)
                for n in range(k1, geometry.n + 1)
            )
            for (k1, _), count in counts.joint.items()
        )
        completion = sum(
            count / samples
            * sum(
                math.comb(geometry.n, n) * p**n * (1.0 - p) ** (geometry.n - n)
                for n in range(k2, geometry.n + 1)
            )
            for (_, k2), count in counts.joint.items()
        )
        return onset + completion - 1.0

    left, right = 0.0, 1.0
    for _ in range(100):
        middle = (left + right) / 2.0
        if balance(middle) < 0.0:
            left = middle
        else:
            right = middle
    root = (left + right) / 2.0

    return {
        "schema": "matching-one.two-activation-rank-mixture.v1",
        "issue": 28,
        "theorem_dependency": "digital Alexander rank identity from Issue 269",
        "oracle": "existing threshold_rank_nz exact axis L=2 histogram",
        "N": geometry.n,
        "permutations": samples,
        "joint_K1_K2_counts": [
            {"K1": k1, "K2": k2, "count": count}
            for (k1, k2), count in sorted(counts.joint.items())
        ],
        "joint_midpoint_gap_counts": [
            {"C": fraction_text(center), "G": gap, "count": count}
            for (center, gap), count in sorted(midpoint_gap.items())
        ],
        "microcanonical_q_by_n": [fraction_text(value) for value in micro_from_thresholds],
        "matching_power_coefficients_ascending": [fraction_text(value) for value in matching_power],
        "density_power_coefficients_ascending": [fraction_text(value) for value in density_power],
        "beta_mixture": "(1/3) Beta(2,3) + (2/3) Beta(3,2)",
        "exact_rational_probes": probe_rows,
        "root_two_CDF_balance": {
            "equation": "E[H_K1(p)]+E[H_K2(p)]=1",
            "exact_root": "sqrt(1-1/sqrt(2))",
            "numerical_root": root,
            "balance_residual": balance(root),
        },
        "midpoint_gap_semantics": {
            "mean_C": fraction_text(mean_midpoint),
            "mixture_mean_EC_over_N_plus_1": fraction_text(mean_midpoint / (geometry.n + 1)),
            "mean_G": fraction_text(mean_gap),
            "integrated_rank1_probability_EG_over_N_plus_1": fraction_text(mean_gap / (geometry.n + 1)),
            "G_definition": "number of n with K1<=n<K2, equivalently r_n=1",
        },
        "checks": {
            "all_rank_traces_monotone": True,
            "rank_trace_equals_two_activation_indicators": True,
            "trace_joint_equals_existing_exact_histogram": True,
            "direct_microcanonical_equals_threshold_profile": True,
            "Bernoulli_convolution_equals_matching_polynomial": True,
            "Beta_mixture_equals_half_matching_derivative": True,
        },
        "scope": "tiny short-period oracle is a machine certificate; the honest-cell proof scope is supplied by Issue 269",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("analysis/two_activation_rank_mixture_exact.json"),
    )
    args = parser.parse_args()
    payload = exact_axis_L2_oracle()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
