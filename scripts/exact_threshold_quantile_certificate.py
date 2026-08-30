#!/usr/bin/env python3
"""Exact rational brackets for threshold-profile quantiles."""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path
from typing import Sequence

from exact_polynomial_root_certificate import evaluate, fraction_text, isolate_roots
from threshold_histogram_profile import (
    density_coefficients,
    integrate_density,
    mixture_weights,
    parse_histogram,
)


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONTRACT = ROOT / "analysis" / "threshold_histogram_profile_contract.json"


def quantile_bracket(
    cdf: Sequence[Fraction], target: Fraction, bits: int = 24
) -> tuple[Fraction, Fraction]:
    target = Fraction(target)
    if not Fraction(0) < target < Fraction(1):
        raise ValueError("quantile target must lie strictly between zero and one")
    equation = [Fraction(value) for value in cdf]
    if not equation:
        raise ValueError("CDF polynomial must not be empty")
    equation[0] -= target
    roots = isolate_roots(equation, bits=bits)
    if len(roots) != 1:
        raise ValueError(f"quantile equation must have exactly one root, found {len(roots)}")
    left, right = roots[0]
    if evaluate(cdf, left) > target or evaluate(cdf, right) < target:
        raise ValueError("quantile bracket does not have the required endpoint signs")
    return left, right


def serialize_bracket(
    cdf: Sequence[Fraction], target: Fraction, bracket: tuple[Fraction, Fraction]
) -> dict[str, str | bool]:
    left, right = bracket
    return {
        "target": fraction_text(target),
        "left": fraction_text(left),
        "right": fraction_text(right),
        "width": fraction_text(right - left),
        "left_cdf": fraction_text(evaluate(cdf, left)),
        "right_cdf": fraction_text(evaluate(cdf, right)),
        "endpoint_signs_certified": evaluate(cdf, left) <= target <= evaluate(cdf, right),
    }


def build_artifact(contract_path: Path = DEFAULT_CONTRACT, bits: int = 24) -> dict:
    contract = json.loads(contract_path.read_text(encoding="utf-8"))
    n = contract["N"]
    minus = parse_histogram(contract["K_minus_counts"], n, "K_minus")
    plus = parse_histogram(contract["K_plus_counts"], n, "K_plus")
    weights = mixture_weights(minus, plus, n)
    cdf = integrate_density(density_coefficients(weights))
    targets = (Fraction(1, 4), Fraction(1, 2), Fraction(3, 4))
    brackets = {target: quantile_bracket(cdf, target, bits) for target in targets}
    lower, median, upper = (brackets[target] for target in targets)
    reflection = lower[0] + upper[1] == 1 and lower[1] + upper[0] == 1
    if median != (Fraction(1, 2), Fraction(1, 2)) or not reflection:
        raise ValueError("frozen quantile reflection certificate drifted")
    iqr = (upper[0] - lower[1], upper[1] - lower[0])
    return {
        "schema": "matching-one/exact-threshold-quantile-certificate/v1",
        "issue": 28,
        "data_class": "exact synthetic rational-polynomial brackets",
        "isolation_bits": bits,
        "maximum_bracket_width": fraction_text(Fraction(1, 1 << bits)),
        "quantiles": [
            serialize_bracket(cdf, target, brackets[target]) for target in targets
        ],
        "quartile_reflection_certified": reflection,
        "median_exact": "1/2",
        "interquartile_range_bracket": {
            "left": fraction_text(iqr[0]),
            "right": fraction_text(iqr[1]),
            "width": fraction_text(iqr[1] - iqr[0]),
        },
        "boundary": (
            "Exact synthetic brackets only: no production histogram, empirical quantile, "
            "bootstrap uncertainty, tail fit, or universality claim."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, default=DEFAULT_CONTRACT)
    parser.add_argument("--bits", type=int, default=24)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    rendered = json.dumps(build_artifact(args.contract, args.bits), indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

