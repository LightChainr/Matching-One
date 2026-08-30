#!/usr/bin/env python3
"""Exact eighth-power fingerprints for the four normalized-P4 channels."""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path
from typing import Dict, Mapping


ROOT = Path(__file__).resolve().parents[1]
CONTRACT = ROOT / "analysis" / "gaussian_channel_fingerprint_contract.json"

EXPONENT_EIGHTHS = {
    "P4_S": 8,
    "P4_Dprime": 5,
    "P4_D": 13,
    "P4_Sprime": 10,
}
NORMS = (2, 5)
NORM5_H12_OVER_H4 = Fraction(-1679, 625)


def fraction_text(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def h4_fingerprint(norm: int, exponent_eighths: int) -> Fraction:
    """Return the eighth power of Q^(-exponent_eighths/8)."""

    if norm <= 1 or exponent_eighths <= 0:
        raise ValueError("norm and exponent numerator must be positive")
    return Fraction(1, norm**exponent_eighths)


def h12_fingerprint(norm: int, exponent_eighths: int) -> Fraction:
    """Return the eighth-power magnitude after the normalized H12/H4 factor."""

    return abs(NORM5_H12_OVER_H4) ** 8 * h4_fingerprint(
        norm, exponent_eighths
    )


def fingerprints_at(norm: int) -> Dict[str, Fraction]:
    return {
        channel: h4_fingerprint(norm, exponent)
        for channel, exponent in EXPONENT_EIGHTHS.items()
    }


def assert_non_alias(fingerprints: Mapping[str, Fraction]) -> None:
    if len(set(fingerprints.values())) != len(fingerprints):
        raise AssertionError("two normalized-P4 channels have aliased fingerprints")


def build_contract() -> Dict[str, object]:
    fingerprints: Dict[str, Dict[str, str]] = {}
    q_cubed_relations: Dict[str, Dict[str, str]] = {}
    for norm in NORMS:
        exact = fingerprints_at(norm)
        assert_non_alias(exact)
        fingerprints[str(norm)] = {
            channel: fraction_text(value) for channel, value in exact.items()
        }
        q_cubed_relations[str(norm)] = {
            "P4_Dprime_over_P4_S": fraction_text(
                exact["P4_Dprime"] / exact["P4_S"]
            ),
            "P4_Sprime_over_P4_D": fraction_text(
                exact["P4_Sprime"] / exact["P4_D"]
            ),
        }

    h12 = {
        channel: {
            "sign": "negative",
            "eighth_power_magnitude": fraction_text(
                h12_fingerprint(5, EXPONENT_EIGHTHS[channel])
            ),
        }
        for channel in ("P4_D", "P4_Sprime")
    }
    return {
        "schema": "matching-one/gaussian-channel-fingerprint/v1",
        "status": "valid_exact_normalized_p4_channel_fingerprint",
        "parent_issue": "remain open",
        "exponent_eighths": dict(EXPONENT_EIGHTHS),
        "h4_eighth_power_fingerprints": fingerprints,
        "q_cubed_relations": q_cubed_relations,
        "norm5_h12_over_h4": fraction_text(NORM5_H12_OVER_H4),
        "norm5_same_radial_h12": h12,
        "uses_production_or_covariance_data": False,
        "selects_a_physical_channel": False,
    }


def load_contract(path: Path = CONTRACT) -> Dict[str, object]:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def validate_contract(path: Path = CONTRACT) -> Dict[str, object]:
    expected = load_contract(path)
    actual = build_contract()
    if actual != expected:
        raise AssertionError("checked-in Gaussian channel fingerprint contract drifted")
    return actual


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, default=CONTRACT)
    args = parser.parse_args()
    print(json.dumps(validate_contract(args.contract), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
