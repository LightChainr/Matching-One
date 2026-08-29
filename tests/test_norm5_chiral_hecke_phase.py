#!/usr/bin/env python3

from __future__ import annotations

from fractions import Fraction
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from norm5_chiral_hecke_phase import (  # noqa: E402
    ONE,
    ZERO,
    build_artifact,
    cadd,
    cconjugate,
    cmul,
    exact_marked_response,
    gaussian_ratio_power,
    zeta_power,
)


def test_exact_hecke_phase_targets() -> None:
    assert gaussian_ratio_power(4) == (-527, -336, 625)
    assert gaussian_ratio_power(8) == (164833, 354144, 390625)
    assert gaussian_ratio_power(12) == (32125393, -242017776, 244140625)
    for spin in (4, 8, 12):
        real, imag, denominator = gaussian_ratio_power(spin)
        assert real * real + imag * imag == denominator * denominator


def test_exact_zeta5_arithmetic_and_conjugation() -> None:
    assert cmul(zeta_power(2), zeta_power(3)) == ONE
    assert cconjugate(zeta_power(1)) == zeta_power(-1)
    root_sum = ZERO
    for power in range(5):
        root_sum = cadd(root_sum, zeta_power(power))
    assert root_sum == ZERO


def test_tiny_opposite_character_response_and_reflection_transport() -> None:
    plus, plus_derivative = exact_marked_response(2, 1)
    minus, minus_derivative = exact_marked_response(2, -1)
    expected = (Fraction(-46, 25), Fraction(0), Fraction(0), Fraction(0))
    assert plus == plus_derivative == expected
    assert minus == minus_derivative == expected
    assert minus == cconjugate(plus)


def test_single_handed_pair_has_three_distinct_targets() -> None:
    artifact = build_artifact()
    predictions = artifact["hecke_eigenfield_predictions"]
    assert predictions["all_three_exactly_distinct"] is True
    assert predictions["minimum_separation_degrees"] > 60.0
    assert "conditional" in predictions["single_handed_pair_discriminates"]
