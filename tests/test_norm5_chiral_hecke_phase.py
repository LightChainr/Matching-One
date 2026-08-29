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
from norm5_chiral_fixedp_mc import (  # noqa: E402
    CHILD_MINUS,
    CHILD_PLUS,
    PRODUCTION_RUN,
    mapping_gate,
    run,
    validate_gate,
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


def test_n65_to_n325_cover_mapping_and_reflection_gate() -> None:
    gate = mapping_gate()
    assert CHILD_PLUS == ((15, -10), (10, 15))
    assert CHILD_MINUS == ((17, 6), (-6, 17))
    assert gate["unique_labels_per_child"] == 325
    assert gate["deck_action"] == "k -> k+1 mod 5"
    assert gate["same_parent_hands_are_reflections"] is False
    assert gate["passed"] is True


def test_fixedp_stream_is_worker_invariant_and_reflection_null_is_exact() -> None:
    rows_one, analysis_one = run(40, 4, 1, 0.5, 2265, 1)
    rows_two, analysis_two = run(40, 4, 2, 0.5, 2265, 1)
    assert rows_one == rows_two
    assert analysis_one == analysis_two
    assert analysis_one["true_reflection_conjugacy_null"]["point_re_im"] == [0.0, 0.0]
    assert analysis_one["true_reflection_conjugacy_null"]["covariance_of_mean"] == [
        [0.0, 0.0],
        [0.0, 0.0],
    ]
    covariance = analysis_one["primary_covariance_of_mean"]
    assert len(covariance) == 4 and all(len(row) == 4 for row in covariance)


def test_production_gate_accepts_only_the_frozen_200k_command() -> None:
    class Args:
        samples = PRODUCTION_RUN["samples"]
        batches = PRODUCTION_RUN["batches"]
        workers = PRODUCTION_RUN["workers"]
        p = PRODUCTION_RUN["p"]
        seed = PRODUCTION_RUN["seed"]
        radius = PRODUCTION_RUN["radius"]
        output = Path(PRODUCTION_RUN["output"])
        production_manifest = (
            ROOT / "experiments" / "p226_norm5_chiral_fixedp_production_20260829.json"
        )

    assert (
        validate_gate(Args(), Path(PRODUCTION_RUN["batches_output"]))
        == "production_under_frozen_manifest"
    )
    Args.seed += 1
    try:
        validate_gate(Args(), Path(PRODUCTION_RUN["batches_output"]))
    except ValueError as error:
        assert "CLI differs" in str(error)
    else:
        raise AssertionError("changed production command passed the gate")
