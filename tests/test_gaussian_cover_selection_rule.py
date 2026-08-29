#!/usr/bin/env python3

from __future__ import annotations

from fractions import Fraction
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from gaussian_cover_selection_rule import (  # noqa: E402
    DETAIL_SIGNS,
    NON_EQUIVARIANT_SIGNS,
    anchored_deck_even_occupancy,
    bernoulli_hessian_score,
    build_oracle,
    deck_translate,
    exact_linear_response,
    mask_from_integer,
    marked_pivotal_detail,
    matching_odd_cross,
    raw_centered_score,
)


def test_deck_characters_transform_exactly_on_every_configuration() -> None:
    for mask in range(1 << 10):
        active = mask_from_integer(mask)
        translated = deck_translate(active)
        assert deck_translate(translated) == active
        assert matching_odd_cross(translated) == matching_odd_cross(active)
        assert anchored_deck_even_occupancy(translated) == anchored_deck_even_occupancy(active)
        assert raw_centered_score(translated, DETAIL_SIGNS) == -raw_centered_score(
            active, DETAIL_SIGNS
        )
        assert marked_pivotal_detail(translated) == -marked_pivotal_detail(active)


def test_non_equivariant_registry_is_detected_not_silently_zeroed() -> None:
    assert exact_linear_response(anchored_deck_even_occupancy, DETAIL_SIGNS) == 0
    assert exact_linear_response(
        anchored_deck_even_occupancy, NON_EQUIVARIANT_SIGNS
    ) == Fraction(2)


def test_oracle_locks_allowed_linear_and_quadratic_channels() -> None:
    result = build_oracle()
    assert result["first_order_selection_rule"]["score_response"]["exact"] == "0"
    assert (
        result["opposite_character_marked_row"]["score_response"]["exact"]
        == "-10944/390625"
    )
    quadratic = result["invariant_second_order_response"]
    assert quadratic["full_bernoulli_hessian"]["exact"] == "109056/78125"
    assert quadratic["full_bernoulli_hessian"] == quadratic["direct_symbolic_second_derivative"]
    assert quadratic["raw_score_product_term"] != quadratic["full_bernoulli_hessian"]


def test_bernoulli_hessian_contains_diagonal_likelihood_correction() -> None:
    active = [True, False] * 5
    full, product, diagonal = bernoulli_hessian_score(active, DETAIL_SIGNS)
    assert full == product - diagonal
    assert diagonal != 0
