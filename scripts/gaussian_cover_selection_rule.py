#!/usr/bin/env python3
"""Exact norm-two deck-character selection-rule oracle for Issue #244."""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path
import sys
from typing import Callable, Sequence

sys.path.insert(0, str(Path(__file__).resolve().parent))

from integer_period_torus import classify_configuration, gaussian_integer_torus  # noqa: E402


Rational = Fraction
Mask = Sequence[bool]
Observable = Callable[[Mask], int]

PARENT = (2, 1)
CHILD = (1, 3)
PARENT_ORDER = 5
CHILD_ORDER = 10
P = Rational(2, 5)
DETAIL_SIGNS = (1,) * 5 + (-1,) * 5
# Balanced but not a character: fibers 0 and 4 have equal, rather than
# opposite, signs.  This mimics a non-equivariant fiber registry.
NON_EQUIVARIANT_SIGNS = (1, 1, 1, 1, -1, 1, -1, -1, -1, -1)
GEOMETRY = gaussian_integer_torus(*CHILD)


def fraction_record(value: Rational) -> dict[str, str | float]:
    return {"exact": str(value), "decimal": float(value)}


def mask_from_integer(mask: int, n: int = CHILD_ORDER) -> list[bool]:
    return [bool(mask & (1 << site)) for site in range(n)]


def deck_translate(active: Mask) -> list[bool]:
    half = len(active) // 2
    return list(active[half:]) + list(active[:half])


def matching_odd_cross(active: Mask) -> int:
    primal, _ = classify_configuration(GEOMETRY, active)
    matching, _ = classify_configuration(
        GEOMETRY, [not value for value in active], matching=True
    )
    return int(primal.cross) - int(matching.cross)


def anchored_deck_even_occupancy(active: Mask) -> int:
    return int(active[0]) + int(active[PARENT_ORDER])


def pivotal_delta(active: Mask, site: int, observable: Observable) -> int:
    absent = list(active)
    absent[site] = False
    present = list(active)
    present[site] = True
    return observable(present) - observable(absent)


def marked_pivotal_detail(active: Mask) -> int:
    """Opposite-character marked row at the two sites over parent fiber 0."""

    return pivotal_delta(active, 0, matching_odd_cross) - pivotal_delta(
        active, PARENT_ORDER, matching_odd_cross
    )


def configuration_weight(active: Mask, p: Rational = P) -> Rational:
    occupied = sum(active)
    return p**occupied * (1 - p) ** (len(active) - occupied)


def raw_centered_score(active: Mask, signs: Sequence[int], p: Rational = P) -> Rational:
    return sum(
        Rational(sign) * (int(state) - p) for sign, state in zip(signs, active)
    )


def linear_likelihood_score(
    active: Mask, signs: Sequence[int], p: Rational = P
) -> Rational:
    return raw_centered_score(active, signs, p) / (p * (1 - p))


def bernoulli_hessian_score(
    active: Mask, signs: Sequence[int], p: Rational = P
) -> tuple[Rational, Rational, Rational]:
    """Return full Hessian, raw score product, and diagonal correction.

    For ``p_i(epsilon)=p+epsilon*h_i``, the second likelihood-ratio
    derivative is ``L'^2 + L''``.  The returned diagonal term is ``-L''``,
    so the full score is ``product - diagonal``.
    """

    first = linear_likelihood_score(active, signs, p)
    product = first * first
    q = 1 - p
    diagonal = sum(
        Rational(sign * sign) * (1 / (p * p) if state else 1 / (q * q))
        for sign, state in zip(signs, active)
    )
    return product - diagonal, product, diagonal


def expectation_polynomial(
    observable: Observable, signs: Sequence[int], p: Rational = P
) -> list[Rational]:
    """Exact coefficients of E_{p+epsilon*h}[O]."""

    coefficients = [Rational(0)] * (len(signs) + 1)
    q = 1 - p
    for mask in range(1 << len(signs)):
        active = mask_from_integer(mask, len(signs))
        probability = [Rational(1)]
        for state, sign in zip(active, signs):
            constant = p if state else q
            slope = Rational(sign if state else -sign)
            updated = [Rational(0)] * (len(probability) + 1)
            for degree, coefficient in enumerate(probability):
                updated[degree] += constant * coefficient
                updated[degree + 1] += slope * coefficient
            probability = updated
        value = observable(active)
        for degree, coefficient in enumerate(probability):
            coefficients[degree] += value * coefficient
    return coefficients


def exact_linear_response(observable: Observable, signs: Sequence[int]) -> Rational:
    return sum(
        observable(active)
        * linear_likelihood_score(active, signs)
        * configuration_weight(active)
        for active in (mask_from_integer(mask) for mask in range(1 << CHILD_ORDER))
    )


def exact_hessian_response(
    observable: Observable, signs: Sequence[int]
) -> tuple[Rational, Rational, Rational]:
    full = Rational(0)
    product = Rational(0)
    diagonal = Rational(0)
    for mask in range(1 << CHILD_ORDER):
        active = mask_from_integer(mask)
        weight = configuration_weight(active)
        value = observable(active)
        score, score_product, correction = bernoulli_hessian_score(active, signs)
        full += value * score * weight
        product += value * score_product * weight
        diagonal += value * correction * weight
    return full, product, diagonal


def verify_group_action() -> dict:
    fixed_points = 0
    paired_orbits = 0
    visited: set[int] = set()
    for mask in range(1 << CHILD_ORDER):
        active = mask_from_integer(mask)
        translated = deck_translate(active)
        translated_mask = sum(int(value) << site for site, value in enumerate(translated))
        assert deck_translate(translated) == active
        assert matching_odd_cross(translated) == matching_odd_cross(active)
        assert anchored_deck_even_occupancy(translated) == anchored_deck_even_occupancy(active)
        assert raw_centered_score(translated, DETAIL_SIGNS) == -raw_centered_score(
            active, DETAIL_SIGNS
        )
        assert marked_pivotal_detail(translated) == -marked_pivotal_detail(active)
        if translated_mask == mask:
            fixed_points += 1
            assert raw_centered_score(active, DETAIL_SIGNS) == 0
        elif mask not in visited:
            paired_orbits += 1
        visited.add(mask)
        visited.add(translated_mask)
    return {
        "configurations": 1 << CHILD_ORDER,
        "fixed_configurations": fixed_points,
        "two_element_orbits": paired_orbits,
        "observable_character": {
            "matching_odd_cross": "+1",
            "anchored_deck_even_occupancy": "+1",
            "marked_pivotal_detail": "-1",
            "detail_score": "-1",
        },
    }


def build_oracle() -> dict:
    group_action = verify_group_action()

    invariant_linear = exact_linear_response(matching_odd_cross, DETAIL_SIGNS)
    invariant_polynomial = expectation_polynomial(matching_odd_cross, DETAIL_SIGNS)
    if invariant_linear != 0 or invariant_polynomial[1] != 0:
        raise AssertionError("deck-invariant linear detail response did not vanish")

    anchored_null = exact_linear_response(anchored_deck_even_occupancy, DETAIL_SIGNS)
    wrong_label_response = exact_linear_response(
        anchored_deck_even_occupancy, NON_EQUIVARIANT_SIGNS
    )
    wrong_polynomial = expectation_polynomial(
        anchored_deck_even_occupancy, NON_EQUIVARIANT_SIGNS
    )
    if anchored_null != 0 or wrong_label_response == 0:
        raise AssertionError("non-equivariant-label regression failed")
    if wrong_label_response != wrong_polynomial[1]:
        raise AssertionError("wrong-label score and symbolic derivative disagree")

    marked_response = exact_linear_response(marked_pivotal_detail, DETAIL_SIGNS)
    marked_polynomial = expectation_polynomial(marked_pivotal_detail, DETAIL_SIGNS)
    if marked_response == 0 or marked_response != marked_polynomial[1]:
        raise AssertionError("allowed marked response failed its symbolic check")

    hessian, raw_product, diagonal = exact_hessian_response(
        matching_odd_cross, DETAIL_SIGNS
    )
    invariant_polynomial = expectation_polynomial(matching_odd_cross, DETAIL_SIGNS)
    direct_second_derivative = 2 * invariant_polynomial[2]
    if hessian == 0 or hessian != direct_second_derivative:
        raise AssertionError("Bernoulli Hessian failed its symbolic check")
    if hessian != raw_product - diagonal:
        raise AssertionError("Bernoulli diagonal correction was not applied")

    return {
        "schema": "matching-one.gaussian-cover-selection-rule.v1",
        "issue": 244,
        "status": "exact_finite_volume_oracle",
        "geometry": {
            "parent_gaussian": list(PARENT),
            "child_gaussian": list(CHILD),
            "parent_order": PARENT_ORDER,
            "child_order": CHILD_ORDER,
            "deck_group": "Z/2",
            "deck_translation": "j -> j+5 mod 10",
            "p": str(P),
        },
        "group_action": group_action,
        "first_order_selection_rule": {
            "observable": "matching_odd_cross",
            "score_character": "detail (-1)",
            "score_response": fraction_record(invariant_linear),
            "direct_symbolic_derivative": fraction_record(invariant_polynomial[1]),
            "exact_zero": True,
        },
        "non_equivariant_negative_control": {
            "observable": "X_0+X_5 (deck even)",
            "correct_detail_response": fraction_record(anchored_null),
            "wrong_balanced_signs": list(NON_EQUIVARIANT_SIGNS),
            "wrong_signs_transform_as_character": False,
            "wrong_label_score_response": fraction_record(wrong_label_response),
            "direct_symbolic_derivative": fraction_record(wrong_polynomial[1]),
        },
        "opposite_character_marked_row": {
            "observable": "pivotal_0(matching_odd_cross)-pivotal_5(matching_odd_cross)",
            "observable_character": "detail (-1)",
            "score_character": "detail (-1)",
            "tensor_product": "trivial (+1)",
            "score_response": fraction_record(marked_response),
            "direct_symbolic_derivative": fraction_record(marked_polynomial[1]),
        },
        "invariant_second_order_response": {
            "observable": "matching_odd_cross",
            "perturbation": "p_j(epsilon)=p+epsilon for j<5; p-epsilon for j>=5",
            "full_bernoulli_hessian": fraction_record(hessian),
            "raw_score_product_term": fraction_record(raw_product),
            "diagonal_likelihood_correction_subtracted": fraction_record(diagonal),
            "direct_symbolic_second_derivative": fraction_record(
                direct_second_derivative
            ),
            "interpretation": "quadratic composite susceptibility, not a linear RG tangent",
        },
        "enumeration": "all 2^10 Bernoulli configurations",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = build_oracle()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
