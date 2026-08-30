from __future__ import annotations

import copy
from fractions import Fraction
import json
from pathlib import Path
import sys
import unittest


SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))

from branching_survival_counterexample import (  # noqa: E402
    A_COORDINATES,
    B_COORDINATES,
    DEFAULT_OUTPUT,
    branch_success_direct,
    branch_success_from_h2,
    build_artifact,
    complete_survival_counts,
    complete_survival_probabilities,
    mask_from_coordinates,
    row_major_mask,
    successor_h2_distribution,
    validate_artifact,
)
from integer_period_torus import axis_integer_torus  # noqa: E402
from rank_one_survival_certificate import RankCache  # noqa: E402


class BranchingSurvivalCounterexampleTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.geometry = axis_integer_torus(4)
        cls.cache = RankCache(cls.geometry)
        cls.mask_a = mask_from_coordinates(cls.geometry, A_COORDINATES)
        cls.mask_b = mask_from_coordinates(cls.geometry, B_COORDINATES)

    def test_checked_artifact_reproduces(self) -> None:
        checked = json.loads(DEFAULT_OUTPUT.read_text(encoding="utf-8"))
        self.assertEqual(checked, build_artifact())
        self.assertEqual(
            validate_artifact(checked)["status"],
            "valid_exact_branching_survival_counterexample",
        )

    def test_coordinate_witnesses_match_declared_masks(self) -> None:
        self.assertEqual(row_major_mask(A_COORDINATES, 4), 12463)
        self.assertEqual(row_major_mask(B_COORDINATES, 4), 4343)

    def test_present_rank_and_line_are_equal(self) -> None:
        self.assertEqual(self.cache.rank_and_line(self.mask_a), (1, (1, 0)))
        self.assertEqual(self.cache.rank_and_line(self.mask_b), (1, (1, 0)))

    def test_complete_survival_laws_are_equal(self) -> None:
        expected_counts = (1, 7, 18, 20, 8, 0, 0, 0, 0)
        expected_probabilities = (
            Fraction(1), Fraction(7, 8), Fraction(9, 14),
            Fraction(5, 14), Fraction(4, 35),
            Fraction(0), Fraction(0), Fraction(0), Fraction(0),
        )
        for mask in (self.mask_a, self.mask_b):
            self.assertEqual(complete_survival_counts(self.cache, mask), expected_counts)
            self.assertEqual(
                complete_survival_probabilities(self.cache, mask),
                expected_probabilities,
            )

    def test_successor_h2_distributions_differ(self) -> None:
        self.assertEqual(successor_h2_distribution(self.cache, self.mask_a), (1, {1: 3, 2: 2, 3: 2}))
        self.assertEqual(successor_h2_distribution(self.cache, self.mask_b), (1, {1: 1, 2: 6}))

    def test_direct_branch_enumeration_is_exact(self) -> None:
        self.assertEqual(branch_success_direct(self.cache, self.mask_a), (190, 392))
        self.assertEqual(branch_success_direct(self.cache, self.mask_b), (186, 392))

    def test_h2_formula_matches_direct_branching(self) -> None:
        probability_a = branch_success_from_h2(self.cache, self.mask_a)
        probability_b = branch_success_from_h2(self.cache, self.mask_b)
        self.assertEqual(probability_a, Fraction(95, 196))
        self.assertEqual(probability_b, Fraction(93, 196))
        self.assertEqual(probability_a - probability_b, Fraction(1, 98))

    def test_immediate_cloning_does_not_distinguish(self) -> None:
        first_a = complete_survival_probabilities(self.cache, self.mask_a)[1]
        first_b = complete_survival_probabilities(self.cache, self.mask_b)[1]
        self.assertEqual(first_a * first_a, Fraction(49, 64))
        self.assertEqual(first_b * first_b, Fraction(49, 64))

    def test_invalid_inputs_fail_closed(self) -> None:
        with self.assertRaisesRegex(ValueError, "integer pairs"):
            mask_from_coordinates(self.geometry, ((0, "x"),))  # type: ignore[arg-type]
        with self.assertRaisesRegex(ValueError, "duplicate"):
            mask_from_coordinates(self.geometry, ((0, 0), (4, 0)))
        with self.assertRaisesRegex(ValueError, "outside"):
            row_major_mask(((4, 0),), 4)
        with self.assertRaisesRegex(ValueError, "rank-one"):
            complete_survival_counts(self.cache, 0)

    def test_artifact_tampering_fails_closed(self) -> None:
        tampered = copy.deepcopy(build_artifact())
        tampered["comparison"]["branch_success_gap"] = "0"
        with self.assertRaisesRegex(ValueError, "does not exactly reproduce"):
            validate_artifact(tampered)


if __name__ == "__main__":
    unittest.main()
