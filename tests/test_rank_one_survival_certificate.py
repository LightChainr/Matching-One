
from __future__ import annotations
from fractions import Fraction
import json
from pathlib import Path
import sys
import unittest

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))

from integer_period_torus import gaussian_integer_torus  # noqa: E402
from rank_one_survival_certificate import (  # noqa: E402
    DEFAULT_OUTPUT,
    RankCache,
    active_from_mask,
    build_artifact,
    killed_kernel_survival,
    mask_from_labels,
    permutation_exit_counts,
    subset_survival,
    trigger_layers,
    two_step_identity,
    validate_artifact,
)


class RankOneSurvivalCertificateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.geometry = gaussian_integer_torus(3, 1)
        cls.cache = RankCache(cls.geometry)
        cls.mask_a = mask_from_labels(cls.geometry, (0, 1, 2, 3, 4))
        cls.mask_b = mask_from_labels(cls.geometry, (0, 1, 2, 3, 5))

    def test_checked_artifact_reproduces(self) -> None:
        checked = json.loads(DEFAULT_OUTPUT.read_text(encoding="utf-8"))
        self.assertEqual(checked, build_artifact())
        self.assertEqual(
            validate_artifact(checked)["status"],
            "valid_exact_rank_one_survival_certificate",
        )

    def test_witnesses_have_equal_declared_present_state(self) -> None:
        self.assertEqual(self.cache.rank_and_line(self.mask_a), (1, (0, 1)))
        self.assertEqual(self.cache.rank_and_line(self.mask_b), (1, (0, 1)))
        singleton_a, _ = trigger_layers(self.cache, self.mask_a)
        singleton_b, _ = trigger_layers(self.cache, self.mask_b)
        self.assertEqual((len(singleton_a), len(singleton_b)), (1, 1))

    def test_minimal_trigger_pairs_are_exact(self) -> None:
        singleton_a, pairs_a = trigger_layers(self.cache, self.mask_a)
        singleton_b, pairs_b = trigger_layers(self.cache, self.mask_b)
        vertex = self.geometry.vertex
        self.assertEqual(singleton_a, (vertex((0, 7)),))
        self.assertEqual(singleton_b, (vertex((0, 8)),))
        self.assertEqual(
            set(pairs_a),
            {
                tuple(sorted((vertex((0, 5)), vertex((0, 8))))),
                tuple(sorted((vertex((0, 6)), vertex((0, 9))))),
            },
        )
        self.assertEqual(
            set(pairs_b),
            {
                tuple(sorted((vertex((0, 4)), vertex((0, 7))))),
                tuple(sorted((vertex((0, 6)), vertex((0, 7))))),
                tuple(sorted((vertex((0, 6)), vertex((0, 9))))),
            },
        )

    def test_subset_survival_separates_witnesses_at_two_steps(self) -> None:
        self.assertEqual(subset_survival(self.cache, self.mask_a, 1), Fraction(4, 5))
        self.assertEqual(subset_survival(self.cache, self.mask_b, 1), Fraction(4, 5))
        self.assertEqual(subset_survival(self.cache, self.mask_a, 2), Fraction(2, 5))
        self.assertEqual(subset_survival(self.cache, self.mask_b, 2), Fraction(3, 10))

    def test_future_permutation_exit_counts_are_independent_controls(self) -> None:
        self.assertEqual(permutation_exit_counts(self.cache, self.mask_a), {1: 24, 2: 48, 3: 48})
        self.assertEqual(permutation_exit_counts(self.cache, self.mask_b), {1: 24, 2: 60, 3: 36})

    def test_killed_kernel_matches_subset_survival(self) -> None:
        for mask in (self.mask_a, self.mask_b):
            for horizon in range(6):
                self.assertEqual(
                    killed_kernel_survival(self.cache, mask, horizon),
                    subset_survival(self.cache, mask, horizon),
                )

    def test_two_step_identity_is_exhaustive_on_n10_rank_one_states(self) -> None:
        checked = 0
        for mask in range(1 << self.geometry.n):
            if self.cache.rank(mask) != 1:
                continue
            q = self.geometry.n - bin(mask).count("1")
            if q >= 2:
                self.assertEqual(
                    two_step_identity(self.cache, mask),
                    subset_survival(self.cache, mask, 2),
                )
                checked += 1
        self.assertEqual(checked, 310)

    def test_invalid_inputs_fail_closed(self) -> None:
        with self.assertRaisesRegex(ValueError, "outside"):
            active_from_mask(self.geometry.n, -1)
        with self.assertRaisesRegex(ValueError, "outside"):
            mask_from_labels(self.geometry, (10,))
        with self.assertRaisesRegex(ValueError, "duplicate"):
            mask_from_labels(self.geometry, (0, 0))
        rank_zero = mask_from_labels(self.geometry, ())
        with self.assertRaisesRegex(ValueError, "rank-one"):
            trigger_layers(self.cache, rank_zero)
        with self.assertRaisesRegex(ValueError, "horizon"):
            subset_survival(self.cache, self.mask_a, 6)


if __name__ == "__main__":
    unittest.main()
