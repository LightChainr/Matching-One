from __future__ import annotations

from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from score_p334_age_iota_site_control import (  # noqa: E402
    add_sufficient,
    jackknife_covariance,
    slope,
)


class P334AgeIotaSiteControlTests(unittest.TestCase):
    def test_iota_noop_relabels_identical_sufficient_statistics(self) -> None:
        primary = {}
        indexed = {}
        rows = [(0, 10, False), (1, 20, True), (2, 15, False)]
        for age, count, event in rows:
            add_sufficient(primary, (1, 0), age, count, event)
            add_sufficient(indexed, (1, 0, 1, 1), age, count, event)
        self.assertEqual(
            slope(primary, 100)["beta_age_per_density"],
            slope(indexed, 100)["beta_age_per_density"],
        )
        self.assertEqual(
            slope(primary, 100)["within_stratum_age_denominator_steps2"],
            slope(indexed, 100)["within_stratum_age_denominator_steps2"],
        )

    def test_fixed_site_class_can_only_remove_between_class_age_variation(self) -> None:
        primary = {}
        controlled = {}
        rows = [(0, 10, False, 1), (1, 20, True, 1), (2, 15, False, 5), (3, 12, True, 5)]
        for age, count, event, site_class in rows:
            add_sufficient(primary, (1, 0), age, count, event)
            add_sufficient(controlled, (1, 0, 1, 1, site_class), age, count, event)
        self.assertLessEqual(
            slope(controlled, 100)["within_stratum_age_denominator_steps2"],
            slope(primary, 100)["within_stratum_age_denominator_steps2"],
        )

    def test_identical_delete_one_columns_have_exact_zero_difference_variance(self) -> None:
        covariance = jackknife_covariance([[1.0, 1.0], [0.9, 0.9], [1.1, 1.1]])
        difference_variance = covariance[1, 1] + covariance[0, 0] - 2 * covariance[1, 0]
        self.assertEqual(difference_variance, 0.0)


if __name__ == "__main__":
    unittest.main()
