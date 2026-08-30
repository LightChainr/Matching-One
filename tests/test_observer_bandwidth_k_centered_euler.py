from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from observer_bandwidth_k_centered_euler import (  # noqa: E402
    build_report,
    conditional_checks,
    conditional_mean,
    degree_one_projections,
    euler_observer_values,
    k_center,
    square_torus_incidence,
)


class KCenteredEulerBandwidthTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.values, cls.incidence = euler_observer_values(3)
        cls.centered = k_center(cls.values, 9)

    def test_periodic_incidence_has_declared_distinct_counts(self) -> None:
        edges, faces = square_torus_incidence(3)
        self.assertEqual(len(edges), 18)
        self.assertEqual(len(faces), 9)
        self.assertTrue(all(len(set(face)) == 4 for face in faces))

    def test_conditional_mean_formula_matches_every_slice(self) -> None:
        rows = conditional_checks(self.values, self.centered, 9)
        self.assertEqual(len(rows), 10)
        self.assertEqual(rows[4]["formula_mean"], str(conditional_mean(9, 4)))
        self.assertTrue(all(row["centered_sum"] == "0" for row in rows))

    def test_centered_degree_one_projection_is_exactly_zero(self) -> None:
        for p in (Fraction(1, 3), Fraction(2, 5), Fraction(1, 2)):
            self.assertEqual(degree_one_projections(self.centered, 9, p), [0] * 9)

    def test_checked_report_is_exactly_reproducible(self) -> None:
        manifest = json.loads(
            (ROOT / "analysis/observer_bandwidth_k_centered_euler_manifest.json").read_text(
                encoding="utf-8"
            )
        )
        checked = json.loads(
            (ROOT / "results/observer-bandwidth-k-centered-euler/latest.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(checked, build_report(manifest))
        for row in checked["p_checks"]:
            self.assertEqual(row["centered_site_projections"], ["0"] * 9)
            self.assertNotIn(1, row["active_degrees"])
            self.assertLessEqual(max(row["active_degrees"]), 4)


if __name__ == "__main__":
    unittest.main()
