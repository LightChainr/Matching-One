from __future__ import annotations

import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from local_matching_zeros import (  # noqa: E402
    ALL_REAL,
    PHYSICAL_ROOTS,
    analyze_catalog,
    complex_zero_route_closed,
    max_min_ratio,
    named_diagnostics,
    payload,
)


class LocalMatchingZeroTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.rows = analyze_catalog()
        cls.payload = payload()

    def test_physical_roots_match_the_committed_pilot(self) -> None:
        observed = {(row.geometry, row.L): row.physical_root for row in self.rows}
        self.assertEqual(set(observed), set(PHYSICAL_ROOTS))
        for key, value in PHYSICAL_ROOTS.items():
            self.assertAlmostEqual(observed[key], value, places=12, msg=key)

    def test_tiny_polynomials_are_all_real(self) -> None:
        for row in self.rows:
            key = (row.geometry, row.L)
            if key in ALL_REAL:
                self.assertIsNone(row.nearest_nonreal, key)
                self.assertEqual(row.n_nonreal, 0, key)
            else:
                self.assertIsNotNone(row.nearest_nonreal, key)
                self.assertGreater(row.n_nonreal, 0, key)

    def test_nearest_nonreal_usually_lies_outside_the_unit_interval(self) -> None:
        inside = [
            (row.geometry, row.L)
            for row in self.rows
            if row.re_in_unit_interval is True
        ]
        self.assertEqual(inside, [("diamond", 3)])

    def test_named_diagnostics_are_not_stable(self) -> None:
        diagnostics = named_diagnostics(self.rows)
        imag = diagnostics["imag_times_L_to_3_over_4"]
        dist = diagnostics["complex_distance_times_L_to_4"]
        self.assertEqual(len(imag), 4)
        self.assertAlmostEqual(imag[0], 0.4288, places=3)  # axis L=3
        self.assertAlmostEqual(imag[1], 0.7300, places=3)  # axis L=4
        self.assertAlmostEqual(imag[2], 0.6853, places=3)  # diamond L=2
        self.assertAlmostEqual(imag[3], 0.7037, places=3)  # diamond L=3
        self.assertGreater(max_min_ratio(imag), 1.5)
        self.assertGreater(max_min_ratio(dist), 5.0)
        self.assertTrue(complex_zero_route_closed(self.rows))
        self.assertEqual(self.payload["complex_zero_scaling_route"], "closed")

    def test_matching_partner_is_one_minus_physical_root(self) -> None:
        for row in self.rows:
            self.assertAlmostEqual(
                row.matching_partner_of_physical,
                1.0 - row.physical_root,
                places=15,
            )
            self.assertAlmostEqual(
                row.physical_self_matching_gap,
                abs(2.0 * row.physical_root - 1.0),
                places=15,
            )

    def test_metric_freeze_file_exists_and_forbids_cft_language(self) -> None:
        freeze = (ROOT / "predictions" / "local_matching_zero_metrics_20260829.yaml").read_text(
            encoding="utf-8"
        )
        self.assertIn("imag_times_L_to_3_over_4", freeze)
        self.assertIn("complex_distance_times_L_to_4", freeze)
        self.assertIn("not_fisher_zeros", freeze)
        self.assertIn("not_lee_yang_zeros", freeze)
        self.assertIn("do_not_fit_a_power", freeze)


if __name__ == "__main__":
    unittest.main()
