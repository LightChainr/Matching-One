from __future__ import annotations

import json
import math
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import mpmath as mp


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
sys.path.insert(0, str(SCRIPTS))

from analyze_threshold_ranks import matching_value  # noqa: E402
from analyze_two_activation_h4 import (  # noqa: E402
    DECISION_METRICS,
    FULL_METRICS,
    SCHEMA,
    _joint_coordinates,
    analyze_archive,
    activation_components,
    jackknife_covariance,
)
from integer_period_torus import axis_integer_torus, gaussian_integer_torus  # noqa: E402
from threshold_rank_nz import enumerate_exact  # noqa: E402


class ExactActivationConventionTests(unittest.TestCase):
    def test_tiny_exact_tori_reconstruct_matching_curve(self) -> None:
        mp.mp.dps = 60
        for geometry in (axis_integer_torus(2), gaussian_integer_torus(2, 1)):
            counts = enumerate_exact(geometry)
            self.assertTrue(all(k1 <= k2 for k1, k2 in counts.joint))
            for p in (0.31, 0.5927460507921, 0.73):
                f1, f2, _d1, _d2 = activation_components(
                    geometry.n,
                    counts.sample_count,
                    counts.kminus,
                    counts.kplus,
                    p,
                )
                canonical = matching_value(
                    geometry.n,
                    counts.sample_count,
                    counts.kminus,
                    counts.kplus,
                    mp.mpf(str(p)),
                )
                self.assertAlmostEqual(f1 + f2 - 1.0, float(canonical), places=13)

    def test_joint_moments_transform_to_midpoint_and_gap(self) -> None:
        pairs = ((1, 3), (2, 5), (4, 4), (3, 6))
        totals = {
            "samples": len(pairs),
            "sum_kminus": sum(left for left, _right in pairs),
            "sum_kplus": sum(right for _left, right in pairs),
            "sum_kminus2": sum(left * left for left, _right in pairs),
            "sum_kplus2": sum(right * right for _left, right in pairs),
            "sum_product": sum(left * right for left, right in pairs),
            "sum_gap": sum(right - left for left, right in pairs),
            "sum_gap2": sum((right - left) ** 2 for left, right in pairs),
        }
        result = _joint_coordinates(totals)
        centers = [(left + right) / 2 for left, right in pairs]
        gaps = [right - left for left, right in pairs]
        center_mean = sum(centers) / len(centers)
        gap_mean = sum(gaps) / len(gaps)
        self.assertAlmostEqual(result["C_mean"], center_mean)
        self.assertAlmostEqual(result["G_mean"], gap_mean)
        self.assertAlmostEqual(
            result["C_variance"],
            sum((value - center_mean) ** 2 for value in centers) / len(centers),
        )
        self.assertAlmostEqual(
            result["G_variance"],
            sum((value - gap_mean) ** 2 for value in gaps) / len(gaps),
        )
        self.assertAlmostEqual(
            result["C_G_covariance"],
            sum(
                (center - center_mean) * (gap - gap_mean)
                for center, gap in zip(centers, gaps)
            )
            / len(pairs),
        )

    def test_aligned_jackknife_covariance_preserves_cross_metric_sign(self) -> None:
        left = [1.0, 2.0, 4.0, 8.0]
        right = [-2.0, -4.0, -8.0, -16.0]
        covariance = jackknife_covariance(left, right)
        self.assertLess(covariance, 0.0)
        self.assertAlmostEqual(covariance, jackknife_covariance(right, left))

    def test_insufficient_archive_is_marked_not_scoreable(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            result = analyze_archive(
                temporary,
                {
                    "N": 7,
                    "histogram": "missing.hist.csv",
                    "moments": "missing.moments.csv",
                    "metadata": "missing.metadata.json",
                    "dependency_group": "missing",
                    "expected_first": [2, 1],
                    "expected_second": [1, 2],
                },
            )
        self.assertEqual(result["public"]["status"], "not_scoreable")
        self.assertIn("missing", result["public"]["reason"])


class ProductionArchiveRegressionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.temporary = tempfile.TemporaryDirectory()
        output = Path(cls.temporary.name)
        cls.output_json = output / "latest.json"
        cls.output_md = output / "latest.md"
        subprocess.run(
            [
                sys.executable,
                str(SCRIPTS / "analyze_two_activation_h4.py"),
                "--manifest",
                str(ROOT / "analysis/two_activation_h4_manifest.yaml"),
                "--output-json",
                str(cls.output_json),
                "--output-md",
                str(cls.output_md),
                "--workers",
                "2",
            ],
            check=True,
            capture_output=True,
            text=True,
        )
        cls.payload = json.loads(cls.output_json.read_text(encoding="utf-8"))

    @classmethod
    def tearDownClass(cls) -> None:
        cls.temporary.cleanup()

    def test_all_declared_sizes_are_scoreable(self) -> None:
        self.assertEqual(self.payload["schema"], SCHEMA)
        self.assertEqual(
            self.payload["scoreable_sizes"],
            [65, 85, 130, 145, 170, 185, 265, 290, 325, 425],
        )
        self.assertEqual(self.payload["not_scoreable_sizes"], [])

    def test_curve_and_root_shift_identities_hold(self) -> None:
        for n in self.payload["scoreable_sizes"]:
            row = self.payload["by_N"][str(n)]
            curves = row["activation_curves_at_p_bar"]
            delta = curves["delta_first_minus_second"]
            shifts = row["root_shift_decomposition"]
            self.assertAlmostEqual(
                delta["delta_M"], delta["delta_F1"] + delta["delta_F2"], places=14
            )
            self.assertAlmostEqual(
                shifts["linearized_root_gap"],
                shifts["delta_p1"] + shifts["delta_p2"],
                places=14,
            )
            self.assertLess(abs(row["identity_audit"]["pooled_M_at_p_bar"]), 3e-14)

    def test_full_covariance_and_dependency_blocks_are_explicit(self) -> None:
        for n in self.payload["scoreable_sizes"]:
            row = self.payload["by_N"][str(n)]
            self.assertEqual(row["estimate_vector_order"], list(FULL_METRICS))
            self.assertEqual(len(row["jackknife_covariance"]), len(FULL_METRICS))
            self.assertTrue(all(len(line) == len(FULL_METRICS) for line in row["jackknife_covariance"]))
            for index in range(len(FULL_METRICS)):
                self.assertGreaterEqual(row["jackknife_covariance"][index][index], -1e-24)

        decision = self.payload["decision_covariance"]
        width = len(self.payload["scoreable_sizes"]) * len(DECISION_METRICS)
        self.assertEqual(len(decision["jackknife_covariance"]), width)
        index = {
            (row["N"], row["metric"]): position
            for position, row in enumerate(decision["metric_order_with_N"])
        }
        shared = decision["jackknife_covariance"][index[(65, "angular_delta_F1")]][
            index[(85, "angular_delta_F1")]
        ]
        independent = decision["jackknife_covariance"][index[(65, "angular_delta_F1")]][
            index[(145, "angular_delta_F1")]
        ]
        self.assertNotEqual(shared, 0.0)
        self.assertEqual(independent, 0.0)

    def test_report_is_a_no_fit_retrospective(self) -> None:
        report = self.output_md.read_text(encoding="utf-8")
        self.assertIn("K1=K_minus", report)
        self.assertIn("generates no Monte Carlo samples and fits no exponent", report)
        self.assertNotIn("not scoreable", report.lower())


if __name__ == "__main__":
    unittest.main()
