from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
sys.path.insert(0, str(SCRIPTS))

from analyze_activation_curve_nodes import (  # noqa: E402
    DECISION_METRICS,
    SCHEMA,
    activation_area,
    bernstein_derivative,
    bernstein_value,
)


class BernsteinIdentityTests(unittest.TestCase):
    def test_value_derivative_and_integral(self) -> None:
        # B_2[0,1,1] = 2p-p^2.
        coefficients = [0.0, 1.0, 1.0]
        self.assertAlmostEqual(bernstein_value(coefficients, 0.3), 0.51)
        self.assertAlmostEqual(bernstein_derivative(coefficients, 0.3), 1.4)
        self.assertAlmostEqual(activation_area(coefficients), 2.0 / 3.0)


class ProductionCurveRegressionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        output = ROOT / "results/activation-curve-nodes"
        cls.output_json = output / "latest.json"
        cls.output_md = output / "latest.md"
        cls.payload = json.loads(cls.output_json.read_text(encoding="utf-8"))

    def test_all_archives_and_curve_identities(self) -> None:
        self.assertEqual(self.payload["schema"], SCHEMA)
        self.assertEqual(
            self.payload["scoreable_sizes"],
            [65, 85, 130, 145, 170, 185, 265, 290, 325, 425],
        )
        self.assertEqual(self.payload["not_scoreable_sizes"], [])
        for n in self.payload["scoreable_sizes"]:
            row = self.payload["by_N"][str(n)]
            self.assertLess(
                max(abs(value) for value in row["identity_audit"].values()), 3e-13
            )
            coefficients = row["full_curve_bernstein"]["D2_coefficients"]
            self.assertAlmostEqual(
                bernstein_value(coefficients, row["p_bar"]),
                row["curve_values_at_p_bar"]["D2"],
                places=14,
            )
            self.assertAlmostEqual(
                bernstein_derivative(coefficients, row["p_bar"]),
                row["curve_values_at_p_bar"]["D2_prime"],
                places=12,
            )

    def test_A2_sign_and_local_negative_nodes(self) -> None:
        findings = self.payload["descriptive_findings"]
        self.assertTrue(findings["A2_has_one_nonzero_sign_across_all_scoreable_sizes"])
        self.assertEqual(findings["A2_sign"], "positive")
        self.assertEqual(findings["negative_D2_at_p_bar_sizes"], [265, 325, 425])
        self.assertEqual(
            findings["negative_points_explained_by_scoreable_nearby_upper_node"],
            [265, 325, 425],
        )
        self.assertTrue(
            findings["upper_node_statement_uses_point_estimate_not_significant_ordering"]
        )
        for n in (265, 325, 425):
            row = self.payload["by_N"][str(n)]
            context = row["local_K2_context"]
            spectrum = row["node_spectra"]["D2"]["critical_window"]
            self.assertEqual(spectrum["status"], "scoreable")
            self.assertEqual(spectrum["point_estimate_branch_count"], 1)
            self.assertGreater(context["nearest_scoreable_K2_node"], row["p_bar"])
            self.assertGreater(row["integrated_areas"]["A2"], 0.0)
            self.assertLess(row["curve_values_at_p_bar"]["D2"], 0.0)

    def test_cross_N_dependency_covariance_is_aligned(self) -> None:
        decision = self.payload["decision_covariance"]
        self.assertTrue(decision["views_must_not_be_added_as_independent_evidence"])
        width = len(self.payload["scoreable_sizes"]) * len(DECISION_METRICS)
        self.assertEqual(len(decision["jackknife_covariance"]), width)
        self.assertTrue(all(len(row) == width for row in decision["jackknife_covariance"]))
        index = {
            (item["N"], item["metric"]): position
            for position, item in enumerate(decision["metric_order_with_N"])
        }
        shared = decision["jackknife_covariance"][index[(65, "A2")]][
            index[(85, "D2_at_p_bar")]
        ]
        independent = decision["jackknife_covariance"][index[(65, "A2")]][
            index[(145, "D2_at_p_bar")]
        ]
        self.assertNotEqual(shared, 0.0)
        self.assertEqual(independent, 0.0)

    def test_report_states_retrospective_boundary(self) -> None:
        report = self.output_md.read_text(encoding="utf-8")
        self.assertIn("generates no Monte Carlo samples and fits no exponent", report)
        self.assertIn("K2/(N+1)=C+W/2", report)


if __name__ == "__main__":
    unittest.main()
