from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from exact_gauge_boundary_oracle import DEFAULT_OUTPUT, build_result, validate_result  # noqa: E402


class ExactGaugeBoundaryOracleTests(unittest.TestCase):
    def test_reachable_source_chart_is_open_but_close_to_boundary(self) -> None:
        chart = build_result()["reachable_source_chart"]
        self.assertEqual(chart["source_minor"], "1/1024")
        self.assertIs(chart["chart_open"], True)
        self.assertEqual(chart["boundary_value"], "0")
        self.assertEqual(chart["largest_normalized_entry"], "1024")

    def test_gauge_normalization_preserves_all_declared_observables(self) -> None:
        checks = build_result()["exact_invariant_checks"]
        self.assertIs(checks["responses_identical"], True)
        self.assertEqual(checks["original_response_rows"], checks["normalized_response_rows"])
        self.assertEqual(checks["trace_before_after"], ["3", "3"])
        self.assertEqual(checks["determinant_before_after"], ["2", "2"])

    def test_certificate_keeps_the_boundary_uncovered(self) -> None:
        result = build_result()
        self.assertEqual(result["boundary_diagnostic"]["at_epsilon_zero"], "reachable-source chart undefined")
        self.assertIn("coverage", result["claim_boundary"]["excluded"])


if __name__ == "__main__":
    unittest.main()
