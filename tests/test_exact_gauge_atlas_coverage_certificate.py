import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from exact_gauge_atlas_coverage_certificate import DEFAULT_OUTPUT, frozen_descriptor, validate_result, verify_atlas  # noqa: E402


class ExactGaugeAtlasCoverageCertificateTests(unittest.TestCase):
    def test_frozen_atlas_covers_boundary_exactly(self) -> None:
        result = verify_atlas(frozen_descriptor())
        self.assertTrue(result["complete_for_declared_domain"])
        self.assertEqual(result["uncovered_components"], [])

    def test_missing_boundary_chart_reports_point_gap(self) -> None:
        descriptor = frozen_descriptor()
        descriptor["charts"].pop(1)
        result = verify_atlas(descriptor)
        self.assertFalse(result["complete_for_declared_domain"])
        self.assertEqual(result["uncovered_components"][0]["lower"], "0")
        self.assertEqual(result["uncovered_components"][0]["upper"], "0")

    def test_checked_in_result_reproduces(self) -> None:
        summary = validate_result(json.loads(DEFAULT_OUTPUT.read_text(encoding="utf-8")))
        self.assertEqual(summary["status"], "valid_exact_one_dimensional_gauge_atlas")


if __name__ == "__main__":
    unittest.main()
