import copy
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from exact_gauge_chart_certificate import (  # noqa: E402
    DEFAULT_OUTPUT,
    build_result,
    validate_result,
    verify_chart,
)


class ExactGaugeChartCertificateTests(unittest.TestCase):
    def test_frozen_chart_is_verified(self) -> None:
        result = build_result()
        verification = result["verification"]
        self.assertEqual(verification["source_minor"], "1/1024")
        self.assertEqual(verification["amplification_factor"], "1024")
        self.assertTrue(verification["responses_identical"])

    def test_singular_source_fails_closed(self) -> None:
        descriptor = copy.deepcopy(build_result()["descriptor"])
        descriptor["source_matrix"][0][0] = "0"
        with self.assertRaisesRegex(ValueError, "singular"):
            verify_chart(descriptor)

    def test_wrong_similarity_fails_closed(self) -> None:
        descriptor = copy.deepcopy(build_result()["descriptor"])
        descriptor["normalizing_similarity"][0][0] = "1023"
        with self.assertRaisesRegex(ValueError, "exact source inverse"):
            verify_chart(descriptor)

    def test_boundary_must_remain_explicit(self) -> None:
        descriptor = copy.deepcopy(build_result()["descriptor"])
        descriptor["boundary_minor_value"] = "1"
        with self.assertRaisesRegex(ValueError, "zero-minor locus"):
            verify_chart(descriptor)

    def test_checked_in_result_reproduces(self) -> None:
        result = json.loads(DEFAULT_OUTPUT.read_text(encoding="utf-8"))
        summary = validate_result(result)
        self.assertEqual(summary["status"], "valid_exact_gauge_chart_certificate")
        self.assertEqual(summary["coverage"], "nonzero_source_minor_only")


if __name__ == "__main__":
    unittest.main()
