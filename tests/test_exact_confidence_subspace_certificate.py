import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from exact_confidence_subspace_certificate import (  # noqa: E402
    DEFAULT_OUTPUT,
    build_result,
    frozen_descriptor,
    validate_result,
    verify_discrepancy,
)


class ExactConfidenceSubspaceCertificateTests(unittest.TestCase):
    def test_frozen_singular_control_is_exact(self) -> None:
        verification = build_result()["verification"]
        self.assertEqual(verification["covariance_rank"], 2)
        self.assertEqual(verification["quadratic_discrepancy"], "5/4")
        self.assertEqual(verification["relation_to_supplied_outer_set"], "inside_or_on")

    def test_residual_outside_covariance_range_fails(self) -> None:
        descriptor = frozen_descriptor()
        descriptor["predicted"][2] = "3"
        with self.assertRaisesRegex(ValueError, "outside covariance range"):
            verify_discrepancy(descriptor)

    def test_invalid_pseudoinverse_fails(self) -> None:
        descriptor = frozen_descriptor()
        descriptor["covariance_pseudoinverse"][1][1] = "1/5"
        with self.assertRaisesRegex(ValueError, "identity failed"):
            verify_discrepancy(descriptor)

    def test_indefinite_covariance_fails(self) -> None:
        descriptor = frozen_descriptor()
        descriptor["covariance"][1][1] = "-4"
        descriptor["covariance_pseudoinverse"][1][1] = "-1/4"
        with self.assertRaisesRegex(ValueError, "positive semidefinite"):
            verify_discrepancy(descriptor)

    def test_exact_outside_relation_is_reported(self) -> None:
        descriptor = frozen_descriptor()
        descriptor["cutoff"] = "1"
        self.assertEqual(verify_discrepancy(descriptor)["relation_to_supplied_outer_set"], "outside")

    def test_checked_in_result_reproduces(self) -> None:
        result = json.loads(DEFAULT_OUTPUT.read_text(encoding="utf-8"))
        summary = validate_result(result)
        self.assertEqual(summary["status"], "valid_exact_confidence_subspace_certificate")
        self.assertEqual(summary["discrepancy"], "5/4")


if __name__ == "__main__":
    unittest.main()
