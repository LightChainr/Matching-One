import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from exact_multivariate_ideal_certificate import DEFAULT_OUTPUT, frozen_descriptor, parse_polynomial, validate_result, verify_witness  # noqa: E402


class ExactMultivariateIdealCertificateTests(unittest.TestCase):
    def test_frozen_bezout_identity_produces_one(self) -> None:
        result = verify_witness(frozen_descriptor())
        self.assertEqual(result["result"], [{"coefficient": "1", "powers": {}}])

    def test_tampered_multiplier_fails_closed(self) -> None:
        descriptor = frozen_descriptor()
        descriptor["multipliers"][0][0]["coefficient"] = "2"
        with self.assertRaisesRegex(ValueError, "does not produce one"):
            verify_witness(descriptor)

    def test_unknown_variable_and_negative_power_fail(self) -> None:
        with self.assertRaisesRegex(ValueError, "unknown"):
            parse_polynomial([{"coefficient": "1", "powers": {"z": 1}}], ["x", "y"])
        with self.assertRaisesRegex(ValueError, "nonnegative"):
            parse_polynomial([{"coefficient": "1", "powers": {"x": -1}}], ["x", "y"])

    def test_count_mismatch_fails_closed(self) -> None:
        descriptor = frozen_descriptor()
        descriptor["multipliers"].pop()
        with self.assertRaisesRegex(ValueError, "count mismatch"):
            verify_witness(descriptor)

    def test_checked_in_result_reproduces(self) -> None:
        summary = validate_result(json.loads(DEFAULT_OUTPUT.read_text(encoding="utf-8")))
        self.assertEqual(summary["status"], "valid_exact_multivariate_ideal_certificate")


if __name__ == "__main__":
    unittest.main()
