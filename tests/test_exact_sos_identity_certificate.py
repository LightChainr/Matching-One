import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from exact_sos_identity_certificate import DEFAULT_OUTPUT, frozen_descriptor, validate_result, verify_identity  # noqa: E402


class ExactSosIdentityCertificateTests(unittest.TestCase):
    def test_frozen_localizing_identity_is_exact(self) -> None:
        result = verify_identity(frozen_descriptor())
        self.assertEqual(result["inequality_count"], 1)
        self.assertEqual(result["status"], "exact_sos_localizing_identity_verified")

    def test_tampered_target_fails_closed(self) -> None:
        descriptor = frozen_descriptor()
        descriptor["target"][0]["coefficient"] = "2"
        with self.assertRaisesRegex(ValueError, "does not match"):
            verify_identity(descriptor)

    def test_count_mismatch_fails_closed(self) -> None:
        descriptor = frozen_descriptor()
        descriptor["localizing_squares"] = []
        with self.assertRaisesRegex(ValueError, "count mismatch"):
            verify_identity(descriptor)

    def test_unknown_variable_fails_closed(self) -> None:
        descriptor = frozen_descriptor()
        descriptor["target"].append({"coefficient": "1", "powers": {"z": 1}})
        with self.assertRaisesRegex(ValueError, "unknown"):
            verify_identity(descriptor)

    def test_checked_in_result_reproduces(self) -> None:
        summary = validate_result(json.loads(DEFAULT_OUTPUT.read_text(encoding="utf-8")))
        self.assertEqual(summary["status"], "valid_exact_sos_localizing_identity")


if __name__ == "__main__":
    unittest.main()
