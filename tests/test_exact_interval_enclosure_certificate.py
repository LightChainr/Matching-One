import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from exact_interval_enclosure_certificate import DEFAULT_OUTPUT, frozen_descriptor, validate_result, verify_enclosure  # noqa: E402


class ExactIntervalEnclosureCertificateTests(unittest.TestCase):
    def test_frozen_operations_are_exact(self) -> None:
        result = verify_enclosure(frozen_descriptor())
        self.assertEqual(result["sum"], ["1/12", "7/6"])
        self.assertEqual(result["product"], ["-1/8", "1/3"])

    def test_sign_crossing_product_uses_all_endpoints(self) -> None:
        descriptor = frozen_descriptor()
        descriptor["declared_product"] = ["-1/12", "1/3"]
        with self.assertRaisesRegex(ValueError, "product"):
            verify_enclosure(descriptor)

    def test_bad_square_root_enclosure_fails(self) -> None:
        descriptor = frozen_descriptor()
        descriptor["sqrt_interval"] = ["3/2", "2"]
        with self.assertRaisesRegex(ValueError, "does not enclose"):
            verify_enclosure(descriptor)

    def test_checked_in_result_reproduces(self) -> None:
        summary = validate_result(json.loads(DEFAULT_OUTPUT.read_text(encoding="utf-8")))
        self.assertEqual(summary["status"], "valid_exact_interval_enclosure")


if __name__ == "__main__":
    unittest.main()
