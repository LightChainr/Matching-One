import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from exact_parameter_box_certificate import DEFAULT_OUTPUT, frozen_parameters, validate_result, verify_parameter_box  # noqa: E402


class ExactParameterBoxCertificateTests(unittest.TestCase):
    def test_frozen_box_reports_active_faces(self) -> None:
        result = verify_parameter_box(frozen_parameters())
        self.assertEqual(result["active_boundaries"], ["k:lower", "lambda:upper"])

    def test_outside_value_fails_closed(self) -> None:
        parameters = frozen_parameters()
        parameters[0]["value"] = "2"
        with self.assertRaisesRegex(ValueError, "outside bounds"):
            verify_parameter_box(parameters)

    def test_reversed_bounds_fail_closed(self) -> None:
        parameters = frozen_parameters()
        parameters[0]["lower"] = "2"
        with self.assertRaisesRegex(ValueError, "reversed"):
            verify_parameter_box(parameters)

    def test_missing_provenance_and_duplicate_fail(self) -> None:
        parameters = frozen_parameters()
        parameters[0]["provenance"] = ""
        with self.assertRaisesRegex(ValueError, "provenance"):
            verify_parameter_box(parameters)
        duplicate = frozen_parameters()
        duplicate[1]["name"] = duplicate[0]["name"]
        with self.assertRaisesRegex(ValueError, "unique"):
            verify_parameter_box(duplicate)

    def test_checked_in_result_reproduces(self) -> None:
        summary = validate_result(json.loads(DEFAULT_OUTPUT.read_text(encoding="utf-8")))
        self.assertEqual(summary["status"], "valid_exact_parameter_box")


if __name__ == "__main__":
    unittest.main()
