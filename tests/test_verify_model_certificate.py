import copy
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from verify_model_certificate import (  # noqa: E402
    DEFAULT_OUTPUT,
    build_result,
    frozen_entries,
    validate_result,
    verify_bundle,
)


class VerifyModelCertificateTests(unittest.TestCase):
    def test_frozen_bundle_dispatches_all_types(self) -> None:
        result = build_result()
        self.assertEqual(result["verification"]["certificate_count"], 4)
        self.assertTrue(result["verification"]["all_type_specific_verifiers_pass"])

    def test_noncanonical_or_duplicate_paths_fail_closed(self) -> None:
        entries = frozen_entries()
        noncanonical = copy.deepcopy(entries)
        noncanonical[0]["path"] = "results/model-certificates/framework/../framework/linear-ideal/latest.json"
        with self.assertRaisesRegex(ValueError, "not canonical"):
            verify_bundle(noncanonical)
        duplicate = [entries[0], entries[0]]
        with self.assertRaisesRegex(ValueError, "duplicate"):
            verify_bundle(duplicate)

    def test_checked_in_result_reproduces(self) -> None:
        result = json.loads(DEFAULT_OUTPUT.read_text(encoding="utf-8"))
        summary = validate_result(result)
        self.assertEqual(summary["status"], "valid_fail_closed_model_certificate_dispatcher")
        self.assertEqual(summary["certificate_count"], 4)


if __name__ == "__main__":
    unittest.main()
