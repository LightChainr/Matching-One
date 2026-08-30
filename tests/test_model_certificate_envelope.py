import copy
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from model_certificate_envelope import load_manifest, validate_manifest  # noqa: E402


class ModelCertificateEnvelopeTests(unittest.TestCase):
    def setUp(self) -> None:
        self.manifest = load_manifest()

    def test_checked_in_envelope_is_valid_and_complete(self) -> None:
        summary = validate_manifest(self.manifest)
        self.assertEqual(summary["status"], "valid_fail_closed_envelope")
        self.assertEqual(summary["claim_level"], "exact")
        self.assertEqual(summary["input_count"], 1)
        self.assertEqual(summary["gauge_coverage"], "complete_for_declared_class")

    def test_digest_tampering_fails_closed(self) -> None:
        tampered = copy.deepcopy(self.manifest)
        tampered["inputs"][0]["sha256"] = "0" * 64
        with self.assertRaisesRegex(ValueError, "digest mismatch"):
            validate_manifest(tampered)

    def test_unknown_or_missing_fields_fail_closed(self) -> None:
        unknown = copy.deepcopy(self.manifest)
        unknown["unexpected"] = True
        with self.assertRaisesRegex(ValueError, "top-level fields drift"):
            validate_manifest(unknown)
        missing = copy.deepcopy(self.manifest)
        del missing["confidence_set"]
        with self.assertRaisesRegex(ValueError, "top-level fields drift"):
            validate_manifest(missing)

    def test_chart_only_gauge_requires_an_uncovered_set(self) -> None:
        tampered = copy.deepcopy(self.manifest)
        tampered["gauge"] = {"kind": "reachable_source", "coverage": "chart_only", "uncovered_set": ""}
        with self.assertRaisesRegex(ValueError, "uncovered set"):
            validate_manifest(tampered)

    def test_parent_issue_must_remain_open(self) -> None:
        tampered = copy.deepcopy(self.manifest)
        tampered["claim_boundary"]["parent_issue"] = "closed"
        with self.assertRaisesRegex(ValueError, "parent issue boundary drift"):
            validate_manifest(tampered)


if __name__ == "__main__":
    unittest.main()
