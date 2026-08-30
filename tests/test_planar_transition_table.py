import copy
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from planar_transition_table import (  # noqa: E402
    DEFAULT_OUTPUT,
    build_result,
    build_width_table,
    validate_result,
)


class PlanarTransitionTableTests(unittest.TestCase):
    def test_width_two_records_are_deterministic(self):
        summary, serialized = build_width_table(2)
        self.assertEqual(summary["states"], 2)
        self.assertEqual(summary["cases_per_operation"], 4)
        self.assertEqual(summary["serialized_records"], 8)
        self.assertEqual(hashlib_sha256(serialized), summary["canonical_jsonl_sha256"])

    def test_checked_in_contract_reproduces_exactly(self):
        checked_in = json.loads(DEFAULT_OUTPUT.read_text(encoding="utf-8"))
        expected = build_result(8)
        self.assertEqual(checked_in, expected)
        summary = validate_result(checked_in)
        self.assertEqual(summary["maximum_width"], 8)
        self.assertEqual(summary["serialized_records"], 31042)

    def test_counts_match_prior_operation_coverage(self):
        result = build_result(8)
        self.assertEqual(result["totals"]["states"], 2055)
        self.assertEqual(result["totals"]["cases_per_operation"], 15521)
        self.assertEqual(result["totals"]["serialized_records"], 31042)
        for row in result["widths"]:
            self.assertEqual(
                sum(int(count) for count in row["row_unique_target_degree_histogram"].values()),
                row["states"],
            )

    def test_hash_and_boundary_tampering_fail_closed(self):
        result = build_result(4)
        altered = copy.deepcopy(result)
        altered["totals"]["canonical_jsonl_sha256"] = "0" * 64
        with self.assertRaisesRegex(ValueError, "does not exactly reproduce"):
            validate_result(altered)
        boundary = copy.deepcopy(result)
        boundary["claim_boundary"]["parent_issue"] = "closed"
        with self.assertRaisesRegex(ValueError, "does not exactly reproduce"):
            validate_result(boundary)


def hashlib_sha256(value: bytes) -> str:
    import hashlib

    return hashlib.sha256(value).hexdigest()


if __name__ == "__main__":
    unittest.main()
