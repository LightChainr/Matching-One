import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

from p334_n9_mark_fiber_theorem import (  # noqa: E402
    fiber_signature,
    image_for_row,
)
from p334_n9_reservoir_obstruction import candidate_rows  # noqa: E402


RESULT = ROOT / "results/p334-n9-mark-fiber-theorem/latest.json"


class TestP334N9MarkFiberTheorem(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = json.loads(RESULT.read_text(encoding="utf-8"))

    def test_committed_exact_arithmetic(self):
        derived = self.result["derived_identities"]
        self.assertEqual(derived["source_demand_in_M_units"], 16)
        self.assertEqual(derived["old_image_in_M_units"], 11)
        self.assertEqual(derived["deficiency_in_M_units"], 5)
        self.assertEqual(derived["two_mark_image_in_M_units"], 48)
        self.assertEqual(derived["two_mark_gain_in_M_units"], 37)

    def test_all_six_rows_share_the_frozen_signature(self):
        rows = self.result["rows"]
        self.assertEqual([row["row_index"] for row in rows], [1, 3, 6, 9, 15, 24])
        self.assertEqual(len({row["fiber_signature"]["sha256"] for row in rows}), 1)
        self.assertTrue(all(row["set_relations"]["two_mark_equals_full_MM"] for row in rows))

    def test_primary_histogram_refutes_complete_old_fibers(self):
        signature = self.result["common_fiber_signature"]
        self.assertEqual(
            signature["old_anchor_partner_count_histogram"],
            {"0": 12, "67": 8, "70": 4, "164": 24},
        )
        self.assertNotIn("432", signature["old_anchor_partner_count_histogram"])
        self.assertEqual(signature["two_mark_anchor_partner_count_histogram"], {"432": 48})

    def test_first_row_recomputes_without_flow(self):
        _index, row = candidate_rows()[0]
        old, two_mark, full = image_for_row(row)
        signature = fiber_signature(row, old, two_mark, full)
        self.assertEqual(len(old), 4752)
        self.assertEqual(len(two_mark), 20736)
        self.assertEqual(two_mark, full)
        self.assertEqual(signature["sha256"], self.result["common_fiber_signature"]["sha256"])


if __name__ == "__main__":
    unittest.main()
