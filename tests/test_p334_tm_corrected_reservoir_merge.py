from __future__ import annotations

from copy import deepcopy
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from p334_tm_corrected_reservoir_merge import (  # noqa: E402
    MERGED_BOUNDARY,
    MERGED_SCHEMA,
    merge_shards,
    render_markdown,
)
from p334_tm_corrected_reservoir_scan import scan_order  # noqa: E402


class P334CorrectedReservoirMergeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.shards = [
            scan_order(6, shard_index=index, shard_count=4)
            for index in range(4)
        ]

    def test_complete_shards_merge_in_deterministic_order(self) -> None:
        result = merge_shards(list(reversed(self.shards)))
        self.assertEqual(result["schema"], MERGED_SCHEMA)
        self.assertEqual(result["order"], 6)
        self.assertEqual(result["shard_count"], 4)
        self.assertEqual(result["merged_shard_indices"], [0, 1, 2, 3])
        self.assertEqual(result["complete_order_row_count"], 4)
        self.assertEqual([row["row_index"] for row in result["rows"]], [0, 1, 2, 3])
        self.assertEqual(
            result["summary"],
            {
                "saturated_rows": 4,
                "failed_rows": 0,
                "minimum_Hall_deficiency": 0,
                "maximum_Hall_deficiency": 0,
                "status": "all_rows_saturate",
            },
        )
        self.assertEqual(result["scientific_boundary"], MERGED_BOUNDARY)
        rendered = render_markdown(result)
        self.assertIn("all_rows_saturate", rendered)
        self.assertIn("4` deterministic shards", rendered)

    def test_duplicate_or_incomplete_shard_indices_are_rejected(self) -> None:
        duplicate = deepcopy(self.shards)
        duplicate[-1]["shard"]["index"] = 2
        with self.assertRaisesRegex(ValueError, "duplicate shard index"):
            merge_shards(duplicate)
        with self.assertRaisesRegex(ValueError, "incomplete"):
            merge_shards(self.shards[:-1])

    def test_schema_order_count_and_boundary_drift_are_rejected(self) -> None:
        mutations = (
            ("schema", lambda rows: rows[0].__setitem__("schema", "wrong")),
            ("orders", lambda rows: rows[0].__setitem__("order", 7)),
            ("counts", lambda rows: rows[0]["shard"].__setitem__("count", 5)),
            (
                "boundary",
                lambda rows: rows[0].__setitem__("scientific_boundary", "drifted"),
            ),
        )
        for message, mutate in mutations:
            with self.subTest(message=message):
                shards = deepcopy(self.shards)
                mutate(shards)
                with self.assertRaisesRegex(ValueError, message):
                    merge_shards(shards)

    def test_row_residue_and_complete_coverage_are_rejected(self) -> None:
        wrong_residue = deepcopy(self.shards)
        wrong_residue[0]["rows"][0]["row_index"] = 1
        with self.assertRaisesRegex(ValueError, "does not belong"):
            merge_shards(wrong_residue)

        missing_row = deepcopy(self.shards)
        missing_row[0]["rows"] = []
        missing_row[0]["selected_row_count"] = 0
        missing_row[0]["summary"] = {
            "saturated_rows": 0,
            "failed_rows": 0,
            "minimum_Hall_deficiency": None,
            "maximum_Hall_deficiency": None,
            "status": "all_selected_rows_saturate",
        }
        with self.assertRaisesRegex(ValueError, "residue class"):
            merge_shards(missing_row)

    def test_selected_count_and_summary_drift_are_rejected(self) -> None:
        selected = deepcopy(self.shards)
        selected[0]["selected_row_count"] = 99
        with self.assertRaisesRegex(ValueError, "selected row count"):
            merge_shards(selected)

        summary = deepcopy(self.shards)
        summary[0]["summary"]["failed_rows"] = 1
        with self.assertRaisesRegex(ValueError, "summary drifted"):
            merge_shards(summary)


if __name__ == "__main__":
    unittest.main()
