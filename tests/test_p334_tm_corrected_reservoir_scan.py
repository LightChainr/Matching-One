from __future__ import annotations

from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from p334_tm_corrected_reservoir_scan import (  # noqa: E402
    render_markdown,
    rows_for_order,
    scan_order,
)


class P334CorrectedReservoirScanTests(unittest.TestCase):
    def test_complete_minimal_order_saturates(self) -> None:
        result = scan_order(6)
        self.assertEqual(result["complete_order_row_count"], 4)
        self.assertEqual(result["selected_row_count"], 4)
        self.assertEqual(result["summary"]["failed_rows"], 0)
        self.assertEqual(result["summary"]["saturated_rows"], 4)
        for row in result["rows"]:
            audit = row["corrected_reservoir_orbit_graph"]
            self.assertEqual(audit["maximum_matching"], 192)
            self.assertEqual(audit["Hall_deficiency"], 0)
            self.assertEqual(audit["compression_factor"], 6)

    def test_shards_partition_deterministic_rows(self) -> None:
        rows = rows_for_order(6)
        left = scan_order(6, shard_index=0, shard_count=2)
        right = scan_order(6, shard_index=1, shard_count=2)
        indices = sorted(
            [row["row_index"] for row in left["rows"]]
            + [row["row_index"] for row in right["rows"]]
        )
        self.assertEqual(indices, list(range(len(rows))))
        self.assertIn("bounded exact result", left["scientific_boundary"])

    def test_empty_order_is_explicit(self) -> None:
        result = scan_order(7)
        self.assertEqual(result["selected_row_count"], 0)
        self.assertEqual(result["summary"]["status"], "all_selected_rows_saturate")
        self.assertIn("0` of `0", render_markdown(result))

    def test_invalid_shard_rejected(self) -> None:
        with self.assertRaises(ValueError):
            scan_order(6, shard_index=2, shard_count=2)


if __name__ == "__main__":
    unittest.main()
