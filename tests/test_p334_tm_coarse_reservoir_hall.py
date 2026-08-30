from __future__ import annotations

from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from p334_tm_coarse_reservoir_hall import (  # noqa: E402
    build_result,
    capacitated_hall,
)


class P334CoarseReservoirHallTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.result = build_result()

    def test_small_capacitated_failure_returns_exact_cut(self) -> None:
        keys = (
            (0, (0, 0, 1), (0, 0, 1)),
            (1, (0, 0, 1), (0, 0, 1)),
        )
        neighborhoods = (
            frozenset({("YN", 0), ("YN", 1)}),
            frozenset({("YN", 1), ("YN", 2)}),
        )
        audit = capacitated_hall(keys, neighborhoods, demand=2, channel="YN")
        self.assertEqual(audit["maximum_flow"], 3)
        self.assertEqual(audit["Hall_deficiency"], 1)
        cut = audit["minimum_cut_certificate"]
        self.assertEqual(cut["demand"], 4)
        self.assertEqual(cut["neighbor_target_count"], 3)
        self.assertEqual(cut["Hall_deficiency"], 1)

    def test_all_N6_rows_have_exact_second_compression(self) -> None:
        rows = self.result["N6_rows"]
        self.assertEqual(len(rows), 4)
        for row in rows:
            compression = row["source_compression"]
            self.assertEqual(compression["raw_sources"], 1152)
            self.assertEqual(compression["translation_orbit_sources"], 192)
            self.assertEqual(compression["coarse_twin_classes"], 32)
            self.assertEqual(compression["twin_class_size"], 6)
            self.assertEqual(compression["raw_to_coarse_factor"], 36)
            self.assertEqual(compression["neighborhood_equality"], "exhaustively_verified")
            self.assertEqual(compression["extra_twins_compared"], 160)

    def test_N6_is_YN_only_and_MM_cuts_are_exact(self) -> None:
        rows = self.result["N6_rows"]
        self.assertEqual(
            [row["channel_flows"]["MM"]["maximum_flow"] for row in rows],
            [152, 150, 152, 150],
        )
        self.assertEqual(
            [row["channel_flows"]["MM"]["Hall_deficiency"] for row in rows],
            [40, 42, 40, 42],
        )
        for row in rows:
            self.assertEqual(row["pure_channel_classification"], "YN")
            self.assertEqual(row["channel_flows"]["combined"]["maximum_flow"], 192)
            self.assertEqual(row["channel_flows"]["YN"]["maximum_flow"], 192)
            cut = row["channel_flows"]["MM"]["minimum_cut_certificate"]
            self.assertEqual(cut["demand"] - cut["neighbor_target_count"], cut["Hall_deficiency"])

    def test_N8_Smith_gate_is_MM_only_with_YN_min_cut(self) -> None:
        row = self.result["N8_Smith_2_4_gate"]
        compression = row["source_compression"]
        self.assertEqual(compression["raw_sources"], 46080)
        self.assertEqual(compression["translation_orbit_sources"], 5760)
        self.assertEqual(compression["coarse_twin_classes"], 720)
        self.assertEqual(compression["twin_class_size"], 8)
        self.assertEqual(compression["raw_to_coarse_factor"], 64)
        self.assertEqual(row["pure_channel_classification"], "MM")
        self.assertEqual(row["channel_flows"]["combined"]["maximum_flow"], 5760)
        self.assertEqual(row["channel_flows"]["MM"]["maximum_flow"], 5760)
        self.assertEqual(row["channel_flows"]["YN"]["maximum_flow"], 4800)
        cut = row["channel_flows"]["YN"]["minimum_cut_certificate"]
        self.assertEqual(cut["class_count"], 576)
        self.assertEqual(cut["demand"], 4608)
        self.assertEqual(cut["neighbor_target_count"], 3648)
        self.assertEqual(cut["Hall_deficiency"], 960)
        self.assertEqual(cut["replica_histogram"], {0: 144, 1: 144, 2: 144, 3: 144})

    def test_general_theorem_and_boundary_are_explicit(self) -> None:
        theorem = self.result["theorem"]
        self.assertIn("Q_relative", theorem["source_bijection"])
        self.assertIn("N|A|", theorem["capacitated_Hall"])
        self.assertIn("does not prove", self.result["scientific_boundary"])


if __name__ == "__main__":
    unittest.main()
