#!/usr/bin/env python3
"""Regression tests for canonical literature provenance artifacts."""

from __future__ import annotations

import csv
from decimal import Decimal
import json
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "data" / "literature_threshold_sources.json"
MERTENS = ROOT / "data" / "mertens_2022_square_site_estimators.csv"
JACOBSEN = ROOT / "data" / "jacobsen_2015_square_site_cylinder.csv"


class LiteratureProvenanceTest(unittest.TestCase):
    def setUp(self) -> None:
        self.manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
        self.sources = {entry["id"]: entry for entry in self.manifest["sources"]}

    def test_mertens_exact_tables_are_complete(self) -> None:
        with MERTENS.open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual([int(row["n"]) for row in rows], list(range(1, 25)))
        self.assertEqual(rows[0]["p_med"], "0.500000000000000000000000000000")
        self.assertEqual(rows[0]["p_cell"], "")
        self.assertEqual(rows[-1]["p_med"], "0.591276289864951685617852112360")
        self.assertEqual(rows[-1]["p_cell"], "0.594703812696743490456949711289")
        self.assertTrue(all(row["p_med_source"] == "Mertens 2022 Table 4" for row in rows))
        self.assertTrue(all(
            row["p_cell_source"] == "Mertens 2022 Table 5"
            for row in rows[1:]
        ))
        self.assertTrue(all(row["decimal_status"] == "exact as printed" for row in rows))
        p_med = [Decimal(row["p_med"]) for row in rows]
        self.assertTrue(all(a < b for a, b in zip(p_med, p_med[1:])))

    def test_jacobsen_sequence_contract_is_preserved(self) -> None:
        with JACOBSEN.open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual([int(row["n"]) for row in rows], list(range(1, 22)))
        self.assertEqual(rows[-1]["value"], "0.592744551481371482002735520463")
        self.assertTrue(all("eigenvalue identity" in row["method"] for row in rows))

    def test_primary_and_pending_2024_claims_are_not_conflated(self) -> None:
        yang = self.sources["yang_zhou_2024_comment"]
        self.assertEqual(yang["source_status"], "primary_abstract_verified_table_pending")
        self.assertIsNone(yang["data_file"])
        self.assertEqual(
            yang["quoted_estimates"][0]["value_text"],
            "0.5927460507896(1)",
        )

        reply = self.sources["jacobsen_2024_reply"]
        self.assertEqual(reply["source_status"], "primary_page1_verified_full_text_pending")
        self.assertNotIn("quoted_estimates", reply)
        self.assertIn("Secondary compilations", reply["secondary_index_note"])

    def test_parenthetical_uncertainties_are_strings_not_combined_intervals(self) -> None:
        for source in self.manifest["sources"]:
            for estimate in source.get("quoted_estimates", []):
                self.assertIsInstance(estimate["central_value"], str)
                self.assertIsInstance(estimate["quoted_uncertainty"], str)
                Decimal(estimate["central_value"])
                Decimal(estimate["quoted_uncertainty"])


if __name__ == "__main__":
    unittest.main()
