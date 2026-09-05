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

    def jacobsen_rows(self) -> dict[int, str]:
        with JACOBSEN.open(newline="", encoding="utf-8") as handle:
            return {int(row["n"]): row["value"] for row in csv.DictReader(handle)}

    def test_jacobsen_sequence_contract_is_preserved(self) -> None:
        with JACOBSEN.open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual([int(row["n"]) for row in rows], list(range(1, 22)))
        self.assertEqual(rows[-1]["value"], "0.592744551481371482002735520463")
        self.assertTrue(all("eigenvalue identity" in row["method"] for row in rows))
        values = [Decimal(row["value"]) for row in rows]
        self.assertTrue(all(a < b for a, b in zip(values, values[1:])))

    def test_jacobsen_rows_keep_the_printed_decimal_count(self) -> None:
        """Table 2 prints 40 decimals for n=1..20 and 30 for n=21.

        Padding or truncating a row still satisfies the file digest and still
        parses to the right number, so nothing else in the suite would notice.
        This is the one row-shape invariant available without the paper in hand.
        """
        source = self.sources["jacobsen_2015_eigenvalue_identity"]
        printed = source["transcription_verification"]["printed_decimal_places"]
        self.assertEqual(printed, {"n_1_to_20": 40, "n_21": 30})
        for n, value in self.jacobsen_rows().items():
            expected = printed["n_21"] if n == 21 else printed["n_1_to_20"]
            self.assertEqual(len(value.split(".")[1]), expected, f"n={n}")

    def test_jacobsen_transcription_corrections_are_applied_and_recorded(self) -> None:
        """The committed rows must be the corrected ones, not the superseded ones.

        A SHA-256 pin only detects drift after a file is committed; it cannot
        detect a digit that was wrong when the file was first written, which is
        how n=4 stayed wrong under a green digest.  These rows were checked
        against the primary source, so the checked values are locked here.
        """
        source = self.sources["jacobsen_2015_eigenvalue_identity"]
        corrections = source["transcription_corrections"]
        self.assertEqual({record["n"] for record in corrections}, {1, 4})
        rows = self.jacobsen_rows()
        for record in corrections:
            self.assertEqual(rows[record["n"]], record["printed"], record["n"])
            self.assertNotEqual(record["previously_committed"], record["printed"])
        quartic = next(record for record in corrections if record["n"] == 4)
        self.assertEqual(quartic["decimal_position"], 20)
        self.assertEqual(
            quartic["printed"],
            "0.5914171708531384817988341017359231779642",
        )
        self.assertNotEqual(
            source["correction_impact"]["content_sha256_before"],
            source["content_sha256"],
        )

    def test_every_committed_table_records_a_digit_level_check(self) -> None:
        """A digest pin says a file has not changed, not that it was ever right.

        Every canonical table therefore has to carry the record of a comparison
        against its primary source, so that adding a table without checking its
        digits fails here rather than years later.
        """
        for source in self.manifest["sources"]:
            if source.get("data_file") is None:
                continue
            verification = source["transcription_verification"]
            self.assertIn("checked", verification, source["id"])
            self.assertIn("coverage", verification, source["id"])
            self.assertIn("method", verification, source["id"])
            self.assertIn("result", verification, source["id"])

    def test_mertens_verification_covers_every_non_empty_cell(self) -> None:
        with MERTENS.open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
        filled = sum(1 for row in rows for key in ("p_med", "p_cell") if row[key])
        self.assertEqual(filled, 47)
        coverage = self.sources["mertens_2022_exact_spanning"]["transcription_verification"]["coverage"]
        self.assertIn(str(filled), coverage)

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
