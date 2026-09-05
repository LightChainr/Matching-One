
from __future__ import annotations
import csv
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from square_bond_primitive_pilot import (  # noqa: E402
    CATEGORIES,
    PILOT_DESIGNS,
    analyze_design,
    exact_oracle,
    run_pilot_batches,
)


class SquareBondPrimitiveExactTests(unittest.TestCase):
    def test_n4_pell_exact_oracle(self) -> None:
        result = exact_oracle()
        self.assertTrue(result["passed"])
        self.assertEqual(result["N_vertices"], 4)
        self.assertEqual(result["N_bonds"], 8)
        self.assertEqual(
            result["counts"],
            {
                "rank0": 75,
                "l0": 57,
                "l1": 24,
                "l2": 24,
                "rank1_other": 1,
                "rank2": 75,
                "invariant_failure": 0,
            },
        )
        self.assertEqual(
            result["rank1_other_lines"],
            [{"engine_winding": [1, -2], "count": 1}],
        )


class SquareBondPrimitivePilotTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.rows = run_pilot_batches(
            samples=400,
            batches=4,
            seed=17,
            workers=1,
        )

    def test_batches_are_integer_and_exhaustive(self) -> None:
        self.assertEqual(len(self.rows), 8)
        for row in self.rows:
            self.assertEqual(set(row.counts), set(CATEGORIES))
            self.assertEqual(sum(row.counts.values()), row.samples)
            self.assertEqual(row.counts["invariant_failure"], 0)

    def test_analysis_uses_declared_sign_map_and_reflection_baseline(self) -> None:
        for identifier, matrix in PILOT_DESIGNS:
            result = analyze_design(identifier, matrix, self.rows, dps=60)
            self.assertEqual(result["engine_windings"], [[1, 0], [0, 1], [1, -1]])
            self.assertEqual(result["paper_types"], [[1, 0], [0, 1], [1, 1]])
            self.assertAlmostEqual(
                result["continuum_baselines"][1],
                result["continuum_baselines"][2],
                places=15,
            )
            self.assertEqual(
                set(result["contrasts"]),
                {"C_nontrivial_real", "Q_reflection_null", "S_scalar"},
            )
            covariance = result["contrast_covariance_of_mean"]
            for first in range(3):
                for second in range(3):
                    self.assertAlmostEqual(
                        covariance[first][second], covariance[second][first]
                    )


class SquareBondPrimitiveArchiveTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.result_dir = (
            ROOT / "results" / "local-20260829" / "P156-square-bond-primitive-pilot"
        )
        cls.payload = json.loads(
            (cls.result_dir / "result.json").read_text(encoding="utf-8")
        )

    def test_archive_has_full_declared_sample_counts(self) -> None:
        self.assertEqual(self.payload["samples_per_design"], 200_000)
        self.assertEqual(self.payload["batches_per_design"], 100)
        self.assertEqual(self.payload["seed"], 20260829)
        self.assertTrue(self.payload["exact_oracle"]["passed"])
        for design in self.payload["pilot"]:
            self.assertEqual(sum(design["category_counts"].values()), 200_000)
            self.assertEqual(design["category_counts"]["invariant_failure"], 0)
            self.assertEqual(len(design["continuum_baselines_50dps"]), 3)

    def test_archive_counts_equal_batch_sufficient_statistics(self) -> None:
        totals = {}
        batch_counts = {}
        with (self.result_dir / "result.batches.csv").open(
            newline="", encoding="utf-8"
        ) as handle:
            for row in csv.DictReader(handle):
                identifier = row["design"]
                totals[identifier] = totals.get(identifier, 0) + int(row["samples"])
                counts = batch_counts.setdefault(
                    identifier, {category: 0 for category in CATEGORIES}
                )
                for category in CATEGORIES:
                    counts[category] += int(row[category])
        for design in self.payload["pilot"]:
            identifier = design["design"]
            self.assertEqual(totals[identifier], 200_000)
            self.assertEqual(batch_counts[identifier], design["category_counts"])

    def test_archive_reports_one_nontrivial_mode_and_reflection_null(self) -> None:
        for design in self.payload["pilot"]:
            contrasts = design["contrasts"]
            self.assertGreater(contrasts["C_nontrivial_real"]["value"], 0)
            self.assertLess(abs(contrasts["Q_reflection_null"]["z"]), 2)
            self.assertIn(
                "conjugate",
                " ".join(self.payload["interpretation_boundary"]),
            )


if __name__ == "__main__":
    unittest.main()
