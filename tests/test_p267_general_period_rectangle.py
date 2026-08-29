from __future__ import annotations

import csv
import importlib.util
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

import numpy as np
import yaml


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "p267_score", ROOT / "scripts/score_p267_gaussian_annulus_rectangle.py"
)
SCORE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(SCORE)


class GeneralPeriodRectangleTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        compiler = shutil.which("c++") or shutil.which("g++")
        if compiler is None:
            raise unittest.SkipTest("C++ compiler unavailable")
        cls.temporary = tempfile.TemporaryDirectory()
        cls.directory = Path(cls.temporary.name)
        cls.binary = cls.directory / "general-period-pivotal"
        subprocess.run(
            [compiler, "-std=c++17", "-O1",
             str(ROOT / "src/general_period_multiradius_pivotal_mc.cpp"),
             "-o", str(cls.binary)],
            check=True,
        )

    @classmethod
    def tearDownClass(cls) -> None:
        cls.temporary.cleanup()

    def test_exact_primitive_nonprimitive_and_two_orbit_gate(self) -> None:
        completed = subprocess.run(
            [str(self.binary), "--self-test"], capture_output=True, text=True
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertIn("nonprimitive N=8 Smith=(2,4) configurations=256", completed.stdout)
        self.assertIn("axis+diagonal scalar/spin4 rank=2", completed.stdout)

    def test_all_frozen_R2_landings_are_injective_and_two_orbit(self) -> None:
        completed = subprocess.run(
            [str(self.binary), "--validate-only"], capture_output=True, text=True
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertIn("validated 16 designs at radii 2", completed.stdout)

    def test_common_field_digest_matches_within_equal_N_pairs(self) -> None:
        prefix = self.directory / "smoke"
        completed = subprocess.run(
            [
                str(self.binary), "--samples", "200", "--batches", "10",
                "--radii", "2", "--seed", "26725360829",
                "--replica-offset", "26725300000", "--binary-sha256", "smoke",
                "--design", "p_first,8,-1,1,8",
                "--design", "p_second,7,-4,4,7",
                "--design", "np_first,16,-2,2,16",
                "--design", "np_second,14,-8,8,14",
                "--output-prefix", str(prefix),
            ],
            capture_output=True, text=True,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        with Path(f"{prefix}.batches.csv").open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
        indexed = {(row["label"], int(row["batch"])): row for row in rows}
        for batch in range(10):
            self.assertEqual(
                indexed["p_first", batch]["common_field_digest"],
                indexed["p_second", batch]["common_field_digest"],
            )
            self.assertEqual(
                indexed["np_first", batch]["common_field_digest"],
                indexed["np_second", batch]["common_field_digest"],
            )

    def test_frozen_residuals_annihilate_each_candidate_basis(self) -> None:
        gaussian_coordinate = np.arange(4, dtype=float)
        annulus_coordinate = np.log2(np.asarray(SCORE.RADII) / SCORE.RADII[0])
        for lam in SCORE.LAMBDAS:
            gaussian_basis = np.stack([
                SCORE.radial_basis(lam, value) for value in gaussian_coordinate
            ])
            annulus_basis = np.stack([
                SCORE.radial_basis(lam, value) for value in annulus_coordinate
            ])
            self.assertTrue(np.allclose(SCORE.gaussian_row(lam) @ gaussian_basis, 0))
            self.assertTrue(np.allclose(SCORE.annulus_row(lam) @ annulus_basis, 0))

    def test_synthetic_shared_and_context_enriched_recovery(self) -> None:
        result = SCORE.synthetic_recovery()
        self.assertEqual(result["shared_truth"]["recovered_shared"], "1/2")
        self.assertEqual(
            result["context_enriched_truth"]["recovered_pair"],
            {"Gaussian": "1/2", "annulus": "1"},
        )
        self.assertGreater(result["context_enriched_truth"]["delta"], 100)

    def test_manifest_records_preregistered_production_binary(self) -> None:
        manifest = yaml.safe_load(
            (ROOT / "analysis/p267_gaussian_annulus_missing_cells_20260829.yaml")
            .read_text(encoding="utf-8")
        )
        self.assertTrue(manifest["production_authorized"])
        self.assertEqual(len(manifest["authorization_gate"]), 4)
        self.assertEqual(
            manifest["production_freeze"]["ARM64_binary_sha256"],
            "f273763dea4736db894f0074a125c52debe78a2eb1c6aa4ecef53481f096fdbb",
        )


if __name__ == "__main__":
    unittest.main()
