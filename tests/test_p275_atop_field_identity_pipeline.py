from __future__ import annotations

import csv
import json
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

import numpy as np
import yaml

from score_p275_atop_field_identity import _design_matrix, fit_models


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "src/atop_field_identity_microcanonical_mc.cpp"
PREDICTION = ROOT / "predictions/p275_atop_q4_field_identity_20260829.yaml"


class P275FieldIdentityPipelineTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        compiler = shutil.which("c++") or shutil.which("g++")
        if compiler is None:
            raise unittest.SkipTest("C++ compiler unavailable")
        cls.temp = tempfile.TemporaryDirectory()
        cls.binary = Path(cls.temp.name) / "atop_mc"
        subprocess.run(
            [compiler, "-std=c++20", "-O2", str(SOURCE), "-o", str(cls.binary)],
            check=True,
        )

    @classmethod
    def tearDownClass(cls):
        cls.temp.cleanup()

    def test_exact_078bd61_and_modulus_sensitive_self_test(self):
        completed = subprocess.run(
            [str(self.binary), "--self-test"], check=True, text=True,
            stdout=subprocess.PIPE,
        )
        self.assertIn("exact 078bd61 covariance", completed.stdout)
        self.assertIn("modulus-sensitive physical H4", completed.stdout)

    def _small_run(self, modulus: str, matrix: list[int], z: list[int], prefix: Path):
        subprocess.run(
            [
                str(self.binary), "--matrix", *(str(value) for value in matrix),
                "--modulus", modulus, "--z", *(str(value) for value in z),
                "--samples", "40", "--batches", "4", "--seed", "2026275050",
                "--replica-offset", "9500000000", "--git-commit", "preflight",
                "--binary-sha256", "preflight-binary", "--output-prefix", str(prefix),
            ], check=True, stdout=subprocess.PIPE, text=True,
        )

    def test_same_N_priority_field_is_byte_identical_across_moduli(self):
        directory = Path(self.temp.name)
        first, second = directory/"i", directory/"2i"
        self._small_run("i", [7,-1,1,7], [7,1], first)
        self._small_run("2i", [4,-6,3,8], [4,3], second)
        digests = []
        for path in (first, second):
            with Path(str(path)+".batches.csv").open(newline="") as handle:
                rows = list(csv.DictReader(handle))
            digests.append({int(row["batch"]): row["priority_field_digest"] for row in rows})
            metadata = json.loads(Path(str(path)+".metadata.json").read_text())
            self.assertEqual(metadata["smith_invariants"], [1,50])
            self.assertEqual(metadata["root_estimator"], "k*last_active_mark+(N-k)*next_inactive_mark")
        self.assertEqual(digests[0], digests[1])

    def test_synthetic_ordinary_Q4_regression(self):
        prediction = yaml.safe_load(PREDICTION.read_text())
        observation = _design_matrix("Q4_epsilon_ordinary", prediction) @ np.array([0.7,-0.2])
        result = fit_models(observation, np.eye(18)*1.0e-4, prediction)
        self.assertEqual(result["selected"], "Q4_epsilon_ordinary")
        self.assertLess(result["scores"]["Q4_epsilon_ordinary"]["chi_square"], 1.0e-15)

    def test_synthetic_Jordan_regression(self):
        prediction = yaml.safe_load(PREDICTION.read_text())
        beta = np.array([0.3,0.12,-0.2, -0.1,0.08,0.15])
        observation = _design_matrix("Q4_energy_Jordan", prediction) @ beta
        result = fit_models(observation, np.eye(18)*1.0e-7, prediction)
        self.assertEqual(result["selected"], "Q4_energy_Jordan")
        self.assertLess(result["scores"]["Q4_energy_Jordan"]["chi_square"], 1.0e-12)
        self.assertLess(result["scores"]["Q4_epsilon_ordinary"]["survival_p"], 1.0e-8)

    def test_phase_contract_has_no_fixed_or_posthoc_p(self):
        prediction = yaml.safe_load(PREDICTION.read_text())
        phase = prediction["phase1_microcanonical_matching_root"]
        self.assertEqual(phase["evaluation_p"], "finite_matching_root_inside_each_delete_one")
        self.assertEqual(phase["microcanonical_levels"], "all_k_from_0_through_N")
        self.assertTrue(prediction["production_authorized"])
        self.assertTrue(phase["production_authorized"])
        self.assertEqual(phase["runner_commit"], "cb83673fb5f221616a47d53f564635c11e7d0680")


if __name__ == "__main__":
    unittest.main()
