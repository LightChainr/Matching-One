from __future__ import annotations

import shutil
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import mpmath as mp


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"scripts"))

from score_q2_local_d4_uv import load_run, size_estimate, solve_training  # noqa: E402


class Q2LocalD4UVTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        mp.mp.dps = 60
        compiler = shutil.which("c++") or shutil.which("g++")
        if compiler is None:
            raise unittest.SkipTest("C++ compiler unavailable")
        cls.temp = tempfile.TemporaryDirectory()
        cls.binary = Path(cls.temp.name)/"threshold"
        subprocess.run([
            compiler, "-O2", "-DNDEBUG", "-std=c++17",
            str(ROOT/"src/threshold_rank_integer_period_mc.cpp"),
            "-o", str(cls.binary),
        ], check=True)

    @classmethod
    def tearDownClass(cls):
        cls.temp.cleanup()

    def test_exact_local_complement_and_q2_injectivity(self):
        completed = subprocess.run(
            [str(self.binary), "--self-test"], check=True, text=True,
            stdout=subprocess.PIPE,
        )
        self.assertIn("R2 exhaustive local complement", completed.stdout)
        self.assertIn("q2 R2/R4 injectivity", completed.stdout)

    def test_small_stream_closes_both_contact_identities(self):
        prefix = Path(self.temp.name)/"n65"
        subprocess.run([
            str(self.binary), "--marked-births", "--n", "65",
            "--samples", "400", "--batches", "4", "--seed", "27520260829",
            "--replica-offset", "16000000000", "--threads", "1",
            "--git-commit", "preflight", "--binary-sha256", "preflight",
            "--output-prefix", str(prefix),
        ], check=True, stdout=subprocess.PIPE, text=True)
        estimate = size_estimate(load_run(prefix))
        self.assertLess(estimate["contact_max_residual"], mp.mpf("1e-40"))
        self.assertEqual(len(estimate["base_vector"]), 20)

    def test_synthetic_parent_training_predicts_frozen_child(self):
        scale = mp.power(2, -mp.mpf(13)/8)
        beta, amplitude = mp.mpf("0.3"), mp.mpf("-0.7")
        parent_cos = {"first": mp.mpf("0.8"), "second": mp.mpf("-0.4")}
        child_cos = {"first": mp.mpf("-0.8"), "second": mp.mpf("0.4")}
        parent_s = {"first": mp.mpf("0.2"), "second": mp.mpf("0.5")}
        child_s = {"first": mp.mpf("-0.1"), "second": mp.mpf("0.6")}
        parent = {
            "cos4": parent_cos,
            "orientations": {
                side: {
                    "shell_S": parent_s[side],
                    "shell_D": beta*parent_s[side]+amplitude*parent_cos[side],
                } for side in ("first", "second")
            },
        }
        child = {
            "cos4": child_cos,
            "orientations": {
                side: {
                    "shell_S": child_s[side],
                    "shell_D": beta*child_s[side]-scale*amplitude*child_cos[side],
                } for side in ("first", "second")
            },
        }
        result = solve_training(parent, child)
        self.assertAlmostEqual(float(result["thermal_beta"]), float(beta))
        self.assertAlmostEqual(float(result["H4_amplitude_N65"]), float(amplitude))
        self.assertTrue(all(abs(value) < mp.mpf("1e-45") for value in result["child_residual"]))


if __name__ == "__main__":
    unittest.main()
