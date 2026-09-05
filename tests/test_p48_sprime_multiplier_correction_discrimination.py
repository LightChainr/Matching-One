
from __future__ import annotations
from pathlib import Path
import unittest

import mpmath as mp
import yaml


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "predictions" / "p48_sprime_multiplier_correction_discrimination_20260828.yaml"


class P48SprimeMultiplierCorrectionTests(unittest.TestCase):
    def setUp(self) -> None:
        mp.mp.dps = 70
        with ARTIFACT.open(encoding="utf-8") as handle:
            self.frozen = yaml.safe_load(handle)

    def model_ratio(self, name: str, N: int, Q: int) -> mp.mpf:
        model = self.frozen["source_models"][name]
        A = mp.mpf(model["A"])
        B = mp.mpf(model["B"])
        if name == "q2":
            y_parent = A + B / N
            y_child = A + B / (Q * N)
        else:
            y_parent = A + B * mp.log(N)
            y_child = A + B * mp.log(Q * N)
        return mp.power(Q, -mp.mpf(5) / 4) * y_child / y_parent

    def test_pure_normalized_ratios_are_positive(self) -> None:
        self.assertEqual(self.frozen["version"], 2)
        for section, q in (("norm2", 2), ("norm5", 5)):
            reported = mp.mpf(self.frozen[section]["normalized_pure_ratio"])
            expected = mp.power(q, -mp.mpf(5) / 4)
            self.assertGreater(reported, 0)
            self.assertLess(abs(reported - expected), mp.mpf("1e-48"))

    def test_q2_and_jordan_transfers_recompute_without_angular_factor(self) -> None:
        for section, q in (("norm2", 2), ("norm5", 5)):
            for lineage in self.frozen[section]["lineages"].values():
                N = int(lineage["parent_N"])
                for model_name, field in (("q2", "q2_ratio"), ("jordan_log", "jordan_log_ratio")):
                    expected = self.model_ratio(model_name, N, q)
                    reported = mp.mpf(lineage[field])
                    self.assertGreater(reported, 0)
                    self.assertLess(abs(reported - expected), mp.mpf("1e-29"))

    def test_raw_angular_sign_is_documented_but_not_applied_to_p4(self) -> None:
        self.assertEqual(self.frozen["norm2"]["raw_H4_angular_ratio"], "-1")
        self.assertEqual(self.frozen["norm5"]["raw_H4_angular_ratio"], "-14/25")
        self.assertIn(
            "No H4 angular factor",
            self.frozen["observable"]["normalized_multiplier_rule"],
        )


if __name__ == "__main__":
    unittest.main()
