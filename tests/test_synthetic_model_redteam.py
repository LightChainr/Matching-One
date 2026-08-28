#!/usr/bin/env python3
from pathlib import Path
import importlib.util
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "synthetic_model_redteam", ROOT / "scripts" / "synthetic_model_redteam.py"
)
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class SyntheticModelRedteamTests(unittest.TestCase):
    def setUp(self):
        self.config, self.rows = MODULE.load_config(
            ROOT / "experiments" / "synthetic_model_redteam_current.json"
        )

    def test_harmonic_recurrence_and_declared_norms(self):
        for row in self.rows:
            for value in (row.first, row.second):
                c4 = MODULE.cosine_harmonic(value, 4)
                c12 = MODULE.cosine_harmonic(value, 12)
                self.assertEqual(c12, 4 * c4**3 - 3 * c4)
            self.assertGreater(row.se, 0)

    def test_truths_cover_declared_models(self):
        truths, coefficients = MODULE.truth_means(self.rows, 100.0, 0.7885, 0.5)
        self.assertEqual(tuple(truths), MODULE.MODEL_ORDER)
        self.assertEqual(tuple(coefficients), MODULE.MODEL_ORDER)
        self.assertTrue(all(len(values) == len(self.rows) for values in truths.values()))
        self.assertEqual(len(coefficients["H4"]), 1)
        self.assertTrue(all(len(coefficients[name]) == 2 for name in MODULE.MODEL_ORDER[2:]))

    def test_simulation_is_reproducible_and_counts_close(self):
        first = MODULE.simulate(self.config, self.rows, 12, 1234, 1.0, 0.5)
        second = MODULE.simulate(self.config, self.rows, 12, 1234, 1.0, 0.5)
        self.assertEqual(first, second)
        for criterion in ("heldout_chi_square", "predictive_deviance"):
            for truth in MODULE.MODEL_ORDER:
                self.assertEqual(sum(first["confusion_counts"][criterion][truth].values()), 12)
                self.assertAlmostEqual(sum(first["confusion_rates"][criterion][truth].values()), 1.0)

    def test_low_noise_identifies_pure_angular_models(self):
        result = MODULE.simulate(self.config, self.rows, 40, 9876, 0.02, 0.5)
        rates = result["correct_selection_rate"]["predictive_deviance"]
        # With only two held-out observations an over-parameterized family can
        # still win on chance residuals, even when the noise is tiny.
        self.assertGreater(rates["H4"], 0.7)
        self.assertGreater(rates["H12"], 0.7)


if __name__ == "__main__":
    unittest.main()
