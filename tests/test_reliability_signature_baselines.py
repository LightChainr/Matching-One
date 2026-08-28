import importlib.util
import pathlib
import unittest
from fractions import Fraction


ROOT = pathlib.Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "reliability_signature_baselines",
    ROOT / "scripts" / "reliability_signature_baselines.py",
)
MOD = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MOD)


class ReliabilitySignatureBaselineTests(unittest.TestCase):
    def test_majority_core_minimum_layers(self):
        self.assertEqual(MOD.beta_majority_layer_count(10, 3, 3), 10)
        self.assertEqual(MOD.beta_majority_layer_count(26, 5, 5), 126)
        self.assertEqual(MOD.beta_majority_layer_count(26, 7, 5), 0)

    def test_beta_signature_is_mirror_symmetric(self):
        d = MOD.beta_majority_domination(26, 5)
        q = MOD.activation_signature_from_domination(d)
        for k in range(27):
            self.assertEqual(d[k] + d[26 - k], 1)
        # q[0] is the rank-zero mass; activation ranks 1..N mirror as r<->N+1-r.
        for r in range(1, 27):
            self.assertEqual(q[r], q[27 - r])

    def test_beta_kappas(self):
        self.assertEqual(MOD.beta_kappas(3), (
            Fraction(-256, 225), Fraction(32768, 16875)
        ))
        self.assertEqual(MOD.beta_kappas(5)[0], Fraction(-131072, 99225))
        self.assertEqual(MOD.beta_kappas(7)[0], Fraction(-4194304, 3006003))

    def test_n10_tangent_interaction_levels(self):
        got = MOD.signed_fourier_level_sums_from_lambda_coefficients({
            1: Fraction(5, 4),
            3: Fraction(0),
            5: Fraction(-4),
        })
        self.assertEqual(got, {
            1: Fraction(5, 8),
            3: Fraction(0),
            5: Fraction(-1, 8),
        })


if __name__ == "__main__":
    unittest.main()
