import importlib.util
import pathlib
import unittest
from fractions import Fraction


ROOT = pathlib.Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "coalescence", ROOT / "scripts" / "verify_norm5_conjugate_coalescence.py"
)
MOD = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MOD)


class Norm5ConjugateCoalescenceTests(unittest.TestCase):
    def test_exact_payload(self):
        payload = MOD.verified_payload()
        self.assertEqual(payload[325]["C_smith"], (5, 65))
        self.assertEqual(payload[425]["C_smith"], (5, 85))
        self.assertEqual(payload[325]["parent_sine_ratio"], "6/11")
        self.assertEqual(payload[425]["parent_sine_ratio"], "33/13")

    def test_h4_weights(self):
        self.assertEqual(
            MOD.interpolation_weights((17, 6), (18, 1), (15, 10), 4),
            (Fraction(11, 5), Fraction(-6, 5)),
        )
        self.assertEqual(
            MOD.interpolation_weights((16, 13), (19, 8), (20, 5), 4),
            (Fraction(-13, 20), Fraction(33, 20)),
        )

    def test_conjugate_branches_coalesce(self):
        self.assertEqual(
            MOD.d4_canonical(MOD.mul((8, 1), (2, 1))),
            MOD.d4_canonical(MOD.mul((7, 4), (2, 1))),
        )
        self.assertEqual(
            MOD.d4_canonical(MOD.mul((9, 2), (2, -1))),
            MOD.d4_canonical(MOD.mul((7, 6), (2, -1))),
        )

    def test_all_frozen_harmonic_weights(self):
        expected_325 = {
            8: (Fraction(22517, 44795), Fraction(22278, 44795)),
            12: (Fraction(363263, 7144145), Fraction(6780882, 7144145)),
        }
        expected_425 = {
            8: (Fraction(89531, 242420), Fraction(152889, 242420)),
            12: (Fraction(181189, 68620), Fraction(-112569, 68620)),
        }
        for spin, expected in expected_325.items():
            self.assertEqual(
                MOD.interpolation_weights((17, 6), (18, 1), (15, 10), spin), expected
            )
        for spin, expected in expected_425.items():
            self.assertEqual(
                MOD.interpolation_weights((16, 13), (19, 8), (20, 5), spin), expected
            )


if __name__ == "__main__":
    unittest.main()
