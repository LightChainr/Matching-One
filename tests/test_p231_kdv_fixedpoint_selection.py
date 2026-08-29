from __future__ import annotations

import importlib.util
import json
from fractions import Fraction
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
SCRIPT = ROOT / "scripts" / "p231_kdv_fixedpoint_selection.py"
SPEC = importlib.util.spec_from_file_location("p231_fixedpoint", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class VacuumKdVFixedPointSelectionTests(unittest.TestCase):
    def test_exact_cyclotomic_character(self) -> None:
        vector = [MODULE.ONE, MODULE.OMEGA, MODULE.OMEGA2]
        self.assertEqual(MODULE.cycle(vector), MODULE.vector_scale(MODULE.OMEGA, vector))
        self.assertEqual(MODULE.projector_numerator(MODULE.ONE, vector), [MODULE.ZERO] * 3)
        self.assertEqual(MODULE.projector_numerator(MODULE.OMEGA2, vector), [MODULE.ZERO] * 3)
        self.assertEqual(
            MODULE.projector_numerator(MODULE.OMEGA, vector),
            MODULE.vector_scale(MODULE.C3(Fraction(3), Fraction(0)), vector),
        )

    def test_reflection_even_direction_is_two_minus_one_minus_one(self) -> None:
        vector = [MODULE.ONE, MODULE.OMEGA, MODULE.OMEGA2]
        self.assertEqual(
            MODULE.reflection_even(vector),
            [MODULE.C3(Fraction(2), Fraction(0)), MODULE.C3(Fraction(-1), Fraction(0)), MODULE.C3(Fraction(-1), Fraction(0))],
        )

    def test_numeric_pinson_oracle_has_exact_character(self) -> None:
        oracle = MODULE.numerical_rho_oracle(80)
        for value in oracle["errors"].values():
            self.assertLess(float(value), 1e-68)

    def test_artifact_is_reproducible(self) -> None:
        expected = json.loads((ROOT / "predictions" / "p231_kdv_fixedpoint_selection_20260829.json").read_text())
        self.assertEqual(MODULE.analyze(90), expected)


if __name__ == "__main__":
    unittest.main()
