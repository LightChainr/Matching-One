from __future__ import annotations

import importlib.util
import json
from fractions import Fraction
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "p263_linearized_crossing_oracle.py"
SPEC = importlib.util.spec_from_file_location("p263_crossing", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class LinearizedCrossingOracleTests(unittest.TestCase):
    def test_rank2_and_rank3_jet_dimensions(self) -> None:
        columns = MODULE.toy_jet_columns([-2, -1, 1, 2])
        m2 = [[columns["block"][i], columns["d_delta_block"][i]] for i in range(4)]
        m3 = [row + [columns["half_d2_delta_block"][i]] for i, row in enumerate(m2)]
        self.assertEqual(MODULE.rank(m2), 2)
        self.assertEqual(MODULE.rank(m3), 3)

    def test_frozen_rank2_nulls(self) -> None:
        columns = MODULE.toy_jet_columns([-2, -1, 1, 2])
        nulls = [
            [Fraction(-2), Fraction(3), Fraction(-1), Fraction(0)],
            [Fraction(0), Fraction(-1), Fraction(3), Fraction(-2)],
        ]
        for covector in nulls:
            self.assertEqual(MODULE.dot(covector, columns["block"]), 0)
            self.assertEqual(MODULE.dot(covector, columns["d_delta_block"]), 0)

    def test_rank3_has_equal_nonzero_quadratic_residues(self) -> None:
        columns = MODULE.toy_jet_columns([-2, -1, 1, 2])
        left = [Fraction(-2), Fraction(3), Fraction(-1), Fraction(0)]
        right = [Fraction(0), Fraction(-1), Fraction(3), Fraction(-2)]
        self.assertEqual(MODULE.dot(left, columns["half_d2_delta_block"]), -3)
        self.assertEqual(MODULE.dot(right, columns["half_d2_delta_block"]), -3)
        difference = [a - b for a, b in zip(left, right)]
        self.assertEqual(MODULE.dot(difference, columns["half_d2_delta_block"]), 0)

    def test_four_points_reject_a_cubic_fourth_direction(self) -> None:
        t = [-2, -1, 1, 2]
        cubic = [Fraction(value**3) for value in t]
        rank3_null = [Fraction(-2), Fraction(4), Fraction(-4), Fraction(2)]
        self.assertNotEqual(MODULE.dot(rank3_null, cubic), 0)

    def test_artifact_is_reproducible(self) -> None:
        expected = json.loads((ROOT / "predictions" / "p263_linearized_crossing_toy_20260829.json").read_text())
        self.assertEqual(MODULE.analyze(), expected)


if __name__ == "__main__":
    unittest.main()
