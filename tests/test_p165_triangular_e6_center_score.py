from __future__ import annotations

import importlib.util
import json
from fractions import Fraction
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
SCRIPT = ROOT / "scripts" / "p165_triangular_e6_center_score.py"
SPEC = importlib.util.spec_from_file_location("p165_e6", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class TriangularE6CenterScoreTests(unittest.TestCase):
    def test_Eisenstein_matrix_and_multiplier(self) -> None:
        m = MODULE.Eisenstein(1, 1)
        matrix = m.multiplication_matrix()
        self.assertEqual(matrix, [[1, -1], [1, 2]])
        self.assertEqual(matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0], m.norm())
        self.assertEqual(m.norm(), 3)
        self.assertEqual(MODULE.inverse_sixth_record(m)["text"], "(-27+0*omega)/729")

    def test_N91_pair_has_exact_nonzero_sextic_contrast(self) -> None:
        first = MODULE.normalized_sixth_phase(MODULE.Eisenstein(1, 9))
        second = MODULE.normalized_sixth_phase(MODULE.Eisenstein(5, 6))
        self.assertEqual(Fraction(first["real"]), Fraction(644221, 753571))
        self.assertEqual(Fraction(second["real"]), Fraction(-716579, 753571))
        self.assertEqual(Fraction(first["real"]) - Fraction(second["real"]), Fraction(1360800, 753571))

    def test_Q6_Ward_coefficient(self) -> None:
        row = MODULE.q6_coefficients()
        self.assertEqual(row["Q6_over_primary_g3"], "-3975/224")
        self.assertEqual(row["Q6_over_primary_pi6_E6"], "-1325/252")
        self.assertEqual(row["Ward_map_rank"], 1)

    def test_existing_center_score_is_nonzero(self) -> None:
        rows = MODULE.exact_center_derivatives(ROOT)
        self.assertEqual([row["M_at_half"] for row in rows], ["0/1", "0/1", "0/1"])
        self.assertEqual([row["Mprime_at_half"] for row in rows], ["3/1", "261/64", "5147/1024"])

    def test_E6_modulus_oracle(self) -> None:
        for value in MODULE.numerical_e6_oracle(80)["errors"].values():
            self.assertLess(float(value), 1e-68)

    def test_artifact_is_reproducible(self) -> None:
        expected = json.loads((ROOT / "predictions" / "p165_triangular_e6_center_score_20260829.json").read_text())
        self.assertEqual(MODULE.analyze(ROOT, 90), expected)


if __name__ == "__main__":
    unittest.main()
