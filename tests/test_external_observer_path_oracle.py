from __future__ import annotations

from fractions import Fraction
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from external_observer_path_oracle import build_artifact  # noqa: E402


class ExternalObserverPathOracleTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.artifact = build_artifact()

    def test_euler_residue_escapes_q_sigma_algebra(self) -> None:
        self.assertEqual(
            self.artifact["selected_observer"]["definition"],
            "O_ext=C_black_NN-C_white_matching-q=V-E+F0",
        )
        witnesses = [
            geometry["q_sigma_algebra_escape_witness"]
            for geometry in self.artifact["geometries"]
            if geometry["q_sigma_algebra_escape_witness"] is not None
        ]
        self.assertTrue(witnesses)
        for witness in witnesses:
            self.assertIn(witness["q"], (-1, 0, 1))
            self.assertNotEqual(witness["first_O_ext"], witness["second_O_ext"])

    def test_exact_path_products_and_index9_backend(self) -> None:
        self.assertIn("axis-L3-index9", [row["id"] for row in self.artifact["geometries"]])
        for geometry in self.artifact["geometries"]:
            for row in geometry["microcanonical_rows"]:
                self.assertTrue(
                    all(Fraction(value) == 0 for value in row["residual"].values())
                )

    def test_qj_is_only_a_contact_control(self) -> None:
        self.assertIn("q_times_J_D4_contact_control", self.artifact["stored_path_statistics"])
        self.assertIn("N325/N425", self.artifact["production_gate"])
        authorization = self.artifact["observer_authorization_gate"]
        self.assertIn("explicitly scalar", authorization["status"])
        self.assertIn("determinant -2", authorization["future_local_H4"])


if __name__ == "__main__":
    unittest.main()
