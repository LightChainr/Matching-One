from __future__ import annotations

from fractions import Fraction
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from marked_birth_path_oracle import build_artifact  # noqa: E402


class MarkedBirthPathOracleTests(unittest.TestCase):
    def test_exact_horvitz_full_source_and_frames(self) -> None:
        artifact = build_artifact()
        self.assertIn("S_active", artifact["exact_boundaries"]["full_source"])
        for geometry in artifact["geometries"]:
            for row in geometry["microcanonical_rows"]:
                self.assertTrue(
                    all(Fraction(value) == 0 for value in row["residual"].values())
                )
                exact = row["exact_absent_site_sum"]
                self.assertEqual(Fraction(exact["active_S"]), Fraction(exact["inactive_S"]))
                self.assertEqual(Fraction(exact["active_D"]), -Fraction(exact["inactive_D"]))
                self.assertEqual(Fraction(exact["site_S"]), Fraction(exact["active_S"]))
                self.assertEqual(Fraction(exact["site_D"]), Fraction(exact["active_D"]))
        axis, gaussian = artifact["geometries"]
        self.assertGreater(Fraction(axis["direct_0_to_2_absent_site_mass_summed_over_k"]), 0)
        self.assertIn(
            {"physical_vector": "2,1", "cos4": "-7/25", "sin4": "24/25"},
            gaussian["lifted_Euclidean_chi4_examples"],
        )


if __name__ == "__main__":
    unittest.main()
