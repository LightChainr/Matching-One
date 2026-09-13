import sys
from fractions import Fraction as F
from pathlib import Path
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from p267_third_geometry_certificate import (certificate, complete_family,
                                             reduced_shape, sigma3, shape_coordinate)


class ThirdGeometryTests(unittest.TestCase):
    def test_finite_family_is_complete_and_minimal(self):
        counts = [len({row["shape"] for row in complete_family(n)})
                  for n in (25, 50, 75, 100, 125)]
        self.assertEqual(counts, [1, 2, 2, 3, 3])
        self.assertEqual({row["shape"] for row in complete_family(50)},
                         {(F(0), F(1)), (F(0), F(2))})

    def test_selected_cells_and_runner_geometry(self):
        data = certificate()
        self.assertEqual(len(data["cells"]), 10)
        self.assertTrue(all(not x["spin4_projector_degenerate"] for x in data["cells"]))
        self.assertEqual([x["delta_chi4"] for x in data["cells"] if x["N"] == 125],
                         ["16128/15625", "-1152/625", "1152/625"])
        self.assertEqual({tuple(x) for x in data["families"]["100"]["shapes"]},
                         {("0", "2"), ("0", "4"), ("1/2", "1")})

    def test_hecke_coefficient_identity_and_display(self):
        for p in (2, 5):
            for n in range(1, 41):
                old = sigma3(n//p) if n % p == 0 else 0
                self.assertEqual(sigma3(p*n)+p**3*old, (1+p**3)*sigma3(n))
        for row in certificate()["exact_conditional_E4_relations"]:
            self.assertLess(abs(row["display_residual"]), 1e-12)

    def test_no_unqualified_real_shape_for_general_tau(self):
        with self.assertRaises(ValueError):
            shape_coordinate((F(1, 3), F(2)))
        self.assertEqual(reduced_shape((F(-1, 2), F(1))), (F(1, 2), F(1)))


if __name__ == "__main__":
    unittest.main()
