
from __future__ import annotations
import math
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from discover_matrix_semigroup import _covariance, _design, _gls  # noqa: E402


class MatrixSemigroupDiscoveryTests(unittest.TestCase):
    def test_affine_cross_either_involution(self) -> None:
        matrix = ((1.0, 0.0), (1.0, -1.0))
        square = [[sum(matrix[i][k] * matrix[k][j] for k in range(2))
                   for j in range(2)] for i in range(2)]
        self.assertEqual(square, [[1.0, 0.0], [0.0, 1.0]])
        s = 0.37
        self.assertAlmostEqual(matrix[1][0] + matrix[1][1] * s, 1.0 - s)

    def test_three_point_rank1_and_jordan_nulls(self) -> None:
        n, a, b = 65.0, 3.2, -91.0
        ordinary = [a + b / n, a + b / (2 * n), a + b / (5 * n)]
        self.assertAlmostEqual(3 * ordinary[0] - 8 * ordinary[1] + 5 * ordinary[2], 0.0)
        jordan = [a + b * math.log(n), a + b * math.log(2 * n), a + b * math.log(5 * n)]
        contrast = 3 * jordan[0] - 8 * jordan[1] + 5 * jordan[2]
        self.assertAlmostEqual(contrast, b * (-8 * math.log(2) + 5 * math.log(5)))

    def test_feature_matrices_form_semigroups(self) -> None:
        q1, q2 = 2.0, 5.0
        ordinary_product = (q1 ** -1) * (q2 ** -1)
        self.assertAlmostEqual(ordinary_product, (q1 * q2) ** -1)
        self.assertAlmostEqual(math.log(q1) + math.log(q2), math.log(q1 * q2))

    def test_gls_recovers_ordinary_feature(self) -> None:
        labels = [(65, "T_Su"), (130, "T_Su"), (325, "T_Su")]
        design, names = _design(labels, "ordinary_Su")
        values = [3.0 - 90.0 / n for n, _metric in labels]
        fit = _gls(values, [[0.01 if i == j else 0.0 for j in range(3)] for i in range(3)], design, names)
        self.assertAlmostEqual(fit["parameters"]["T_Su:constant"], 3.0)
        self.assertAlmostEqual(fit["parameters"]["T_Su:inverse_N"], -90.0)
        self.assertAlmostEqual(fit["chi_square"], 0.0, places=10)

    def test_delete_one_covariance(self) -> None:
        covariance = _covariance([[1.0, 2.0], [2.0, 4.0], [3.0, 6.0]])
        self.assertAlmostEqual(covariance[0][1], 2 * covariance[0][0])
        self.assertAlmostEqual(covariance[1][1], 4 * covariance[0][0])


if __name__ == "__main__":
    unittest.main()
