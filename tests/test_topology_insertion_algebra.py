from __future__ import annotations

from fractions import Fraction
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from topology_insertion_algebra import build_artifact  # noqa: E402


def matmul(left, right):
    return [
        [sum(left[i][k] * right[k][j] for k in range(len(right))) for j in range(len(right[0]))]
        for i in range(len(left))
    ]


def add(*matrices):
    return [
        [sum(matrix[i][j] for matrix in matrices) for j in range(len(matrices[0][0]))]
        for i in range(len(matrices[0]))
    ]


def scale(value, matrix):
    return [[value * item for item in row] for row in matrix]


def identity(n):
    return [[Fraction(i == j) for j in range(n)] for i in range(n)]


def numeric(matrix):
    return [[Fraction(item) for item in row] for row in matrix]


class TopologyInsertionAlgebraTests(unittest.TestCase):
    def test_minimal_polynomials_and_exact_oracle(self) -> None:
        artifact = build_artifact()
        algebra = artifact["finite_algebra"]
        q = numeric(algebra["Q_multiplication"])
        s = numeric(algebra["S_multiplication"])
        d = numeric(algebra["D_multiplication"])
        q2, q3 = matmul(q, q), matmul(matmul(q, q), q)
        self.assertEqual(q3, q)
        s2, s3 = matmul(s, s), matmul(matmul(s, s), s)
        self.assertEqual(add(s3, scale(-3, s2), scale(2, s)), [[Fraction(0)] * 4 for _ in range(4)])
        d3 = matmul(matmul(d, d), d)
        self.assertEqual(d3, d)
        self.assertEqual(matmul(s, d), d)
        self.assertEqual(matmul(d, d), add(scale(2, s), scale(-1, s2)))
        self.assertEqual(sum(row["exact_pointwise_checks"] for row in artifact["tiny_exact_oracle"]), 8960)
        self.assertTrue(all(not row["failures"] for row in artifact["tiny_exact_oracle"]))
        self.assertIn("G(0,0)=0", artifact["theorem"]["marked_source"])


if __name__ == "__main__":
    unittest.main()
