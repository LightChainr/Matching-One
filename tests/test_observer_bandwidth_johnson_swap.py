from __future__ import annotations

from fractions import Fraction
import json
import math
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from observer_bandwidth_johnson_swap import (  # noqa: E402
    build_report,
    expected_multiplicity,
    expected_swap_eigenvalue,
    shifted_nullity,
    swap_matrix,
)


class JohnsonSwapBandwidthTests(unittest.TestCase):
    def test_swap_matrix_is_exactly_stochastic_and_symmetric(self) -> None:
        _, matrix = swap_matrix(6, 3)
        self.assertTrue(all(sum(row) == 1 for row in matrix))
        self.assertEqual(matrix, [list(row) for row in zip(*matrix)])

    def test_small_johnson_spectrum_has_exact_multiplicities(self) -> None:
        _, matrix = swap_matrix(6, 3)
        for degree in range(4):
            eigenvalue = expected_swap_eigenvalue(6, 3, degree)
            self.assertEqual(
                shifted_nullity(matrix, eigenvalue),
                expected_multiplicity(6, degree),
            )

    def test_endpoint_slices_do_not_divide_by_zero(self) -> None:
        self.assertEqual(swap_matrix(7, 0)[1], [[Fraction(1)]])
        self.assertEqual(swap_matrix(7, 7)[1], [[Fraction(1)]])
        self.assertEqual(expected_swap_eigenvalue(7, 0, 0), 1)

    def test_checked_report_is_exactly_reproducible(self) -> None:
        manifest = json.loads(
            (ROOT / "analysis/observer_bandwidth_johnson_swap_manifest.json").read_text(
                encoding="utf-8"
            )
        )
        checked = json.loads(
            (ROOT / "results/observer-bandwidth-johnson-swap/latest.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(checked, build_report(manifest))
        self.assertEqual(sum(row["exact_nullity"] for row in checked["spectrum"]), 70)
        self.assertEqual(
            [row["evaluation_rank"] for row in checked["degree_spaces"]],
            [math.comb(8, degree) for degree in range(5)],
        )


if __name__ == "__main__":
    unittest.main()
