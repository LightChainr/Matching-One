from __future__ import annotations

import sys
from pathlib import Path
import math
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from score_p234_projective_jordan import flow_rates, gate_vector, pencil  # noqa: E402


class ProjectiveJordanPencilTests(unittest.TestCase):
    def test_canonical_log_pair_is_lower_unipotent(self) -> None:
        first = [0.0, 2.0, 6.0]
        second = [0.0, 5.0, 25.0]
        transfer = pencil(first, second)
        self.assertEqual(transfer, [[1.0, 0.0], [2.0, 1.0]])
        self.assertEqual(gate_vector(first, second), [0.0, 0.0, 0.0])

    def test_field_rescalings_cancel(self) -> None:
        first = [0.0, 2.0, 6.0]
        second = [0.0, 5.0, 25.0]
        first_raw = [0.0, 2.0 / (3.0 * 7.0), 6.0 / 7.0**2]
        second_raw = [0.0, 5.0 / (11.0 * 13.0), 25.0 / 13.0**2]
        recovered = gate_vector(first_raw, second_raw, 3.0 / 7.0, 11.0 / 13.0)
        for value in recovered:
            self.assertAlmostEqual(value, 0.0)

    def test_ordinary_two_direction_change_fails(self) -> None:
        gate = gate_vector([2.0, 1.0, 3.0], [4.0, 1.0, 9.0])
        self.assertTrue(any(abs(value) > 1e-9 for value in gate))

    def test_two_size_gate_is_blind_to_size_dependent_gauge(self) -> None:
        first_raw = [0.0, 2.0 / (3.0 * 7.0), 6.0 / 7.0**2]
        second_raw = [0.0, 5.0 / (11.0 * 13.0), 25.0 / 13.0**2]
        self.assertEqual(gate_vector(first_raw, second_raw), [0.0, 0.0, 0.0])

    def test_three_sizes_recover_constant_projective_flow(self) -> None:
        sizes = [1.0, 2.0, 4.0, 8.0]
        vectors = [[0.0, 1.0, 3.0 + 2.5 * math.log(size)] for size in sizes]
        rates = flow_rates(vectors, sizes, [1.0] * 4)
        for rate in rates:
            self.assertAlmostEqual(rate, 2.5)


if __name__ == "__main__":
    unittest.main()
