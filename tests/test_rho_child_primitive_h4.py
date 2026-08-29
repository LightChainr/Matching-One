from __future__ import annotations

import sys
from pathlib import Path
import unittest

import mpmath as mp


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from rho_child_primitive_h4_mc import (  # noqa: E402
    CHILD_DESIGNS,
    child_gate,
    counter_mask,
    physical_phase,
    run,
    tiny_oracle,
)
from score_rho_child_primitive_h4 import complex_design, gls, real_vector  # noqa: E402


class RhoChildPrimitiveH4Tests(unittest.TestCase):
    def test_child_gate_and_taus(self):
        gate = child_gate()
        self.assertTrue(gate["passed"])
        self.assertTrue(gate["direction_alias_gate"]["all_rank_two"])
        self.assertTrue(all(row["rank"] == 2 for row in gate["direction_alias_gate"]["children"]))
        expected = [1 + 1.75j, 0.25 + 0.4375j, 0.75 + 0.4375j]
        for (_, matrix), target in zip(CHILD_DESIGNS, expected):
            (a, b), (c, d) = matrix
            self.assertAlmostEqual(complex(b, d) / complex(a, c), target)

    def test_counter_and_phase_contract(self):
        self.assertEqual(counter_mask(267156112, 91, 224), counter_mask(267156112, 91, 224))
        self.assertNotEqual(counter_mask(267156112, 91, 224), counter_mask(267156112, 92, 224))
        matrix = CHILD_DESIGNS[0][1]
        self.assertAlmostEqual(physical_phase(matrix, (1, -2)), physical_phase(matrix, (-1, 2)))

    def test_tiny_exact_oracle(self):
        oracle = tiny_oracle()
        self.assertTrue(oracle["passed"])
        self.assertEqual(oracle["invalid_counts"], [0, 0, 0])
        self.assertEqual(len(oracle["digest_sha256"]), 64)

    def test_single_child_stream_has_one_complex_covariance_block(self):
        rows, summary = run(
            20, 2, 1, 2671562001,
            child="2omega", replica_offset=17000000000,
        )
        self.assertEqual(summary["primary_order"], ["2omega_re", "2omega_im"])
        self.assertEqual(len(summary["full_common_field_covariance_6x6"]), 2)
        self.assertEqual(rows[0]["replica_first"], 17000000000)

    def test_synthetic_gls_selects_named_column(self):
        mp.mp.dps = 60
        zeta = mp.exp(2 * mp.pi * mp.j / 3)
        shapes = [
            [mp.mpc(1), zeta, zeta**2],
            [mp.mpc(1), mp.mpc(1), mp.mpc(1)],
            [mp.mpc(1), zeta**2, zeta],
        ]
        covariance = mp.eye(6) * mp.mpf("1e-4")
        for index, shape in enumerate(shapes):
            y = real_vector([mp.mpc("0.2", "-0.1") * value for value in shape])
            scores = [gls(y, covariance, complex_design([candidate])) for candidate in shapes]
            self.assertLess(mp.mpf(scores[index]["chi_square"]), mp.mpf("1e-40"))
            self.assertTrue(all(mp.mpf(scores[j]["chi_square"]) > 1 for j in range(3) if j != index))


if __name__ == "__main__":
    unittest.main()
