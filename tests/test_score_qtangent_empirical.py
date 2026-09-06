#!/usr/bin/env python3
"""Lock the #581 empirical control to its committed artifact and its pure machinery."""

from __future__ import annotations

import json
import random
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import score_qtangent_empirical as emp  # noqa: E402
from square_bond_kappa3 import square_bond_pairs  # noqa: E402

#: Exact L=3 numbers from the enumerated gate (results/qtangent-scale-decomposition).
#: ``Cov(wrap, X)`` and ``Cov(wrap, B_even)`` in natural units, from the
#: covariance_split fractions doubled: ambient 18865/131072 * 2, betti -116247365/8589934592 * 2.
EXACT_L3_COV_WRAP_X = 18865 / 131072 * 2
EXACT_L3_COV_WRAP_B_EVEN = -116247365 / 8589934592 * 2
#: X's degree-one Walsh coefficient at L=3 (every bond equals this one value).
EXACT_L3_X_DEGREE_ONE = -8721 / 65536


class PureMachinery(unittest.TestCase):
    def test_coarse_scales_is_nested_and_exhaustive(self) -> None:
        """The filtration must be strictly nested and end at every bond."""
        for length in (3, 4, 8):
            scales = emp.coarse_scales(length, emp.N_SCALES)
            bonds = 2 * length * length
            self.assertEqual(scales[-1][0], (1 << bonds) - 1)
            for (a, _), (b, _) in zip(scales, scales[1:]):
                self.assertTrue(a < b, f"L={length}: not strictly nested")

    def test_identity_t_minus_t_star_is_x(self) -> None:
        """``T - T* = X`` must hold configuration by configuration, exactly."""
        rng = random.Random(1)
        for length in (3, 4, 8):
            pairs = square_bond_pairs(length)
            primal_index = {pair.primal: i for i, pair in enumerate(pairs)}
            bonds = 2 * length * length
            for _ in range(200):
                mask = rng.getrandbits(bonds)
                p = emp._primal_statistics(length, mask, pairs)
                d = emp._primal_statistics(length, emp._dual_mask(mask, pairs, primal_index), pairs)
                twice_t = 2 * p["components"] + p["open_edges"]
                twice_t_star = 2 * d["components"] + d["open_edges"]
                x = p["span_rank"] - 1
                self.assertEqual(twice_t - twice_t_star, 2 * x,
                                 f"L={length}: T - T* != X")

    def test_paired_estimator_telescopes(self) -> None:
        """The per-scale increments must sum to the total covariance."""
        length = 4
        pairs = square_bond_pairs(length)
        primal_index = {pair.primal: i for i, pair in enumerate(pairs)}
        scales = emp.coarse_scales(length, emp.N_SCALES)
        est = emp.estimate_gammas(length, scales, 2000, 3, pairs, primal_index)
        residual = max(abs(x) for row in est["telescoping_residual"] for x in row)
        self.assertLess(residual, 1e-9)

    def test_paired_estimator_vanishes_without_structure(self) -> None:
        """A constant observable has no covariance at any scale."""
        length = 3
        pairs = square_bond_pairs(length)
        primal_index = {pair.primal: i for i, pair in enumerate(pairs)}
        scales = emp.coarse_scales(length, emp.N_SCALES)
        original = emp.evaluate
        try:
            # Monkey-patch evaluate to a constant joint vector.
            emp.evaluate = lambda *a, **k: (1, 0, 0)
            est = emp.estimate_gammas(length, scales, 2000, 5, pairs, primal_index)
        finally:
            emp.evaluate = original
        for gamma in est["gammas"]:
            self.assertAlmostEqual(gamma[0][1], 0.0, places=9)
            self.assertAlmostEqual(gamma[0][2], 0.0, places=9)


class ArtifactLock(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = json.loads(emp.DEFAULT_OUTPUT.read_text(encoding="utf-8"))
        cls.by_length = {row["length"]: row for row in cls.payload["sizes"]}

    def test_l3_x_channel_reproduces_the_enumerated_gate(self) -> None:
        row = self.by_length[3]
        self.assertAlmostEqual(row["total_cov_wrap_x"], EXACT_L3_COV_WRAP_X, delta=0.012)

    def test_l3_betti_channel_reproduces_the_enumerated_gate(self) -> None:
        row = self.by_length[3]
        self.assertAlmostEqual(row["total_cov_wrap_b_even"], EXACT_L3_COV_WRAP_B_EVEN, delta=0.012)

    def test_l3_degree_one_x_matches_the_exact_value(self) -> None:
        row = self.by_length[3]
        self.assertAlmostEqual(row["degree_one"]["x_degree_one_value"],
                               EXACT_L3_X_DEGREE_ONE, delta=0.012)

    def test_betti_degree_one_stays_small_at_every_size(self) -> None:
        for row in self.payload["sizes"]:
            self.assertLess(row["degree_one"]["betti_even_degree_one_max_abs"], 0.05,
                            f"L={row['length']}: B_even degree-one did not stay near zero")

    def test_profiles_stay_distinguishable_at_every_size(self) -> None:
        for row in self.payload["sizes"]:
            x = abs(row["total_cov_wrap_x"])
            b = abs(row["total_cov_wrap_b_even"])
            self.assertGreater(x, 3 * b,
                               f"L={row['length']}: X coupling did not dominate B_even")
            self.assertTrue(row["x_profile_single_signed"],
                            f"L={row['length']}: X profile lost its single sign")
            self.assertTrue(row["betti_profile_alternates"],
                            f"L={row['length']}: B_even profile stopped alternating")


if __name__ == "__main__":
    unittest.main()
