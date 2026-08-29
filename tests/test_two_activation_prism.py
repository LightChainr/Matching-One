from __future__ import annotations

import math
import sys
import unittest
from fractions import Fraction
from pathlib import Path

import mpmath as mp


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from analyze_two_activation_h4 import Archive, HistogramBatch  # noqa: E402
from analyze_two_activation_prism import (  # noqa: E402
    block_covariance,
    component_model_vector,
    fit_component_pair,
    fixed_p_estimate,
    score_archive,
)


def _batch(orientation: str, batch: int, k1_rank: int, k2_rank: int) -> HistogramBatch:
    k1 = [0, 0, 0, 0, 0]
    k2 = [0, 0, 0, 0, 0]
    k1[k1_rank] = 10
    k2[k2_rank] = 10
    return HistogramBatch(
        n=4,
        a=2 if orientation == "first" else 1,
        b=0 if orientation == "first" else 1,
        orientation=orientation,
        batch=batch,
        samples=10,
        k1=tuple(k1),
        k2=tuple(k2),
    )


def _synthetic_archive() -> Archive:
    first = tuple(
        _batch("first", batch, k1_rank, k2_rank)
        for batch, (k1_rank, k2_rank) in enumerate(((1, 3), (2, 3), (1, 4), (2, 4)))
    )
    second = tuple(
        _batch("second", batch, k1_rank, k2_rank)
        for batch, (k1_rank, k2_rank) in enumerate(((2, 4), (2, 4), (3, 4), (2, 3)))
    )
    return Archive(
        n=4,
        dependency_group="synthetic-aligned",
        histograms={"first": first, "second": second},
        moments={},
        metadata={},
        paths={},
    )


class FixedProbabilityActivationTests(unittest.TestCase):
    def test_components_exactly_reconstruct_matching_contrast(self) -> None:
        point = fixed_p_estimate(_synthetic_archive(), 0.4)
        self.assertAlmostEqual(point["delta_M"], point["delta_F1"] + point["delta_F2"])
        self.assertAlmostEqual(point["reconstruction_residual"], 0.0)

    def test_delete_one_uses_same_batch_for_both_components(self) -> None:
        archive = _synthetic_archive()
        scored = score_archive(archive, 0.4)
        self.assertEqual(scored["batch_ids"], [0, 1, 2, 3])
        direct = [fixed_p_estimate(archive, 0.4, batch) for batch in range(4)]
        self.assertEqual(scored["deleted"], direct)
        expected_cross = (3 / 4) * sum(
            (row["delta_F1"] - sum(item["delta_F1"] for item in direct) / 4)
            * (row["delta_F2"] - sum(item["delta_F2"] for item in direct) / 4)
            for row in direct
        )
        self.assertAlmostEqual(scored["covariance"][0][1], expected_cross)

    def test_block_covariance_keeps_cross_N_zero(self) -> None:
        by_size = {
            5: {"covariance": [[2.0, 0.5], [0.5, 3.0]]},
            7: {"covariance": [[4.0, -0.25], [-0.25, 5.0]]},
        }
        result = block_covariance(by_size, [5, 7])
        self.assertEqual(result[0][1], 0.5)
        self.assertEqual(result[2][3], -0.25)
        self.assertEqual(result[0][2], 0.0)
        self.assertEqual(result[1][3], 0.0)


class ComponentPairGLSTests(unittest.TestCase):
    def test_joint_gls_recovers_two_amplitudes_with_correlated_blocks(self) -> None:
        first = [mp.mpf("1.0"), mp.mpf("0.5"), mp.mpf("0.2")]
        second = [mp.mpf("0.7"), mp.mpf("-0.3"), mp.mpf("0.4")]
        truth = (mp.mpf("0.8"), mp.mpf("0.3"))
        observed = []
        for left, right in zip(first, second):
            observed.extend([float(truth[0] * left), float(truth[1] * right)])
        covariance = [[0.0 for _ in range(6)] for _ in range(6)]
        for offset in (0, 2, 4):
            covariance[offset][offset] = 0.04
            covariance[offset + 1][offset + 1] = 0.09
            covariance[offset][offset + 1] = covariance[offset + 1][offset] = 0.01
        result = fit_component_pair(
            observed, covariance, first, second, "H4", "H12"
        )
        self.assertAlmostEqual(result["amplitudes"]["A1"], 0.8)
        self.assertAlmostEqual(result["amplitudes"]["A2"], 0.3)
        self.assertLess(result["chi_square"], 1e-28)
        self.assertEqual(result["degrees_of_freedom"], 4)
        self.assertAlmostEqual(result["K2_fitted_fraction"], 3 / 11)

    def test_character_vector_keeps_frozen_exponent_and_exact_sign(self) -> None:
        vector = component_model_vector(
            [25, 50, 125],
            {25: Fraction(1), 50: Fraction(-2), 125: Fraction(3)},
            Fraction(13, 8),
        )
        self.assertGreater(vector[0], 0)
        self.assertLess(vector[1], 0)
        self.assertGreater(vector[2], 0)
        self.assertAlmostEqual(float(vector[0]), 25 ** (-13 / 8))


if __name__ == "__main__":
    unittest.main()
