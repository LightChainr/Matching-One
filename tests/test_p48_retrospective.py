import sys
import unittest
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from analyze_p48_retrospective import (  # noqa: E402
    Histogram,
    covariance_of_mean,
    gls_score,
    pseudovalues,
    tail_and_derivative,
    validate_alignment,
)


class P48RetrospectiveTests(unittest.TestCase):
    def test_threshold_tail_and_derivative(self):
        # A deterministic K=2 threshold at N=3 has tail 3p^2(1-p)+p^3.
        p = 0.37
        value, derivative = tail_and_derivative([0, 0, 1, 0], 1, p)
        self.assertAlmostEqual(value, 3 * p * p * (1 - p) + p ** 3, places=14)
        self.assertAlmostEqual(derivative, 6 * p * (1 - p), places=14)

    def test_pseudovalue_covariance_normalization(self):
        self.assertEqual(pseudovalues(2.0, [1.5, 2.0, 2.5]), [3.0, 2.0, 1.0])
        covariance = covariance_of_mean([[3.0, 6.0], [2.0, 4.0], [1.0, 2.0]])
        self.assertAlmostEqual(covariance[0][0], 1.0 / 3.0)
        self.assertAlmostEqual(covariance[0][1], 2.0 / 3.0)
        self.assertAlmostEqual(covariance[1][1], 4.0 / 3.0)

    def test_gls_retains_train_heldout_cross_covariance(self):
        values = [1.0, 1.2, 1.1]
        covariance = [
            [0.04, 0.00, 0.01],
            [0.00, 0.09, 0.02],
            [0.01, 0.02, 0.16],
        ]
        score = gls_score(values, covariance, [0, 1], [2])
        expected_amplitude = (1.0 / 0.04 * 1.0 + 1.0 / 0.09 * 1.2) / (1.0 / 0.04 + 1.0 / 0.09)
        self.assertAlmostEqual(score["amplitude"], expected_amplitude)
        influence = score["influence"]
        cross = covariance[0][2] * influence[0] + covariance[1][2] * influence[1]
        expected_variance = 0.16 + score["amplitude_se"] ** 2 - 2 * cross
        self.assertAlmostEqual(score["residual_covariance"][0][0], expected_variance)

    def test_alignment_rejects_unsynchronized_batches(self):
        records = {}
        for orientation in ("first", "second"):
            for batch in (0, 1):
                n = 2
                records[(n, orientation, batch)] = Histogram(
                    n, 1, 1, orientation, batch, 1, [0, 1, 0], [0, 1, 0]
                )
        # A second N is missing batch 1, so its counter blocks cannot be
        # deleted synchronously with the first N.
        for orientation in ("first", "second"):
            records[(5, orientation, 0)] = Histogram(
                5, 2, 1, orientation, 0, 1, [0, 1, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0]
            )
        with self.assertRaisesRegex(ValueError, "same contiguous batch ids"):
            validate_alignment(records)


if __name__ == "__main__":
    unittest.main()
