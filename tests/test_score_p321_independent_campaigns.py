from __future__ import annotations

from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from score_p321_independent_campaigns import combine_scored_campaigns  # noqa: E402


def campaign(mean: float, variance: float, samples: int) -> dict:
    return {
        "N": 144,
        "roots": [mean + index for index in range(5)],
        "root_covariance": [
            [variance if row == column else 0.0 for column in range(5)]
            for row in range(5)
        ],
        "batches": 10,
        "samples_per_shape": samples,
        "square_histograms_byte_identical": True,
        "square_moments_byte_identical": True,
        "elapsed_seconds_all_four_pairs": 1.0,
    }


class IndependentP321CampaignTests(unittest.TestCase):
    def test_inverse_variance_pooling(self) -> None:
        combined = combine_scored_campaigns(
            [campaign(1.0, 4.0, 100), campaign(3.0, 1.0, 400)],
            ["first", "second"],
        )
        self.assertAlmostEqual(combined["roots"][0], 2.6)
        self.assertAlmostEqual(combined["root_covariance"][0][0], 0.8)
        self.assertEqual(combined["samples_per_shape"], 500)
        self.assertEqual(combined["batches"], 20)
        for observed, expected in zip(
            combined["contrasts_to_square"], [1.0, 2.0, 3.0, 4.0]
        ):
            self.assertAlmostEqual(observed, expected)

    def test_mismatched_N_rejected(self) -> None:
        other = campaign(2.0, 1.0, 10)
        other["N"] = 576
        with self.assertRaisesRegex(ValueError, "common N"):
            combine_scored_campaigns([campaign(1.0, 1.0, 10), other], ["a", "b"])


if __name__ == "__main__":
    unittest.main()
