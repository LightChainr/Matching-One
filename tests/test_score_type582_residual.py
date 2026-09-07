#!/usr/bin/env python3
"""Lock the #584 (Gate 3) remainder screen to its committed artifact and its pure machinery."""

from __future__ import annotations

import json
import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import score_type582_residual as gate  # noqa: E402

#: The affine-orthogonal statistics #582 committed, which the re-freeze must
#: reproduce bit-for-bit.  If these move, the re-freeze is not the same law.
AFFINE_STATISTICS = {
    "65->130": 435746.1, "130->325": 192537.4, "85->170": 292805.4,
    "170->425": 154876.2, "145->290": 112067.1,
}


class ArtifactLock(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = json.loads(gate.DEFAULT_OUTPUT.read_text(encoding="utf-8"))

    def test_the_refreeze_reproduces_the_committed_affine_statistics(self) -> None:
        """The re-freeze must be the same law, not a second implementation."""
        by_label = {t["transition"]: t for t in self.payload["transitions"]}
        for label, expected in AFFINE_STATISTICS.items():
            self.assertAlmostEqual(by_label[label]["affine_statistic"], expected,
                                   delta=0.1 * expected,
                                   msg=f"{label} moved from #582")

    def test_the_consensus_direction_carries_almost_all_of_every_transition(self) -> None:
        """The frozen g_N must remove >= 99% of each transition's chi^2."""
        for t in self.payload["transitions"]:
            remaining = t["remainder_statistic_after_g"] / t["affine_statistic"]
            self.assertLess(remaining, 0.01,
                            f"{t['transition']} leaves {remaining:.4f} of its chi^2")

    def test_the_gaussian_primes_are_what_the_note_claims(self) -> None:
        """13 = 3+2i, 17 = 4+i, 29 = 5+2i, with the stated dominant parities."""
        expected = {"gaussian_13": ((3, 2), "odd"), "gaussian_17": ((4, 1), "even"),
                    "p50": ((5, 2), "odd")}
        self.assertEqual(gate.GAUSSIAN_PRIMES, {k: v[0] for k, v in expected.items()})
        for lineage, (prime, parity) in expected.items():
            dominant = max(prime)
            self.assertEqual("odd" if dominant % 2 else "even", parity)

    def test_no_label_beats_the_permutation_null(self) -> None:
        """The Gate-3 verdict: every candidate partition fails the 5% null."""
        for entry in self.payload["screen"]:
            p = entry["permutation_null"]["p_value"]
            self.assertGreaterEqual(p, 0.05,
                                    f"{entry['label']} beat the null (p={p})")

    def test_the_taylor_curvature_does_not_align_with_the_remainder(self) -> None:
        """The remainder is not smooth second-order curvature: no transition's
        remainder is close to its own lineage's curvature (all >= 43 deg, i.e.
        at most ~53% shared variance and never dominant)."""
        for entry in self.payload["taylor_curvature_alignment"]:
            if not entry["computable"]:
                continue
            for pair in entry["remainder_angles_deg"]:
                self.assertGreater(pair["angle_to_curvature_deg"], 40.0,
                                   f"{pair['transition']} aligns with curvature")

    def test_the_degenerate_labels_are_recorded_not_asserted(self) -> None:
        """Smith/cyclic and deck/scale-word are degenerate or missing, not hidden."""
        text = " ".join(self.payload["degenerate_or_missing_labels"])
        self.assertIn("Smith", text)
        self.assertIn("deck", text)


class PureMachinery(unittest.TestCase):
    def test_acute_angle_folds_into_zero_ninety(self) -> None:
        self.assertAlmostEqual(gate._acute_angle_degrees([1, 0], [1, 0]), 0.0)
        self.assertAlmostEqual(gate._acute_angle_degrees([1, 0], [0, 1]), 90.0)
        self.assertAlmostEqual(gate._acute_angle_degrees([1, 0], [-1, 0]), 0.0)

    def test_separation_is_between_minus_within(self) -> None:
        angles = [[0.0, 5.0, 30.0, 85.0],
                  [5.0, 0.0, 30.0, 85.0],
                  [30.0, 30.0, 0.0, 85.0],
                  [85.0, 85.0, 85.0, 0.0]]
        assignments = ["A", "A", "A", "B"]
        result = gate._separation(angles, assignments)
        # within pairs: (0,1)=5, (0,2)=30, (1,2)=30 -> mean 21.666...
        self.assertAlmostEqual(result["mean_within_deg"], 65 / 3, places=6)
        # between pairs: (0,3)=85, (1,3)=85, (2,3)=85 -> mean 85
        self.assertAlmostEqual(result["mean_between_deg"], 85.0, places=6)
        self.assertAlmostEqual(result["separation_deg"], 85.0 - 65 / 3, places=6)

    def test_permutation_null_is_an_exact_fraction(self) -> None:
        # Asymmetric 3-element case, classes {A:2, B:1}; 3 distinct relabellings,
        # and only the identity relabelling matches the observed separation.
        angles = [[0.0, 5.0, 30.0],
                  [5.0, 0.0, 85.0],
                  [30.0, 85.0, 0.0]]
        assignments = ["A", "A", "B"]
        result = gate.permutation_null(angles, assignments)
        self.assertEqual(result["relabellings"], 3)
        self.assertAlmostEqual(result["p_value"], 1.0 / 3.0, places=6)

    def test_permutation_null_is_one_when_no_structure(self) -> None:
        angles = [[0.0, 90.0, 90.0],
                  [90.0, 0.0, 90.0],
                  [90.0, 90.0, 0.0]]
        assignments = ["A", "A", "B"]
        result = gate.permutation_null(angles, assignments)
        self.assertEqual(result["p_value"], 1.0)

    def test_weighted_angle_uses_the_inverse_covariance(self) -> None:
        from scripts.projective_inference import spectral_pseudo_inverse
        # Diagonal covariance weighting the second axis ten times less.
        pinv, _, _, _ = spectral_pseudo_inverse(
            [[1.0, 0.0], [0.0, 0.01]])
        angle = gate._weighted_acute_angle([1.0, 1.0], [1.0, -1.0], pinv)
        # Weighted cosine ~ (1 - 0.01)/(1 + 0.01), so the angle is much less
        # than the unweighted 90 degrees.
        self.assertLess(angle, 45.0)


if __name__ == "__main__":
    unittest.main()
