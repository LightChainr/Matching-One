from __future__ import annotations

import hashlib
import json
from pathlib import Path
import unittest

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
RESULT = ROOT / "results/local-20260830/P250-projective-leg-augmented-map"


class P250AugmentedMapArtifactTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.score = json.loads((RESULT / "score.json").read_text())
        cls.influences = np.load(RESULT / "influences.npz")

    def test_influence_hash_and_shapes(self) -> None:
        digest = hashlib.sha256((RESULT / "influences.npz").read_bytes()).hexdigest()
        self.assertEqual(digest, self.score["influences"]["sha256"])
        self.assertEqual(
            self.influences["old_centered_jackknife_influences"].shape, (5, 400, 70),
        )
        self.assertEqual(
            self.influences["fresh_centered_jackknife_influences"].shape, (5, 400, 70),
        )

    def test_each_covariance_reconstructs_from_independent_sources(self) -> None:
        old = self.influences["old_centered_jackknife_influences"]
        fresh = self.influences["fresh_centered_jackknife_influences"]
        for index, name in enumerate(self.score["candidate_order"]):
            observed = np.asarray(self.score["candidate_maps"][name]["score"]["covariance"])
            expected = old[index].T @ old[index] + fresh[index].T @ fresh[index]
            np.testing.assert_allclose(observed, expected, rtol=1e-13, atol=1e-18)

    def test_candidates_are_unique_and_all_rejected(self) -> None:
        names = list(self.score["candidate_order"])
        self.assertEqual(len(names), len(set(names)))
        self.assertEqual(names, list(self.influences["candidate_names"]))
        self.assertEqual(self.score["surviving_maps"], [])
        self.assertEqual(self.score["alexander_union_decision"], "rejected")
        self.assertFalse(
            self.score["deduplication"]["R2_R3_old_and_fresh_gates_treated_as_independent"]
        )


if __name__ == "__main__":
    unittest.main()
