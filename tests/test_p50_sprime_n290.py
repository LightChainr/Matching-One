#!/usr/bin/env python3
from __future__ import annotations

import math
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from score_p50_sprime_n290 import scalar_prediction_score  # noqa: E402


class P50SprimeN290Tests(unittest.TestCase):
    def test_scalar_score_combines_independent_source_and_target_variance(self) -> None:
        score = scalar_prediction_score(3.0, 0.3, 2.0, 0.4)
        self.assertAlmostEqual(score["variance"], 0.25)
        self.assertAlmostEqual(score["signed_z"], 2.0)
        self.assertAlmostEqual(score["two_sided_p"], math.erfc(math.sqrt(2.0)))


if __name__ == "__main__":
    unittest.main()
