from __future__ import annotations

import importlib.util
import json
import math
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
SPEC = importlib.util.spec_from_file_location(
    "score_angular_root_amplitude", ROOT / "scripts" / "score_angular_root_amplitude.py"
)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)
RESULTS = ROOT / "results" / "server-20260828" / "P45-root-amplitude"


class AngularRootAmplitudeScoreTests(unittest.TestCase):
    def test_committed_primary_score_and_algebra(self) -> None:
        payload = json.loads((RESULTS / "score.json").read_text(encoding="utf-8"))
        self.assertTrue(
            math.isclose(payload["frozen_prediction"]["chi_square"], 2.426668731892274)
        )
        self.assertTrue(
            math.isclose(payload["free_common_amplitude"]["value"], 0.41300725114754655)
        )
        for n in ("65", "85"):
            estimate = payload["by_size"][n]["estimate"]
            self.assertTrue(
                math.isclose(
                    estimate["A_p"],
                    estimate["closure_C"] * estimate["A_p_predicted_from_A_M_over_B"],
                    rel_tol=2e-13,
                )
            )

    def test_frozen_chi_square_recomputes(self) -> None:
        payload = json.loads((RESULTS / "score.json").read_text(encoding="utf-8"))
        frozen = payload["frozen_prediction"]
        inverse = MODULE.invert_2x2(frozen["residual_covariance"])
        self.assertTrue(
            math.isclose(
                MODULE.quadratic(frozen["residuals"], inverse),
                frozen["chi_square"],
                rel_tol=1e-14,
            )
        )

    def test_thread_reproducibility_artifacts_are_identical(self) -> None:
        self.assertEqual(
            (RESULTS / "thread1.hist.csv").read_bytes(),
            (RESULTS / "thread4.hist.csv").read_bytes(),
        )
        self.assertEqual(
            (RESULTS / "thread1.moments.csv").read_bytes(),
            (RESULTS / "thread4.moments.csv").read_bytes(),
        )


if __name__ == "__main__":
    unittest.main()
