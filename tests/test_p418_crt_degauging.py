import importlib.util
import json
from pathlib import Path
import sys
import unittest

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))
spec = importlib.util.spec_from_file_location(
    "p418", ROOT / "scripts" / "score_p418_crt_degauging.py"
)
p418 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(p418)

RESULT = ROOT / "results/huawei-20260830/P418-crt-degauging/score.json"


class TestP418CrtDegauging(unittest.TestCase):
    def test_exact_crt_section_and_masks(self):
        exact = p418.exact_section_and_masks()
        self.assertTrue(exact["passed"])
        self.assertEqual(exact["arithmetic"]["405_mod_101"], 1)
        self.assertEqual(exact["arithmetic"]["405_mod_5"], 0)
        for hand in ("plus", "minus"):
            self.assertEqual(exact["hands"][hand]["gates"]["homomorphism_failures"], 0)
            self.assertEqual(exact["hands"][hand]["residual_phase_counts"], [21, 20, 20, 20, 20])
            for charge in ("r1", "r2"):
                mask = exact["hands"][hand]["masks"][charge]
                self.assertEqual(mask["zero_count"], 0)
                self.assertEqual(mask["A0_residual"], 0.0)
                self.assertGreater(mask["mask_fourier_min_real"], -1e-12)
                self.assertLess(mask["mask_fourier_max_abs_imag"], 1e-12)

    def test_masked_design_is_exact_forward_multiplication(self):
        exact = p418.exact_section_and_masks()
        coordinates = [(0, 0), (1, 0), (0, 1)]
        raw = p418.p406.design(coordinates)
        masked = p418.masked_design(exact, "plus", 1, coordinates)
        weights = np.zeros(101)
        weights[[1, 7, 23]] = [2.0, 0.4, 1.2]
        for index, coordinate in enumerate(coordinates):
            value = p418.mask_complex(exact, "plus", 1, coordinate)
            observed = complex(masked[2 * index] @ weights, masked[2 * index + 1] @ weights)
            latent = complex(raw[2 * index] @ weights, raw[2 * index + 1] @ weights)
            self.assertLess(abs(observed - value * latent), 1e-12)

    def test_published_result_replays_raw_and_rejects_masked(self):
        result = json.loads(RESULT.read_text())
        self.assertTrue(result["raw_P406_replay"]["passed"])
        self.assertEqual(result["raw_P406_replay"]["maximum_distance_error"], 0.0)
        self.assertEqual(result["raw_P406_replay"]["maximum_p_error"], 0.0)
        self.assertEqual(result["decision"], "mask_times_positive_fourier_cone_rejected_in_all_channels")
        for row in result["channels"].values():
            self.assertEqual(row["masked_cone"]["bootstrap_p"], 1 / 251)
            self.assertGreater(row["increment_over_raw_cone"], 500)
            self.assertEqual(row["design"]["rank"], 69)
            self.assertEqual(row["design"]["linear_nullity"], 32)
            self.assertFalse(row["prediction_envelope"]["accepted_model_envelope"])

    def test_full_archive_score_recomputes(self):
        freeze = json.loads((ROOT / "analysis/p418_crt_degauging_freeze.json").read_text())
        paths = [
            (ROOT / row["path"], row["sha256"])
            for row in freeze["inputs"]
        ]
        observed = json.loads(RESULT.read_text())
        self.assertEqual(p418.build_result(paths, 250, 40610120260830), observed)


if __name__ == "__main__":
    unittest.main()
