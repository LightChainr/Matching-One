from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from etop_clock_plane_elimination import build_report  # noqa: E402


class EtopClockPlaneEliminationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.report = build_report()

    def test_replays_existing_global_ray(self):
        ray = self.report["primary"]["models"]["A_ray"]["production"]
        self.assertAlmostEqual(ray["chi2"], 15.5147250211756, places=8)
        self.assertAlmostEqual(ray["beta"]["P4_A_top"], -0.59183498, places=6)

    def test_p205_reconstruction_replays_AE(self):
        rows = self.report["primary"]["P205_rows"]
        self.assertEqual([row["N"] for row in rows], [25, 50, 125])
        self.assertLess(max(row["replay_AE_max_abs_error"] for row in rows), 1e-14)

    def test_clock_plane_changes_model_ranking(self):
        models = self.report["primary"]["models"]
        self.assertLess(models["A_ray"]["joint"]["p"], 0.01)
        self.assertGreater(models["A_plus_C"]["joint"]["p"], 0.01)
        self.assertGreater(models["A_plus_C"]["P205_profiled_prediction"]["p"], 0.01)
        self.assertGreater(
            models["A_plus_C_plus_W"]["joint"]["improvement_over_A_plus_C"]["p"],
            0.01,
        )


if __name__ == "__main__":
    unittest.main()
