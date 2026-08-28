import json
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from audit_issue43_even_channel_contract import audit  # noqa: E402


class Issue43EvenChannelContractTests(unittest.TestCase):
    def test_real_artifacts_expose_exact_channel_sign_transport(self):
        result_root = ROOT / "results" / "server-20260828" / "P43-heldout-fullcurve-500m"
        result = audit(
            ROOT / "results" / "server-20260828" / "P31"
            / "p31_confirmation_seed2026093001.analysis.csv",
            ROOT / "results" / "server-20260828" / "P31"
            / "p31_confirmation_seed2026093001.batches.csv",
            result_root / "analysis" / "primary_score.json",
            [result_root / "raw" / "n185.metadata.json", result_root / "raw" / "n265.metadata.json"],
        )
        self.assertEqual(
            result["classification"],
            "frozen prediction generator/protocol channel mismatch",
        )
        fits = result["source"]["constant_amplitude_fits"]
        self.assertAlmostEqual(fits["cross"]["mean"], -fits["either"]["mean"], places=15)
        self.assertAlmostEqual(fits["cross"]["mean"], -0.010603216462677735, places=15)
        transport = result["exact_transport"]
        self.assertAlmostEqual(transport["transported_chi_square"], 0.5700315435551193, places=12)
        self.assertAlmostEqual(transport["literal_frozen_positive_chi_square"], 240.2472113766935, places=10)
        self.assertEqual(
            transport["transported_cross_mean"],
            [-value for value in transport["frozen_either_mean"]],
        )
        self.assertEqual(
            result["source"]["configuration_level_batch_check"],
            {"size_count": 5, "batch_count": 500, "integer_identity_max_abs_residual": 0},
        )

    def test_rejects_non_cross_production_metadata(self):
        result_root = ROOT / "results" / "server-20260828" / "P43-heldout-fullcurve-500m"
        with tempfile.TemporaryDirectory() as directory:
            wrong = Path(directory) / "n185.metadata.json"
            metadata = json.loads((result_root / "raw" / "n185.metadata.json").read_text())
            metadata["channel"] = "either wrapping"
            wrong.write_text(json.dumps(metadata))
            with self.assertRaisesRegex(ValueError, "expected cross channel"):
                audit(
                    ROOT / "results" / "server-20260828" / "P31"
                    / "p31_confirmation_seed2026093001.analysis.csv",
                    ROOT / "results" / "server-20260828" / "P31"
                    / "p31_confirmation_seed2026093001.batches.csv",
                    result_root / "analysis" / "primary_score.json",
                    [wrong, result_root / "raw" / "n265.metadata.json"],
                )


if __name__ == "__main__":
    unittest.main()
