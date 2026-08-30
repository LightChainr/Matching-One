from __future__ import annotations

import copy
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from verify_etop_production_model_certificate import verify_certificate  # noqa: E402


class EtopProductionModelCertificateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.manifest = json.loads(
            (ROOT / "analysis/etop_production_model_certificate_manifest.json").read_text()
        )
        cls.certificate = json.loads(
            (ROOT / "results/etop-production-model-certificate/latest.json").read_text()
        )

    def test_checked_certificate_replays_without_remote_source(self) -> None:
        replay = verify_certificate(self.certificate, self.manifest)
        self.assertTrue(replay["verified"])
        self.assertFalse(replay["source_verified"])
        self.assertEqual(
            replay["decisions"]["eliminated"],
            ["M_ETOP_ZERO", "M_F2_ZERO", "M_F1_ZERO"],
        )
        self.assertEqual(replay["decisions"]["not_eliminated"], ["M_COMMON_STATE_LINE"])

    def test_corrupted_model_decision_fails_closed(self) -> None:
        corrupted = copy.deepcopy(self.certificate)
        corrupted["fixed_ratio_models"][0]["decision"] = "not_eliminated"
        with self.assertRaisesRegex(ValueError, "arithmetic mismatch"):
            verify_certificate(corrupted, self.manifest)


if __name__ == "__main__":
    unittest.main()
