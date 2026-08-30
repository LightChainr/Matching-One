import copy
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from certify_etop_production_rational_confidence import build, validate  # noqa: E402
from model_certificate_envelope import load_manifest, validate_manifest  # noqa: E402


class EtopProductionRationalConfidenceTests(unittest.TestCase):
    def test_production_model_ordering(self):
        result = build()
        self.assertEqual(result["decision_summary"]["eliminated"],
                         ["M_ETOP_ZERO", "M_F2_ZERO", "M_F1_ZERO"])
        self.assertEqual(result["free_ratio_model"]["feasible_witness_r"], "-2/3")
        self.assertTrue(all(row["inside_outer_band"]
                            for row in result["free_ratio_model"]["rows"]))

    def test_checked_in_certificate_reproduces(self):
        path = ROOT / "results/model-certificates/production/etop-rational-confidence/latest.json"
        summary = validate(json.loads(path.read_text(encoding="utf-8")))
        self.assertEqual(summary["production_rows"], 8)

    def test_margin_tampering_fails_closed(self):
        result = build()
        tampered = copy.deepcopy(result)
        tampered["fixed_ratio_models"][0]["rows"][0]["inside_outer_band"] = True
        with self.assertRaisesRegex(ValueError, "does not exactly reproduce"):
            validate(tampered)

    def test_production_envelope_validates(self):
        path = ROOT / "analysis/etop_production_rational_confidence_envelope.yaml"
        summary = validate_manifest(load_manifest(path))
        self.assertEqual(summary["claim_level"], "robust_statistical")
        self.assertEqual(summary["gauge_coverage"], "complete_for_declared_class")


if __name__ == "__main__":
    unittest.main()
