import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from projective_direct02_decomposition import (  # noqa: E402
    analyze_batch,
    build_artifact,
    synthetic_batch,
    validate_artifact,
)


class ProjectiveDirect02DecompositionTest(unittest.TestCase):
    def test_committed_certificate_reproduces_exactly(self):
        path = ROOT / "analysis" / "projective_direct02_decomposition_certificate.json"
        artifact = json.loads(path.read_text(encoding="utf-8"))
        self.assertEqual(artifact, build_artifact())
        self.assertEqual(validate_artifact(artifact)["orientation_count"], 2)

    def test_exact_total_and_sector_contrasts(self):
        report = analyze_batch(synthetic_batch())
        self.assertEqual(
            report["H4_contrasts"],
            {
                "A_M_total": "-1/12",
                "A_M_DIRECT_RANK2_contribution": "-1/6",
                "A_M_plateau_contribution": "1/12",
            },
        )
        self.assertTrue(all(report["exact_checks"].values()))

    def test_first_orientation_exact_mechanism_values(self):
        first = analyze_batch(synthetic_batch())["orientations"]["first"]
        self.assertEqual(first["P_direct02"], "1/3")
        self.assertEqual(first["M_with_direct02"], "1/3")
        self.assertEqual(first["M_without_direct02"], "0")
        self.assertEqual(first["additive_M_contributions"], {"DIRECT_RANK2": "1/3", "plateau": "0"})

    def test_second_orientation_exact_mechanism_values(self):
        second = analyze_batch(synthetic_batch())["orientations"]["second"]
        self.assertEqual(second["P_direct02"], "1/6")
        self.assertEqual(second["M_with_direct02"], "1/6")
        self.assertEqual(second["M_without_direct02"], "1/4")
        self.assertEqual(second["additive_M_contributions"], {"DIRECT_RANK2": "0", "plateau": "1/6"})

    def test_conditioned_value_is_distinct_from_additive_contribution(self):
        second = analyze_batch(synthetic_batch())["orientations"]["second"]
        self.assertNotEqual(second["M_without_direct02"], second["additive_M_contributions"]["plateau"])

    def test_equal_covectors_fail_closed(self):
        batch = synthetic_batch()
        batch["orientations"][1]["covector"] = "-1"
        with self.assertRaises(ValueError):
            analyze_batch(batch)

    def test_unequal_orientation_counts_fail_closed(self):
        batch = synthetic_batch()
        batch["orientations"][1]["rows"].pop()
        with self.assertRaises(ValueError):
            analyze_batch(batch)

    def test_float_threshold_and_malformed_row_fail_closed(self):
        batch = synthetic_batch()
        batch["threshold"] = 0.5
        with self.assertRaises(TypeError):
            analyze_batch(batch)
        batch = synthetic_batch()
        del batch["orientations"][0]["rows"][0]["line"]
        with self.assertRaises(ValueError):
            analyze_batch(batch)


if __name__ == "__main__":
    unittest.main()
