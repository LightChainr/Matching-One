from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys
import unittest

try:
    import numpy as np
    import yaml
except ModuleNotFoundError as error:
    raise unittest.SkipTest(f"optional numerical analysis dependencies unavailable: {error}")


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from score_etop_production_elimination import (  # noqa: E402
    OUTPUT_SCHEMA,
    build_exact_transform,
    build_report,
    load_manifest,
)


class EtopExactTransformTests(unittest.TestCase):
    def test_integer_transform_has_declared_state_coordinates(self) -> None:
        order = [
            {"N": 7, "metric": "angular_delta_F1"},
            {"N": 7, "metric": "angular_delta_F2"},
            {"N": 7, "metric": "angular_delta_M"},
        ]
        transform, output, _index = build_exact_transform(order, [7])
        observed = transform @ np.asarray([2.0, 5.0, 7.0])
        self.assertEqual([item["coordinate"] for item in output], ["A_top", "E_top", "F1", "F2"])
        self.assertTrue(np.array_equal(observed, np.asarray([7.0, 3.0, 2.0, 5.0])))
        self.assertTrue(np.array_equal(transform, transform.astype(int)))


class EtopProductionRegressionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.input_path = ROOT / "results/two-activation-h4/latest.json"
        cls.manifest_path = ROOT / "analysis/etop_production_elimination_manifest.yaml"
        cls.manifest = load_manifest(cls.manifest_path)
        cls.payload = json.loads(cls.input_path.read_text(encoding="utf-8"))
        cls.report = build_report(
            cls.payload,
            cls.manifest,
            input_path=cls.input_path,
            manifest_path=cls.manifest_path,
        )

    def test_input_is_content_and_schema_pinned(self) -> None:
        digest = hashlib.sha256(self.input_path.read_bytes()).hexdigest()
        self.assertEqual(digest, self.manifest["input"]["sha256"])
        self.assertEqual(self.payload["schema"], self.manifest["input"]["schema"])
        self.assertEqual(self.report["schema"], OUTPUT_SCHEMA)

    def test_full_transformed_covariance_and_dependencies_are_preserved(self) -> None:
        block = self.report["transformed_production_block"]
        self.assertEqual(block["dimension"], 40)
        self.assertEqual(len(block["covariance"]), 40)
        self.assertTrue(all(len(row) == 40 for row in block["covariance"]))
        self.assertEqual(block["spectrum"]["rank"], 20)
        self.assertEqual(block["spectrum"]["discarded_modes"], 20)
        order = self.report["exact_transform"]["output_coordinate_order"]
        index = {(row["N"], row["coordinate"]): position for position, row in enumerate(order)}
        covariance = block["covariance"]
        self.assertNotEqual(covariance[index[(65, "E_top")]][index[(85, "E_top")]], 0.0)
        self.assertEqual(covariance[index[(65, "E_top")]][index[(145, "E_top")]], 0.0)
        self.assertTrue(self.report["exact_transform"]["identity_audit"]["passed"])

    def test_frozen_models_reproduce_production_scores(self) -> None:
        models = self.report["models"]
        expected = {
            "M0_PURE_ALEXANDER_ODD": (445.6184112524901, 10),
            "M1_SECOND_ACTIVATION_DIRECTIONAL_RESPONSE_ZERO": (182.9045183223509, 10),
            "M2_FIRST_ACTIVATION_DIRECTIONAL_RESPONSE_ZERO": (1041.0486164716365, 10),
            "M3_COMMON_PROJECTIVE_RANK_PLANE_LINE": (28.593006303007932, 9),
            "M4_SINGLE_FIXED_H4_POWER": (37.48203247166376, 9),
        }
        for model_id, (chi_square, degrees) in expected.items():
            score = models[model_id]["score"]
            self.assertAlmostEqual(score["mahalanobis_chi_square"], chi_square, places=8)
            self.assertEqual(score["degrees_of_freedom"], degrees)
            self.assertTrue(score["decision"]["model_image_excluded"])
        self.assertAlmostEqual(
            models["M3_COMMON_PROJECTIVE_RANK_PLANE_LINE"]["fitted_parameters"]["lambda"],
            -0.4327100379427825,
            places=8,
        )
        self.assertAlmostEqual(
            models["M4_SINGLE_FIXED_H4_POWER"]["fitted_parameters"]["c"],
            -0.3119371144635731,
            places=8,
        )

    def test_model_boundaries_do_not_overclaim(self) -> None:
        m3 = self.report["models"]["M3_COMMON_PROJECTIVE_RANK_PLANE_LINE"]["contract"]
        m4 = self.report["models"]["M4_SINGLE_FIXED_H4_POWER"]["contract"]
        self.assertIn("one_common", m3["scope"])
        self.assertIn("without_corrections", m4["scope"])
        self.assertIn(
            "not_claims_about_whether_K1_or_K2_exists",
            self.report["exact_transform"]["semantic_boundary"]["F1_F2"],
        )
        self.assertIn("SOS_infeasibility_certificate", self.report["certificate_boundary"]["not_claimed"])


if __name__ == "__main__":
    unittest.main()
