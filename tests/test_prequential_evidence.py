import importlib.util
import pathlib
import unittest


ROOT = pathlib.Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "score_prequential_evidence", ROOT / "scripts" / "score_prequential_evidence.py"
)
MOD = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MOD)


class PrequentialEvidenceTest(unittest.TestCase):
    def test_correlated_joint_score_is_not_double_counted_marginals(self):
        residual = [1.0, 1.0]
        joint = MOD.gaussian_score(residual, [[1.0, 0.8], [0.8, 1.0]])
        marginals = sum(MOD.gaussian_score([x], [[1.0]])["nlpd"] for x in residual)
        self.assertNotAlmostEqual(joint["nlpd"], marginals)
        self.assertAlmostEqual(joint["chi_square"], 1.1111111111111112)

    def test_duplicate_primary_raw_group_is_rejected(self):
        block = {
            "raw_data_group": "same", "role": "primary", "status": "SCORED",
            "channel": {"source": "odd", "target": "odd"},
            "observation": [0.0], "models": {"m": {"mean": [0.0], "covariance": [[1.0]]}}
        }
        manifest = {"schema_version": 1, "governance": {}, "blocks": [dict(block, id="a"), dict(block, id="b")]}
        with self.assertRaisesRegex(ValueError, "multiple additive primary views"):
            MOD.score_manifest(manifest)

    def test_channel_mismatch_cannot_be_scored_without_exact_map(self):
        manifest = {"schema_version": 1, "governance": {}, "blocks": [{
            "id": "bad", "raw_data_group": "g", "role": "primary", "status": "SCORED",
            "channel": {"source": "either", "target": "cross", "exact_map": None},
            "observation": [0.0], "models": {}
        }]}
        with self.assertRaisesRegex(ValueError, "channel mismatch"):
            MOD.score_manifest(manifest)

    def test_issue95_regression_nlpds(self):
        analytic = MOD.gaussian_score(
            [-0.09729390621321743, 0.03778282994882565],
            [[0.07075921446106298, 0.07145989058713698], [0.07145989058713698, 0.09757393747678815]],
        )
        jordan = MOD.gaussian_score(
            [-0.3686919350486826, -0.4825156636056427],
            [[0.11443450284008791, 0.1428792399110507], [0.1428792399110507, 0.2096667139352909]],
        )
        self.assertAlmostEqual(analytic["nlpd"], -0.891628740, places=8)
        self.assertAlmostEqual(jordan["nlpd"], -0.376709058, places=8)
        self.assertAlmostEqual(jordan["nlpd"] - analytic["nlpd"], 0.514919682, places=8)


if __name__ == "__main__":
    unittest.main()
