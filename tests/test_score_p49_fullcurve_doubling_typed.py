from __future__ import annotations

import copy
import json
import math
from pathlib import Path
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import score_p49_fullcurve_doubling_typed as typed  # noqa: E402


def assert_nested_close(test: unittest.TestCase, actual, expected, path="root") -> None:
    test.assertIs(type(actual), type(expected), path)
    if isinstance(actual, dict):
        test.assertEqual(set(actual), set(expected), path)
        for key in actual:
            assert_nested_close(test, actual[key], expected[key], f"{path}/{key}")
    elif isinstance(actual, list):
        test.assertEqual(len(actual), len(expected), path)
        for index, (left, right) in enumerate(zip(actual, expected)):
            assert_nested_close(test, left, right, f"{path}/{index}")
    elif isinstance(actual, float):
        test.assertTrue(
            math.isclose(actual, expected, rel_tol=2e-12, abs_tol=2e-15),
            f"{path}: {actual!r} != {expected!r}",
        )
    else:
        test.assertEqual(actual, expected, path)


class P49FullcurveDoublingTypedTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.gate = json.loads((ROOT / typed.SEMANTIC_GATE).read_text(encoding="utf-8"))
        cls.histograms = [Path(row["path"]) for row in cls.gate["canonical_inputs"]]
        cls.typed_result = typed.score_typed(ROOT, cls.histograms)
        cls.frozen_result = copy.deepcopy(cls.typed_result)
        cls.semantics = cls.frozen_result.pop("observable_semantics")

    def test_registered_sign_and_normalized_projector_maps(self) -> None:
        _, validated = typed.load_semantic_gate(ROOT)
        matching = validated["matching_map"]
        self.assertEqual((matching.scale, matching.offset), (-1.0, 0.0))
        self.assertTrue(all(
            (row["transform"].scale, row["transform"].offset) == (1.0, 0.0)
            for row in validated["projectors"].values()
        ))

    def test_canonical_replay_matches_committed_score(self) -> None:
        committed = json.loads(
            (ROOT / self.gate["committed_score"]["path"]).read_text(encoding="utf-8")
        )
        replay = json.loads(json.dumps(self.frozen_result))
        assert_nested_close(self, replay, committed)
        self.assertEqual(
            self.semantics["validation_order"],
            "semantic_maps_and_canonical_inputs_before_frozen_fullcurve_score",
        )

    def test_noncanonical_input_fails_before_runner(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            altered = Path(directory) / "n65.hist.csv"
            altered.write_bytes(self.histograms[0].read_bytes() + b"\n")
            inputs = [altered, *self.histograms[1:]]
            called = False

            def runner(_paths) -> dict:
                nonlocal called
                called = True
                return {}

            with self.assertRaisesRegex(ValueError, "canonical histogram 0 git blob"):
                typed.score_typed(ROOT, inputs, runner=runner)
            self.assertFalse(called)

    def test_joint_order_drift_fails_closed(self) -> None:
        bad = copy.deepcopy(self.frozen_result)
        name, value = bad["joint_scores"].popitem()
        bad["joint_scores"] = {name: value, **bad["joint_scores"]}
        with self.assertRaisesRegex(ValueError, "joint score order"):
            typed.score_typed(ROOT, self.histograms, runner=lambda _: bad)

    def test_refit_boundary_drift_fails_closed(self) -> None:
        bad = copy.deepcopy(self.frozen_result)
        bad["P48_Sprime_fresh_seed_replication"]["classification"] = "refit"
        with self.assertRaisesRegex(ValueError, "refit boundary"):
            typed.score_typed(ROOT, self.histograms, runner=lambda _: bad)

    def test_response_coordinates_are_not_promoted_to_topology(self) -> None:
        self.assertIn("response/model coordinates", self.gate["semantic_boundary"])
        self.assertIn("not independent evidence rows", self.gate["evidence_boundary"])


if __name__ == "__main__":
    unittest.main()
