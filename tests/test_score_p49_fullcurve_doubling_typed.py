
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

    def test_response_coordinates_are_not_promoted_to_topology(self) -> None:
        self.assertIn("response/model coordinates", self.gate["semantic_boundary"])
        self.assertIn("not independent evidence rows", self.gate["evidence_boundary"])


if __name__ == "__main__":
    unittest.main()
