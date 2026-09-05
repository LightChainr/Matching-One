
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


def assert_nested_close(
    test: unittest.TestCase, actual, expected,
    rel_tol: float, abs_tol: float, path: str = "root",
) -> None:
    """Compare two decoded JSON trees: structure exactly, floats to a tolerance.

    The tolerances are arguments rather than defaults because the only caller
    reads them from the semantic gate, where they are recorded next to the
    measurement that justifies them.
    """
    test.assertIs(type(actual), type(expected), path)
    if isinstance(actual, dict):
        test.assertEqual(set(actual), set(expected), path)
        for key in actual:
            assert_nested_close(
                test, actual[key], expected[key], rel_tol, abs_tol, f"{path}/{key}"
            )
    elif isinstance(actual, list):
        test.assertEqual(len(actual), len(expected), path)
        for index, (left, right) in enumerate(zip(actual, expected)):
            assert_nested_close(
                test, left, right, rel_tol, abs_tol, f"{path}/{index}"
            )
    elif isinstance(actual, float):
        test.assertTrue(
            math.isclose(actual, expected, rel_tol=rel_tol, abs_tol=abs_tol),
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
        """Stops us believing a score.json the frozen kernel no longer produces.

        The wrong number this test would catch is a changed projector, a
        flipped orientation, a dropped batch, or a changed jackknife
        normalization -- any of which move these numbers by a percent or more.

        It is deliberately not a bitwise test. The committed score is a dated
        reveal, pinned by blob in analysis/evidence_ledger_manifest.yaml against
        a prediction frozen before it, so it is not regenerated when a shared
        numerical routine improves. Commit 844173c re-anchored the binomial tail
        recurrence at the distribution mode; that reassociation moved 565 of the
        598 float leaves here, by at most 6.8e-9 relative, and no chi-square by
        more than 2.4e-10. The tolerance and that measurement live together in
        the gate's replay_reproducibility block.
        """
        committed = json.loads(
            (ROOT / self.gate["committed_score"]["path"]).read_text(encoding="utf-8")
        )
        replay = json.loads(json.dumps(self.frozen_result))
        reproducibility = self.gate["replay_reproducibility"]
        assert_nested_close(
            self, replay, committed,
            rel_tol=reproducibility["replay_rel_tol"],
            abs_tol=reproducibility["replay_abs_tol"],
        )
        self.assertEqual(
            self.semantics["validation_order"],
            "semantic_maps_and_canonical_inputs_before_frozen_fullcurve_score",
        )

    def test_the_replay_tolerance_stays_tied_to_its_measurement(self) -> None:
        """Stops the replay tolerance from being widened on its own.

        The wrong number here is a tolerance raised to make a real regression
        pass. Loosening it requires re-measuring the drift and writing the new
        measurement into the gate beside it, where a reader sees both.
        """
        reproducibility = self.gate["replay_reproducibility"]
        self.assertEqual(reproducibility["contract"], "semantic_not_bitwise")
        self.assertIn("844173c", reproducibility["why_not_bitwise"])
        self.assertLess(
            reproducibility["largest_relative_move_outside_covariance"],
            reproducibility["replay_rel_tol"],
        )
        self.assertLess(
            reproducibility["largest_chi_square_relative_move"],
            reproducibility["replay_rel_tol"],
        )
        self.assertEqual(reproducibility["structural_differences"], 0)
        # A tolerance this side of 1e-6 is still orders of magnitude below the
        # smallest drift any of the failures above would produce.
        self.assertLessEqual(reproducibility["replay_rel_tol"], 1e-6)

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
