from __future__ import annotations

import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock

import mpmath as mp


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import score_angular_root_amplitude_typed as typed  # noqa: E402


def fixture_records() -> dict:
    records = {}
    for size, first, second in ((65, (8, 1), (7, 4)), (85, (9, 2), (7, 6))):
        records[(size, "first", 0)] = {"a": first[0], "b": first[1]}
        records[(size, "second", 0)] = {"a": second[0], "b": second[1]}
    return records


def fixture_result() -> dict:
    return {
        "schema": "frozen angular-normalized root amplitude score v1",
        "p_ref": "0.59274605079",
        "sizes": [65, 85],
        "batch_count": 100,
        "by_size": {"65": {"estimate": {"A_p": 0.42}}, "85": {"estimate": {"A_p": 0.39}}},
        "A_p_cross_size_jackknife_covariance": [[1.0, 0.2], [0.2, 2.0]],
        "frozen_prediction": {
            "value": 0.4510066187069702,
            "source_standard_error": 0.02013371335254959,
            "chi_square": 2.426668731892274,
        },
    }


class AngularRootAmplitudeTypedTests(unittest.TestCase):
    def test_gate_applies_size_specific_raw_to_angular_maps(self) -> None:
        gate, source, target, transforms = typed.load_semantic_gate(ROOT)
        self.assertEqual(source.normalization.value, "raw")
        self.assertEqual(target.normalization.value, "angular_normalized")
        for size in gate["primary_sizes"]:
            delta = gate["designs"][str(size)]["signed_delta_cos4"]
            self.assertAlmostEqual(transforms[size].scale * delta, 1.0, places=15)

    def test_typed_score_preserves_every_frozen_number(self) -> None:
        frozen = fixture_result()
        scorer = mock.Mock(return_value=copy.deepcopy(frozen))
        result = typed.score_typed(
            ROOT,
            fixture_records(),
            mp.mpf("0.592746050790"),
            0.4510066187069702,
            0.02013371335254959,
            scorer=scorer,
        )
        semantics = result.pop("observable_semantics")
        self.assertEqual(result, frozen)
        self.assertEqual(semantics["response_coordinate"], "implicit matching-root location")
        scorer.assert_called_once()

    def test_bad_angular_factor_fails_before_frozen_score(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            destination = root / typed.SEMANTIC_GATE
            destination.parent.mkdir(parents=True)
            gate = json.loads((ROOT / typed.SEMANTIC_GATE).read_text(encoding="utf-8"))
            gate["designs"]["65"]["signed_delta_cos4"] = 2.0
            destination.write_text(json.dumps(gate), encoding="utf-8")
            scorer = mock.Mock()
            with self.assertRaisesRegex(ValueError, "registered map"):
                typed.score_typed(
                    root, fixture_records(), mp.mpf("0.592746050790"),
                    0.4510066187069702, 0.02013371335254959, scorer=scorer,
                )
            scorer.assert_not_called()

    def test_orientation_geometry_is_signed_and_frozen(self) -> None:
        records = fixture_records()
        records[(65, "first", 0)]["a"] = 7
        scorer = mock.Mock()
        with self.assertRaisesRegex(ValueError, "orientation geometry"):
            typed.score_typed(
                ROOT, records, mp.mpf("0.592746050790"),
                0.4510066187069702, 0.02013371335254959, scorer=scorer,
            )
        scorer.assert_not_called()


if __name__ == "__main__":
    unittest.main()
