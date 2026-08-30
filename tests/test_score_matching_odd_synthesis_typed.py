from __future__ import annotations

import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import score_matching_odd_synthesis as frozen_kernel  # noqa: E402
from score_matching_odd_synthesis_typed import (  # noqa: E402
    SEMANTIC_MANIFEST,
    load_semantic_gate,
    synthesize_typed,
)


def score(chi_square: float, nlpd: float) -> dict:
    return {"chi_square": chi_square, "dimension": 2, "nlpd": nlpd}


def fixture_ledger() -> dict:
    return {
        "schema_version": "v1",
        "manifest_version": "v1",
        "blocks": [
            {
                "id": "issue43_n185_n265_deltaM",
                "role": "primary",
                "status": "SCORED",
                "channel": {"source": "matching_odd", "target": "matching_odd"},
                "raw_data_group": "issue43_n185_n265_500m_histograms",
                "scores": {
                    "zero_effect": score(2.0, 10.0),
                    "H4_x21_over_4": score(1.0, 4.0),
                },
            },
            {
                "id": "issue57_norm5",
                "role": "primary",
                "status": "SCORED",
                "channel": {"source": "matching_odd", "target": "matching_odd"},
                "raw_data_group": "issue57_norm5_production",
                "scores": {
                    "zero_effect": score(3.0, 5.0),
                    "H4_norm5": score(2.0, 2.0),
                },
            },
        ],
    }


class TypedMatchingOddSynthesisTests(unittest.TestCase):
    def test_both_registered_maps_are_exact_identity(self) -> None:
        _, blocks = load_semantic_gate(ROOT)
        self.assertEqual([block["id"] for block in blocks], [
            "issue43_n185_n265_deltaM", "issue57_norm5"
        ])
        for block in blocks:
            self.assertEqual(
                (block["applied_transform"]["scale"], block["applied_transform"]["offset"]),
                (1.0, 0.0),
            )

    def test_typed_synthesis_preserves_frozen_scores(self) -> None:
        ledger = fixture_ledger()
        expected = frozen_kernel.synthesize(ledger, source_sha256="a" * 64)
        result = synthesize_typed(ROOT, ledger, source_sha256="a" * 64)
        semantics = result.pop("observable_semantics")
        self.assertEqual(result, expected)
        self.assertEqual(
            semantics["validation_order"],
            "all_registered_maps_before_frozen_synthesis",
        )
        self.assertEqual(len(semantics["blocks"]), 2)
        self.assertEqual(result["predictive_comparison"]["delta_nlpd_fixed_H4_minus_zero_effect"], -9.0)

    def test_descriptor_drift_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            target = root / SEMANTIC_MANIFEST
            target.parent.mkdir(parents=True)
            payload = json.loads((ROOT / SEMANTIC_MANIFEST).read_text(encoding="utf-8"))
            payload["blocks"][0]["target_descriptor"]["orientation_order"] = "second_minus_first"
            target.write_text(json.dumps(payload), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "differs"):
                load_semantic_gate(root)

    def test_legacy_channel_drift_still_fails(self) -> None:
        ledger = copy.deepcopy(fixture_ledger())
        ledger["blocks"][0]["channel"]["target"] = "matching_even"
        with self.assertRaisesRegex(ValueError, "matching_odd"):
            synthesize_typed(ROOT, ledger)

    def test_semantic_block_id_drift_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            target = root / SEMANTIC_MANIFEST
            target.parent.mkdir(parents=True)
            payload = json.loads((ROOT / SEMANTIC_MANIFEST).read_text(encoding="utf-8"))
            payload["blocks"][0]["id"] = "other"
            target.write_text(json.dumps(payload), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "block IDs"):
                load_semantic_gate(root)


if __name__ == "__main__":
    unittest.main()
