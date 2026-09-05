
from __future__ import annotations
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from score_c4_self_matching_n26_typed import (  # noqa: E402
    SEMANTIC_MANIFEST,
    load_semantic_gate,
    score_typed,
)


RESULT = ROOT / "results" / "local-20260828" / "P115-n26-self-matching-exact"
PREDICTION = ROOT / "predictions" / "c4_self_matching_n26_beta_targets_20260828.json"


class TypedC4SelfMatchingN26Tests(unittest.TestCase):
    def test_registered_map_is_exact_identity(self) -> None:
        _, source, target, transform = load_semantic_gate(ROOT)
        self.assertEqual(source.to_dict(), target.to_dict())
        self.assertEqual((transform.scale, transform.offset), (1.0, 0.0))

    def test_typed_replay_preserves_frozen_result(self) -> None:
        result = score_typed(
            ROOT,
            PREDICTION,
            RESULT / "raw" / "n26_threads10.json",
            RESULT / "raw" / "n26_threads1.json",
        )
        self.assertEqual(
            result["protocol_conclusion"],
            "BOTH_FROZEN_LAWS_FAILED_STOP_NO_GENERALIZED_BETA_FIT",
        )
        semantics = result["observable_semantics"]
        self.assertEqual(semantics["applied_transform"]["scale"], 1.0)
        self.assertEqual(semantics["validation_order"], "semantic_map_before_frozen_kernel_score")


if __name__ == "__main__":
    unittest.main()
