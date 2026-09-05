
from __future__ import annotations
import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import score_axis_pair_annihilator_stable_typed as stable_typed  # noqa: E402
import score_axis_pair_annihilator_typed as typed  # noqa: E402


class AxisPairAnnihilatorTypedTests(unittest.TestCase):
    def test_gate_registers_exact_matching_odd_identity(self) -> None:
        gate, source, target, transform = typed.load_semantic_gate(ROOT)
        self.assertEqual(source.to_dict(), target.to_dict())
        self.assertEqual(source.channel.value, "cross")
        self.assertEqual(source.combination.value, "odd")
        self.assertEqual((transform.scale, transform.offset), (1.0, 0.0))
        self.assertEqual(gate["candidate_q_default_order"], [2.0, 3.0, 4.0, 6.0])

    def test_typed_score_preserves_frozen_numbers(self) -> None:
        frozen = {
            "format_version": 1, "p_ref": 0.592746050790,
            "train_max_L": 14, "pairs": [{"F_p_ref": 1.25}],
            "F_shape_models": [{"heldout_chi_square": 2.5}],
            "accelerated_root_models": [{"heldout_chi_square": 3.5}],
            "primary_candidate_order": [2.0, 3.0, 4.0, 6.0],
        }
        calculator = mock.Mock(return_value=copy.deepcopy(frozen))
        result = typed.calculate_typed(
            ROOT, [Path("hist.csv")], 0.592746050790, 14,
            [2.0, 3.0, 4.0, 6.0], calculator=calculator,
        )
        semantics = result.pop("observable_semantics")
        self.assertEqual(result, frozen)
        self.assertEqual(semantics["frozen_kernel"], "base")
        calculator.assert_called_once()

    def test_invalid_map_fails_before_frozen_kernel(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            destination = root / typed.SEMANTIC_GATE
            destination.parent.mkdir(parents=True)
            gate = json.loads((ROOT / typed.SEMANTIC_GATE).read_text(encoding="utf-8"))
            gate["target_descriptor"]["channel"] = "direction_0"
            destination.write_text(json.dumps(gate), encoding="utf-8")
            calculator = mock.Mock()
            with self.assertRaisesRegex(ValueError, "no exact topology map"):
                typed.calculate_typed(root, [], 0.592746050790, 14, [2.0], calculator=calculator)
            calculator.assert_not_called()

    def test_stable_entrypoint_selects_stable_reader_kernel(self) -> None:
        frozen = {"p_ref": 0.592746050790, "train_max_L": 14,
                  "primary_candidate_order": [3.0]}
        with mock.patch.object(stable_typed.stable, "calculate", return_value=frozen.copy()) as kernel:
            result = stable_typed.calculate_typed(
                ROOT, [Path("hist.csv")], 0.592746050790, 14, [3.0]
            )
        self.assertEqual(result["observable_semantics"]["frozen_kernel"], "stable")
        kernel.assert_called_once()


if __name__ == "__main__":
    unittest.main()
