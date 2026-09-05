
from __future__ import annotations
import copy
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import score_p48_new_geometry_channels_typed as typed  # noqa: E402


def inputs() -> tuple[dict, dict[int, dict]]:
    gate, _ = typed.load_semantic_gate(ROOT)
    canonical = gate["canonical_inputs"]
    return typed.load_canonical_inputs(
        gate,
        ROOT / canonical["source"]["path"],
        ROOT / canonical["185"]["path"],
        ROOT / canonical["265"]["path"],
    )


class TypedP48NewGeometryChannelsTests(unittest.TestCase):
    def test_four_projector_maps_and_powers_are_frozen(self) -> None:
        gate, validated = typed.load_semantic_gate(ROOT)
        self.assertEqual(gate["projector_order"], list(typed.frozen.CHANNELS))
        for projector, values in validated.items():
            self.assertEqual(
                (values["transform"].scale, values["transform"].offset),
                (1.0, 0.0),
                projector,
            )
        self.assertEqual(
            validated["P4_S"]["source_descriptor"].combination.value, "even"
        )
        self.assertEqual(
            validated["P4_D"]["source_descriptor"].combination.value, "odd"
        )
        self.assertEqual(
            validated["P4_S_prime"]["normalization_power_in_N"],
            {"numerator": 5, "denominator": 4},
        )

    def test_typed_score_preserves_every_frozen_number(self) -> None:
        source, targets = inputs()
        expected = typed.frozen.score(copy.deepcopy(source), copy.deepcopy(targets))
        result = typed.score_typed(ROOT, source, targets)
        semantics = result.pop("observable_semantics")
        self.assertEqual(result, expected)
        self.assertEqual(
            semantics["validation_order"],
            "semantic_maps_before_frozen_four_projector_score",
        )
        self.assertEqual(
            semantics["projectors"]["P4_D_prime"]["response_coordinate"],
            "first_p_derivative",
        )


if __name__ == "__main__":
    unittest.main()
