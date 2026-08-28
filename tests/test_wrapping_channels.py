from __future__ import annotations

import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from wrapping_channels import (  # noqa: E402
    Combination,
    Normalization,
    ObservableDescriptor,
    ObservableMappingError,
    OrientationOrder,
    ProbabilityCoordinate,
    Quantity,
    TopologyChannel,
    map_observable,
)


def descriptor(
    channel: TopologyChannel,
    combination: Combination,
    *,
    coordinate: ProbabilityCoordinate = ProbabilityCoordinate.P,
    order: OrientationOrder = OrientationOrder.FIRST_MINUS_SECOND,
    normalization: Normalization = Normalization.RAW,
    quantity: Quantity = Quantity.ORIENTATION_CONTRAST,
) -> ObservableDescriptor:
    return ObservableDescriptor(
        channel=channel,
        combination=combination,
        coordinate=coordinate,
        orientation_order=order,
        normalization=normalization,
        quantity=quantity,
    )


class WrappingChannelMapTests(unittest.TestCase):
    def test_issue43_even_contrast_map_negates_without_offset(self) -> None:
        transform = map_observable(
            descriptor(TopologyChannel.EITHER, Combination.EVEN),
            descriptor(TopologyChannel.CROSS, Combination.EVEN),
        )
        self.assertEqual((transform.scale, transform.offset), (-1.0, 0.0))
        self.assertAlmostEqual(
            transform.apply(0.010603216462677735), -0.010603216462677735
        )

    def test_even_values_are_affine_but_odd_channels_are_identical(self) -> None:
        even_value = map_observable(
            descriptor(
                TopologyChannel.EITHER,
                Combination.EVEN,
                order=OrientationOrder.NONE,
                quantity=Quantity.VALUE,
            ),
            descriptor(
                TopologyChannel.CROSS,
                Combination.EVEN,
                order=OrientationOrder.NONE,
                quantity=Quantity.VALUE,
            ),
        )
        self.assertAlmostEqual(even_value.apply(0.37), 0.63)

        odd = map_observable(
            descriptor(TopologyChannel.EITHER, Combination.ODD),
            descriptor(TopologyChannel.CROSS, Combination.ODD),
        )
        self.assertEqual((odd.scale, odd.offset), (1.0, 0.0))

    def test_primal_matching_complement_identity(self) -> None:
        transform = map_observable(
            descriptor(
                TopologyChannel.EITHER,
                Combination.PRIMAL,
                coordinate=ProbabilityCoordinate.P,
                order=OrientationOrder.NONE,
                quantity=Quantity.VALUE,
            ),
            descriptor(
                TopologyChannel.CROSS,
                Combination.MATCHING,
                coordinate=ProbabilityCoordinate.COMPLEMENT,
                order=OrientationOrder.NONE,
                quantity=Quantity.VALUE,
            ),
        )
        self.assertAlmostEqual(transform.apply(0.81), 0.19)

    def test_raw_normalized_boundary_requires_angular_factor(self) -> None:
        source = descriptor(TopologyChannel.CROSS, Combination.ODD)
        target = descriptor(
            TopologyChannel.CROSS,
            Combination.ODD,
            normalization=Normalization.ANGULAR_NORMALIZED,
        )
        with self.assertRaises(ObservableMappingError):
            map_observable(source, target)
        transform = map_observable(source, target, source_angular_factor=-0.5)
        self.assertEqual(transform.scale, -2.0)

    def test_unsupported_channel_change_fails_closed(self) -> None:
        with self.assertRaises(ObservableMappingError):
            map_observable(
                descriptor(TopologyChannel.BOTH, Combination.EVEN),
                descriptor(TopologyChannel.EITHER, Combination.EVEN),
            )

    def test_machine_readable_audit_maps_validate(self) -> None:
        path = ROOT / "predictions" / "wrapping_channel_audit_20260828.json"
        payload = json.loads(path.read_text(encoding="utf-8"))
        registered = 0
        for row in payload["records"]:
            if row["status"] != "registered":
                continue
            registered += 1
            if "source_descriptor" in row and "target_descriptor" in row:
                transform = map_observable(
                    ObservableDescriptor.from_dict(row["source_descriptor"]),
                    ObservableDescriptor.from_dict(row["target_descriptor"]),
                )
                self.assertEqual(transform.scale, row["expected_transform"]["scale"])
                self.assertEqual(transform.offset, row["expected_transform"]["offset"])
                continue

            primitives = row.get("primitive_descriptors")
            if not isinstance(primitives, dict) or not primitives:
                self.fail(
                    f"registered audit row {row.get('id')} lacks a supported descriptor schema"
                )
            expected = row["expected_transform"]
            for name, primitive in primitives.items():
                with self.subTest(audit=row.get("id"), primitive=name):
                    descriptor_value = ObservableDescriptor.from_dict(primitive)
                    transform = map_observable(descriptor_value, descriptor_value)
                    self.assertEqual(transform.scale, expected["scale"])
                    self.assertEqual(transform.offset, expected["offset"])
        self.assertGreaterEqual(registered, 4)


if __name__ == "__main__":
    unittest.main()
