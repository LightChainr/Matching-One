#!/usr/bin/env python3
"""Type-safe semantics for torus wrapping observables.

The module deliberately registers only exact topology identities.  Unsupported
channel, coordinate, combination, or normalization changes raise
``ObservableMappingError`` instead of silently treating similarly named
statistics as interchangeable.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Any, Dict, Mapping, Optional, Tuple


class ObservableMappingError(ValueError):
    """Raised when no exact registered observable map exists."""


class TopologyChannel(str, Enum):
    CROSS = "cross"
    EITHER = "either"
    BOTH = "both"
    DIRECTION_0 = "direction_0"
    DIRECTION_1 = "direction_1"


class ModularBehavior(str, Enum):
    """Behavior of a topology label under an SL(2,Z) homology-basis change."""

    SCALAR = "modular_scalar"
    BASIS_DEPENDENT = "basis_dependent"


class Combination(str, Enum):
    PRIMAL = "primal"
    MATCHING = "matching"
    EVEN = "even"
    ODD = "odd"


class ProbabilityCoordinate(str, Enum):
    P = "p"
    COMPLEMENT = "1-p"


class OrientationOrder(str, Enum):
    NONE = "none"
    FIRST_MINUS_SECOND = "first_minus_second"
    SECOND_MINUS_FIRST = "second_minus_first"


class Normalization(str, Enum):
    RAW = "raw"
    ANGULAR_NORMALIZED = "angular_normalized"


class Quantity(str, Enum):
    VALUE = "value"
    ORIENTATION_CONTRAST = "orientation_contrast"


@dataclass(frozen=True)
class ObservableDescriptor:
    channel: TopologyChannel
    combination: Combination
    coordinate: ProbabilityCoordinate
    orientation_order: OrientationOrder
    normalization: Normalization
    quantity: Quantity

    def __post_init__(self) -> None:
        if self.quantity is Quantity.VALUE:
            if self.orientation_order is not OrientationOrder.NONE:
                raise ValueError("a scalar value must use orientation_order=none")
            if self.normalization is not Normalization.RAW:
                raise ValueError("angular normalization requires an orientation contrast")
        elif self.orientation_order is OrientationOrder.NONE:
            raise ValueError("an orientation contrast must declare its signed order")

    def to_dict(self) -> Dict[str, str]:
        return {
            "channel": self.channel.value,
            "combination": self.combination.value,
            "coordinate": self.coordinate.value,
            "orientation_order": self.orientation_order.value,
            "normalization": self.normalization.value,
            "quantity": self.quantity.value,
        }

    @property
    def modular_behavior(self) -> ModularBehavior:
        """Classify only the topology label, not a physical response field."""

        return modular_behavior(self.channel)

    @classmethod
    def from_dict(cls, payload: Mapping[str, Any]) -> "ObservableDescriptor":
        required = {
            "channel",
            "combination",
            "coordinate",
            "orientation_order",
            "normalization",
            "quantity",
        }
        missing = sorted(required - set(payload))
        if missing:
            raise ValueError("observable descriptor lacks: " + ", ".join(missing))
        return cls(
            channel=TopologyChannel(str(payload["channel"])),
            combination=Combination(str(payload["combination"])),
            coordinate=ProbabilityCoordinate(str(payload["coordinate"])),
            orientation_order=OrientationOrder(str(payload["orientation_order"])),
            normalization=Normalization(str(payload["normalization"])),
            quantity=Quantity(str(payload["quantity"])),
        )


@dataclass(frozen=True)
class AffineTransform:
    scale: float
    offset: float
    exact_identity: str

    def apply(self, value: float) -> float:
        return self.scale * float(value) + self.offset

    def apply_standard_error(self, standard_error: float) -> float:
        return abs(self.scale) * float(standard_error)

    def to_dict(self) -> Dict[str, object]:
        return {
            "scale": self.scale,
            "offset": self.offset,
            "exact_identity": self.exact_identity,
        }


def modular_behavior(channel: TopologyChannel) -> ModularBehavior:
    """Return the exact homology-basis behavior of a registered channel.

    ``either`` is the event ``rank > 0`` and ``cross`` is ``rank == 2``;
    rational subgroup rank is invariant under SL(2,Z).  The other registered
    channels refer to the selected quotient generators.  In particular, a
    rank-one spiral can have ``both=True`` and shear to a single generator.

    This is a topology-label contract only.  It does not assert that an
    orientation response is a homogeneous CFT spin field.
    """

    if channel in (TopologyChannel.CROSS, TopologyChannel.EITHER):
        return ModularBehavior.SCALAR
    if channel in (
        TopologyChannel.BOTH,
        TopologyChannel.DIRECTION_0,
        TopologyChannel.DIRECTION_1,
    ):
        return ModularBehavior.BASIS_DEPENDENT
    raise ObservableMappingError(f"unregistered topology channel {channel!r}")


def require_modular_scalar_topology(
    descriptor: ObservableDescriptor,
) -> ModularBehavior:
    """Fail closed unless ``descriptor`` uses a modular-scalar topology label."""

    behavior = descriptor.modular_behavior
    if behavior is not ModularBehavior.SCALAR:
        raise ObservableMappingError(
            f"channel {descriptor.channel.value!r} is {behavior.value}, not modular_scalar"
        )
    return behavior


def _opposite_coordinates(
    source: ProbabilityCoordinate, target: ProbabilityCoordinate
) -> bool:
    return {source, target} == {
        ProbabilityCoordinate.P,
        ProbabilityCoordinate.COMPLEMENT,
    }


def _topology_transform(
    source: ObservableDescriptor, target: ObservableDescriptor
) -> Tuple[float, float, str]:
    if (
        source.channel is target.channel
        and source.combination is target.combination
        and source.coordinate is target.coordinate
    ):
        return 1.0, 0.0, "identical channel, combination, and coordinate"

    cross_either = {source.channel, target.channel} == {
        TopologyChannel.CROSS,
        TopologyChannel.EITHER,
    }

    if (
        cross_either
        and source.combination is target.combination
        and source.coordinate is target.coordinate
    ):
        if source.combination is Combination.EVEN:
            offset = 1.0 if source.quantity is Quantity.VALUE else 0.0
            return -1.0, offset, "S_either = 1 - S_cross"
        if source.combination is Combination.ODD:
            return 1.0, 0.0, "D_either = D_cross"

    primal_matching = {source.combination, target.combination} == {
        Combination.PRIMAL,
        Combination.MATCHING,
    }
    if cross_either and primal_matching and _opposite_coordinates(
        source.coordinate, target.coordinate
    ):
        offset = 1.0 if source.quantity is Quantity.VALUE else 0.0
        return -1.0, offset, "R_primal,either(p) = 1 - R_matching,cross(1-p)"

    raise ObservableMappingError(
        "no exact topology map from "
        f"{source.to_dict()} to {target.to_dict()}"
    )


def map_observable(
    source: ObservableDescriptor,
    target: ObservableDescriptor,
    *,
    source_angular_factor: Optional[float] = None,
    target_angular_factor: Optional[float] = None,
) -> AffineTransform:
    """Return the exact affine map from ``source`` to ``target``.

    Signed angular factors are required only when crossing the raw/normalized
    boundary.  They must use the orientation order declared by the associated
    raw descriptor.
    """

    if source.quantity is not target.quantity:
        raise ObservableMappingError("cannot map a scalar value to an orientation contrast")

    scale, offset, identity = _topology_transform(source, target)

    if source.quantity is Quantity.VALUE:
        return AffineTransform(scale, offset, identity)

    if source.normalization is Normalization.RAW and target.normalization is Normalization.RAW:
        if source.orientation_order is not target.orientation_order:
            scale = -scale
            identity += "; reversed orientation order"
    elif (
        source.normalization is Normalization.ANGULAR_NORMALIZED
        and target.normalization is Normalization.ANGULAR_NORMALIZED
    ):
        identity += "; normalized contrast is invariant under order reversal"
    elif (
        source.normalization is Normalization.RAW
        and target.normalization is Normalization.ANGULAR_NORMALIZED
    ):
        if source_angular_factor is None or source_angular_factor == 0:
            raise ObservableMappingError(
                "raw-to-normalized mapping requires a nonzero source_angular_factor"
            )
        scale /= float(source_angular_factor)
        identity += "; divide by the source signed angular factor"
    else:
        if target_angular_factor is None or target_angular_factor == 0:
            raise ObservableMappingError(
                "normalized-to-raw mapping requires a nonzero target_angular_factor"
            )
        scale *= float(target_angular_factor)
        identity += "; multiply by the target signed angular factor"

    if offset != 0.0:
        raise AssertionError("orientation contrasts cannot retain an affine offset")
    return AffineTransform(scale, 0.0, identity)
