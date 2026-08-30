#!/usr/bin/env python3
"""Exhaustive finite-state audit of the canonical observable descriptor map."""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from functools import lru_cache
from pathlib import Path
from typing import Any

from wrapping_channels import (
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


ANGULAR_FACTORS = {
    OrientationOrder.FIRST_MINUS_SECOND: 2.0,
    OrientationOrder.SECOND_MINUS_FIRST: -2.0,
}


def enumerate_descriptors() -> tuple[ObservableDescriptor, ...]:
    result = []
    for channel in TopologyChannel:
        for combination in Combination:
            for coordinate in ProbabilityCoordinate:
                result.append(ObservableDescriptor(
                    channel=channel,
                    combination=combination,
                    coordinate=coordinate,
                    orientation_order=OrientationOrder.NONE,
                    normalization=Normalization.RAW,
                    quantity=Quantity.VALUE,
                ))
                for order in (
                    OrientationOrder.FIRST_MINUS_SECOND,
                    OrientationOrder.SECOND_MINUS_FIRST,
                ):
                    for normalization in Normalization:
                        result.append(ObservableDescriptor(
                            channel=channel,
                            combination=combination,
                            coordinate=coordinate,
                            orientation_order=order,
                            normalization=normalization,
                            quantity=Quantity.ORIENTATION_CONTRAST,
                        ))
    return tuple(result)


def descriptor_id(value: ObservableDescriptor) -> str:
    return "/".join(value.to_dict()[key] for key in (
        "channel", "combination", "coordinate", "quantity", "orientation_order", "normalization"
    ))


def angular_factor(value: ObservableDescriptor) -> float | None:
    return ANGULAR_FACTORS.get(value.orientation_order)


def registered_map(source: ObservableDescriptor, target: ObservableDescriptor):
    return map_observable(
        source,
        target,
        source_angular_factor=angular_factor(source),
        target_angular_factor=angular_factor(target),
    )


def connected_components(
    descriptors: tuple[ObservableDescriptor, ...],
    maps: dict[tuple[int, int], Any],
) -> list[list[int]]:
    unseen = set(range(len(descriptors)))
    components = []
    while unseen:
        root = min(unseen)
        unseen.remove(root)
        stack = [root]
        component = []
        while stack:
            current = stack.pop()
            component.append(current)
            neighbors = [
                candidate for candidate in sorted(unseen)
                if maps.get((current, candidate)) is not None
                or maps.get((candidate, current)) is not None
            ]
            for candidate in neighbors:
                unseen.remove(candidate)
                stack.append(candidate)
        components.append(sorted(component))
    return components


@lru_cache(maxsize=1)
def build_artifact() -> dict[str, Any]:
    descriptors = enumerate_descriptors()
    if len({descriptor_id(value) for value in descriptors}) != len(descriptors):
        raise AssertionError("descriptor enumeration is not unique")

    maps: dict[tuple[int, int], Any] = {}
    blocked_reasons: Counter[str] = Counter()
    canonical_edges = []
    quantity_successes: Counter[str] = Counter()
    for i, source in enumerate(descriptors):
        for j, target in enumerate(descriptors):
            try:
                transform = registered_map(source, target)
            except ObservableMappingError as exc:
                maps[(i, j)] = None
                reason = str(exc)
                if reason.startswith("no exact topology map"):
                    reason = "no exact topology map"
                blocked_reasons[reason] += 1
                continue
            maps[(i, j)] = transform
            quantity_successes[source.quantity.value] += 1
            canonical_edges.append(json.dumps({
                "source": descriptor_id(source),
                "target": descriptor_id(target),
                "scale": transform.scale,
                "offset": transform.offset,
                "identity": transform.exact_identity,
            }, sort_keys=True, separators=(",", ":")))

    inverse_failures = []
    for (i, j), transform in maps.items():
        if transform is None:
            continue
        reverse = maps.get((j, i))
        if reverse is None or (
            reverse.scale * transform.scale,
            reverse.scale * transform.offset + reverse.offset,
        ) != (1.0, 0.0):
            inverse_failures.append([descriptor_id(descriptors[i]), descriptor_id(descriptors[j])])

    composition_paths = 0
    composition_failures = []
    for (i, j), first in maps.items():
        if first is None:
            continue
        for k in range(len(descriptors)):
            second = maps.get((j, k))
            if second is None:
                continue
            composition_paths += 1
            direct = maps.get((i, k))
            composed = (
                second.scale * first.scale,
                second.scale * first.offset + second.offset,
            )
            if direct is None or (direct.scale, direct.offset) != composed:
                composition_failures.append([
                    descriptor_id(descriptors[i]),
                    descriptor_id(descriptors[j]),
                    descriptor_id(descriptors[k]),
                ])

    components = connected_components(descriptors, maps)
    component_records = [
        {
            "id": index,
            "size": len(component),
            "members": [descriptors[item].to_dict() for item in component],
        }
        for index, component in enumerate(components)
    ]
    size_histogram = Counter(len(component) for component in components)
    successful = len(canonical_edges)
    ordered_pairs = len(descriptors) ** 2
    artifact = {
        "schema": "matching-one/observable-descriptor-map-audit/v1",
        "issue": 146,
        "data_class": "exhaustive finite registry audit",
        "angular_factor_fixture": {
            "first_minus_second": 2,
            "second_minus_first": -2,
            "purpose": "exercise raw/normalized maps; not a physical geometry claim",
        },
        "valid_descriptors": len(descriptors),
        "ordered_descriptor_pairs": ordered_pairs,
        "registered_maps": successful,
        "blocked_pairs": ordered_pairs - successful,
        "registered_by_source_quantity": dict(sorted(quantity_successes.items())),
        "blocked_reasons": dict(sorted(blocked_reasons.items())),
        "registered_edge_sha256": hashlib.sha256(("\n".join(canonical_edges) + "\n").encode()).hexdigest(),
        "inverse_checks": {"checked": successful, "failures": inverse_failures},
        "composition_checks": {"checked_paths": composition_paths, "failures": composition_failures},
        "connected_components": component_records,
        "component_size_histogram": {str(key): value for key, value in sorted(size_histogram.items())},
        "boundary": (
            "This exhausts the current finite descriptor registry only. It does not scan all scorers or "
            "artifacts, register a new topology identity, or prove repository-wide descriptor adoption."
        ),
    }
    assert successful == 952
    assert ordered_pairs - successful == 39048
    assert not inverse_failures and not composition_failures
    assert composition_paths == 5720
    assert size_histogram == Counter({1: 24, 2: 8, 4: 24, 8: 8})
    return artifact


def render_markdown(artifact: dict[str, Any]) -> str:
    return "\n".join([
        "# Observable descriptor map audit", "",
        "This is a complete audit of the finite map registry, not a scan of every repository artifact.", "",
        "| item | exact count |",
        "|---|---:|",
        f"| valid descriptors | {artifact['valid_descriptors']} |",
        f"| ordered descriptor pairs | {artifact['ordered_descriptor_pairs']} |",
        f"| registered exact maps | {artifact['registered_maps']} |",
        f"| fail-closed pairs | {artifact['blocked_pairs']} |",
        f"| inverse checks | {artifact['inverse_checks']['checked']} |",
        f"| composable paths | {artifact['composition_checks']['checked_paths']} |",
        f"| inverse/composition failures | {len(artifact['inverse_checks']['failures']) + len(artifact['composition_checks']['failures'])} |", "",
        f"Connected-component size histogram: `{artifact['component_size_histogram']}`.", "",
        f"Blocked reasons: `{artifact['blocked_reasons']}`.", "",
        f"Registered-edge SHA-256: `{artifact['registered_edge_sha256']}`.", "",
        "## Interpretation boundary", "", artifact["boundary"], "",
    ])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--format", choices=("json", "markdown"), default="json")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    artifact = build_artifact()
    rendered = json.dumps(artifact, indent=2, sort_keys=True) + "\n" if args.format == "json" else render_markdown(artifact)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
