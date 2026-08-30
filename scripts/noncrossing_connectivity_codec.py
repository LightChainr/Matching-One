#!/usr/bin/env python3
"""Exact small-width codec for circular noncrossing connectivity states."""

from __future__ import annotations

import argparse
from functools import lru_cache
import json
import math
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONTRACT = ROOT / "analysis" / "noncrossing_connectivity_codec_contract.json"
EXPECTED_SCHEMA = "matching-one/noncrossing-connectivity-codec/v1"


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def canonical_rgs(labels: Iterable[int]) -> tuple[int, ...]:
    """Canonicalize arbitrary integer block labels by first occurrence."""

    mapping: dict[int, int] = {}
    output: list[int] = []
    for label in labels:
        _require(isinstance(label, int), "block labels must be integers")
        if label not in mapping:
            mapping[label] = len(mapping)
        output.append(mapping[label])
    return tuple(output)


def validate_rgs(state: Sequence[int]) -> tuple[int, ...]:
    state = tuple(state)
    _require(bool(state), "connectivity state must not be empty")
    _require(all(isinstance(label, int) for label in state), "block labels must be integers")
    _require(state[0] == 0, "restricted-growth string must start at zero")
    maximum = 0
    for index, label in enumerate(state[1:], start=1):
        _require(0 <= label <= maximum + 1, f"invalid restricted-growth label at index {index}")
        maximum = max(maximum, label)
    _require(canonical_rgs(state) == state, "restricted-growth string is not canonical")
    return state


def is_noncrossing_rgs(state: Sequence[int]) -> bool:
    """Return false exactly when two blocks contain an alternating a<b<c<d."""

    state = validate_rgs(state)
    width = len(state)
    for a in range(width):
        for b in range(a + 1, width):
            if state[a] == state[b]:
                continue
            for c in range(b + 1, width):
                if state[c] != state[a]:
                    continue
                for d in range(c + 1, width):
                    if state[d] == state[b]:
                        return False
    return True


def generate_rgs_partitions(width: int) -> tuple[tuple[int, ...], ...]:
    """Generate set partitions directly as restricted-growth strings."""

    _require(width >= 1, "width must be positive")
    output: list[tuple[int, ...]] = []

    def visit(prefix: tuple[int, ...], maximum: int) -> None:
        if len(prefix) == width:
            output.append(prefix)
            return
        for label in range(maximum + 2):
            visit(prefix + (label,), max(maximum, label))

    visit((0,), 0)
    return tuple(output)


def generate_block_partitions(width: int) -> tuple[tuple[frozenset[int], ...], ...]:
    """Independently generate unlabelled set partitions by block insertion."""

    _require(width >= 1, "width must be positive")
    partitions: list[tuple[frozenset[int], ...]] = [(frozenset({0}),)]
    for point in range(1, width):
        following: list[tuple[frozenset[int], ...]] = []
        for partition in partitions:
            following.append(partition + (frozenset({point}),))
            for index in range(len(partition)):
                blocks = list(partition)
                blocks[index] = blocks[index] | {point}
                following.append(tuple(blocks))
        partitions = following
    return tuple(partitions)


def blocks_to_rgs(blocks: Sequence[frozenset[int]], width: int) -> tuple[int, ...]:
    _require(width >= 1, "width must be positive")
    labels = [-1] * width
    for block_index, block in enumerate(blocks):
        _require(bool(block), "blocks must not be empty")
        for point in block:
            _require(0 <= point < width, "block point is out of range")
            _require(labels[point] == -1, "blocks must be disjoint")
            labels[point] = block_index
    _require(all(label >= 0 for label in labels), "blocks must cover every point")
    return canonical_rgs(labels)


def is_noncrossing_blocks(blocks: Sequence[frozenset[int]]) -> bool:
    """Independent block-pair crossing check used by the second enumerator."""

    ordered = [tuple(sorted(block)) for block in blocks]
    for left_index, left in enumerate(ordered):
        for right in ordered[left_index + 1 :]:
            merged = sorted([(point, 0) for point in left] + [(point, 1) for point in right])
            colors = [color for _point, color in merged]
            for start in range(len(colors) - 3):
                for b in range(start + 1, len(colors) - 2):
                    for c in range(b + 1, len(colors) - 1):
                        for d in range(c + 1, len(colors)):
                            if colors[start] == colors[c] and colors[b] == colors[d] != colors[start]:
                                return False
    return True


@lru_cache(maxsize=None)
def noncrossing_states(width: int) -> tuple[tuple[int, ...], ...]:
    return tuple(state for state in generate_rgs_partitions(width) if is_noncrossing_rgs(state))


def independent_noncrossing_states(width: int) -> tuple[tuple[int, ...], ...]:
    states = [
        blocks_to_rgs(blocks, width)
        for blocks in generate_block_partitions(width)
        if is_noncrossing_blocks(blocks)
    ]
    return tuple(sorted(states))


def rank_state(state: Sequence[int]) -> int:
    state = validate_rgs(state)
    _require(is_noncrossing_rgs(state), "connectivity state is crossing")
    states = noncrossing_states(len(state))
    return states.index(state)


def unrank_state(width: int, rank: int) -> tuple[int, ...]:
    _require(isinstance(rank, int), "rank must be an integer")
    states = noncrossing_states(width)
    _require(0 <= rank < len(states), "rank is out of range")
    return states[rank]


def catalan(width: int) -> int:
    _require(width >= 0, "Catalan index must be nonnegative")
    return math.comb(2 * width, width) // (width + 1)


def validate_contract(contract: Mapping[str, Any]) -> dict[str, Any]:
    _require(contract.get("schema") == EXPECTED_SCHEMA, "schema drifted")
    _require(contract.get("issue") == 11, "issue must be 11")
    _require(contract.get("status") == "exact_small_width_oracle_only", "status drifted")
    maximum = contract.get("maximum_checked_width")
    _require(isinstance(maximum, int) and 1 <= maximum <= 9, "checked width is invalid")
    expected = contract.get("expected_counts")
    _require(isinstance(expected, dict), "expected counts must be an object")

    rows = []
    for width in range(1, maximum + 1):
        primary = noncrossing_states(width)
        independent = independent_noncrossing_states(width)
        declared = expected.get(str(width))
        _require(declared == catalan(width), f"declared Catalan count drifted at width {width}")
        _require(primary == independent, f"independent enumerators disagree at width {width}")
        _require(len(primary) == declared, f"state count drifted at width {width}")
        _require(
            all(unrank_state(width, rank_state(state)) == state for state in primary),
            f"rank round trip failed at width {width}",
        )
        rows.append({"width": width, "states": len(primary), "catalan": declared})

    return {
        "schema": EXPECTED_SCHEMA,
        "status": "valid_exact_small_width_codec",
        "maximum_checked_width": maximum,
        "counts": rows,
        "independent_enumerators_agree": True,
        "all_rank_round_trips_exact": True,
        "contains_transfer_matrix_result": False,
        "parent_issue": "remain open",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, default=DEFAULT_CONTRACT)
    args = parser.parse_args()
    contract = json.loads(args.contract.read_text(encoding="utf-8"))
    print(json.dumps(validate_contract(contract), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
