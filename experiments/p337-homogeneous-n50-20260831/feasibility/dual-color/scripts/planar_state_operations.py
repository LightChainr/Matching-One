#!/usr/bin/env python3
"""Exact detach and adjacent-join operations on noncrossing boundary states."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Mapping, Sequence

from noncrossing_connectivity_codec import (
    blocks_to_rgs,
    canonical_rgs,
    is_noncrossing_rgs,
    noncrossing_states,
    validate_rgs,
)


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONTRACT = ROOT / "analysis" / "planar_state_operations_contract.json"
EXPECTED_SCHEMA = "matching-one/planar-state-operations/v1"


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _checked_state_and_point(state: Sequence[int], point: int) -> tuple[tuple[int, ...], int]:
    state = validate_rgs(state)
    _require(is_noncrossing_rgs(state), "connectivity state is crossing")
    _require(isinstance(point, int) and not isinstance(point, bool), "point must be an integer")
    _require(0 <= point < len(state), "point is out of range")
    return state, point


def detach_rgs(state: Sequence[int], point: int) -> tuple[int, ...]:
    """Detach one boundary point into its own singleton block."""

    state, point = _checked_state_and_point(state, point)
    labels = list(state)
    labels[point] = max(state) + 1
    output = canonical_rgs(labels)
    _require(is_noncrossing_rgs(output), "detach produced a crossing state")
    return output


def join_cyclic_adjacent_rgs(state: Sequence[int], point: int) -> tuple[int, ...]:
    """Join ``point`` with its next cyclic neighbor."""

    state, point = _checked_state_and_point(state, point)
    neighbor = (point + 1) % len(state)
    source = state[neighbor]
    target = state[point]
    labels = [target if label == source else label for label in state]
    output = canonical_rgs(labels)
    _require(is_noncrossing_rgs(output), "cyclic-adjacent join produced a crossing state")
    return output


def rgs_to_blocks(state: Sequence[int]) -> tuple[frozenset[int], ...]:
    state = validate_rgs(state)
    return tuple(
        frozenset(index for index, label in enumerate(state) if label == block)
        for block in range(max(state) + 1)
    )


def detach_blocks(state: Sequence[int], point: int) -> tuple[int, ...]:
    """Independent detach implementation on explicit block sets."""

    state, point = _checked_state_and_point(state, point)
    blocks = [set(block) for block in rgs_to_blocks(state)]
    containing = next(index for index, block in enumerate(blocks) if point in block)
    blocks[containing].remove(point)
    blocks = [block for block in blocks if block]
    blocks.append({point})
    return blocks_to_rgs(tuple(frozenset(block) for block in blocks), len(state))


def join_cyclic_adjacent_blocks(state: Sequence[int], point: int) -> tuple[int, ...]:
    """Independent adjacent join implementation on explicit block sets."""

    state, point = _checked_state_and_point(state, point)
    neighbor = (point + 1) % len(state)
    blocks = [set(block) for block in rgs_to_blocks(state)]
    left = next(index for index, block in enumerate(blocks) if point in block)
    right = next(index for index, block in enumerate(blocks) if neighbor in block)
    if left != right:
        blocks[left] |= blocks[right]
        del blocks[right]
    return blocks_to_rgs(tuple(frozenset(block) for block in blocks), len(state))


def validate_contract(contract: Mapping[str, Any]) -> dict[str, Any]:
    _require(contract.get("schema") == EXPECTED_SCHEMA, "schema drifted")
    _require(contract.get("issue") == 11, "issue must be 11")
    _require(contract.get("status") == "exact_small_width_operations_only", "status drifted")
    maximum = contract.get("maximum_checked_width")
    _require(isinstance(maximum, int) and 1 <= maximum <= 9, "checked width is invalid")

    states_checked = 0
    cases = 0
    changed_detaches = 0
    changed_joins = 0
    for width in range(1, maximum + 1):
        for state in noncrossing_states(width):
            states_checked += 1
            for point in range(width):
                cases += 1
                detached = detach_rgs(state, point)
                joined = join_cyclic_adjacent_rgs(state, point)
                _require(detached == detach_blocks(state, point), "detach implementations disagree")
                _require(
                    joined == join_cyclic_adjacent_blocks(state, point),
                    "join implementations disagree",
                )
                _require(detach_rgs(detached, point) == detached, "detach is not idempotent")
                _require(join_cyclic_adjacent_rgs(joined, point) == joined, "join is not idempotent")
                _require(detached in noncrossing_states(width), "detach left the state space")
                _require(joined in noncrossing_states(width), "join left the state space")
                changed_detaches += detached != state
                changed_joins += joined != state

    _require(states_checked == contract.get("expected_states_checked"), "state coverage drifted")
    _require(cases == contract.get("expected_cases_per_operation"), "operation coverage drifted")
    return {
        "schema": EXPECTED_SCHEMA,
        "status": "valid_exact_planar_state_operations",
        "maximum_checked_width": maximum,
        "states_checked": states_checked,
        "cases_per_operation": cases,
        "combined_operation_cases": 2 * cases,
        "changed_detaches": changed_detaches,
        "changed_joins": changed_joins,
        "independent_implementations_agree": True,
        "all_outputs_canonical_noncrossing": True,
        "all_operations_idempotent": True,
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
