#!/usr/bin/env python3
"""Enumerate the 16 local square certificates for the 4/8 proof bridge."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


CORNERS = ("SW", "SE", "NE", "NW")
NN_EDGES = ((0, 1), (1, 2), (2, 3), (3, 0))
DIAGONALS = ((0, 2), (1, 3))


def canonical_edge(edge: tuple[int, int]) -> tuple[int, int]:
    return tuple(sorted(edge))  # type: ignore[return-value]


def connected(vertices: list[int], edges: list[tuple[int, int]], source: int, target: int) -> bool:
    seen = {source}
    frontier = [source]
    adjacency = {vertex: [] for vertex in vertices}
    for left, right in edges:
        adjacency[left].append(right)
        adjacency[right].append(left)
    while frontier:
        vertex = frontier.pop()
        for neighbour in adjacency[vertex]:
            if neighbour not in seen:
                seen.add(neighbour)
                frontier.append(neighbour)
    return target in seen


def classify(mask: int) -> dict[str, object]:
    white = [index for index in range(4) if mask & (1 << index)]
    white_set = set(white)
    nn = [canonical_edge(edge) for edge in NN_EDGES if set(edge) <= white_set]
    diagonals = [canonical_edge(edge) for edge in DIAGONALS if set(edge) <= white_set]
    full_edges = sorted(set(nn + diagonals))
    kept_diagonals = diagonals if len(white) == 2 and len(diagonals) == 1 else []
    spine_edges = sorted(set(nn + kept_diagonals))
    removed = sorted(set(diagonals) - set(kept_diagonals))
    replacements = []
    for left, right in removed:
        if not connected(white, spine_edges, left, right):
            raise AssertionError(f"mask {mask}: removed diagonal has no NN replacement")
        replacements.append([CORNERS[left], CORNERS[right]])

    if not white:
        category = "all_black"
    elif len(white) == 1:
        category = "one_white"
    elif len(white) == 2 and diagonals:
        category = "two_diagonal_white"
    elif len(white) == 2:
        category = "two_adjacent_white"
    elif len(white) == 3:
        category = "three_white"
    else:
        category = "all_white"

    return {
        "mask": format(mask, "04b"),
        "white_corners": [CORNERS[index] for index in white],
        "category": category,
        "matching_edges": [[CORNERS[a], CORNERS[b]] for a, b in full_edges],
        "incidence_spokes": [[CORNERS[index], "c_f"] for index in white],
        "spine_edges": [[CORNERS[a], CORNERS[b]] for a, b in spine_edges],
        "removed_diagonals_with_NN_replacement": replacements,
        "nerve_simplex_dimension": len(white) - 1 if white else None,
        "local_two_cell": len(white) == 4,
    }


def build_certificate() -> dict[str, object]:
    rows = [classify(mask) for mask in range(16)]
    counts: dict[str, int] = {}
    for row in rows:
        category = str(row["category"])
        counts[category] = counts.get(category, 0) + 1
    expected = {
        "all_black": 1,
        "one_white": 4,
        "two_adjacent_white": 4,
        "two_diagonal_white": 2,
        "three_white": 4,
        "all_white": 1,
    }
    if counts != expected:
        raise AssertionError((counts, expected))
    return {
        "schema": "matching-one/digital-alexander-local-certificate/v1",
        "issue": 269,
        "status": "exact_16_pattern_local_certificate",
        "corner_order": list(CORNERS),
        "category_counts": counts,
        "checks": {
            "all_16_masks_present": len(rows) == 16,
            "only_two_diagonal_white_keeps_a_diagonal": all(
                (
                    any(set(edge) in ({"SW", "NE"}, {"SE", "NW"}) for edge in row["spine_edges"])
                )
                == (row["category"] == "two_diagonal_white")
                for row in rows
            ),
            "every_removed_diagonal_has_NN_replacement": True,
        },
        "patterns": rows,
        "boundary": "local regular-cell certificate; small nonregular period quotients remain covered by the finite exact oracle",
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("analysis/digital_alexander_local_certificate.json"),
    )
    args = parser.parse_args()
    payload = build_certificate()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
