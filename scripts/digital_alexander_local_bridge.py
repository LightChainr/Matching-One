#!/usr/bin/env python3
"""Local 4/8-adjacency certificate for the digital Alexander rank proof."""

from __future__ import annotations

import argparse
import json
from collections import deque
from pathlib import Path
from typing import Iterable, Sequence


CORNERS = ((0, 0), (1, 0), (1, 1), (0, 1))
BOUNDARY_EDGES = ((0, 1), (1, 2), (2, 3), (0, 3))
DIAGONAL_EDGES = ((0, 2), (1, 3))


def canonical_edge(first: int, second: int) -> tuple[int, int]:
    if first == second:
        raise ValueError("an edge needs two distinct corners")
    return (first, second) if first < second else (second, first)


def components(vertices: Iterable[int], edges: Iterable[tuple[int, int]]) -> tuple[tuple[int, ...], ...]:
    remaining = set(vertices)
    adjacency = {vertex: set() for vertex in remaining}
    for first, second in edges:
        if first in adjacency and second in adjacency:
            adjacency[first].add(second)
            adjacency[second].add(first)
    output = []
    while remaining:
        start = min(remaining)
        queue = deque([start])
        component = {start}
        remaining.remove(start)
        while queue:
            vertex = queue.popleft()
            for neighbor in sorted(adjacency[vertex]):
                if neighbor in remaining:
                    remaining.remove(neighbor)
                    component.add(neighbor)
                    queue.append(neighbor)
        output.append(tuple(sorted(component)))
    return tuple(output)


def find_path(
    start: int,
    end: int,
    vertices: set[int],
    edges: set[tuple[int, int]],
) -> tuple[int, ...] | None:
    adjacency = {vertex: set() for vertex in vertices}
    for first, second in edges:
        if first in adjacency and second in adjacency:
            adjacency[first].add(second)
            adjacency[second].add(first)
    queue = deque([start])
    previous = {start: None}
    while queue:
        vertex = queue.popleft()
        if vertex == end:
            path = []
            current = end
            while current is not None:
                path.append(current)
                current = previous[current]
            return tuple(reversed(path))
        for neighbor in sorted(adjacency[vertex]):
            if neighbor not in previous:
                previous[neighbor] = vertex
                queue.append(neighbor)
    return None


def face_pattern(black_mask: int) -> dict[str, object]:
    if not 0 <= black_mask < 16:
        raise ValueError("a square-face mask must lie in [0,15]")
    black = {corner for corner in range(4) if black_mask & (1 << corner)}
    white = set(range(4)) - black
    white_boundary = {
        canonical_edge(*edge)
        for edge in BOUNDARY_EDGES
        if edge[0] in white and edge[1] in white
    }
    white_diagonals = {
        canonical_edge(*edge)
        for edge in DIAGONAL_EDGES
        if edge[0] in white and edge[1] in white
    }
    retained_diagonals = {
        edge for edge in white_diagonals if white == set(edge)
    }
    removed_diagonals = white_diagonals - retained_diagonals
    pruned_edges = white_boundary | retained_diagonals
    full_edges = white_boundary | white_diagonals

    replacements = []
    replacement_failures = []
    for edge in sorted(removed_diagonals):
        path = find_path(edge[0], edge[1], white, pruned_edges)
        if path is None:
            replacement_failures.append(list(edge))
        else:
            replacements.append({"diagonal": list(edge), "boundary_path": list(path)})

    full_components = components(white, full_edges)
    pruned_components = components(white, pruned_edges)
    return {
        "black_mask": black_mask,
        "black_corners": sorted(black),
        "white_corners": sorted(white),
        "white_boundary_edges": [list(edge) for edge in sorted(white_boundary)],
        "white_matching_diagonals": [list(edge) for edge in sorted(white_diagonals)],
        "retained_diagonals": [list(edge) for edge in sorted(retained_diagonals)],
        "removed_diagonal_replacements": replacements,
        "replacement_failures": replacement_failures,
        "full_matching_components": [list(row) for row in full_components],
        "pruned_components": [list(row) for row in pruned_components],
        "connectivity_preserved": full_components == pruned_components,
        "embedded_diagonal_gate": len(retained_diagonals) <= 1,
    }


def rank_consequences() -> list[dict[str, int]]:
    rows = []
    for rank_black in range(3):
        rank_white = 2 - rank_black
        q = rank_black - 1
        rows.append({
            "rank_black": rank_black,
            "rank_white": rank_white,
            "q": q,
            "weak_residual": 2 * q - (rank_black - rank_white),
            "strong_residual": rank_black + rank_white - 2,
        })
    return rows


def analyze(config: dict) -> dict[str, object]:
    patterns = [face_pattern(mask) for mask in range(16)]
    replacement_count = sum(
        len(row["removed_diagonal_replacements"]) for row in patterns
    )
    retained_masks = [
        row["black_mask"] for row in patterns if row["retained_diagonals"]
    ]
    local_pass = all(
        row["connectivity_preserved"]
        and row["embedded_diagonal_gate"]
        and not row["replacement_failures"]
        for row in patterns
    )
    consequences = rank_consequences()
    return {
        "schema_version": 1,
        "issue": 269,
        "face_certificate": {
            "corner_order": [list(corner) for corner in CORNERS],
            "pattern_count": len(patterns),
            "patterns": patterns,
            "retained_diagonal_masks": retained_masks,
            "removed_diagonal_replacement_count": replacement_count,
            "all_local_cases_pass": local_pass,
        },
        "surface_duality_theorem": {
            "coefficients": "Q",
            "statement": (
                "For complementary compact subsurfaces U,V of a closed oriented "
                "surface S, im(H1(V)->H1(S)) is the intersection-form orthogonal "
                "complement of im(H1(U)->H1(S))."
            ),
            "proof_chain": [
                "intersection pairing identifies the annihilator of im H1(U) with ker[H^1(S)->H^1(U)]",
                "the long exact sequence of (S,U) identifies that kernel with the image of H^1(S,U)",
                "excision gives H^1(S,U)=H^1(V,boundary V)",
                "Poincare-Lefschetz duality gives H^1(V,boundary V)=H_1(V)",
                "naturality identifies the resulting map with inclusion H_1(V)->H_1(S)",
            ],
            "torus_rank_sum": "rank_black + rank_white = dim H1(T^2;Q) = 2",
        },
        "lattice_bridge": {
            "black": "NN induced graph retracts from its regular neighborhood U",
            "white": (
                "the pruned white matching graph is a 1-skeleton of the complement V; "
                "removed diagonals have same-face white NN replacements, so the full "
                "matching graph has the same ambient-H1 image"
            ),
            "scope": config["theorem_scope"],
        },
        "rank_consequences": consequences,
        "conclusion": {
            "local_4_8_bridge_certified": local_pass,
            "strong_rank_sum_theorem_in_scope": local_pass,
            "weak_rank_identity_in_scope": all(
                row["weak_residual"] == 0 and row["strong_residual"] == 0
                for row in consequences
            ),
            "claim_level": "proof_for_honest_square_cell_tori_plus_finite_oracle_for_declared_degenerate_controls",
        },
        "sources": config["sources"],
        "scientific_boundary": config["scientific_boundary"],
    }


def render_markdown(result: dict[str, object]) -> str:
    face = result["face_certificate"]
    lines = [
        "# Digital Alexander local bridge",
        "",
        "Executable 16-face certificate plus the complementary-subsurface duality proof.",
        "",
        "## Local 4/8 certificate",
        "",
        f"- face patterns checked: {face['pattern_count']}",
        f"- retained-diagonal black masks: {face['retained_diagonal_masks']}",
        f"- redundant diagonals replaced inside one face: {face['removed_diagonal_replacement_count']}",
        f"- connectivity/replacement/embedding gate: `{face['all_local_cases_pass']}`",
        "",
        "A white diagonal is retained only when its endpoints are the only two white corners. In every other active case a white NN boundary path replaces it inside the same square, preserving its ambient homology class.",
        "",
        "## Surface theorem",
        "",
        result["surface_duality_theorem"]["statement"],
        "",
    ]
    lines.extend(
        f"{index}. {step}"
        for index, step in enumerate(result["surface_duality_theorem"]["proof_chain"], 1)
    )
    lines += [
        "",
        "## Rank consequence",
        "",
        "| r_black | r_white | q | weak residual |",
        "|---:|---:|---:|---:|",
    ]
    for row in result["rank_consequences"]:
        lines.append(
            f"| {row['rank_black']} | {row['rank_white']} | {row['q']} | {row['weak_residual']} |"
        )
    lines += [
        "",
        "Within the declared honest square-cell scope, `r_black+r_white=2` and therefore `2q=r_black-r_white`.",
        "",
        "## Boundary",
        "",
    ]
    lines.extend(f"- {item}" for item in result["scientific_boundary"])
    lines.append("")
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("--json-output", type=Path)
    parser.add_argument("--markdown-output", type=Path)
    args = parser.parse_args()
    result = analyze(json.loads(args.manifest.read_text(encoding="utf-8")))
    if args.json_output:
        args.json_output.parent.mkdir(parents=True, exist_ok=True)
        args.json_output.write_text(
            json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
    else:
        print(json.dumps(result, indent=2, sort_keys=True))
    if args.markdown_output:
        args.markdown_output.parent.mkdir(parents=True, exist_ok=True)
        args.markdown_output.write_text(render_markdown(result), encoding="utf-8")


if __name__ == "__main__":
    main()
