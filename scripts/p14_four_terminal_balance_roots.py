#!/usr/bin/env python3
"""Score the frozen four-terminal corpus by an exact necessary balance root."""

from __future__ import annotations

import argparse
from fractions import Fraction
from itertools import combinations, permutations
import json
from math import comb, gcd
from pathlib import Path
from typing import Iterable, Sequence


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = ROOT / "results" / "terminal-reliability" / "bounded-four-terminal-corpus.json"
REFERENCE_PC = Fraction(5927460507896, 10**13)


def trim(poly: Sequence[Fraction]) -> list[Fraction]:
    out = list(poly)
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def derivative(poly: Sequence[Fraction]) -> list[Fraction]:
    return trim([Fraction(i) * poly[i] for i in range(1, len(poly))] or [Fraction(0)])


def divmod_poly(a: Sequence[Fraction], b: Sequence[Fraction]) -> tuple[list[Fraction], list[Fraction]]:
    dividend, divisor = trim(a), trim(b)
    if divisor == [0]:
        raise ZeroDivisionError("zero polynomial")
    quotient = [Fraction(0)] * max(1, len(dividend) - len(divisor) + 1)
    remainder = dividend[:]
    while remainder != [0] and len(remainder) >= len(divisor):
        shift = len(remainder) - len(divisor)
        factor = remainder[-1] / divisor[-1]
        quotient[shift] += factor
        for index, value in enumerate(divisor):
            remainder[index + shift] -= factor * value
        remainder = trim(remainder)
    return trim(quotient), remainder


def monic(poly: Sequence[Fraction]) -> list[Fraction]:
    out = trim(poly)
    if out == [0]:
        return out
    return [value / out[-1] for value in out]


def gcd_poly(a: Sequence[Fraction], b: Sequence[Fraction]) -> list[Fraction]:
    left, right = trim(a), trim(b)
    while right != [0]:
        _, remainder = divmod_poly(left, right)
        left, right = right, remainder
    return monic(left)


def square_free(poly: Sequence[Fraction]) -> list[Fraction]:
    out = trim(poly)
    common = gcd_poly(out, derivative(out))
    quotient, remainder = divmod_poly(out, common)
    if remainder != [0]:
        raise AssertionError("square-free division failed")
    return trim(quotient)


def eval_poly(poly: Sequence[Fraction], x: Fraction) -> Fraction:
    value = Fraction(0)
    for coefficient in reversed(poly):
        value = value * x + coefficient
    return value


def sturm_sequence(poly: Sequence[Fraction]) -> list[list[Fraction]]:
    first = square_free(poly)
    second = derivative(first)
    sequence = [first, second]
    while second != [0]:
        _, remainder = divmod_poly(first, second)
        if remainder == [0]:
            break
        remainder = [-value for value in remainder]
        sequence.append(remainder)
        first, second = second, remainder
    return sequence


def variations(sequence: Sequence[Sequence[Fraction]], x: Fraction) -> int:
    signs = []
    for poly in sequence:
        value = eval_poly(poly, x)
        if value:
            signs.append(1 if value > 0 else -1)
    return sum(left != right for left, right in zip(signs, signs[1:]))


def isolate_roots(poly: Sequence[Fraction], bits: int = 96) -> list[tuple[Fraction, Fraction]]:
    sequence = sturm_sequence(poly)
    lo, hi = Fraction(0), Fraction(1)
    if eval_poly(sequence[0], lo) == 0 or eval_poly(sequence[0], hi) == 0:
        raise ValueError("open-unit root isolator requires nonroot endpoints")
    total = variations(sequence, lo) - variations(sequence, hi)
    result: list[tuple[Fraction, Fraction]] = []

    def recurse(left: Fraction, right: Fraction, count: int) -> None:
        if count == 0:
            return
        if count == 1 and right - left <= Fraction(1, 1 << bits):
            result.append((left, right))
            return
        middle = (left + right) / 2
        vmid = variations(sequence, middle)
        # Sturm's V(left)-V(right) convention counts roots in (left,right].
        # Keeping that convention also handles an exact dyadic midpoint without
        # passing an ambiguous open endpoint into the next recursion.
        left_count = variations(sequence, left) - vmid
        right_count = vmid - variations(sequence, right)
        if min(left_count, right_count) < 0 or left_count + right_count != count:
            raise AssertionError("Sturm accounting failed")
        recurse(left, middle, left_count)
        recurse(middle, right, right_count)

    recurse(lo, hi, total)
    return sorted(result)


def primitive_integer(poly: Sequence[Fraction]) -> list[int]:
    den = 1
    for value in poly:
        den = den * value.denominator // gcd(den, value.denominator)
    integers = [int(value * den) for value in poly]
    common = 0
    for value in integers:
        common = gcd(common, abs(value))
    integers = [value // common for value in integers]
    if integers[-1] < 0:
        integers = [-value for value in integers]
    return integers


def bernstein_to_power(counts: Sequence[int]) -> list[Fraction]:
    degree = len(counts) - 1
    result = [Fraction(0)] * (degree + 1)
    for open_edges, count in enumerate(counts):
        for tail in range(degree - open_edges + 1):
            result[open_edges + tail] += count * ((-1) ** tail) * comb(degree - open_edges, tail)
    return trim(result)


def decode_edges(encoding: str) -> tuple[int, int, tuple[tuple[int, int], ...]]:
    terminal_count, vertex_count, bits = encoding.split(":")
    terminal_count, vertex_count = int(terminal_count), int(vertex_count)
    slots = tuple(combinations(range(vertex_count), 2))
    if len(bits) != len(slots):
        raise ValueError("bad graph encoding")
    return terminal_count, vertex_count, tuple(edge for edge, bit in zip(slots, bits) if bit == "1")


def component_count(vertices: Iterable[int], edges: Sequence[tuple[int, int]]) -> int:
    remaining = set(vertices)
    count = 0
    while remaining:
        count += 1
        stack = [remaining.pop()]
        while stack:
            current = stack.pop()
            neighbors = {v if u == current else u for u, v in edges if u == current or v == current}
            for neighbor in neighbors & remaining:
                remaining.remove(neighbor)
                stack.append(neighbor)
    return count


def graph_structure(encoding: str) -> dict[str, object]:
    terminal_count, vertex_count, edges = decode_edges(encoding)
    vertices = tuple(range(vertex_count))
    base_components = component_count(vertices, edges)
    bridges = []
    for edge in edges:
        reduced = tuple(item for item in edges if item != edge)
        if component_count(vertices, reduced) > base_components:
            bridges.append(list(edge))
    articulations = []
    for vertex in vertices:
        kept = tuple(item for item in vertices if item != vertex)
        reduced = tuple((u, v) for u, v in edges if u != vertex and v != vertex)
        if len(kept) > 1 and component_count(kept, reduced) > base_components:
            articulations.append(vertex)

    edge_set = set(edges)
    automorphisms = []
    for perm in permutations(range(terminal_count)):
        mapping = tuple(perm) + tuple(range(terminal_count, vertex_count))
        moved = {tuple(sorted((mapping[u], mapping[v]))) for u, v in edges}
        if moved == edge_set:
            automorphisms.append(perm)
    orbit_zero = sorted({perm[0] for perm in automorphisms})
    return {
        "bridges": bridges,
        "articulations": articulations,
        "terminal_automorphism_order": len(automorphisms),
        "terminal_zero_orbit": orbit_zero,
        "terminal_transitive": len(orbit_zero) == terminal_count,
    }


def partition_key(values: Sequence[int]) -> str:
    groups: dict[int, list[str]] = {}
    for index, value in enumerate(values):
        groups.setdefault(value, []).append(str(index))
    return "|".join("".join(group) for group in groups.values())


def candidate_score(candidate: dict[str, object]) -> dict[str, object]:
    rows = {
        partition_key(row["terminal_partition_rgs"]): row["bernstein_counts_by_open_edge_count"]
        for row in candidate["terminal_partition_polynomials"]
    }
    degree = int(candidate["edge_count"])
    zeros = [0] * (degree + 1)
    all_connected = bernstein_to_power(rows.get("0123", zeros))
    all_separate = bernstein_to_power(rows.get("0|1|2|3", zeros))
    length = max(len(all_connected), len(all_separate))
    all_connected += [Fraction(0)] * (length - len(all_connected))
    all_separate += [Fraction(0)] * (length - len(all_separate))
    balance = trim([left - right for left, right in zip(all_connected, all_separate)])
    squarefree = primitive_integer(square_free(balance))
    intervals = isolate_roots(balance)

    pair_keys = ("01|23", "02|13", "03|12")
    zero_crossings = [key for key in pair_keys if not any(rows.get(key, zeros))]
    structure = graph_structure(str(candidate["canonical_graph_encoding"]))
    roots = []
    for lower, upper in intervals:
        midpoint = (lower + upper) / 2
        roots.append(
            {
                "lower": str(lower),
                "upper": str(upper),
                "midpoint": format(float(midpoint), ".16g"),
                "distance_to_reference": format(float(abs(midpoint - REFERENCE_PC)), ".16g"),
                "exact_rational": None,
            }
        )
    roots.sort(key=lambda row: float(row["distance_to_reference"]))
    gate = (
        not structure["bridges"]
        and not structure["articulations"]
        and structure["terminal_transitive"]
        and bool(zero_crossings)
    )
    return {
        "canonical_graph_encoding": candidate["canonical_graph_encoding"],
        "edge_count": degree,
        "internal_degree": candidate["internal_degree"],
        "balance_power_coefficients_low_to_high": primitive_integer(balance),
        "balance_squarefree_coefficients_low_to_high": squarefree,
        "open_unit_roots": roots,
        "zero_crossing_pair_partitions": zero_crossings,
        "structure": structure,
        "passes_primary_structural_gate": gate,
    }


def build_result(source: Path) -> dict[str, object]:
    payload = json.loads(source.read_text(encoding="utf-8"))
    candidates = [candidate_score(row) for row in payload["candidates"]]
    candidates.sort(
        key=lambda row: min(
            (float(root["distance_to_reference"]) for root in row["open_unit_roots"]),
            default=float("inf"),
        )
    )
    gated = [row for row in candidates if row["passes_primary_structural_gate"]]
    return {
        "schema": "matching-one/p14-four-terminal-balance-roots/v1",
        "input_schema": payload["schema"],
        "input_corpus_sha256": payload["enumeration"]["corpus_sha256"],
        "reference_pc_descriptive": str(REFERENCE_PC),
        "balance_definition": "P(0123)-P(0|1|2|3)",
        "candidate_count": len(candidates),
        "candidate_with_open_unit_root_count": sum(bool(row["open_unit_roots"]) for row in candidates),
        "primary_structural_gate_survivor_count": len(gated),
        "candidates": candidates,
        "claim_boundary": "necessary scalar screen only; no embedding, duality theorem, critical point, bound, or square-site identification",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = build_result(args.input)
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
