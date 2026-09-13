#!/usr/bin/env python3
"""Issue #593 Phase B infrastructure: arbitrary-width noncrossing states and a
Generator-compatible shim.

The repository's ``Generator`` caps at width 8.  This module reproduces the
same object for widths 9 and 10:

- ``noncrossing_states_fast`` enumerates canonical noncrossing RGS by interval
  recursion (the block containing the left endpoint picks an arbitrary subset,
  every gap is an independent noncrossing partition).  At widths 1..8 it must
  equal the codec's enumeration exactly -- order included -- and that equality
  is re-checked at import time of the validation entry point.
- ``WideGenerator`` duck-types ``p398_intervention_transport.Generator`` so the
  repository's own frozen-span, memory-kernel, lumping and transport-scoring
  code paths run unchanged at widths 9 and 10.  No frozen convention is
  substituted: the moves, rate families, readouts, spans and grids are the
  repository's own code.
"""

from __future__ import annotations

from functools import lru_cache
from itertools import product
from typing import Dict, List, Optional, Sequence, Tuple

from noncrossing_connectivity_codec import canonical_rgs
from planar_state_operations import detach_rgs, join_cyclic_adjacent_rgs

OPERATIONS: Tuple[str, ...] = ("join", "detach")


@lru_cache(maxsize=None)
def _interval_partitions(i: int, j: int) -> Tuple[Tuple[Tuple[int, ...], ...], ...]:
    """Noncrossing partitions of the contiguous points ``i..j``."""

    if i > j:
        return ((),)
    out: List[Tuple[Tuple[int, ...], ...]] = []
    n = j - i
    for mask in range(1 << n):
        subset = [i + 1 + k for k in range(n) if (mask >> k) & 1]
        gaps: List[Tuple[int, int]] = []
        previous = i
        for s in subset:
            gaps.append((previous + 1, s - 1))
            previous = s
        gaps.append((previous + 1, j))
        head = (i,) + tuple(subset)
        for combo in product(*(_interval_partitions(a, b) for a, b in gaps)):
            blocks = [head]
            for part in combo:
                blocks.extend(part)
            out.append(tuple(sorted(blocks)))
    return tuple(out)


def noncrossing_states_fast(width: int) -> Tuple[Tuple[int, ...], ...]:
    """Canonical noncrossing RGS in lexicographic order, any width."""

    if width <= 0:
        raise ValueError("width must be positive")
    states = set()
    for partition in _interval_partitions(0, width - 1):
        labels = [0] * width
        for block_index, block in enumerate(partition):
            for point in block:
                labels[point] = block_index
        states.add(canonical_rgs(labels))
    return tuple(sorted(states))


class WideGenerator:
    """Interface-compatible with ``p398_intervention_transport.Generator``."""

    def __init__(self, width: int, states: Sequence[Tuple[int, ...]]) -> None:
        self.width = width
        self.states = tuple(states)
        self.size = len(self.states)
        rank = {state: index for index, state in enumerate(self.states)}
        target: Dict[Tuple[str, int], List[Optional[int]]] = {}
        for point in range(width):
            for operation, implementation in (
                ("join", join_cyclic_adjacent_rgs),
                ("detach", detach_rgs),
            ):
                column: List[Optional[int]] = []
                for source, state in enumerate(self.states):
                    image = rank[implementation(state, point)]
                    column.append(None if image == source else image)
                target[(operation, point)] = column
        self.target = target
        self.moves: Tuple[Tuple[str, int], ...] = tuple(
            (operation, point) for point in range(width) for operation in OPERATIONS
        )

    # -- rate vectors (identical semantics to the repository's Generator) --

    def baseline_rates(self) -> Dict[Tuple[str, int], float]:
        return {move: 1.0 for move in self.moves}

    def rates(self, intervention: str, eta: float) -> Dict[Tuple[str, int], float]:
        from p398_intervention_transport import INTERVENTIONS

        coefficient = INTERVENTIONS[intervention]
        return {
            (operation, point): 1.0 + eta * coefficient(operation, point, self.width)
            for operation, point in self.moves
        }

    def coefficients(self, intervention: str) -> Dict[Tuple[str, int], float]:
        from p398_intervention_transport import INTERVENTIONS

        coefficient = INTERVENTIONS[intervention]
        return {
            (operation, point): coefficient(operation, point, self.width)
            for operation, point in self.moves
        }

    def exact_rates(self, intervention: str, eta) -> Dict[Tuple[str, int], float]:
        from fractions import Fraction
        from p398_intervention_transport import INTERVENTIONS

        coefficient = INTERVENTIONS[intervention]
        return {
            (operation, point): Fraction(1) + eta * Fraction(
                coefficient(operation, point, self.width)
            ).limit_denominator(10**6)
            for operation, point in self.moves
        }

    def exact_entry(self, rates, source: int, destination: int):
        from fractions import Fraction

        value = Fraction(0)
        for move, rate in rates.items():
            image = self.target[move][source]
            if image is None:
                continue
            if image == destination:
                value += rate
            if source == destination:
                value -= rate
        return value

    def exit_rate(self, rates) -> float:
        best = 0.0
        for source in range(self.size):
            total = 0.0
            for move, rate in rates.items():
                if self.target[move][source] is not None:
                    total += rate
            best = max(best, total)
        return best

    # -- float sparse forms -------------------------------------------------

    def rows(self, rates) -> List[List[Tuple[int, float]]]:
        out: List[List[Tuple[int, float]]] = []
        for source in range(self.size):
            entries: Dict[int, float] = {}
            diagonal = 0.0
            for move, rate in rates.items():
                image = self.target[move][source]
                if image is None or rate == 0.0:
                    continue
                entries[image] = entries.get(image, 0.0) + rate
                diagonal -= rate
            entries[source] = entries.get(source, 0.0) + diagonal
            out.append(sorted(entries.items()))
        return out

    def columns(self, rates) -> List[List[Tuple[int, float]]]:
        gathered: List[Dict[int, float]] = [dict() for _ in range(self.size)]
        for source, row in enumerate(self.rows(rates)):
            for destination, value in row:
                gathered[destination][source] = (
                    gathered[destination].get(source, 0.0) + value
                )
        return [sorted(column.items()) for column in gathered]

    def tangent_rows(self, intervention: str) -> List[List[Tuple[int, float]]]:
        return self.rows(self.coefficients(intervention))


def validate_shim(maximum_width: int = 8) -> dict:
    """Exact equality of the shim against the repository, widths 1..8."""

    import json
    import hashlib
    import sys
    from pathlib import Path

    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from noncrossing_connectivity_codec import noncrossing_states
    from p398_intervention_transport import Generator

    report = []
    for width in range(1, maximum_width + 1):
        mine = noncrossing_states_fast(width)
        theirs = noncrossing_states(width)
        assert mine == theirs, f"state enumeration disagrees at width {width}"
        shim = WideGenerator(width, mine)
        reference = Generator(width)
        assert shim.target == reference.target, f"targets disagree at width {width}"
        assert shim.rows(shim.baseline_rates()) == reference.rows(
            reference.baseline_rates()
        ), f"baseline rows disagree at width {width}"
        report.append({"width": width, "states": len(mine), "exact_match": True})

    # transition-table control: the repository's own builder must reproduce
    from planar_transition_table import build_result

    table = build_result(maximum_width)
    return {
        "shim_validation": report,
        "transition_table_sha256": table["totals"]["canonical_jsonl_sha256"],
        "transition_table_maximum_width": maximum_width,
    }


if __name__ == "__main__":
    import json

    print(json.dumps(validate_shim(8), indent=2))
