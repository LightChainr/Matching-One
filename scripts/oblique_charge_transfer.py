#!/usr/bin/env python3
"""Generic primitive-direction safe charge transfer for square-site percolation.

Choose a primitive integer circumference direction u=(a,b).  Complete it to
an SL(2,Z) basis (u,v).  In coordinates x=s*u+t*v, quotient s modulo n and
advance in t.  Physical NN/matching edges become a fixed finite-range edge set
(ds,dt); a frontier memory of R=max positive dt is exact.

The safe transfer rejects any occupied cycle with nonzero lifted-s deck gain.
The charge root solves equality of the NN Perron root at p and the matching
Perron root at 1-p.

This script is intended for small directional controls.  State spaces grow
rapidly with both n and the row memory R.
"""
from __future__ import annotations

import argparse
import json
import math
from collections import Counter, deque
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigs


@dataclass(frozen=True)
class State:
    labels: tuple[int, ...]
    gains: tuple[int, ...]


class DSU:
    def __init__(self, size: int):
        self.parent = list(range(size))
        self.delta = [0] * size
        self.bad = False

    def find(self, a: int) -> tuple[int, int]:
        if self.parent[a] != a:
            root, gain = self.find(self.parent[a])
            self.delta[a] += gain
            self.parent[a] = root
        return self.parent[a], self.delta[a]

    def join(self, a: int, b: int, gain_b_minus_a: int) -> None:
        ra, da = self.find(a)
        rb, db = self.find(b)
        if ra == rb:
            if db - da != gain_b_minus_a:
                self.bad = True
            return
        self.parent[rb] = ra
        self.delta[rb] = gain_b_minus_a + da - db


def extended_gcd(a: int, b: int) -> tuple[int, int, int]:
    if b == 0:
        return abs(a), 1 if a >= 0 else -1, 0
    g, x1, y1 = extended_gcd(b, a % b)
    return g, y1, x1 - (a // b) * y1


def bezout_complement(a: int, b: int) -> tuple[int, int]:
    """Return v=(c,d) with det((a,b),(c,d))=1."""
    g, x, y = extended_gcd(a, b)
    if g != 1:
        raise ValueError("direction must be primitive")
    d = x
    c = -y
    if a * d - b * c != 1:
        raise AssertionError("Bezout orientation failure")
    return c, d


def transformed_edges(direction: tuple[int, int], *, matching: bool):
    a, b = direction
    c, d = bezout_complement(a, b)
    generators = [(1, 0), (0, 1)]
    if matching:
        generators += [(1, 1), (1, -1)]

    edges = set()
    for dx, dy in generators:
        ds = d * dx - c * dy
        dt = -b * dx + a * dy
        if dt < 0 or (dt == 0 and ds < 0):
            ds, dt = -ds, -dt
        edges.add((ds, dt))

    ordered = sorted(edges, key=lambda item: (item[1], item[0]))
    memory = max(dt for _, dt in ordered)
    return (c, d), ordered, memory


def empty_state(width: int, memory: int) -> State:
    return State((-1,) * (width * memory), (0,) * (width * memory))


def step(
    state: State,
    mask: int,
    width: int,
    edges: list[tuple[int, int]],
    memory: int,
) -> State | None:
    old_size = memory * width
    new_base = old_size
    dsu = DSU((memory + 1) * width)
    representatives: dict[int, int] = {}

    for index, label in enumerate(state.labels):
        if label < 0:
            continue
        if label in representatives:
            dsu.join(representatives[label], index, state.gains[index])
        else:
            representatives[label] = index

    for ds, dt in edges:
        if dt == 0:
            for i in range(width):
                if not (mask >> i) & 1:
                    continue
                j = (i + ds) % width
                if (mask >> j) & 1:
                    dsu.join(new_base + i, new_base + j, (i + ds) // width)
            continue

        row_position = memory - dt
        for i in range(width):
            old_index = row_position * width + i
            if state.labels[old_index] < 0:
                continue
            j = (i + ds) % width
            if (mask >> j) & 1:
                dsu.join(old_index, new_base + j, (i + ds) // width)

    if dsu.bad:
        return None

    retained = list(range(width, old_size)) + [new_base + i for i in range(width)]
    labels = [-1] * old_size
    gains = [0] * old_size
    canonical: dict[int, tuple[int, int]] = {}
    next_label = 0

    for output_index, index in enumerate(retained):
        occupied = (
            state.labels[index] >= 0
            if index < old_size
            else bool(mask >> (index - new_base) & 1)
        )
        if not occupied:
            continue
        root, gain = dsu.find(index)
        if root not in canonical:
            canonical[root] = (next_label, gain)
            next_label += 1
        label, base_gain = canonical[root]
        labels[output_index] = label
        gains[output_index] = gain - base_gain

    return State(tuple(labels), tuple(gains))


def build(
    width: int,
    direction: tuple[int, int],
    *,
    matching: bool,
    state_cap: int,
):
    complement, edges, memory = transformed_edges(direction, matching=matching)
    start = empty_state(width, memory)
    states = [start]
    index = {start: 0}
    queue = deque([start])
    transitions = []

    while queue:
        state = queue.popleft()
        row = []
        for mask in range(1 << width):
            nxt = step(state, mask, width, edges, memory)
            if nxt is None:
                row.append(-1)
                continue
            if nxt not in index:
                if len(states) >= state_cap:
                    raise RuntimeError("state cap reached")
                index[nxt] = len(states)
                states.append(nxt)
                queue.append(nxt)
            row.append(index[nxt])
        transitions.append(row)

    return states, transitions, complement, edges, memory


class ObliqueSafeTransfer:
    def __init__(
        self,
        width: int,
        direction: tuple[int, int],
        *,
        matching: bool,
        state_cap: int,
    ):
        self.width = width
        self.direction = direction
        (
            self.states,
            transitions,
            self.complement,
            self.edges,
            self.memory,
        ) = build(
            width,
            direction,
            matching=matching,
            state_cap=state_cap,
        )

        counter: Counter[tuple[int, int, int]] = Counter()
        for source, row in enumerate(transitions):
            for mask, destination in enumerate(row):
                if destination >= 0:
                    counter[(source, destination, mask.bit_count())] += 1
        keys = list(counter)
        self.src = np.array([key[0] for key in keys], dtype=np.int32)
        self.dst = np.array([key[1] for key in keys], dtype=np.int32)
        self.occ = np.array([key[2] for key in keys], dtype=np.int16)
        self.counts = np.array([counter[key] for key in keys], dtype=float)

    def matrix(self, p: float):
        q = 1.0 - p
        data = self.counts * p**self.occ * q ** (self.width - self.occ)
        return coo_matrix(
            (data, (self.src, self.dst)),
            shape=(len(self.states), len(self.states)),
        ).tocsr()

    def lambda0(self, p: float) -> float:
        matrix = self.matrix(p)
        if matrix.shape[0] < 200:
            return float(max(np.linalg.eigvals(matrix.toarray()).real))
        value = eigs(
            matrix,
            k=1,
            which="LM",
            tol=1e-10,
            maxiter=500000,
            return_eigenvectors=False,
        )[0]
        return float(value.real)


def cos4(direction: tuple[int, int]) -> float:
    a, b = direction
    denominator = (a * a + b * b) ** 2
    return (a**4 - 6 * a * a * b * b + b**4) / denominator


def run_width(
    width: int,
    direction: tuple[int, int],
    reference_pc: float,
    state_cap: int,
) -> dict[str, object]:
    g4 = ObliqueSafeTransfer(
        width, direction, matching=False, state_cap=state_cap
    )
    g8 = ObliqueSafeTransfer(
        width, direction, matching=True, state_cap=state_cap
    )

    def equation(p: float) -> float:
        return math.log(g4.lambda0(p)) - math.log(g8.lambda0(1.0 - p))

    root = brentq(equation, 0.5, 0.8, xtol=3e-11)
    a, b = direction
    length = math.hypot(a, b)
    ell = width * length
    angular = cos4(direction)
    scaled_shift = (root - reference_pc) * ell**4
    amplitude_estimate = (
        -scaled_shift / angular if abs(angular) > 1e-12 else None
    )
    return {
        "n": width,
        "period_vector": [width * a, width * b],
        "primitive_direction": [a, b],
        "bezout_complement_G4": list(g4.complement),
        "transformed_edges_G4": [list(edge) for edge in g4.edges],
        "transformed_edges_G8": [list(edge) for edge in g8.edges],
        "row_memory_G4": g4.memory,
        "row_memory_G8": g8.memory,
        "safe_states_G4": len(g4.states),
        "safe_states_G8": len(g8.states),
        "physical_circumference": ell,
        "cos_4theta": angular,
        "charge_root_p": root,
        "root_minus_reference_pc": root - reference_pc,
        "root_shift_times_ell_fourth": scaled_shift,
        "spin4_amplitude_estimate_minus_shift_over_cos4": amplitude_estimate,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--a", type=int, required=True)
    parser.add_argument("--b", type=int, required=True)
    parser.add_argument("--min-n", type=int, default=2)
    parser.add_argument("--max-n", type=int, default=4)
    parser.add_argument("--state-cap", type=int, default=300000)
    parser.add_argument("--reference-pc", type=float, default=0.59274605079)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    direction = (args.a, args.b)
    if math.gcd(abs(args.a), abs(args.b)) != 1:
        raise ValueError("(a,b) must be primitive")

    result = {
        "schema": "oblique-charge-transfer-v1",
        "model": "same square-site NN / complementary matching graph in an SL(2,Z) oblique basis",
        "direction": list(direction),
        "reference_pc": args.reference_pc,
        "spin4_prediction": "p_root-pc ~ -A*cos(4theta)/ell^4",
        "claim_boundary": [
            "charge roots are located from safe Perron equality only",
            "reference_pc is diagnostic only",
            "state complexity grows rapidly with row memory; small n can have larger corrections"
        ],
        "records": [
            run_width(n, direction, args.reference_pc, args.state_cap)
            for n in range(args.min_n, args.max_n + 1)
        ],
    }
    text = json.dumps(result, indent=2, sort_keys=True)
    if args.output:
        args.output.write_text(text + "\n", encoding="utf-8")
    print(text)


if __name__ == "__main__":
    main()
