#!/usr/bin/env python3
"""Safe charge transfer for the square lattice with period n*(1,1).

This is the SAME physical square NN / matching site graph, not a rotated
interaction.  Use integer coordinates

    (x,y) = s*(1,1) + t*(0,1) = (s, s+t),

and quotient s modulo n.  The physical circumference is ell=n*sqrt(2).

In these coordinates NN edges with positive t increment are

    (ds,dt) = (0,1), (-1,1),

while matching adds

    (1,0)      [the e1+e2 diagonal],
    (-1,2)     [the -e1+e2 diagonal].

Therefore two frontier rows are sufficient.  Any cycle with nonzero lifted-s
deck gain is rejected as soon as it forms.

The charge-coexistence root solves equality of the NN safe Perron root at p
and the matching safe Perron root at 1-p.

Small-width deterministic research control; SciPy/NumPy required.
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
    labels: tuple[int, ...]  # oldest row, newest row
    gains: tuple[int, ...]


class DSU:
    def __init__(self, n: int):
        self.parent = list(range(n))
        self.delta = [0] * n
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


def empty_state(width: int) -> State:
    return State((-1,) * (2 * width), (0,) * (2 * width))


def step(state: State, mask: int, width: int, *, matching: bool) -> State | None:
    old_size = 2 * width
    new_base = old_size
    dsu = DSU(3 * width)
    representatives: dict[int, int] = {}

    for index, label in enumerate(state.labels):
        if label < 0:
            continue
        if label in representatives:
            dsu.join(representatives[label], index, state.gains[index])
        else:
            representatives[label] = index

    occupied_new = [i for i in range(width) if (mask >> i) & 1]

    # Matching edge e1+e2 = (ds,dt)=(1,0), within the new row.
    if matching:
        for i in occupied_new:
            j = (i + 1) % width
            if (mask >> j) & 1:
                dsu.join(new_base + i, new_base + j, (i + 1) // width)

    # NN edges with dt=1: e2=(0,1) and -e1=(-1,1).
    for i in range(width):
        old_index = width + i  # newest retained old row
        if state.labels[old_index] < 0:
            continue
        for dx in (0, -1):
            j = (i + dx) % width
            if (mask >> j) & 1:
                dsu.join(old_index, new_base + j, (i + dx) // width)

    # Matching edge -e1+e2=(-1,2), from oldest retained row to new row.
    if matching:
        for i in range(width):
            old_index = i
            if state.labels[old_index] < 0:
                continue
            j = (i - 1) % width
            if (mask >> j) & 1:
                dsu.join(old_index, new_base + j, (i - 1) // width)

    if dsu.bad:
        return None

    retained = list(range(width, 2 * width)) + [new_base + i for i in range(width)]
    labels = [-1] * (2 * width)
    gains = [0] * (2 * width)
    canonical: dict[int, tuple[int, int]] = {}
    next_label = 0

    for output_index, index in enumerate(retained):
        occupied = (
            state.labels[index] >= 0
            if index < 2 * width
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


def build(width: int, *, matching: bool):
    start = empty_state(width)
    states = [start]
    index = {start: 0}
    queue = deque([start])
    transitions: list[list[int]] = []

    while queue:
        state = queue.popleft()
        row = []
        for mask in range(1 << width):
            nxt = step(state, mask, width, matching=matching)
            if nxt is None:
                row.append(-1)
                continue
            if nxt not in index:
                index[nxt] = len(states)
                states.append(nxt)
                queue.append(nxt)
            row.append(index[nxt])
        transitions.append(row)

    if len(states) != len(transitions):
        raise AssertionError("incomplete state discovery")
    return states, transitions


class DiagonalSafeTransfer:
    def __init__(self, width: int, *, matching: bool):
        self.width = width
        self.matching = matching
        self.states, transitions = build(width, matching=matching)
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
            tol=1e-11,
            maxiter=500000,
            return_eigenvectors=False,
        )[0]
        return float(value.real)


def run_width(width: int, reference_pc: float) -> dict[str, object]:
    g4 = DiagonalSafeTransfer(width, matching=False)
    g8 = DiagonalSafeTransfer(width, matching=True)

    def equation(p: float) -> float:
        return math.log(g4.lambda0(p)) - math.log(g8.lambda0(1.0 - p))

    root = brentq(equation, 0.5, 0.8, xtol=2e-12)
    ell = width * math.sqrt(2.0)
    shift = root - reference_pc
    return {
        "period_vector": [width, width],
        "width_parameter_n": width,
        "physical_circumference": ell,
        "safe_states_G4_two_row": len(g4.states),
        "safe_states_G8_two_row": len(g8.states),
        "charge_root_p": root,
        "log_perron_ratio_at_root": equation(root),
        "root_minus_reference_pc": shift,
        "root_shift_times_physical_circumference_fourth": shift * ell**4,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--min-width", type=int, default=2)
    parser.add_argument("--max-width", type=int, default=5)
    parser.add_argument("--reference-pc", type=float, default=0.59274605079)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    result = {
        "schema": "diagonal-charge-transfer-v1",
        "model": "same square-site NN / complementary matching graph; period n*(1,1)",
        "coordinate_map": "(x,y)=s*(1,1)+t*(0,1), s mod n",
        "spin4_control": "axis cos(4theta)=+1; diagonal theta=pi/4 gives cos(4theta)=-1",
        "claim_boundary": [
            "root is located from safe Perron equality only",
            "reference_pc is diagnostic only",
            "n=2 may have short-period degeneracies; asymptotic interpretation should use increasing n",
        ],
        "reference_pc": args.reference_pc,
        "records": [
            run_width(width, args.reference_pc)
            for width in range(args.min_width, args.max_width + 1)
        ],
    }
    text = json.dumps(result, indent=2, sort_keys=True)
    if args.output:
        args.output.write_text(text + "\n", encoding="utf-8")
    print(text)


if __name__ == "__main__":
    main()
