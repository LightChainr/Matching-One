#!/usr/bin/env python3
"""Safe-frontier transfer for fixed-width charge coexistence.

For a cylinder of circumference w, construct the finite frontier automaton that
rejects a row transition as soon as horizontal homology is created. The Perron
root lambda^0_G,w(p) of the safe substochastic transfer gives

    I^0_G,w(p) = -log lambda^0_G,w(p).

The charge-coexistence root solves

    lambda^0_4,w(p) = lambda^0_8,w(1-p).

The implementation keeps lifted horizontal gains, so periodic seam winding is
detected rather than inferred from coarse connectivity flags. Width two
therefore retains the two physically distinct lifted horizontal bonds.

This is a small-width research control, not a production pc estimator.
SciPy and NumPy are required.
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
from scipy.sparse import coo_matrix, csr_matrix
from scipy.sparse.linalg import eigs


@dataclass(frozen=True)
class State:
    labels: tuple[int, ...]
    gains: tuple[int, ...]


class DSU:
    """Potential union-find with integer lifted-x differences."""

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
    return State((-1,) * width, (0,) * width)


def step(state: State, mask: int, *, matching: bool) -> State | None:
    """Advance one row; return None if horizontal winding is created."""
    width = len(state.labels)
    dsu = DSU(2 * width)
    representatives: dict[int, int] = {}

    for i, label in enumerate(state.labels):
        if label < 0:
            continue
        if label in representatives:
            dsu.join(representatives[label], i, state.gains[i])
        else:
            representatives[label] = i

    occupied_new = [i for i in range(width) if (mask >> i) & 1]

    for i in occupied_new:
        j = (i + 1) % width
        if (mask >> j) & 1:
            dsu.join(width + i, width + j, (i + 1) // width)

    displacements = (-1, 0, 1) if matching else (0,)
    for i in occupied_new:
        for dx in displacements:
            j = (i + dx) % width
            if state.labels[j] >= 0:
                dsu.join(width + i, j, (i + dx) // width)

    if dsu.bad:
        return None

    labels = [-1] * width
    gains = [0] * width
    canonical: dict[int, tuple[int, int]] = {}
    next_label = 0
    for i in occupied_new:
        root, gain = dsu.find(width + i)
        if root not in canonical:
            canonical[root] = (next_label, gain)
            next_label += 1
        label, base_gain = canonical[root]
        labels[i] = label
        gains[i] = gain - base_gain

    return State(tuple(labels), tuple(gains))


def build_automaton(width: int, *, matching: bool):
    if width < 2:
        raise ValueError("width must be >=2")
    start = empty_state(width)
    states = [start]
    index = {start: 0}
    queue = deque([start])
    transitions: list[list[int]] = []

    while queue:
        state = queue.popleft()
        row: list[int] = []
        for mask in range(1 << width):
            nxt = step(state, mask, matching=matching)
            if nxt is None:
                row.append(-1)
                continue
            if nxt not in index:
                index[nxt] = len(states)
                states.append(nxt)
                queue.append(nxt)
            row.append(index[nxt])
        transitions.append(row)

    if len(transitions) != len(states):
        raise AssertionError("incomplete BFS")
    return states, transitions


def central_trinomial(width: int) -> int:
    return sum(
        math.comb(width, 2 * r) * math.comb(2 * r, r)
        for r in range(width // 2 + 1)
    )


def aggregate_transitions(transitions: list[list[int]], width: int):
    counter: Counter[tuple[int, int, int]] = Counter()
    for source, row in enumerate(transitions):
        for mask, destination in enumerate(row):
            if destination >= 0:
                counter[(source, destination, mask.bit_count())] += 1
    keys = list(counter)
    counts = np.array([counter[key] for key in keys], dtype=float)
    src = np.array([key[0] for key in keys], dtype=np.int32)
    dst = np.array([key[1] for key in keys], dtype=np.int32)
    occ = np.array([key[2] for key in keys], dtype=np.int16)
    return src, dst, occ, counts


class SafeTransfer:
    def __init__(self, width: int, *, matching: bool):
        self.width = width
        self.matching = matching
        self.states, transitions = build_automaton(width, matching=matching)
        self.src, self.dst, self.occ, self.counts = aggregate_transitions(
            transitions, width
        )

    def matrix(self, p: float) -> csr_matrix:
        q = 1.0 - p
        data = self.counts * (p ** self.occ) * (q ** (self.width - self.occ))
        return coo_matrix(
            (data, (self.src, self.dst)),
            shape=(len(self.states), len(self.states)),
        ).tocsr()

    def derivative_matrix(self, p: float) -> csr_matrix:
        q = 1.0 - p
        weight = self.counts * (p ** self.occ) * (q ** (self.width - self.occ))
        factor = self.occ / p - (self.width - self.occ) / q
        return coo_matrix(
            (weight * factor, (self.src, self.dst)),
            shape=(len(self.states), len(self.states)),
        ).tocsr()

    @staticmethod
    def _dominant_pair(matrix: csr_matrix):
        n = matrix.shape[0]
        if n <= 4:
            dense = matrix.toarray()
            vals, vecs = np.linalg.eig(dense)
            k = int(np.argmax(vals.real))
            lam = float(vals[k].real)
            right = np.asarray(vecs[:, k].real)
            vals_l, vecs_l = np.linalg.eig(dense.T)
            j = int(np.argmin(np.abs(vals_l - vals[k])))
            left = np.asarray(vecs_l[:, j].real)
        else:
            val, vec = eigs(matrix, k=1, which="LM", tol=1e-12, maxiter=200000)
            lam = float(val[0].real)
            right = np.asarray(vec[:, 0].real)
            _, vec_l = eigs(
                matrix.T, k=1, which="LM", tol=1e-12, maxiter=200000
            )
            left = np.asarray(vec_l[:, 0].real)

        if right.sum() < 0:
            right = -right
        if left.sum() < 0:
            left = -left
        inner = float(left @ right)
        if inner <= 0:
            raise ArithmeticError("Perron left/right normalization failed")
        left /= inner
        return lam, left, right

    def perron(self, p: float) -> dict[str, float]:
        matrix = self.matrix(p)
        derivative = self.derivative_matrix(p)
        lam, left, right = self._dominant_pair(matrix)
        dlam = float(left @ (derivative @ right))
        I0 = -math.log(lam)
        I0_prime = -dlam / lam
        mean_occupied = self.width * p - p * (1 - p) * I0_prime
        residual = float(
            np.linalg.norm(matrix @ right - lam * right, ord=np.inf)
            / max(1.0, np.linalg.norm(right, ord=np.inf))
        )
        return {
            "lambda": lam,
            "I0": I0,
            "I0_prime": I0_prime,
            "mean_occupied_safe_row": mean_occupied,
            "right_residual_inf": residual,
        }


def coexistence_record(
    width: int,
    *,
    reference_pc: float,
    root_xtol: float = 2e-14,
) -> dict[str, object]:
    g4 = SafeTransfer(width, matching=False)
    g8 = SafeTransfer(width, matching=True)

    expected = central_trinomial(width)
    if len(g4.states) != expected or len(g8.states) != expected:
        raise AssertionError(
            f"safe-state count mismatch at w={width}: "
            f"{len(g4.states)}, {len(g8.states)}, expected {expected}"
        )
    same_state_set = set(g4.states) == set(g8.states)

    def log_ratio(p: float) -> float:
        a = g4.perron(p)["lambda"]
        b = g8.perron(1.0 - p)["lambda"]
        return math.log(a) - math.log(b)

    root = brentq(log_ratio, 0.5000000001, 0.999999999, xtol=root_xtol)
    p4 = g4.perron(root)
    p8 = g8.perron(1.0 - root)
    theta_slope = p4["I0_prime"] + p8["I0_prime"]

    pc4 = g4.perron(reference_pc)
    pc8 = g8.perron(1.0 - reference_pc)
    theta_pc = pc4["I0"] - pc8["I0"]

    return {
        "width": width,
        "safe_state_count_G4": len(g4.states),
        "safe_state_count_G8": len(g8.states),
        "central_trinomial_count": expected,
        "G4_G8_safe_state_sets_identical": same_state_set,
        "charge_coexistence_root": root,
        "log_perron_ratio_at_root": log_ratio(root),
        "lambda_root_G4": p4["lambda"],
        "lambda_root_G8_complement": p8["lambda"],
        "theta_slope_at_root": theta_slope,
        "theta_slope_times_w_quarter": theta_slope * width ** 0.25,
        "mean_occupied_G4_safe_row_at_root": p4["mean_occupied_safe_row"],
        "mean_occupied_G8_safe_row_at_complement_root": p8[
            "mean_occupied_safe_row"
        ],
        "theta_at_reference_pc": theta_pc,
        "theta_at_reference_pc_times_w_17_over_4": theta_pc
        * width ** (17.0 / 4.0),
        "reference_pc_minus_root": reference_pc - root,
        "root_shift_times_w4": (reference_pc - root) * width**4,
        "perron_residuals": {
            "G4_root": p4["right_residual_inf"],
            "G8_root": p8["right_residual_inf"],
            "G4_reference": pc4["right_residual_inf"],
            "G8_reference": pc8["right_residual_inf"],
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--min-width", type=int, default=2)
    parser.add_argument("--max-width", type=int, default=8)
    parser.add_argument(
        "--reference-pc",
        type=float,
        default=0.59274605079,
        help="diagnostic reference only; not used to locate charge roots",
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    records = [
        coexistence_record(width, reference_pc=args.reference_pc)
        for width in range(args.min_width, args.max_width + 1)
    ]
    result = {
        "schema": "fixed-width-charge-transfer-v1",
        "model": "independent square-site NN / complementary matching site",
        "claim_boundary": [
            "Charge roots are located only from equality of safe Perron roots.",
            "reference_pc is used only for finite-size scaling diagnostics.",
            "This small-width implementation is not a new production pc estimator.",
        ],
        "reference_pc": args.reference_pc,
        "records": records,
    }
    text = json.dumps(result, indent=2, sort_keys=True)
    if args.output:
        args.output.write_text(text + "\n", encoding="utf-8")
    print(text)


if __name__ == "__main__":
    main()
