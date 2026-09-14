#!/usr/bin/env python3
"""Momentum-resolved first descendant of the safe charge transfer.

At each fixed width locate the charge-coexistence root, diagonalize the safe
G4 and complementary-G8 transfer matrices, and inspect the first excited
(twofold-degenerate) eigenspace under one-column rotation.

The expected magnetic level-one signature is

    w log(lambda0/|lambda1|) -> 2*pi,
    rotation eigenvalues -> exp(+/-2*pi*i/w).

This is deterministic transfer spectroscopy, not Monte Carlo.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
from scipy.sparse.linalg import eigs

from fixed_width_charge_transfer import SafeTransfer, State


def rotate_state(state: State) -> State:
    """Translate the frontier one physical column to the right."""
    width = len(state.labels)
    entries = []
    for column, label in enumerate(state.labels):
        if label < 0:
            continue
        new_column = (column + 1) % width
        deck = (column + 1) // width
        entries.append((new_column, label, state.gains[column] + deck))
    entries.sort()

    labels = [-1] * width
    gains = [0] * width
    new_label = {}
    base_gain = {}
    next_label = 0
    for column, label, gain in entries:
        if label not in new_label:
            new_label[label] = next_label
            base_gain[label] = gain
            next_label += 1
        labels[column] = new_label[label]
        gains[column] = gain - base_gain[label]
    return State(tuple(labels), tuple(gains))


def leading_lambda(transfer: SafeTransfer, p: float) -> float:
    matrix = transfer.matrix(p)
    if matrix.shape[0] <= 4:
        return float(max(np.linalg.eigvals(matrix.toarray()).real))
    value = eigs(
        matrix,
        k=1,
        which="LM",
        tol=1e-12,
        maxiter=200000,
        return_eigenvectors=False,
    )[0]
    return float(value.real)


def charge_root(g4: SafeTransfer, g8: SafeTransfer) -> float:
    def equation(p: float) -> float:
        return math.log(leading_lambda(g4, p)) - math.log(
            leading_lambda(g8, 1.0 - p)
        )

    return brentq(equation, 0.5000000001, 0.999999999, xtol=2e-14)


def first_descendant(transfer: SafeTransfer, p: float) -> dict[str, object]:
    matrix = transfer.matrix(p)
    values, vectors = eigs(matrix, k=6, which="LM", tol=1e-11, maxiter=200000)
    order = np.argsort(-np.abs(values))
    values = values[order]
    vectors = vectors[:, order]

    lam0 = values[0]
    lam1 = values[1]
    pair = [i for i in range(1, len(values)) if abs(values[i] - lam1) < 1e-7]
    if len(pair) < 2:
        raise ArithmeticError("first excited eigenvalue is not resolved as a pair")
    pair = pair[:2]
    basis = vectors[:, pair]

    index = {state: i for i, state in enumerate(transfer.states)}
    rotation = np.array([index[rotate_state(state)] for state in transfer.states])

    # Verify exact commutation of the state permutation with the weighted kernel.
    dense_test = matrix[:, rotation][rotation, :] - matrix
    if dense_test.nnz and np.max(np.abs(dense_test.data)) > 1e-13:
        raise ArithmeticError("one-column translation does not commute with transfer")

    rotated_basis = np.empty_like(basis)
    for column in range(basis.shape[1]):
        vector = np.empty(basis.shape[0], dtype=complex)
        vector[rotation] = basis[:, column]
        rotated_basis[:, column] = vector

    representation = np.linalg.lstsq(basis, rotated_basis, rcond=None)[0]
    phases = np.linalg.eigvals(representation)
    phases = sorted(phases, key=lambda z: np.angle(z))

    expected = [
        np.exp(-2j * np.pi / transfer.width),
        np.exp(+2j * np.pi / transfer.width),
    ]
    phase_error = max(abs(a - b) for a, b in zip(phases, expected))

    gap = float(transfer.width * math.log(abs(lam0 / lam1)))
    return {
        "lambda0": float(lam0.real),
        "lambda1": float(lam1.real),
        "first_excited_multiplicity_resolved": 2,
        "scaled_gap_w_log_lambda0_over_lambda1": gap,
        "scaled_gap_over_2pi": gap / (2 * math.pi),
        "rotation_eigenvalues": [
            {"real": float(z.real), "imag": float(z.imag)} for z in phases
        ],
        "expected_momenta": [-1, 1],
        "max_rotation_phase_error": float(phase_error),
    }


def run_width(width: int) -> dict[str, object]:
    g4 = SafeTransfer(width, matching=False)
    g8 = SafeTransfer(width, matching=True)
    root = charge_root(g4, g8)
    return {
        "width": width,
        "charge_root_p": root,
        "safe_states": len(g4.states),
        "G4": first_descendant(g4, root),
        "G8_complement": first_descendant(g8, 1.0 - root),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--min-width", type=int, default=5)
    parser.add_argument("--max-width", type=int, default=9)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    result = {
        "schema": "safe-transfer-momentum-spectrum-v1",
        "claim_boundary": [
            "deterministic small-width transfer spectroscopy",
            "level-one magnetic-descendant interpretation uses the pDTL/CFT sector dictionary",
            "higher raw eigenvalues are not assigned to one Verma tower by this script",
        ],
        "records": [run_width(w) for w in range(args.min_width, args.max_width + 1)],
    }
    text = json.dumps(result, indent=2, sort_keys=True)
    if args.output:
        args.output.write_text(text + "\n", encoding="utf-8")
    print(text)


if __name__ == "__main__":
    main()
