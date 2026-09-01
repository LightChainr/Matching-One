#!/usr/bin/env python3
"""Exact three-fibre certificate for the P537 ordinary-four-arm minor."""

from __future__ import annotations

import importlib.util
import json
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
REF = ROOT / "experiments/p337-thermal-gate-audit-20260901/thermal_gate.py"
SPEC = importlib.util.spec_from_file_location("p537_thermal_gate", REF)
assert SPEC and SPEC.loader
MOD = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MOD)


def mask_from_coordinates(torus, occupied):
    index = {(x % 4, y % 4): i for i, (x, y) in enumerate(torus.reps)}
    return sum(1 << index[(x % 4, y % 4)] for x, y in occupied)


def site_from_coordinate(torus, coordinate):
    target = (coordinate[0] % 4, coordinate[1] % 4)
    return next(i for i, p in enumerate(torus.reps) if (p[0] % 4, p[1] % 4) == target)


FIBRES = {
    "entry": {
        "occupied": ((0, 0), (1, 0), (0, 1), (2, 0), (1, 1), (0, 2)),
        "z": (3, 0),
        "pair": ((2, 1), (1, 2)),
        "expected": (0, 1, 0, 0),
    },
    "kernel": {
        "occupied": ((0, 0), (1, 0), (0, 1), (2, 0), (1, 1), (2, 2)),
        "z": (2, 3),
        "pair": ((2, 1), (1, 2)),
        "expected": (0, 0, 4, 0),
    },
    "completion": {
        "occupied": ((0, 0), (1, 0), (0, 1), (2, 0), (0, 2), (3, 0)),
        "z": (0, 3),
        "pair": ((1, 1), (2, 1)),
        "expected": (1, 2, 4, 4),
    },
}


def main():
    torus = MOD.Torus(4, 0)
    rows = {}
    for name, spec in FIBRES.items():
        mask = mask_from_coordinates(torus, spec["occupied"])
        z = site_from_coordinate(torus, spec["z"])
        x, y = (site_from_coordinate(torus, p) for p in spec["pair"])
        assert not (mask >> z) & 1 and not (mask >> x) & 1 and not (mask >> y) & 1
        incident = tuple((mask >> u) & 1 for u in torus.nb[z])
        assert incident in ((1, 0, 1, 0), (0, 1, 0, 1))
        rank0, rank1 = torus.rank(mask), torus.rank(mask | (1 << z))
        pair_index = torus.pairs.index(tuple(sorted((x, y))))
        g0 = torus.pair_values(mask)[pair_index]
        g1 = torus.pair_values(mask | (1 << z))[pair_index]
        assert (rank0, rank1, g0, g1) == spec["expected"]
        rows[name] = {
            "occupied": spec["occupied"],
            "z": spec["z"],
            "pair": spec["pair"],
            "incident_NESW": incident,
            "rank_before_after": [rank0, rank1],
            "g16_before_after": [g0, g1],
            "delta_g": str(Fraction(g1 - g0, 16)),
        }

    # Columns are kernel-entry and completion-entry.  The source row is
    # independent of the Schur parameters; the readout row contains chi,R.
    # det([[-1/4,0],[chi+R,2chi]]) = -chi/2.
    determinant_coefficient_of_chi = Fraction(-1, 2)
    print(json.dumps({
        "geometry": [4, 0],
        "N": 16,
        "background_occupation_each": 6,
        "fibres": rows,
        "difference_matrix": [
            ["-1/4", "0"],
            ["chi+R", "2*chi"],
        ],
        "determinant": "-chi/2",
        "determinant_coefficient_of_chi": str(determinant_coefficient_of_chi),
    }, indent=2))


if __name__ == "__main__":
    main()
