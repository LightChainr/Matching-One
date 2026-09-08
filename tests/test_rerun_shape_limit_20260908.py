#!/usr/bin/env python3
"""Test for the 2026-09-08 #622 re-run (rerun_shape_limit_20260908.py).

Checks the four machine-checked statements against merged, cited inputs.
Read-only on results: writes nothing.
"""

import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
HERE = ROOT / "scripts" / "probe_invariant_shape"
sys.path.insert(0, str(HERE))

import rerun_shape_limit_20260908 as rr  # noqa: E402


def test_t1_bond_antisymmetry_exact() -> None:
    census = rr.load_census()
    t1 = rr.check_t1(census)
    assert t1["dual_fail"] == 0
    assert t1["pairs_symmetric"] and t1["M_half_is_zero"]
    assert t1["Z_antisymmetry_max_err"] < 1e-12
    assert t1["Z_half_is_half"]


def test_t2_duality_class_does_not_pin_shape() -> None:
    t2 = rr.check_t2()
    # both toy members obey the T1 constraint to integration error
    assert t2["T1_constraint_max_err_A"] < 1e-3
    assert t2["T1_constraint_max_err_B"] < 1e-3
    # yet their shapes differ well above integration error
    assert t2["shape_spread_within_duality_class"] > 0.1


def test_t3_decomposition_consistency() -> None:
    census = rr.load_census()
    n725 = rr.load_n725()
    t3 = rr.check_t3(census, n725)
    # A + S must reproduce Z at u <= 1/2
    for name, z in t3["Z_by_object"].items():
        for i in range(5):
            a = t3["antisymmetric_part_A_u0.1_to_0.5"][name][i]
            s = t3["symmetric_part_S_u0.1_to_0.5"][name][i]
            assert abs((a + s) - z[i]) < 1e-12, (name, i)
    # bond self-dual object has S ≡ 0 exactly (to float noise)
    assert max(abs(x) for x in
               t3["symmetric_part_S_u0.1_to_0.5"]["bond_L3_selfdual"]) < 1e-12
    # spread ordering recorded in the note: A-spread < 0.05
    assert t3["max_spread_A_overall"] < 0.05


def test_n725_source_is_the_corrected_block() -> None:
    n725 = rr.load_n725()
    assert n725["schema"] == (
        "matching-one.probe-invariant-shape.n725-zflow-corrected.v1")
    # the withdrawn object must not be silently used
    blob = json.dumps(n725)
    assert "n725-zflow-corrected" in blob or n725.get("issue") == 633


def test_results_json_reproducible() -> None:
    """The committed rerun JSON (if present) matches a fresh computation."""
    dest = ROOT / "results" / "probe-invariant-shape" / "rerun-20260908.json"
    if not dest.exists():
        return
    committed = json.loads(dest.read_text())
    fresh = {
        "T1": rr.check_t1(rr.load_census()),
        "T2": rr.check_t2(),
    }
    assert committed["T1_duality_symmetry_theorem"]["Z_antisymmetry_max_err"] \
        == fresh["T1"]["Z_antisymmetry_max_err"]
    assert abs(
        committed["T2_duality_does_not_pin_shape"][
            "shape_spread_within_duality_class"]
        - fresh["T2"]["shape_spread_within_duality_class"]) < 1e-9


if __name__ == "__main__":
    test_t1_bond_antisymmetry_exact()
    test_t2_duality_class_does_not_pin_shape()
    test_t3_decomposition_consistency()
    test_n725_source_is_the_corrected_block()
    test_results_json_reproducible()
    print("all rerun-20260908 checks passed")
