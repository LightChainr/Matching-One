#!/usr/bin/env python3
"""#674 tests: spiral cell identification and involution actions, axis L=2,3 (seconds).

Run directly:  python3 tests/test_issue674_spiral_involution.py
(no pytest required)
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from issue674_spiral_involution import build_maps, classify, in_cell  # noqa: E402
from matched_torus_reference import axis_geometry  # noqa: E402
from probe635_sector_decomposition import sector_label, wrap_homology  # noqa: E402


def test_L3_cell_size_is_6() -> None:
    g = axis_geometry(3)
    cell = [m for m in range(1 << 9) if in_cell(m, g)]
    assert len(cell) == 6, cell
    assert all(bin(m).count("1") == 6 for m in cell)


def test_L3_cell_configwise_cancellation() -> None:
    """Every question-cell config has either=1 on both sides => D(C)=0 config-wise."""
    g = axis_geometry(3)
    for m in range(1 << 9):
        if not in_cell(m, g):
            continue
        rb, chb, lb = classify(m, g, matching=False)
        rw, chw, lw = classify((~m) & ((1 << 9) - 1), g, matching=True)
        assert (lb, lw, rb, rw) == ("both-same", "both-same", 1, 1)
        assert chb.either and chw.either


def test_L3_fixed_point_free_reflection_pairings() -> None:
    """Axis reflections are fixed-point-free k-preserving involutions on the cell."""
    g = axis_geometry(3)
    cell = [m for m in range(1 << 9) if in_cell(m, g)]
    idx = {m: j for j, m in enumerate(cell)}
    maps = build_maps(3)
    for name in ("refl-x(0)", "refl-x(1)", "refl-x(2)", "refl-y(0)", "refl-y(1)", "refl-y(2)"):
        perm = maps[name]
        image = [apply(perm, m, g) for m in cell]
        assert all(x in idx for x in image), name
        image_idx = [idx[x] for x in image]
        assert all(image_idx[image_idx[j]] == j for j in range(6)), name
        assert all(a != b for a, b in zip(image_idx, range(6))), name
        assert all(bin(cell[a]).count("1") == bin(cell[b]).count("1")
                   for a, b in zip(image_idx, range(6))), name


def apply(perm, mask, geometry):
    out = 0
    for i in range(geometry.n):
        if (mask >> i) & 1:
            out |= 1 << perm[i]
    return out


def test_colour_flip_never_preserves_rank1_cell() -> None:
    """Complement of a rank-1 both-same x both-same spiral is never in the cell (L=2,3)."""
    for L in (2, 3):
        g = axis_geometry(L)
        n = g.n
        full = (1 << n) - 1
        cell = {m for m in range(1 << n) if in_cell(m, g)}
        for m in cell:
            assert (~m) & full not in cell, (L, m)


def test_both_two_is_empty_L23() -> None:
    """Structural both-two = 0: x-wrap and y-wrap forces both-same, L=2,3, both lattices."""
    for L in (2, 3):
        g = axis_geometry(L)
        n = g.n
        for edges in (g.primal_edges, g.matching_edges):
            for mask in range(1 << n):
                black = [bool((mask >> i) & 1) for i in range(n)]
                b = wrap_homology(black, edges, n)
                if b[0] and b[1]:
                    assert sector_label(*b) == "both-same", (L, mask)


def main() -> None:
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
        print(f"PASS {t.__name__}")
    print(f"{len(tests)} tests, all pass")


if __name__ == "__main__":
    main()
