"""Exact controls for the axial RSW-ring lemma, not a simulation or proof of RSW.

Run: python scripts/axial_ring_gluing.py [--output NEW_JSON]
Only Python's standard library is required. Existing result files are not replaced.
"""
from __future__ import annotations

import argparse
import itertools
import json
from fractions import Fraction
from pathlib import Path


def cell_boundaries(width: int, scale: int) -> list[int]:
    """Balanced integer cells, each of length in [scale, 2*scale]."""
    if scale < 1 or width < 4 * scale:
        raise ValueError("Require scale >= 1 and width >= 4*scale")
    count = width // scale
    base, extra = divmod(width, count)
    lengths = [base + (i < extra) for i in range(count)]
    bounds = [0]
    for length in lengths:
        bounds.append(bounds[-1] + length)
    return bounds


def crossing(mask: int, width: int, height: int,
             left: int, right: int, vertical: bool) -> bool:
    """Free rectangle in the integer lift; only occupancy is read modulo width."""
    if not 0 < right - left < width:
        raise ValueError("Rectangle must inject into the horizontal cylinder")

    def occupied(x: int, y: int) -> bool:
        return bool(mask & (1 << (y * width + x % width)))

    starts = ([(x, 0) for x in range(left, right + 1)] if vertical
              else [(left, y) for y in range(height + 1)])
    seen = {v for v in starts if occupied(*v)}
    stack = list(seen)
    while stack:
        x, y = stack.pop()
        if (vertical and y == height) or (not vertical and x == right):
            return True
        for nx, ny in ((x+1, y), (x-1, y), (x, y+1), (x, y-1)):
            if (left <= nx <= right and 0 <= ny <= height
                    and (nx, ny) not in seen and occupied(nx, ny)):
                seen.add((nx, ny))
                stack.append((nx, ny))
    return False


def ring_events(mask: int, width: int, scale: int) -> list[bool]:
    bounds = cell_boundaries(width, scale)
    count = len(bounds) - 1
    extended = bounds + [width + bounds[1]]
    vertical = [crossing(mask, width, scale, bounds[i], bounds[i+1], True)
                for i in range(count)]
    horizontal = [crossing(mask, width, scale, extended[i], extended[i+2], False)
                  for i in range(count)]
    return vertical + horizontal


def horizontal_winding(mask: int, width: int, height: int) -> bool:
    """Independent graph-potential test; no crossing/gluing routine is called."""
    potential: dict[int, int] = {}
    for vertex in range(width * (height + 1)):
        if not (mask >> vertex) & 1 or vertex in potential:
            continue
        potential[vertex] = vertex % width
        stack = [vertex]
        while stack:
            v = stack.pop()
            x, y = v % width, v // width
            neighbors = [(y*width+(x+1)%width, 1),
                         (y*width+(x-1)%width, -1)]
            if y: neighbors.append((v-width, 0))
            if y < height: neighbors.append((v+width, 0))
            for u, dx in neighbors:
                if not (mask >> u) & 1:
                    continue
                expected = potential[v] + dx
                if u in potential:
                    if potential[u] != expected:
                        assert (expected-potential[u]) % width == 0
                        return True
                else:
                    potential[u] = expected
                    stack.append(u)
    return False


def probability(counts: list[int], p: Fraction) -> Fraction:
    n = len(counts) - 1
    return sum((Fraction(c)*p**k*(1-p)**(n-k)
                for k, c in enumerate(counts)), Fraction(0))


def exact_census(width: int, scale: int) -> dict:
    n = width * (scale + 1)
    if n > 16:
        raise ValueError("This control deliberately caps exhaustive enumeration at 16 sites")
    good_counts = [0] * (n+1)
    winding_counts = [0] * (n+1)
    event_count = 2 * (width // scale)
    marginals = [[0]*(n+1) for _ in range(event_count)]
    for mask in range(1 << n):
        events = ring_events(mask, width, scale)
        good = all(events)
        wraps = horizontal_winding(mask, width, scale)
        if good and not wraps:
            raise AssertionError(f"Gluing counterexample: w={width}, s={scale}, mask={mask}")
        k = mask.bit_count()
        good_counts[k] += good
        winding_counts[k] += wraps
        for row, event in zip(marginals, events):
            row[k] += event
    checks = []
    for p in (Fraction(1, 3), Fraction(1, 2), Fraction(2, 3)):
        joint = probability(good_counts, p)
        product = Fraction(1)
        for row in marginals:
            product *= probability(row, p)
        wrap = probability(winding_counts, p)
        assert product <= joint <= wrap
        checks.append({"p": str(p), "product_of_marginals": str(product),
                       "joint_ring_probability": str(joint),
                       "winding_probability": str(wrap)})
    return {"width": width, "scale": scale, "configurations": 1 << n,
            "good_counts_by_occupation": good_counts,
            "winding_counts_by_occupation": winding_counts,
            "implication_failures": 0, "exact_fkg_checks": checks}


def uneven_cell_controls() -> dict:
    """Nonexhaustive deterministic controls, including nondivisible circumferences."""
    checked = good = 0
    for width, scale in ((9, 2), (11, 2), (13, 3)):
        boundaries = cell_boundaries(width, scale)
        assert boundaries[-1] == width
        assert all(scale <= b-a <= 2*scale for a,b in zip(boundaries, boundaries[1:]))
        n = width*(scale+1)
        full = (1 << n)-1
        masks = {0, full}
        for holes in itertools.chain(itertools.combinations(range(n), 1),
                                     itertools.combinations(range(n), 2)):
            masks.add(full ^ sum(1 << v for v in holes))
        # Add striped and checkerboard masks, not samples from any probability law.
        masks.update(sum(1 << (y*width+x) for x in range(width) for y in range(scale+1)
                         if (x+2*y+shift) % modulus)
                     for modulus in (2, 3, 4, 5) for shift in range(modulus))
        for mask in sorted(masks):
            yes = all(ring_events(mask, width, scale))
            assert not yes or horizontal_winding(mask, width, scale)
            good += yes
            checked += 1
    return {"configurations": checked, "ring_positive": good,
            "implication_failures": 0, "exhaustive": False,
            "role": "deterministic geometry controls, not probability estimates"}


def run() -> dict:
    censuses = [exact_census(w, 1) for w in (4, 5, 6, 7)]
    return {"schema": "matching-one.axial-ring-controls.v1", "date": "2026-09-13",
            "model": "NN independent square sites; horizontal cylinder; vertical free boundary",
            "proof_note": "notes/axial-quantile-geometry-criterion-20260913.md",
            "censuses": censuses, "uneven_cells": uneven_cell_controls(),
            "total_exhaustive_configurations": sum(c["configurations"] for c in censuses),
            "scope": "Checks the finite ring implication and FKG arithmetic. Does not verify RSW, all-size necessity, or an oblique-torus extension.",
            "full_repository_ci_run": False}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = json.dumps(run(), indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with args.output.open("x", encoding="utf-8") as stream:
            stream.write(result)
        print(f"Written {args.output}")
    else:
        print(result, end="")


if __name__ == "__main__":
    main()
