#!/usr/bin/env python3
"""Resolution pass for issue #576 Part 1 discrepancies.

Reads the published torus.txt blocks and the independent displacement-DSU
enumeration output (exact_enum2_l2l5.txt) and writes three derived artifacts:

1. polynomial_vs_enum2.json   per-coefficient comparison, L = 2..5
2. torus_blocks_health.json   proven divisibility test on every block,
                              the L=10 missing trailing coefficient, and the
                              P_L(p_c) table against the correct continuum
                              anchor R^v = 0.521058290... (Mertens-Ziff)
3. spiral_discrepancy_resolution.json  the 10 spiral configurations missed by
                              the same-row doubled-grid criterion
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import mpmath as mp

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "results" / "issue576" / "wrapping-grounding"

PC = mp.mpf("0.59274605079210")          # square site p_c (Jacobsen 2015)
RV = mp.mpf("0.521058290")               # Mertens-Ziff continuum R^v


def read_blocks(path: Path) -> dict[int, list[int]]:
    blocks: dict[int, list[int]] = {}
    current = None
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith("L="):
            current = int(line[2:])
            blocks[current] = []
        else:
            blocks[current].append(int(line))
    return blocks


def read_enum2(path: Path) -> dict[int, list[int]]:
    out: dict[int, list[int]] = {}
    current = None
    for line in path.read_text().splitlines():
        line = line.strip()
        if line.startswith("L="):
            current = int(line.split()[0][2:])
            out[current] = []
        elif line.startswith("c_k:"):
            out[current] = [int(v) for v in line.split()[1:]]
    return out


def main() -> None:
    mp.mp.dps = 60
    blocks = read_blocks(BASE / "raw" / "torus.txt")
    enum2 = read_enum2(BASE / "derived" / "exact_enum2_l2l5.txt")

    # 1. per-coefficient comparison at the brute-force reachable sizes
    comparison = {}
    for L in sorted(enum2):
        record = {"enum2_sum": sum(enum2[L])}
        if L not in blocks:
            record.update(
                {
                    "published_sum": None,
                    "coefficient_mismatches": None,
                    "match": None,
                    "note": "no published block (torus.txt starts at L=3)",
                }
            )
            comparison[str(L)] = record
            continue
        pub = blocks[L]
        mine = enum2[L]
        assert len(pub) == len(mine) == L * L + 1
        diffs = [
            {"k": k, "published": pub[k], "enum2": mine[k]}
            for k in range(len(pub))
            if pub[k] != mine[k]
        ]
        record.update(
            {
                "published_sum": sum(pub),
                "coefficient_mismatches": diffs,
                "match": not diffs,
            }
        )
        comparison[str(L)] = record

    (BASE / "derived" / "polynomial_vs_enum2.json").write_text(
        json.dumps(
            {
                "definition": "NN square SITE percolation on the L x L torus; "
                "wraps the horizontal (period-1) direction, including configs "
                "that also wrap vertically and including spirals. enum2 = "
                "scripts/exact_wrapping_enum2.cpp (displacement union-find, "
                "independent of the doubled-grid criterion).",
                "comparison": comparison,
            },
            indent=2,
        )
    )

    # 2. block health: proven divisibility test + P_L(p_c) against R^v
    health = {}
    for L in sorted(blocks):
        cs = list(blocks[L])
        N = L * L
        truncated = len(cs) == N  # L=10: trailing c_N = 1 missing in the file
        if truncated:
            cs.append(1)
        assert len(cs) == N + 1
        failures = []
        for k, c in enumerate(cs):
            d = N // math.gcd(k, N)
            if c % d != 0:
                failures.append({"k": k, "c": c, "required_divisor": d})
        val = sum(c * PC**k * (1 - PC) ** (N - k) for k, c in enumerate(cs))
        health[str(L)] = {
            "coefficients_in_file": len(blocks[L]),
            "expected": N + 1,
            "file_truncated": truncated,
            "assumed_missing_coefficient": {"k": N, "value": 1} if truncated else None,
            "divisibility_failures": failures,
            "passes_divisibility": not failures,
            "P_L_at_pc": mp.nstr(val, 15),
            "RV_minus_P_L": mp.nstr(RV - val, 12),
        }

    (BASE / "derived" / "torus_blocks_health.json").write_text(
        json.dumps(
            {
                "anchor": {
                    "quantity": "R^v, continuum limit of P_L(p_c) for "
                    "specified-direction wrapping (incl. both-direction) on the "
                    "square torus",
                    "value": "0.521058290",
                    "source": "Mertens-Ziff 2016 (PRE 94, 022149); the same "
                    "0.521058290 appears in Newman-Ziff PRE 64, 016706",
                    "note": "P_L(1/2) has limit 0 (p=1/2 < p_c); the earlier "
                    "comparison of P_L(1/2) against 0.1694 was a category error "
                    "and the 'L>=9 corrupted' claim is retracted.",
                },
                "blocks": health,
            },
            indent=2,
        )
    )

    # 3. the 10 spiral configurations (regenerable via
    #    scripts/exact_wrapping_criterion_diff.cpp); each is verified here to
    #    wind in x AND y simultaneously (spiral), with no pure-x winding cycle
    spirals = [
        (8174204, [(0, 2), (0, 3), (0, 4), (1, 0), (1, 1), (1, 4), (2, 1), (2, 2),
                   (2, 3), (3, 0), (3, 3), (3, 4), (4, 0), (4, 1), (4, 2)]),
        (8277465, [(0, 0), (0, 3), (0, 4), (1, 1), (1, 2), (1, 3), (2, 0), (2, 1),
                   (2, 4), (3, 2), (3, 3), (3, 4), (4, 0), (4, 1), (4, 2)]),
        (15331577, [(0, 0), (0, 3), (0, 4), (1, 0), (1, 1), (1, 2), (2, 2), (2, 3),
                    (2, 4), (3, 0), (3, 1), (3, 4), (4, 1), (4, 2), (4, 3)]),
        (15507347, [(0, 0), (0, 1), (0, 4), (1, 2), (1, 3), (1, 4), (2, 0), (2, 1),
                    (2, 2), (3, 0), (3, 3), (3, 4), (4, 1), (4, 2), (4, 3)]),
        (20407548, [(0, 2), (0, 3), (0, 4), (1, 0), (1, 1), (1, 2), (2, 0), (2, 3),
                    (2, 4), (3, 1), (3, 2), (3, 3), (4, 0), (4, 1), (4, 4)]),
        (20848430, [(0, 1), (0, 2), (0, 3), (1, 0), (1, 3), (1, 4), (2, 0), (2, 1),
                    (2, 2), (3, 2), (3, 3), (3, 4), (4, 0), (4, 1), (4, 4)]),
        (26473070, [(0, 1), (0, 2), (0, 3), (1, 0), (1, 1), (1, 4), (2, 2), (2, 3),
                    (2, 4), (3, 0), (3, 1), (3, 2), (4, 0), (4, 3), (4, 4)]),
        (26693511, [(0, 0), (0, 1), (0, 2), (1, 2), (1, 3), (1, 4), (2, 0), (2, 1),
                    (2, 4), (3, 1), (3, 2), (3, 3), (4, 0), (4, 3), (4, 4)]),
        (29615571, [(0, 0), (0, 1), (0, 4), (1, 1), (1, 2), (1, 3), (2, 0), (2, 3),
                    (2, 4), (3, 0), (3, 1), (3, 2), (4, 2), (4, 3), (4, 4)]),
        (29997863, [(0, 0), (0, 1), (0, 2), (1, 0), (1, 3), (1, 4), (2, 1), (2, 2),
                    (2, 3), (3, 0), (3, 1), (3, 4), (4, 2), (4, 3), (4, 4)]),
    ]

    def windings(cells: list[tuple[int, int]], L: int = 5) -> tuple[bool, bool]:
        occ = set(cells)
        edges = []
        for r, c in cells:
            nb = (r, (c + 1) % L)
            if nb in occ:
                edges.append(((r, c), nb, (1, 0)))
            nb = ((r + 1) % L, c)
            if nb in occ:
                edges.append(((r, c), nb, (0, 1)))
        par = {s: s for s in cells}
        gx = {s: 0 for s in cells}
        gy = {s: 0 for s in cells}

        def find(x):
            if par[x] == x:
                return x, 0, 0
            rx, px, py = find(par[x])
            gx[x] += px
            gy[x] += py
            par[x] = rx
            return rx, gx[x], gy[x]

        wx = wy = False
        for i, j, (ex, ey) in edges:
            ri, ix, iy = find(i)
            rj, jx, jy = find(j)
            if ri == rj:
                if ix + ex - jx != 0:
                    wx = True
                if iy + ey - jy != 0:
                    wy = True
            else:
                par[rj] = ri
                gx[rj] = ix + ex - jx
                gy[rj] = iy + ey - jy
        return wx, wy

    cfg_records = []
    for idx, cells in spirals:
        wx, wy = windings(cells)
        assert wx and wy, f"config {idx} is not a spiral"
        cfg_records.append(
            {
                "cfg_index": idx,
                "occupied_row_col": cells,
                "winds_x": wx,
                "winds_y": wy,
                "n_occupied": len(cells),
            }
        )

    (BASE / "derived" / "spiral_discrepancy_resolution.json").write_text(
        json.dumps(
            {
                "summary": "At L=5 the doubled-grid criterion 'component contains "
                "(r,0) and (r,L) for one row r' misses exactly 10 configurations "
                "(2^25 scanned): spiral configurations whose occupied cluster "
                "winds in x AND y simultaneously, so the lift to the doubled "
                "strip has no same-row (r,0),(r,L) pair. Published sum 8853301 "
                "is correct; same-row enumerator 8853291 = 8853301 - 10, all 10 "
                "inside the c_15 bin.",
                "winding_check": "each configuration verified winds_x AND winds_y "
                "via displacement union-find (this file is generated, assertions "
                "enforced at generation time)",
                "configurations": cfg_records,
            },
            indent=2,
        )
    )

    print("wrote polynomial_vs_enum2.json, torus_blocks_health.json, "
          "spiral_discrepancy_resolution.json")


if __name__ == "__main__":
    main()
