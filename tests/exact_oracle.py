#!/usr/bin/env python3
"""Independent lifted-graph oracle for finite-torus matching identities.

This module must NOT reproduce the C++ displacement DSU. Clusters and
winding are obtained by BFS on the universal cover of the torus: the same
modulo vertex reached at two distinct cell translations is a wrapping
event. Horizontal wrapping is a nonzero x-translation; vertical wrapping
is a nonzero y-translation.

Site index i = y*L + x. Bit i of the configuration integer is 1 if the
site is occupied on the square NN graph G, and 0 if it is occupied on the
matching NN+NNN graph G*.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import deque
from pathlib import Path
from typing import Iterable


# One-sided lattice displacements (undirected graphs). BFS uses both signs.
G_OFFSETS = ((1, 0), (0, 1))
GSTAR_OFFSETS = ((1, 0), (0, 1), (1, 1), (1, -1))


def _both_directions(offsets: Iterable[tuple[int, int]]) -> tuple[tuple[int, int], ...]:
    seen: set[tuple[int, int]] = set()
    out: list[tuple[int, int]] = []
    for dx, dy in offsets:
        for vec in ((dx, dy), (-dx, -dy)):
            if vec not in seen:
                seen.add(vec)
                out.append(vec)
    return tuple(out)


G_NEIGHBORS = _both_directions(G_OFFSETS)
GSTAR_NEIGHBORS = _both_directions(GSTAR_OFFSETS)


def analyze_component_lifted(
    L: int,
    occupied: list[bool],
    start: int,
    neighbors: tuple[tuple[int, int], ...],
    visited: list[bool],
) -> tuple[bool, bool]:
    """BFS one cluster on the covering space. Returns (wrap_H, wrap_V)."""
    n = L * L
    first_tx = [None] * n
    first_ty = [None] * n
    sx = start % L
    sy = start // L
    q: deque[tuple[int, int, int, int]] = deque()
    q.append((sx, sy, 0, 0))
    first_tx[start] = 0
    first_ty[start] = 0
    visited[start] = True
    wrap_h = False
    wrap_v = False

    while q:
        x, y, tx, ty = q.popleft()
        for dx, dy in neighbors:
            nx = x + dx
            ny = y + dy
            ntx = tx
            nty = ty
            # Preserve original lattice displacement; fold into cell translation.
            if nx >= L:
                nx -= L
                ntx += 1
            elif nx < 0:
                nx += L
                ntx -= 1
            if ny >= L:
                ny -= L
                nty += 1
            elif ny < 0:
                ny += L
                nty -= 1
            nid = ny * L + nx
            if not occupied[nid]:
                continue
            if first_tx[nid] is None:
                first_tx[nid] = ntx
                first_ty[nid] = nty
                visited[nid] = True
                q.append((nx, ny, ntx, nty))
            else:
                if ntx != first_tx[nid]:
                    wrap_h = True
                if nty != first_ty[nid]:
                    wrap_v = True
    return wrap_h, wrap_v


def analyze_graph(L: int, occupied: list[bool], neighbors: tuple[tuple[int, int], ...]):
    n = L * L
    visited = [False] * n
    n_clusters = 0
    wrap_h = False
    wrap_v = False
    for start in range(n):
        if not occupied[start] or visited[start]:
            continue
        n_clusters += 1
        h, v = analyze_component_lifted(L, occupied, start, neighbors, visited)
        wrap_h = wrap_h or h
        wrap_v = wrap_v or v
    return n_clusters, wrap_h, wrap_v


def analyze_config(L: int, mask: int) -> dict[str, int]:
    n = L * L
    occ_g = [bool((mask >> i) & 1) for i in range(n)]
    occ_gs = [not v for v in occ_g]
    c_g, h_g, v_g = analyze_graph(L, occ_g, G_NEIGHBORS)
    c_gs, h_gs, v_gs = analyze_graph(L, occ_gs, GSTAR_NEIGHBORS)
    return {
        "mask": mask,
        "k": occ_g.count(True),
        "clusters_G": c_g,
        "clusters_Gstar": c_gs,
        "H_G": int(h_g),
        "V_G": int(v_g),
        "E_G": int(h_g or v_g),
        "B_G": int(h_g and v_g),
        "H_Gstar": int(h_gs),
        "V_Gstar": int(v_gs),
        "E_Gstar": int(h_gs or v_gs),
        "B_Gstar": int(h_gs and v_gs),
    }


def euler_vef0(L: int, mask: int) -> int:
    """V - E + F0 on the occupied NN subgraph of G (black sites)."""
    n = L * L
    occupied = [bool((mask >> i) & 1) for i in range(n)]
    V = sum(occupied)
    E = 0
    for y in range(L):
        for x in range(L):
            i = y * L + x
            if not occupied[i]:
                continue
            for dx, dy in G_OFFSETS:
                nx = (x + dx) % L
                ny = (y + dy) % L
                j = ny * L + nx
                if occupied[j]:
                    E += 1
    F0 = 0
    for y in range(L):
        for x in range(L):
            corners = (
                y * L + x,
                y * L + (x + 1) % L,
                ((y + 1) % L) * L + x,
                ((y + 1) % L) * L + (x + 1) % L,
            )
            if all(occupied[c] for c in corners):
                F0 += 1
    return V - E + F0


def bernstein_to_monomial(A: list[int], n: int) -> list[int]:
    """Expand sum_k A_k p^k (1-p)^{n-k} into monomial coefficients."""
    coeff = [0] * (n + 1)
    for k, ak in enumerate(A):
        if ak == 0:
            continue
        binom = 1
        nk = n - k
        for j in range(0, nk + 1):
            if j > 0:
                binom = binom * (nk - j + 1) // j
            term = ak * binom
            if j & 1:
                term = -term
            coeff[k + j] += term
    return coeff


def matching_chi_coeffs(n: int) -> list[int]:
    """L^2 * (p - 2 p^2 + p^4) as monomial coefficients of degree n."""
    c = [0] * (n + 1)
    c[1] += n
    c[2] += -2 * n
    if n >= 4:
        c[4] += n
    return c


def sub_poly(a: list[int], b: list[int]) -> list[int]:
    m = max(len(a), len(b))
    out = [0] * m
    for i in range(m):
        av = a[i] if i < len(a) else 0
        bv = b[i] if i < len(b) else 0
        out[i] = av - bv
    return out


def first_mismatch(diff: list[int]) -> int | None:
    for i, v in enumerate(diff):
        if v != 0:
            return i
    return None


def enumerate_microcanonical(L: int):
    n = L * L
    total = 1 << n
    count = [0] * (n + 1)
    cG = [0] * (n + 1)
    cGs = [0] * (n + 1)
    hG = [0] * (n + 1)
    vG = [0] * (n + 1)
    eG = [0] * (n + 1)
    bG = [0] * (n + 1)
    hGs = [0] * (n + 1)
    vGs = [0] * (n + 1)
    eGs = [0] * (n + 1)
    bGs = [0] * (n + 1)
    euler_fail = []
    for mask in range(total):
        o = analyze_config(L, mask)
        k = o["k"]
        count[k] += 1
        cG[k] += o["clusters_G"]
        cGs[k] += o["clusters_Gstar"]
        hG[k] += o["H_G"]
        vG[k] += o["V_G"]
        eG[k] += o["E_G"]
        bG[k] += o["B_G"]
        hGs[k] += o["H_Gstar"]
        vGs[k] += o["V_Gstar"]
        eGs[k] += o["E_Gstar"]
        bGs[k] += o["B_Gstar"]
        euler = euler_vef0(L, mask)
        lhs = o["clusters_G"] - o["clusters_Gstar"] - euler
        rhs_h = o["H_G"] - o["H_Gstar"]
        rhs_v = o["V_G"] - o["V_Gstar"]
        rhs_e = o["E_G"] - o["E_Gstar"]
        rhs_b = o["B_G"] - o["B_Gstar"]
        if not (lhs == rhs_h == rhs_v == rhs_e == rhs_b):
            euler_fail.append(
                {
                    "mask": mask,
                    "lhs_cluster_euler": lhs,
                    "H": rhs_h,
                    "V": rhs_v,
                    "E": rhs_e,
                    "B": rhs_b,
                    "obs": o,
                    "euler_V_E_F0": euler,
                }
            )
    return {
        "count": count,
        "cG": cG,
        "cGs": cGs,
        "hG": hG,
        "vG": vG,
        "eG": eG,
        "bG": bG,
        "hGs": hGs,
        "vGs": vGs,
        "eGs": eGs,
        "bGs": bGs,
        "euler_fail": euler_fail,
    }


def identities_from_microcanonical(L: int, mc: dict) -> dict:
    n = L * L
    N_poly = bernstein_to_monomial(mc["cG"], n)
    Nhat_poly = bernstein_to_monomial(mc["cGs"], n)
    chi = matching_chi_coeffs(n)
    M = sub_poly(sub_poly(N_poly, Nhat_poly), chi)

    def wrap_diff(g_key: str, gs_key: str) -> list[int]:
        return sub_poly(bernstein_to_monomial(mc[g_key], n), bernstein_to_monomial(mc[gs_key], n))

    classes = {
        "H": wrap_diff("hG", "hGs"),
        "V": wrap_diff("vG", "vGs"),
        "E": wrap_diff("eG", "eGs"),
        "B": wrap_diff("bG", "bGs"),
    }
    results = {}
    overall = "PASS"
    for name, poly in classes.items():
        diff = sub_poly(M, poly)
        mm = first_mismatch(diff)
        status = "PASS" if mm is None else "FAIL"
        if status == "FAIL":
            overall = "FAIL"
        results[name] = {
            "status": status,
            "first_mismatching_coefficient": mm,
            "difference_coefficients": diff,
            "coefficients": poly,
        }
    max_bits = 0
    for c in M:
        if c != 0:
            max_bits = max(max_bits, c.bit_length())
    return {
        "L": L,
        "N": n,
        "degree": n,
        "coefficients": M,
        "max_coefficient_bit_length": max_bits,
        "identity": overall,
        "wrapping": results,
        "euler_failures": len(mc["euler_fail"]),
    }


def load_cpp_dump(path: Path) -> list[dict[str, int]]:
    rows = []
    with path.open() as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append({k: int(v) for k, v in row.items()})
    return rows


def compare_dump(L: int, dump_path: Path) -> dict:
    rows = load_cpp_dump(dump_path)
    mismatches = []
    expected_total = 1 << (L * L)
    if len(rows) != expected_total:
        return {
            "status": "FAIL",
            "reason": f"dump has {len(rows)} rows, expected {expected_total}",
            "mismatches": [],
        }
    for row in rows:
        mask = row["mask"]
        oracle = analyze_config(L, mask)
        for key in (
            "k",
            "clusters_G",
            "clusters_Gstar",
            "H_G",
            "V_G",
            "E_G",
            "B_G",
            "H_Gstar",
            "V_Gstar",
            "E_Gstar",
            "B_Gstar",
        ):
            if row[key] != oracle[key]:
                mismatches.append({"mask": mask, "field": key, "cpp": row[key], "oracle": oracle[key]})
                break
        if len(mismatches) >= 32:
            break
    return {
        "status": "PASS" if not mismatches else "FAIL",
        "compared": len(rows),
        "mismatches": mismatches,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--L", type=int, nargs="+", default=[2, 3])
    parser.add_argument("--compare", type=Path, nargs="*", default=[])
    parser.add_argument("--json-out", type=Path, default=None)
    parser.add_argument("--self-identities", action="store_true", default=True)
    args = parser.parse_args()

    report: dict = {"oracle": "lifted-graph BFS", "results": {}}
    overall = "PASS"

    for L in args.L:
        if L * L > 18:
            print(f"refusing L={L}: python oracle is for L=2,3 (optionally 4)", file=sys.stderr)
            return 2
        mc = enumerate_microcanonical(L)
        ident = identities_from_microcanonical(L, mc)
        ident["first_euler_counterexample"] = mc["euler_fail"][:3]
        report["results"][str(L)] = ident
        print(
            f"L={L} identity={ident['identity']} euler_failures={ident['euler_failures']} "
            f"M_L={ident['coefficients']}"
        )
        for name, wr in ident["wrapping"].items():
            print(
                f"  wrap_{name}: {wr['status']} first_mismatch={wr['first_mismatching_coefficient']}"
            )
        if ident["identity"] != "PASS" or ident["euler_failures"]:
            overall = "FAIL"

    compare_status = {}
    for path in args.compare:
        name = path.name
        # Infer L from filename L02_configs.csv or similar.
        L = None
        for cand in args.L:
            if f"L{cand:02d}" in path.name or f"l{cand:02d}" in path.name.lower():
                L = cand
                break
        if L is None:
            # try reading first? default from stem
            stem = path.stem
            for cand in range(2, 8):
                if f"{cand:02d}" in stem:
                    L = cand
                    break
        if L is None:
            print(f"cannot infer L from {path}", file=sys.stderr)
            overall = "FAIL"
            continue
        cmp = compare_dump(L, path)
        compare_status[str(path)] = cmp
        print(f"compare {path} L={L}: {cmp['status']} mismatches={len(cmp['mismatches'])}")
        if cmp["status"] != "PASS":
            overall = "FAIL"
            for m in cmp["mismatches"][:5]:
                print(f"  {m}")

    report["compare"] = compare_status
    report["overall"] = overall
    if args.json_out is not None:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(f"ORACLE_OVERALL={overall}")
    return 0 if overall == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
