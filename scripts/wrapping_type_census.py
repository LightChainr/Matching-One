#!/usr/bin/env python3
"""#640 wrapping-type 4x4-by-rank census for tiny axis/diamond tori.

For each occupation k, count configurations by
  (primal wrap type) x (matching wrap type)
where wrap type is neither / dir0 / dir1 / both, using torus_homology
generator channels. Collapsing with
  D = 1{primal either} - 1{matching either}
must reproduce committed Bernstein integers.

L=3,4 are seconds on a laptop. L=5 (2^25) belongs on Huawei.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from exact_matching_polynomial import bernstein_counts  # noqa: E402
from matched_torus_reference import (  # noqa: E402
    axis_geometry,
    diamond_geometry,
    popcount,
)
from torus_homology import (  # noqa: E402
    component_homologies,
    geometry_periods,
    wrapping_channels,
)

WRAP_LABELS = ("neither", "dir0", "dir1", "both")


def wrap_label(channels) -> str:
    if channels.both:
        return "both"
    if channels.direction_0:
        return "dir0"
    if channels.direction_1:
        return "dir1"
    return "neither"


def either(label: str) -> bool:
    return label != "neither"


def census_one(geometry, mask_lo: int, mask_hi: int) -> dict:
    n = geometry.n
    periods = geometry_periods(geometry)
    tables = {
        k: {a: {b: 0 for b in WRAP_LABELS} for a in WRAP_LABELS}
        for k in range(n + 1)
    }
    for mask in range(mask_lo, mask_hi):
        k = popcount(mask)
        black = [bool((mask >> i) & 1) for i in range(n)]
        white = [not value for value in black]
        pb = wrap_label(wrapping_channels(
            component_homologies(black, geometry.primal_edges, periods)))
        pw = wrap_label(wrapping_channels(
            component_homologies(white, geometry.matching_edges, periods)))
        tables[k][pb][pw] += 1
    collapsed = [0] * (n + 1)
    for k in range(n + 1):
        for a in WRAP_LABELS:
            for b in WRAP_LABELS:
                collapsed[k] += tables[k][a][b] * (int(either(a)) - int(either(b)))
    return {"tables": tables, "collapsed": collapsed}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--geometry", choices=("axis", "diamond"), required=True)
    ap.add_argument("--L", type=int, required=True)
    ap.add_argument("--mask-lo", type=int, default=0)
    ap.add_argument("--mask-hi", type=int, default=None)
    args = ap.parse_args()
    geo = axis_geometry(args.L) if args.geometry == "axis" else diamond_geometry(args.L)
    lo = args.mask_lo
    hi = args.mask_hi if args.mask_hi is not None else (1 << geo.n)
    result = census_one(geo, lo, hi)
    if lo == 0 and hi == (1 << geo.n):
        committed = bernstein_counts(geo)
        if result["collapsed"] != committed:
            raise SystemExit(
                f"tripwire failed: collapsed={result['collapsed']} "
                f"committed={committed}"
            )
        result["tripwire"] = "pass"
        result["committed_bernstein"] = committed
    result["geometry"] = args.geometry
    result["L"] = args.L
    result["N"] = geo.n
    result["mask_lo"] = lo
    result["mask_hi"] = hi
    print(json.dumps(result))


if __name__ == "__main__":
    main()
