#!/usr/bin/env python3
"""Local entry: self-check, independent witness verify, optional EXH cache.

Usage (from package root or repo root):
  python3 research/summary_search/run_local.py
  python3 research/summary_search/run_local.py --exh     # also enumerate n<=5
  python3 research/summary_search/verify_witness.py
"""
from __future__ import annotations

import argparse
import pickle
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from bounded_summary_search import exhaustive_n, self_checks  # noqa: E402


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--exh", action="store_true", help="enumerate planar two-terminal n<=5 and pickle")
    args = p.parse_args()

    print("=== self-checks ===", flush=True)
    self_checks()

    print("=== independent witness ===", flush=True)
    import verify_witness

    rc = verify_witness.main()
    if rc != 0:
        return rc

    cache = HERE / "cache_exh.pkl"
    if args.exh or not cache.exists():
        print("=== exhaustive n<=5 ===", flush=True)
        exh = []
        for n in range(1, 6):
            got = exhaustive_n(n)
            print(f"  n={n} {len(got)}", flush=True)
            exh.extend(got)
        cache.write_bytes(pickle.dumps(exh, protocol=pickle.HIGHEST_PROTOCOL))
        print("wrote", cache, "N=", len(exh), flush=True)
    else:
        print("EXH cache present:", cache, flush=True)

    print("OK. Next: python3 research/summary_search/run_full_search.py", flush=True)
    print("  (needs EXH cache; SP cache optional — will generate SP n<=10 if missing)", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
