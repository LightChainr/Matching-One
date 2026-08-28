#!/usr/bin/env python3
"""Compare C++ per-config dumps against the independent lifted-graph oracle."""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import exact_oracle as oracle  # noqa: E402

DUMPS = [
    (2, ROOT / "results/issue-7/L02_configs.csv"),
    (3, ROOT / "results/issue-7/L03_configs.csv"),
]


def main() -> int:
    report = {"oracle": "lifted-graph BFS, not C++ DSU", "compares": []}
    overall = "PASS"
    for L, path in DUMPS:
        if not path.is_file():
            print(f"MISSING {path}")
            overall = "FAIL"
            report["compares"].append({"L": L, "status": "FAIL", "reason": "missing dump"})
            continue
        cmp = oracle.compare_dump(L, path)
        ident = oracle.identities_from_microcanonical(L, oracle.enumerate_microcanonical(L))
        entry = {
            "L": L,
            "dump": str(path),
            "cpp_vs_oracle": cmp["status"],
            "oracle_identities": ident["identity"],
            "mismatches": cmp.get("mismatches", [])[:8],
            "wrapping": {k: v["status"] for k, v in ident["wrapping"].items()},
        }
        report["compares"].append(entry)
        print(
            f"L={L} cpp_vs_oracle={cmp['status']} mismatches={len(cmp.get('mismatches', []))} "
            f"oracle_identity={ident['identity']}"
        )
        if cmp["status"] != "PASS" or ident["identity"] != "PASS":
            overall = "FAIL"
    report["overall"] = overall
    out = ROOT / "results/issue-7/oracle_compare.json"
    out.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"ORACLE_OVERALL={overall}")
    return 0 if overall == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
