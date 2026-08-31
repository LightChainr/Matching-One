#!/usr/bin/env python3
"""Canonical Kreg two-site Q1 kernel on the 4140 Bell8 port connectivities."""
import argparse
import csv
from collections import Counter, defaultdict
from functools import lru_cache
from math import factorial
from pathlib import Path
import json


@lru_cache(None)
def partitions(n):
    def grow(prefix):
        if len(prefix) == n:
            yield prefix
        else:
            for label in range(max(prefix) + 2):
                yield from grow(prefix + (label,))
    return tuple(grow((0,)))


@lru_cache(None)
def four_kreg_q1(p):
    """4*Kreg on an exact colour-equality pattern; integral and phase-fixed."""
    total = 0
    for v in (p, p[1:] + p[:1]):
        a, b, c, d = v
        if a != b and c != d:
            total += (int(a == c and b == d) + int(a == d and b == c)
                      + sum((a == c, a == d, b == c, b == d)) - 4)
    return total


def packed_key(p):
    return sum(label << (3 * i) for i, label in enumerate(p))


def g16(p):
    answer = 0
    for merger in partitions(max(p) + 1):
        k = max(merger) + 1
        if k >= 2:
            actual = tuple(merger[label] for label in p)
            answer += ((-1) ** (k - 2) * factorial(k - 2)
                       * four_kreg_q1(actual[:4]) * four_kreg_q1(actual[4:]))
    return answer


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        default=Path("analysis/regular_pair_spatial_kernel.tsv"))
    args = parser.parse_args()
    rows = []
    classes = defaultdict(list)
    for p in partitions(8):
        value = g16(p)
        shared = len(set(p[:4]) & set(p[4:]))
        classes[shared].append((p, value))
        if value:
            rows.append((packed_key(p), value))
    rows.sort()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(("key", "g16"))
        writer.writerows(rows)
    summary = {
        "port_order": ["xN", "xE", "xS", "xW", "yN", "yE", "yS", "yW"],
        "kernel": "g = g16 / 16; omitted VALID canonical keys are exactly zero",
        "total_Bell8": sum(len(v) for v in classes.values()),
        "nonzero_rows": len(rows),
        "shared_component_classes": {},
    }
    for shared, values in sorted(classes.items()):
        nonzero = [(p, v) for p, v in values if v]
        summary["shared_component_classes"][shared] = {
            "total": len(values), "nonzero": len(nonzero),
            "positive": sum(v > 0 for _, v in values),
            "negative": sum(v < 0 for _, v in values),
            "g16_histogram": dict(sorted(Counter(v for _, v in values).items())),
            "simplest_positive": next(({"labels": list(p), "g16": v}
                for p, v in sorted(nonzero, key=lambda r: (max(r[0]), r[0]))
                if v > 0), None),
            "simplest_negative": next(({"labels": list(p), "g16": v}
                for p, v in sorted(nonzero, key=lambda r: (max(r[0]), r[0]))
                if v < 0), None),
        }
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
