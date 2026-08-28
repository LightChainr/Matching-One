#!/usr/bin/env python3
"""Reconstruct M(p), M'(p), and the finite root from threshold-rank counts."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import List, Optional, Sequence, Tuple

import mpmath as mp


def read_histogram(path: Path, rank_name: str, n: int) -> List[int]:
    values = [0] * (n + 1)
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            rank = int(row[rank_name])
            count = int(row["count"])
            if not 1 <= rank <= n or count < 0:
                raise ValueError(f"invalid histogram row in {path}")
            values[rank] += count
    return values


def read_counts(input_dir: Path) -> Tuple[int, int, List[int], List[int], dict]:
    metadata = json.loads((input_dir / "metadata.json").read_text(encoding="utf-8"))
    n = int(metadata["N"])
    samples = int(metadata["sample_count"])
    minus = read_histogram(input_dir / "kminus_hist.csv", "K_minus", n)
    plus = read_histogram(input_dir / "kplus_hist.csv", "K_plus", n)
    if sum(minus) != samples or sum(plus) != samples:
        raise ValueError("histogram totals do not equal metadata sample_count")
    return n, samples, minus, plus, metadata


def threshold_cdf(n: int, rank: int, p: mp.mpf) -> mp.mpf:
    if p <= 0:
        return mp.mpf(0)
    if p >= 1:
        return mp.mpf(1)
    return mp.betainc(rank, n + 1 - rank, 0, p, regularized=True)


def threshold_pdf(n: int, rank: int, p: mp.mpf) -> mp.mpf:
    if not 0 < p < 1:
        return mp.mpf(0)
    return p ** (rank - 1) * (1 - p) ** (n - rank) / mp.beta(
        rank, n + 1 - rank
    )


def matching_value(
    n: int,
    samples: int,
    minus: Sequence[int],
    plus: Sequence[int],
    p: mp.mpf,
) -> mp.mpf:
    if p <= 0:
        return mp.mpf(-1)
    if p >= 1:
        return mp.mpf(1)
    # Rearrange sum_k h_k P(Bin(N,p)>=k) as a single binomial
    # convolution against cumulative threshold counts.  This avoids N calls
    # to incomplete beta functions for every root iteration.
    q = 1 - p
    probability = q**n
    cumulative = 0
    total = mp.mpf(0)
    for occupied in range(n + 1):
        if occupied:
            cumulative += minus[occupied] + plus[occupied]
        total += cumulative * probability
        if occupied < n:
            probability *= (n - occupied) * p / ((occupied + 1) * q)
    return total / samples - 1


def matching_derivative(
    n: int,
    samples: int,
    minus: Sequence[int],
    plus: Sequence[int],
    p: mp.mpf,
) -> mp.mpf:
    if not 0 < p < 1:
        return mp.mpf(0)
    # BetaPDF(k,N+1-k) = N*C(N-1,k-1)*p^(k-1)*(1-p)^(N-k).
    q = 1 - p
    density = n * q ** (n - 1)
    total = mp.mpf(0)
    for rank in range(1, n + 1):
        total += (minus[rank] + plus[rank]) * density
        if rank < n:
            density *= (n - rank) * p / (rank * q)
    return total / samples


def matching_root(
    n: int,
    samples: int,
    minus: Sequence[int],
    plus: Sequence[int],
    *,
    iterations: Optional[int] = None,
) -> mp.mpf:
    if iterations is None:
        iterations = max(128, math.ceil(mp.mp.dps * math.log2(10)) + 8)
    lower = mp.mpf(0)
    upper = mp.mpf(1)
    for _ in range(iterations):
        midpoint = (lower + upper) / 2
        value = matching_value(n, samples, minus, plus, midpoint)
        if value < 0:
            lower = midpoint
        else:
            upper = midpoint
    return (lower + upper) / 2


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--p", action="append", default=[])
    parser.add_argument("--dps", type=int, default=50)
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()
    if args.dps < 30:
        raise SystemExit("--dps must be at least 30")
    mp.mp.dps = args.dps

    n, samples, minus, plus, metadata = read_counts(args.input_dir)
    probabilities = [mp.mpf(value) for value in args.p]
    if not probabilities:
        probabilities = [mp.mpf("0.5")]
    for p in probabilities:
        if not 0 < p < 1:
            raise SystemExit("each --p must lie strictly between zero and one")

    root = matching_root(n, samples, minus, plus)
    evaluations = [
        {
            "p": mp.nstr(p, args.dps),
            "M": mp.nstr(matching_value(n, samples, minus, plus, p), args.dps),
            "M_prime": mp.nstr(
                matching_derivative(n, samples, minus, plus, p), args.dps
            ),
        }
        for p in probabilities
    ]
    payload = {
        "format_version": 1,
        "source_geometry": metadata["geometry"],
        "N": n,
        "sample_count": samples,
        "root": mp.nstr(root, args.dps),
        "evaluations": evaluations,
    }
    print(json.dumps(payload, indent=2))
    if args.json is not None:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
