#!/usr/bin/env python3
"""Finite controls for the two-birth reduction (not a new width scan).

The universal concentration constant comes from Friedgut--Kalai and is NOT
estimated here. Width-two formulas are an existing exact calibration model.
Only mpmath is required beyond the standard library. No random sampling.
"""
from __future__ import annotations
import argparse
from fractions import Fraction
from itertools import permutations
import json
from math import comb, factorial
from pathlib import Path
import mpmath as mp


def lifted_rank(width: int, length: int, mask: int) -> int:
    """Ambient rank of the occupied NN graph, retaining parallel lifted edges."""
    if width < 2 or length < 2 or not 0 <= mask < (1 << (width * length)):
        raise ValueError("Expected width,length >= 2 and a legal site mask")
    potentials: dict[int, tuple[int, int]] = {}
    generator: tuple[int, int] | None = None
    for root in range(width * length):
        if not (mask >> root & 1) or root in potentials:
            continue
        potentials[root] = (0, 0)
        stack = [root]
        while stack:
            vertex = stack.pop()
            x, y = vertex % width, vertex // width
            px, py = potentials[vertex]
            for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                other = ((y + dy) % length) * width + (x + dx) % width
                if not (mask >> other & 1):
                    continue
                expected = (px + dx, py + dy)
                if other not in potentials:
                    potentials[other] = expected
                    stack.append(other)
                else:
                    old = potentials[other]
                    gain = (expected[0] - old[0], expected[1] - old[1])
                    if gain != (0, 0):
                        if generator is None:
                            generator = gain
                        elif generator[0] * gain[1] != generator[1] * gain[0]:
                            return 2
    return int(generator is not None)


def census(width: int, length: int) -> list[list[int]]:
    n = width * length
    if n > 16:
        raise ValueError("This control deliberately caps exhaustive enumeration at 16 sites")
    counts = [[0] * (n + 1) for _ in range(3)]
    for mask in range(1 << n):
        counts[lifted_rank(width, length, mask)][mask.bit_count()] += 1
    return counts


def beta_integral(counts: list[int]) -> Fraction:
    """Integral of sum_k counts[k] p^k(1-p)^(N-k), exactly."""
    n = len(counts) - 1
    return sum((Fraction(c, (n + 1) * comb(n, k))
                for k, c in enumerate(counts)), Fraction())


def census_probability(counts: list[int], p):
    n = len(counts) - 1
    return sum(c * p**k * (1 - p)**(n - k) for k, c in enumerate(counts))


def sector_probabilities(length: int, p):
    """Exact algebraic width-two expressions, evaluated at current mp precision."""
    if length < 2 or not 0 <= p <= 1:
        raise ValueError("Expected length >= 2 and p in [0,1]")
    p = mp.mpf(p)
    if p == 0:
        return mp.mpf(1), mp.mpf(0), mp.mpf(0)
    if p == 1:
        return mp.mpf(0), mp.mpf(0), mp.mpf(1)
    x, y = p * (1 - p), p * p
    lp = p * (1 + mp.sqrt(1 + 4 * p * (1 - p))) / 2
    lm = -x * y / lp  # stable product relation, instead of subtracting square roots
    p0 = (1 - y)**length - 2 * x**length
    p2 = lp**length + lm**length - x**length
    return p0, 1 - p0 - p2, p2


def cdf(length: int, channel: str, p):
    p0, _, p2 = sector_probabilities(length, p)
    if channel == "first":
        return 1 - p0
    if channel == "second":
        return p2
    if channel == "mixture":
        return (1 - p0 + p2) / 2
    raise ValueError("channel must be first, second, or mixture")


def quantile(length: int, channel: str, u):
    u = mp.mpf(u)
    if not 0 < u < 1:
        raise ValueError("Quantile level must be strictly between 0 and 1")
    lo, hi = mp.mpf(0), mp.mpf(1)
    # Reporting precision is intentionally far below this working precision.
    for _ in range(4 * mp.mp.dps):
        mid = (lo + hi) / 2
        if mid == lo or mid == hi:
            break
        if cdf(length, channel, mid) < u:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def exact_small_controls() -> dict:
    records = []
    for length in (2, 3, 4):
        counts = census(2, length)
        n = 2 * length
        for k in range(n + 1):
            assert sum(row[k] for row in counts) == comb(n, k)
        for rational in (Fraction(1, 3), Fraction(1, 2), Fraction(2, 3)):
            p = mp.mpf(rational.numerator) / rational.denominator
            got = sector_probabilities(length, p)
            for rank in range(3):
                expected = census_probability(counts[rank], rational)
                assert abs(got[rank] - mp.mpf(expected.numerator) / expected.denominator) < mp.mpf("1e-50")
        mu1 = beta_integral(counts[0])
        mu2 = 1 - beta_integral(counts[2])
        gap = beta_integral(counts[1])
        assert mu2 - mu1 == gap
        records.append({"width": 2, "length": length, "configurations": 1 << n,
                        "counts_by_rank_and_occupation": counts,
                        "mean_first_birth": str(mu1), "mean_second_birth": str(mu2),
                        "integral_P1_exact": str(gap)})
    # Independent occupation-order check of the normalized rank-gap identity.
    n, sum_first, sum_second, sum_product = 6, 0, 0, 0
    for order in permutations(range(n)):
        mask, k1, k2 = 0, None, None
        for k, v in enumerate(order, 1):
            mask |= 1 << v
            r = lifted_rank(2, 3, mask)
            if r >= 1 and k1 is None:
                k1 = k
            if r == 2:
                k2 = k
                break
        assert k1 is not None and k2 is not None
        sum_first += k1
        sum_second += k2
        sum_product += k1 * (k2 + 1)
    mu1 = Fraction(sum_first, factorial(n) * (n + 1))
    mu2 = Fraction(sum_second, factorial(n) * (n + 1))
    covariance = Fraction(sum_product, factorial(n) * (n + 1) * (n + 2)) - mu1 * mu2
    assert str(mu2 - mu1) == records[1]["integral_P1_exact"]
    assert covariance > 0  # independent births are not assumed, even in tiny controls
    return {"censuses": records, "total_configurations": sum(r["configurations"] for r in records),
            "permutation_control": {"width": 2, "length": 3, "permutations": factorial(n),
                "normalized_mean_birth_rank_gap": str(mu2 - mu1),
                "birth_time_covariance": str(covariance)}}


def finite_record(length: int, digits: int = 24) -> dict:
    a = quantile(length, "first", mp.mpf("0.5"))
    b = quantile(length, "second", mp.mpf("0.5"))
    q = quantile(length, "mixture", mp.mpf("0.5"))
    q25 = quantile(length, "mixture", mp.mpf("0.25"))
    q75 = quantile(length, "mixture", mp.mpf("0.75"))
    assert q25 <= a <= q <= b <= q75
    f1 = lambda p: cdf(length, "first", p)
    f2 = lambda p: cdf(length, "second", p)
    f = lambda p: cdf(length, "mixture", p)
    mu1 = mp.quad(lambda p: 1 - f1(p), [0, a, b, 1])
    mu2 = mp.quad(lambda p: 1 - f2(p), [0, a, b, 1])
    gap = mp.quad(lambda p: sector_probabilities(length, p)[1], [0, a, b, 1])
    assert abs(mu2 - mu1 - gap) < mp.mpf("1e-45")
    mad1 = mp.quad(f1, [0, a]) + mp.quad(lambda p: 1 - f1(p), [a, b, 1])
    mad2 = mp.quad(f2, [0, a, b]) + mp.quad(lambda p: 1 - f2(p), [b, 1])
    w1 = (mp.quad(f, [0, a]) + mp.quad(lambda p: mp.mpf(".5") - f(p), [a, q])
          + mp.quad(lambda p: f(p) - mp.mpf(".5"), [q, b])
          + mp.quad(lambda p: 1 - f(p), [b, 1]))
    assert w1 <= (mad1 + mad2) / 2 + mp.mpf("1e-45")
    assert abs(gap - (b - a)) <= mad1 + mad2 + mp.mpf("1e-45")
    values = {"first_birth_median": a, "matching_root": q, "second_birth_median": b,
              "median_separation": b - a, "mixture_q25": q25, "mixture_q75": q75,
              "mixture_IQR": q75 - q25, "integral_P1": gap,
              "first_birth_MAD": mad1, "second_birth_MAD": mad2,
              "mixture_W1_to_two_median_atoms": w1,
              "component_coupling_W1_upper": (mad1 + mad2) / 2}
    return {"width": 2, "length": length, "sites": 2 * length,
            **{key: mp.nstr(value, digits) for key, value in values.items()}}


def make_report(dps: int = 80) -> dict:
    if dps < 60:
        raise ValueError("Use at least 60 working decimal digits")
    with mp.workdps(dps):
        return {"schema": "matching-one.two-birth-reduction.v1",
            "scope": "Existing width-two exact oracle; no new width, no Monte Carlo",
            "precision": {"working_decimal_digits": dps, "reported_digits": 24,
                "quadrature_status": "high precision numerical control, not interval certification"},
            "universal_constant": "Friedgut--Kalai constant left unspecified; not fitted or priced",
            "exact_controls": exact_small_controls(),
            "finite_length_controls": [finite_record(m) for m in (2, 4, 8, 16, 32, 128)],
            "interpretation": "Fixed width 2 does not approach the infinite square-site critical point"}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--dps", type=int, default=80)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("Refusing to overwrite an existing result file")
    report = make_report(args.dps)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
