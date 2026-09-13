#!/usr/bin/env python3
"""Exact controls for generator-dependent rank and multi-step birth prediction.

This standard-library reference does not consume Monte Carlo data or identify
continuum fields. The independent topology calculation is deliberately limited
to the N=10 square-NN quotient with period columns (3,1),(-1,3).
"""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from itertools import combinations, permutations
import json
from math import comb, factorial, gcd
from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence, Tuple

Vector = Tuple[int, int]
Matrix = List[List[Fraction]]
Mark = Tuple[int, Optional[Vector]]


def _positive_int(value: int, label: str) -> None:
    if not isinstance(value, int) or isinstance(value, bool) or value < 1:
        raise ValueError(label + " must be a positive integer")


def exact_rank(matrix: Sequence[Sequence[Fraction]]) -> int:
    """Rational row rank, with a rectangular-input gate."""
    if not matrix:
        return 0
    width = len(matrix[0])
    if any(len(row) != width for row in matrix):
        raise ValueError("matrix rows have different widths")
    a = [[Fraction(x) for x in row] for row in matrix]
    rank = 0
    for col in range(width):
        pivot = next((i for i in range(rank, len(a)) if a[i][col]), None)
        if pivot is None:
            continue
        a[rank], a[pivot] = a[pivot], a[rank]
        scale = a[rank][col]
        a[rank] = [x / scale for x in a[rank]]
        for row in range(rank + 1, len(a)):
            scale = a[row][col]
            if scale:
                a[row] = [x - scale*y for x, y in zip(a[row], a[rank])]
        rank += 1
        if rank == len(a):
            break
    return rank


def exact_determinant(matrix: Sequence[Sequence[Fraction]]) -> Fraction:
    """Determinant by exact elimination, independent of the Hilbert product."""
    n = len(matrix)
    if any(len(row) != n for row in matrix):
        raise ValueError("determinant requires a square matrix")
    a = [[Fraction(x) for x in row] for row in matrix]
    result = Fraction(1)
    for col in range(n):
        pivot = next((i for i in range(col, n) if a[i][col]), None)
        if pivot is None:
            return Fraction(0)
        if pivot != col:
            a[col], a[pivot] = a[pivot], a[col]
            result = -result
        value = a[col][col]
        result *= value
        for row in range(col + 1, n):
            scale = a[row][col] / value
            for j in range(col + 1, n):
                a[row][j] -= scale * a[col][j]
    return result


def hankel(sequence: Callable[[int], Fraction], size: int) -> Matrix:
    _positive_int(size, "size")
    return [[Fraction(sequence(i+j)) for j in range(size)] for i in range(size)]


def hilbert_cauchy_determinant(size: int) -> Fraction:
    _positive_int(size, "size")
    numerator = denominator = 1
    for i in range(size):
        for j in range(i + 1, size):
            numerator *= (j-i)**2
        for j in range(size):
            denominator *= i+j+1
    return Fraction(numerator, denominator)


def generator_controls(max_size: int = 10) -> Dict:
    _positive_int(max_size, "max_size")
    rows = []
    for size in range(1, max_size + 1):
        power = hankel(lambda n: Fraction(1, n+1), size)
        determinant = exact_determinant(power)
        if determinant != hilbert_cauchy_determinant(size) or determinant <= 0:
            raise AssertionError("Hilbert determinant certificate failed")
        ranks = {
            "additive_power_rank": exact_rank(power),
            "geometric_power_rank": exact_rank(hankel(lambda n: Fraction(1, 2**n), size)),
            "geometric_log_pair_rank": exact_rank(hankel(lambda n: Fraction(n+1, 2**n), size)),
        }
        if list(ranks.values()) != [size, 1, min(size, 2)]:
            raise AssertionError("generator-dependent rank control failed")
        rows.append({"size": size, **ranks, "hilbert_determinant": str(determinant)})
    # g(r)=r^-1 at r=n+1 versus r=2^n. For the logarithmic control,
    # g(r)=r^-1(1+log_2 r), so no irrational numerical arithmetic is needed.
    g = lambda n: Fraction(1, n+1)
    recurrence = [(n+2)*g(n+1)-(n+1)*g(n) for n in range(2*max_size)]
    if any(recurrence):
        raise AssertionError("variable-coefficient recurrence failed")
    false_rank_one_residual = Fraction(1, 3) - Fraction(1, 2)*Fraction(1, 2)
    return {
        "meaning": "exact kinematic controls, not a production rank decision",
        "power_function": "g(r)=1/r",
        "additive_grid": "r_n=n+1",
        "geometric_grid": "r_n=2^n",
        "log_pair_function": "g(r)=(1+log_2(r))/r",
        "rows": rows,
        "variable_coefficient_recurrence": "(n+2)g_(n+1)-(n+1)g_n=0",
        "recurrence_residuals": [str(x) for x in recurrence],
        "false_constant_rank_one_first_heldout_residual": str(false_rank_one_residual),
    }


def one_clock_covariance(quantiles: Sequence[Fraction]) -> Matrix:
    """Covariance of 1[U<=u_i], U uniform; not a multi-field state count."""
    values = [Fraction(value) for value in quantiles]
    if not values or not 0 < values[0] or not values[-1] < 1:
        raise ValueError("quantiles must lie strictly between zero and one")
    if any(a >= b for a, b in zip(values, values[1:])):
        raise ValueError("quantiles must be strictly increasing")
    return [[min(u, v)-u*v for v in values] for u in values]


def one_clock_controls(max_size: int = 10) -> Dict:
    _positive_int(max_size, "max_size")
    rows = []
    for size in range(1, max_size+1):
        matrix = one_clock_covariance([Fraction(i, size+1) for i in range(1, size+1)])
        determinant = exact_determinant(matrix)
        if exact_rank(matrix) != size or determinant != Fraction(1, (size+1)**(size+1)):
            raise AssertionError("single-clock covariance control failed")
        rows.append({"size": size, "covariance_rank": size, "determinant": str(determinant)})
    return {
        "meaning": "a single threshold clock can have full temporal covariance rank",
        "process": "X(u)=1[U<=u], U uniform on (0,1)",
        "kernel": "min(u,v)-u*v",
        "rows": rows,
        "boundary": "not a replay of the P334 temporal kernel or its eigenvalue fractions",
    }


class N10BirthOracle:
    """Small independent lifted-BFS control, not a general production backend.

    Site j is (0,j) modulo [[3,-1],[1,3]]. Physical unit steps +/-x map
    to j+/-3; +/-y map to j+/-1. Lifted cycle displacements are transformed
    to period coordinates by adj(P)/10. Filled contractible faces do not
    change the ambient H1 image.
    """
    n = 10
    periods = ((3, -1), (1, 3))
    steps = ((3, 1, 0), (-3, -1, 0), (1, 0, 1), (-1, 0, -1))

    def __init__(self) -> None:
        self.marks = tuple(self._mark(mask) for mask in range(1 << self.n))
        for mask in range(1 << self.n):
            rank = self.marks[mask][0]
            for v in self.vacant(mask):
                if self.marks[mask | (1 << v)][0] < rank:
                    raise AssertionError("ambient rank is not monotone")

    def _validate_mask(self, mask: int) -> None:
        if not isinstance(mask, int) or isinstance(mask, bool) or not 0 <= mask < 1 << self.n:
            raise ValueError("mask must encode the ten declared sites")

    def vacant(self, mask: int) -> List[int]:
        self._validate_mask(mask)
        return [v for v in range(self.n) if not (mask >> v) & 1]

    def _mark(self, mask: int) -> Mark:
        lifts: Dict[int, Vector] = {}
        generators = []
        for root in range(self.n):
            if not (mask >> root) & 1 or root in lifts:
                continue
            lifts[root] = (0, 0)
            stack = [root]
            while stack:
                u = stack.pop()
                x, y = lifts[u]
                for offset, dx, dy in self.steps:
                    v = (u+offset) % self.n
                    if not (mask >> v) & 1:
                        continue
                    candidate = (x+dx, y+dy)
                    if v not in lifts:
                        lifts[v] = candidate
                        stack.append(v)
                        continue
                    X, Y = candidate[0]-lifts[v][0], candidate[1]-lifts[v][1]
                    if X or Y:
                        first, second = 3*X+Y, -X+3*Y
                        if first % 10 or second % 10:
                            raise AssertionError("cycle is not in the period lattice")
                        generators.append((first//10, second//10))
        if not generators:
            return 0, None
        x, y = generators[0]
        if any(x*Y-y*X for X, Y in generators[1:]):
            return 2, None
        divisor = gcd(abs(x), abs(y))
        x, y = x//divisor, y//divisor
        if x < 0 or (x == 0 and y < 0):
            x, y = -x, -y
        return 1, (x, y)

    @lru_cache(maxsize=None, typed=True)
    def survival(self, mask: int, steps: int) -> Fraction:
        """Killed-kernel survival from a rank-one state; zero after killing."""
        self._validate_mask(mask)
        if not isinstance(steps, int) or isinstance(steps, bool) or steps < 0:
            raise ValueError("steps must be a nonnegative integer")
        if self.marks[mask][0] != 1:
            return Fraction(0)
        if steps == 0:
            return Fraction(1)
        vacant = self.vacant(mask)
        if steps > len(vacant):
            return Fraction(0)
        return sum((self.survival(mask | (1 << v), steps-1) for v in vacant), Fraction(0))/len(vacant)

    def subset_survival(self, mask: int, steps: int) -> Fraction:
        """Independent unordered-subset evaluation of the same survival law."""
        self._validate_mask(mask)
        if not isinstance(steps, int) or isinstance(steps, bool) or steps < 0:
            raise ValueError("steps must be a nonnegative integer")
        vacant = self.vacant(mask)
        if self.marks[mask][0] != 1 or steps > len(vacant):
            return Fraction(0)
        count = 0
        for chosen in combinations(vacant, steps):
            target = mask
            for v in chosen:
                target |= 1 << v
            count += self.marks[target][0] == 1
        return Fraction(count, comb(len(vacant), steps))

    def triggers(self, mask: int) -> Tuple[List[int], List[Tuple[int, int]]]:
        self._validate_mask(mask)
        if self.marks[mask][0] != 1:
            raise ValueError("trigger counts require a current rank-one state")
        vacant = self.vacant(mask)
        singles = [v for v in vacant if self.marks[mask | (1 << v)][0] == 2]
        safe = [v for v in vacant if v not in singles]
        pairs = [(u, v) for u, v in combinations(safe, 2)
                 if self.marks[mask | (1 << u) | (1 << v)][0] == 2]
        return singles, pairs

    def exit_counts(self, mask: int) -> Dict[int, int]:
        self._validate_mask(mask)
        if self.marks[mask][0] != 1:
            raise ValueError("exit law requires a current rank-one state")
        counts: Counter = Counter()
        for sequence in permutations(self.vacant(mask)):
            target = mask
            for time, vertex in enumerate(sequence, 1):
                target |= 1 << vertex
                if self.marks[target][0] == 2:
                    counts[time] += 1
                    break
            else:
                raise AssertionError("fully occupied state failed to exit")
        return dict(sorted(counts.items()))


def direct_birth_priority_control(oracle: N10BirthOracle) -> Dict:
    """Compare directed-edge Beta weights with an independent permutation DP."""
    edge_counts: Counter = Counter()
    for mask, (rank, _) in enumerate(oracle.marks):
        if rank != 0:
            continue
        k = bin(mask).count("1")
        for vertex in oracle.vacant(mask):
            if oracle.marks[mask | (1 << vertex)][0] == 2:
                edge_counts[k] += 1
    denominator = factorial(oracle.n)
    beta_sum = sum((Fraction(count*factorial(k)*factorial(oracle.n-k-1), denominator)
                    for k, count in edge_counts.items()), Fraction(0))

    @lru_cache(maxsize=None)
    def count_future_paths(mask: int) -> int:
        if oracle.marks[mask][0] != 0:
            return 0
        vacant = oracle.vacant(mask)
        total = 0
        for vertex in vacant:
            child = mask | (1 << vertex)
            rank = oracle.marks[child][0]
            if rank == 2:
                total += factorial(len(vacant)-1)
            elif rank == 0:
                total += count_future_paths(child)
        return total

    paths = count_future_paths(0)
    if beta_sum != Fraction(paths, denominator):
        raise AssertionError("direct-birth priority/Beta identity failed")
    return {
        "meaning": "exact untyped N10 priority identity; no arm-exponent claim",
        "directed_edge_counts_by_k": {str(k): c for k, c in sorted(edge_counts.items())},
        "directed_edge_count": sum(edge_counts.values()),
        "permutation_paths_with_direct_birth": paths,
        "total_permutation_paths": denominator,
        "probability_from_beta_weights": str(beta_sum),
        "probability_from_permutation_dp": str(Fraction(paths, denominator)),
    }


def birth_controls() -> Dict:
    oracle = N10BirthOracle()
    checked = 0
    for mask, (rank, _) in enumerate(oracle.marks):
        if rank != 1:
            continue
        vacant = oracle.vacant(mask)
        q = len(vacant)
        singles, pairs = oracle.triggers(mask)
        if oracle.survival(mask, 1) != 1-Fraction(len(singles), q):
            raise AssertionError("one-step ceiling failed")
        if q >= 2:
            safe = q-len(singles)
            choose_safe_two = comb(safe, 2) if safe >= 2 else 0
            predicted = Fraction(choose_safe_two-len(pairs), comb(q, 2))
            if oracle.survival(mask, 2) != predicted:
                raise AssertionError("two-step trigger identity failed")
        for steps in range(q+1):
            if oracle.survival(mask, steps) != oracle.subset_survival(mask, steps):
                raise AssertionError("killed/subset survival mismatch")
            checked += 1
    witnesses = []
    for mask in (31, 47):
        singles, pairs = oracle.triggers(mask)
        q = len(oracle.vacant(mask))
        counts = oracle.exit_counts(mask)
        for steps in range(q+1):
            survival = Fraction(sum(count for t, count in counts.items() if t > steps), factorial(q))
            if survival != oracle.survival(mask, steps):
                raise AssertionError("future permutation check failed")
        witnesses.append({
            "mask": mask, "occupied_site_labels": [v for v in range(10) if mask >> v & 1],
            "rank": oracle.marks[mask][0], "line": list(oracle.marks[mask][1]),
            "singleton_triggers_H2": singles, "minimal_trigger_pairs": [list(pair) for pair in pairs],
            "one_step_exit_probability": str(Fraction(len(singles), q)),
            "survival_probabilities": [str(oracle.survival(mask, j)) for j in range(q+1)],
            "future_permutation_count": factorial(q),
            "exit_counts": {str(t): count for t, count in counts.items()},
        })
    if witnesses[0]["survival_probabilities"][2] != "2/5" or witnesses[1]["survival_probabilities"][2] != "3/10":
        raise AssertionError("declared N10 witness changed")
    return {
        "meaning": "exact finite counterexample, not a production memory claim",
        "period_matrix": [list(row) for row in oracle.periods],
        "site_label": "j corresponds to (0,j), j=0,...,9",
        "configuration_count": 1 << oracle.n,
        "rank_counts": {str(r): sum(mark[0] == r for mark in oracle.marks) for r in range(3)},
        "killed_vs_subset_checks": checked,
        "direct_birth_priority_control": direct_birth_priority_control(oracle),
        "witnesses": witnesses,
        "two_step_survival_gap": str(oracle.survival(31, 2)-oracle.survival(47, 2)),
    }


def conditional_arm_arithmetic() -> Dict:
    """Arithmetic only: no claim that typed square-site arm assumptions hold."""
    rows = []
    for arms in (6, 8):
        alpha = Fraction(arms*arms-1, 12)
        area_decay = (alpha+Fraction(3, 4)-2)/2
        rows.append({"arms": arms, "arm_exponent": str(alpha), "area_decay_exponent": str(area_decay)})
    return {
        "claim": "conditional scaling arithmetic only; not an arm correspondence theorem",
        "rows": rows,
        "relative_area_decay_exponent": str(Fraction(2)-Fraction(5, 6)),
        "required_inputs": ["typed-arm frequency", "near-critical window", "integrable scaling profile", "square-site universality"],
    }


def build_artifact() -> Dict:
    return {
        "schema": "matching-one/predictive-state-counterexamples/v1",
        "issues": [400, 403, 405],
        "data_class": "independent exact controls; no Monte Carlo inputs",
        "generator_controls": generator_controls(),
        "one_clock_covariance_controls": one_clock_controls(),
        "birth_controls": birth_controls(),
        "conditional_arm_arithmetic": conditional_arm_arithmetic(),
        "boundary": "Does not replay P250/P334/P337 production, infer physical field count, identify Jordan, or prove square-site critical exponents.",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    rendered = json.dumps(build_artifact(), indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
