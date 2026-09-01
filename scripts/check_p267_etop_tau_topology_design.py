#!/usr/bin/env python3
"""Exact certificate for the frozen N=50 tau x topology-map factorial."""

from __future__ import annotations

from fractions import Fraction
import math


P00 = ((7, -1), (1, 7))
P01 = ((5, -5), (5, 5))
P10 = ((3, -8), (4, 6))
P11 = ((0, -10), (5, 0))
O = ((Fraction(4, 5), Fraction(-3, 5)),
     (Fraction(3, 5), Fraction(4, 5)))


def multiply(left, right):
    return tuple(tuple(sum(left[i][k] * right[k][j] for k in range(2))
                       for j in range(2)) for i in range(2))


def determinant(matrix):
    return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]


def gram(matrix):
    return tuple(tuple(sum(matrix[k][i] * matrix[k][j] for k in range(2))
                       for j in range(2)) for i in range(2))


def smith(matrix):
    first = math.gcd(*(abs(value) for row in matrix for value in row))
    return first, abs(determinant(matrix)) // first


def chi4(vector):
    a, b = vector
    norm = a * a + b * b
    return Fraction(a**4 - 6 * a * a * b * b + b**4, norm * norm)


def reduce_tau(real: Fraction, imag: Fraction):
    """Exact PSL2Z reduction, with reflection x -> |x| identified."""
    for _ in range(100):
        nearest = math.floor(real + Fraction(1, 2))
        real -= nearest
        radius = real * real + imag * imag
        if radius < 1 or (radius == 1 and real < 0):
            real, imag = -real / radius, imag / radius
            continue
        return abs(real), imag
    raise RuntimeError("tau reduction did not terminate")


def hnf_classes(n: int):
    output = {}
    for a in range(1, n + 1):
        if n % a:
            continue
        d = n // a
        for b in range(a):
            modular = reduce_tau(Fraction(b, a), Fraction(d, a))
            topology = (math.gcd(a, math.gcd(b, d)), 0)
            topology = (topology[0], n // topology[0])
            output.setdefault(modular, set()).add(topology)
    return output


def has_factorial(n: int) -> bool:
    groups = hnf_classes(n)
    keys = list(groups)
    return any(len(groups[keys[i]] & groups[keys[j]]) >= 2
               for i in range(len(keys)) for j in range(i + 1, len(keys)))


def verify() -> dict:
    assert all(determinant(matrix) == 50 for matrix in (P00, P01, P10, P11))
    assert gram(P00) == gram(P01) == ((50, 0), (0, 50))
    assert gram(P10) == gram(P11) == ((25, 0), (0, 100))
    assert multiply(O, P00) == P01
    assert multiply(O, P10) == P11
    assert smith(P00) == smith(P10) == (1, 50)
    assert smith(P01) == smith(P11) == (5, 10)
    square_delta = chi4((7, 1)) - chi4((5, 5))
    rectangle_delta = chi4((3, 4)) - chi4((0, 5))
    assert square_delta == Fraction(1152, 625)
    assert rectangle_delta == -square_delta
    assert not any(has_factorial(n) for n in range(2, 50))
    assert has_factorial(50)
    return {
        "determinant": 50,
        "smith": {"cyclic": [1, 50], "noncyclic": [5, 10]},
        "gram": {"i": [[50, 0], [0, 50]], "2i": [[25, 0], [0, 100]]},
        "delta_chi4": {"i": str(square_delta), "2i": str(rectangle_delta)},
        "minimal_determinant": 50,
    }


if __name__ == "__main__":
    print(verify())
