#!/usr/bin/env python3
"""Exact harmonic projectors for the four primitive N=1105 Gaussian orientations.

The four orientations are (33,4), (32,9), (31,12), (24,23).  For a square-lattice
observable truncated as

    X(theta) = A0 + A4 cos(4 theta) + A8 cos(8 theta) + A12 cos(12 theta),

this script constructs exact rational weights recovering A0, A4, A8 and A12.
No percolation outcome is read.
"""

from fractions import Fraction
from math import sqrt

N = 1105
ORIENTATIONS = ((33, 4), (32, 9), (31, 12), (24, 23))


def gaussian_mul(z, w):
    a, b = z
    c, d = w
    return (a * c - b * d, a * d + b * c)


def gaussian_pow(z, power):
    result = (1, 0)
    base = z
    while power:
        if power & 1:
            result = gaussian_mul(result, base)
        base = gaussian_mul(base, base)
        power >>= 1
    return result


def cos_4m(a, b, m):
    real, _ = gaussian_pow((a, b), 4 * m)
    return Fraction(real, N ** (2 * m))


def solve_fraction(matrix, rhs):
    a = [list(map(Fraction, row)) + [Fraction(v)] for row, v in zip(matrix, rhs)]
    n = len(a)
    for col in range(n):
        pivot = next(row for row in range(col, n) if a[row][col])
        a[col], a[pivot] = a[pivot], a[col]
        scale = a[col][col]
        a[col] = [v / scale for v in a[col]]
        for row in range(n):
            if row == col:
                continue
            scale = a[row][col]
            if scale:
                a[row] = [x - scale * y for x, y in zip(a[row], a[col])]
    return [a[row][-1] for row in range(n)]


def projector(target_index):
    # Conditions: sum_i w_i basis_j(theta_i) = delta_{j,target_index}
    basis = [
        [Fraction(1) for _ in ORIENTATIONS],
        [cos_4m(*orientation, 1) for orientation in ORIENTATIONS],
        [cos_4m(*orientation, 2) for orientation in ORIENTATIONS],
        [cos_4m(*orientation, 3) for orientation in ORIENTATIONS],
    ]
    rhs = [Fraction(int(j == target_index)) for j in range(4)]
    return solve_fraction(basis, rhs)


def response(weights, harmonic_index):
    if harmonic_index == 0:
        return sum(weights, Fraction(0))
    return sum(
        (weight * cos_4m(*orientation, harmonic_index) for weight, orientation in zip(weights, ORIENTATIONS)),
        Fraction(0),
    )


def norms(weights):
    l1 = sum(abs(w) for w in weights)
    l2 = sqrt(float(sum(w * w for w in weights)))
    return float(l1), l2


def main():
    labels = ("H0", "H4", "H8", "H12")
    for idx, label in enumerate(labels):
        weights = projector(idx)
        l1, l2 = norms(weights)
        print(label)
        for orientation, weight in zip(ORIENTATIONS, weights):
            print(f"  {orientation}: {weight} = {float(weight):+.15f}")
        print(f"  l1={l1:.15f} l2={l2:.15f}")
        for j, check_label in enumerate(labels):
            assert response(weights, j) == Fraction(int(j == idx))
        print()

    h0 = projector(0)
    # Higher harmonics are not killed automatically; report the first few leakages.
    for m in range(4, 9):
        value = response(h0, m)
        print(f"H{4*m} leakage: {value} = {float(value):+.15f}")


if __name__ == "__main__":
    main()
