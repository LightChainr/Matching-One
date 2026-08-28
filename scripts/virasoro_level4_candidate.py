#!/usr/bin/env python3
"""Exact rational check of the c=0, h=5/8 level-4 spin-4 candidate."""

from fractions import Fraction as F

# PBW level-4 basis:
# [L_-4, L_-3 L_-1, L_-2^2, L_-2 L_-1^2, L_-1^4]
# L1 maps it to [L_-3, L_-2 L_-1, L_-1^3].
L1 = [
    [F(5), F(5,4), F(3), F(0), F(0)],
    [F(0), F(4), F(6), F(9,2), F(0)],
    [F(0), F(0), F(0), F(3), F(17)],
]

# Level-2 null state: chi2 = (L_-2 - 2/3 L_-1^2)|h>.
# Its level-4 descendants L_-2 chi2 and L_-1^2 chi2.
NULL4 = [
    [F(0), F(0), F(1), F(-2,3), F(0)],
    [F(2), F(2), F(0), F(1), F(-2,3)],
]

# Integer-normalized quasiprimary candidate.
Q4 = [F(-9), F(-60), F(40), F(0), F(0)]


def matvec(matrix, vector):
    return [sum(row[j] * vector[j] for j in range(len(vector))) for row in matrix]


def in_null_span(vector):
    # Solve using first two components. NULL4[0][0:2]=(0,0),
    # NULL4[1][0:2]=(2,2), so any vector in the span must have v0=v1.
    return vector[0] == vector[1]


def main():
    assert matvec(L1, Q4) == [F(0), F(0), F(0)]
    assert not in_null_span(Q4)

    h = F(5,8)
    x = 2*h + 4
    spin = 4
    assert x == F(21,4)

    print("c=0, h=5/8")
    print("level-2 null: L_-2 - (2/3)L_-1^2")
    print("Q4 = 40 L_-2^2 - 60 L_-3 L_-1 - 9 L_-4")
    print("L1 Q4 = 0; Q4 is not in the level-2 null descendant span")
    print(f"bulk weights: (37/8,5/8) or (5/8,37/8); x={x}, spin=+/-{spin}")


if __name__ == "__main__":
    main()
