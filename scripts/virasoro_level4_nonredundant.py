#!/usr/bin/env python3
"""Check that the h=5/8 level-4 quasiprimary is not null or a total derivative."""

from fractions import Fraction as F

# Level-4 PBW basis:
# [L_-4, L_-3 L_-1, L_-2^2, L_-2 L_-1^2, L_-1^4]
Q4 = [F(-9), F(-60), F(40), F(0), F(0)]

# Descendants of chi2=(L_-2-2/3 L_-1^2)|h> at level 4.
null_columns = [
    [F(0), F(0), F(1), F(-2, 3), F(0)],
    [F(2), F(2), F(0), F(1), F(-2, 3)],
]

# Image of L_-1 acting on level-3 PBW basis
# [L_-3, L_-2 L_-1, L_-1^3].
total_derivative_columns = [
    [F(2), F(1), F(0), F(0), F(0)],
    [F(0), F(1), F(0), F(1), F(0)],
    [F(0), F(0), F(0), F(0), F(1)],
]


def rank(matrix_columns):
    if not matrix_columns:
        return 0
    rows = [list(row) for row in zip(*matrix_columns)]
    m, n = len(rows), len(rows[0])
    r = 0
    for c in range(n):
        pivot = next((i for i in range(r, m) if rows[i][c]), None)
        if pivot is None:
            continue
        rows[r], rows[pivot] = rows[pivot], rows[r]
        pv = rows[r][c]
        rows[r] = [x / pv for x in rows[r]]
        for i in range(m):
            if i != r and rows[i][c]:
                f = rows[i][c]
                rows[i] = [rows[i][j] - f * rows[r][j] for j in range(n)]
        r += 1
        if r == m:
            break
    return r


def in_span(vector, columns):
    return rank(columns) == rank(columns + [vector])


def main():
    assert not in_span(Q4, null_columns)
    assert not in_span(Q4, total_derivative_columns + null_columns)
    print("Q4 is non-null modulo chi2 descendants")
    print("Q4 is not in Im(L_-1) modulo the null submodule")
    print("Therefore the candidate is not eliminated as a pure total derivative.")


if __name__ == "__main__":
    main()
