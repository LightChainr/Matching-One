"""Ordered distinct-quartet products, via the partition-lattice identity."""
from functools import lru_cache
from itertools import permutations, product
from math import factorial
from fractions import Fraction
import numpy as np


@lru_cache(None)
def partitions(m):
    if m == 0:
        return ((),)
    rows = []
    for previous in partitions(m-1):
        rows.append(previous+((m-1,),))
        for i in range(len(previous)):
            rows.append(previous[:i]+(previous[i]+(m-1,),)+previous[i+1:])
    return tuple(rows)


def distinct_product(*arrays):
    """Arrays have (prefix,quartet,...); trailing dimensions may broadcast."""
    variables = np.broadcast_arrays(*arrays)
    m, q = len(variables), variables[0].shape[1]
    if q < m:
        raise ValueError('Fewer quartets than the U-statistic order')
    block_sums = {}
    answer = np.zeros_like(variables[0].sum(axis=1), dtype=float)
    for partition in partitions(m):
        coefficient = (-1)**(m-len(partition))
        term = 1.
        for block in partition:
            coefficient *= factorial(len(block)-1)
            if block not in block_sums:
                value = variables[block[0]]
                for j in block[1:]:
                    value = value*variables[j]
                block_sums[block] = value.sum(axis=1)
            term = term*block_sums[block]
        answer += coefficient*term
    return answer/(factorial(q)//factorial(q-m))


def shape_estimators(b, h):
    """b: P,Q,O,5; h: P,Q,S,O,5; moments C,W,CC,CW,WW."""
    c, w = b[:, :, None, :, 0], b[:, :, None, :, 1]
    hc, hw = h[..., 0], h[..., 1]
    cc, cw, ww = h[..., 2], h[..., 3], h[..., 4]
    cov = np.stack((cc.mean(1)-2*distinct_product(c, hc),
                    cw.mean(1)-distinct_product(c, hw)-distinct_product(w, hc),
                    ww.mean(1)-2*distinct_product(w, hw)), axis=-1)
    square_cc = (distinct_product(cc, cc)-4*distinct_product(cc, c, hc)
                 +4*distinct_product(c, c, hc, hc))
    square_ww = (distinct_product(ww, ww)-4*distinct_product(ww, w, hw)
                 +4*distinct_product(w, w, hw, hw))
    square_cw = (distinct_product(cw, cw)-2*distinct_product(cw, c, hw)
                 -2*distinct_product(cw, w, hc)+distinct_product(c, c, hw, hw)
                 +distinct_product(w, w, hc, hc)+2*distinct_product(c, w, hc, hw))
    squares = np.stack((square_cc, square_cw, square_ww), axis=-1)
    energy = (squares*np.array([1., 2., 1.])).sum(axis=(1, 2, 3))
    return cov, energy, squares


def check_algebra():
    data = np.array([[1., 2., 3., 4.], [2., -1., 4., 1.], [-2., 1., 2., 3.],
                     [3., 0., 1., -1.], [1., -3., 2., 0.], [0., 1., -2., 4.],
                     [2., 2., -1., 1.], [-1., 3., 0., 2.]])
    maximum = 0.
    for q in (4, 8):
        for m in (2, 3, 4):
            direct = sum(np.prod([data[index[j], j] for j in range(m)])
                         for index in permutations(range(q), m))/(factorial(q)//factorial(q-m))
            fast = distinct_product(*(data[None, :q, j] for j in range(m)))[0]
            maximum = max(maximum, abs(fast-direct))
    assert maximum < 1e-12

    # Exhaust all iid 4-quartet outcomes of a two-point finite distribution.
    # Columns: bC,bW,hC,hW,hCC,hCW,hWW; each atom has probability 1/2.
    atoms = ((1, 2, 1, -1, 3, 2, 4), (2, -1, -2, 2, 1, -3, 2))
    means = [Fraction(atoms[0][j]+atoms[1][j], 2) for j in range(7)]
    c, w, hc, hw, cc, cw, ww = means
    true_shapes = (cc-2*c*hc, cw-c*hw-w*hc, ww-2*w*hw)
    expected = sum(v*v*weight for v, weight in zip(true_shapes, (1, 2, 1)))
    rational_total = Fraction(0)
    min_energy, max_difference = float('inf'), 0.
    for outcome in product(range(2), repeat=4):
        rows = [atoms[i] for i in outcome]

        def exact(*columns):
            m = len(columns)
            return sum((Fraction(np.prod([rows[ix[j]][columns[j]] for j in range(m)]))
                        for ix in permutations(range(4), m)), Fraction(0))/(factorial(4)//factorial(4-m))

        ecc = exact(4, 4)-4*exact(4, 0, 2)+4*exact(0, 0, 2, 2)
        eww = exact(6, 6)-4*exact(6, 1, 3)+4*exact(1, 1, 3, 3)
        ecw = exact(5, 5)-2*exact(5, 0, 3)-2*exact(5, 1, 2)+exact(0, 0, 3, 3)+exact(1, 1, 2, 2)+2*exact(0, 1, 2, 3)
        exact_energy = ecc+2*ecw+eww
        rational_total += exact_energy/16
        x = np.array(rows, dtype=float)
        b = np.zeros((1, 4, 1, 5)); h = np.zeros((1, 4, 1, 1, 5))
        b[0, :, 0, :2] = x[:, :2]
        h[0, :, 0, 0] = x[:, 2:]
        observed = shape_estimators(b, h)[1][0]
        max_difference = max(max_difference, abs(observed-float(exact_energy)))
        min_energy = min(min_energy, observed)
    assert rational_total == expected
    assert max_difference < 1e-10
    assert min_energy < 0  # An unbiased estimate of a nonnegative target need not be positive.
    return {'ordered_product_vs_permutations_max_error': maximum,
            'exact_finite_distribution_expected_energy': str(expected),
            'exact_average_estimator_energy': str(rational_total),
            'energy_float_vs_rational_max_error': max_difference,
            'negative_energy_example': min_energy,
            'finite_distribution_outcomes_enumerated': 16}
