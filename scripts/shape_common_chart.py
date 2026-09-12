"""A necessary, coordinate-covariant test of one common symmetrizing chart.

R_i are decreasing involutions of the same interval. Their composition test
requires no fitted polynomial. Zero is necessary, not sufficient.
"""
from typing import Callable, TypeVar
T=TypeVar('T')


def reflection_commutator(r1: Callable[[T], T], r2: Callable[[T], T],
                          r3: Callable[[T], T], p: T):
    # (R1 R2)(R2 R3) = R1 R3 exactly; opposite order has four factors.
    left=r1(r3(p))
    right=r2(r3(r1(r2(p))))
    return {'left':left,'right':right,'difference':left-right}
