#!/usr/bin/env python3
"""dpfloor: shared definitions.

Object under test (identical for every path, by construction):
  * square-site NN percolation on a periodic L1 x L2 (here: n x n) torus,
    sliced along a primitive direction u=(a,b) into a cylinder of `width` sites
    whose physical circumference is  ell = n * |u|.
  * G4 = black NN graph, parameter p (black density)
  * G8 = white graph = NN + both diagonals ("matching"), parameter q = 1-p
  * SAFE transfer R^0_{G}: the row transfer restricted to states whose frontier
    partition has NO component carrying non-zero horizontal (deck) gain
    (i.e. no occupied cluster wraps the cylinder).  Perron root lambda0, and
    I0_G(x) = -log lambda0_G(x).
  * charge-coexistence root  p_root :  I0_G4(p) - I0_G8(1-p) = 0      (Delta_w(p)=0)
  * Omega := -(p_root - p_c) * ell^4
"""
from __future__ import annotations
import math

# The single, common p_c used by EVERY path in this task.
PC = 0.59274605079210

# also carry the OTHER value that appears in the shipped oblique file, so the
# convention split can be quantified.
PC_OBLIQUE_FILE = 0.59274605079

# geometry table: tag, direction u, n, ell=n|u|, and the engine-cylinder width
# (only axis directions have a plain-cylinder realisation in sector802_lib).
GEOMS = [
    # tag              u      n    ell
    ("axis_n2",       (1, 0), 2,  2.0),
    ("axis_n3",       (1, 0), 3,  3.0),
    ("axis_n4",       (1, 0), 4,  4.0),
    ("axis_n5",       (1, 0), 5,  5.0),
    ("axis_n6",       (1, 0), 6,  6.0),
    ("axis_n7",       (1, 0), 7,  7.0),
    ("axis_n8",       (1, 0), 8,  8.0),          # <-- the l=8 case of the 2.1e-12 story
    ("axis_n9",       (1, 0), 9,  9.0),
    ("diag_n4",       (1, 1), 4,  4.0 * math.sqrt(2.0)),
    ("diag_n5",       (1, 1), 5,  5.0 * math.sqrt(2.0)),
    ("slope21_n3",    (2, 1), 3,  3.0 * math.sqrt(5.0)),
    ("slope21_n4",    (2, 1), 4,  4.0 * math.sqrt(5.0)),
    ("slope31_n3",    (3, 1), 3,  3.0 * math.sqrt(10.0)),
    ("slope32_n2",    (3, 2), 2,  2.0 * math.sqrt(13.0)),
    ("slope52_n2",    (5, 2), 2,  2.0 * math.sqrt(29.0)),
]

GEOM_BY_TAG = {g[0]: g for g in GEOMS}


def cos4(direction):
    """cos(4*theta) for the angle theta of the primitive direction."""
    a, b = direction
    den = (a * a + b * b) ** 2
    return (a ** 4 - 6 * a * a * b * b + b ** 4) / den


def ell4(tag):
    return GEOM_BY_TAG[tag][3] ** 4


def omega(p_root, tag, pc=PC):
    """Omega := -(p_root - p_c) * ell^4  for the given geometry tag."""
    return -(p_root - pc) * ell4(tag)


def bisect(f, lo, hi, xtol=1e-16, maxiter=200):
    """Plain bisection with a hard absolute bracket tolerance.

    Returns (root, width, nfev, f(root)).
    """
    flo, fhi = f(lo), f(hi)
    if flo == 0.0:
        return lo, 0.0, 1, 0.0
    if fhi == 0.0:
        return hi, 0.0, 1, 0.0
    if (flo > 0) == (fhi > 0):
        raise ValueError("no sign change on [%r,%r]: f=%r,%r" % (lo, hi, flo, fhi))
    nfev = 2
    a, b = lo, hi
    for _ in range(maxiter):
        if b - a < xtol:
            break
        c = 0.5 * (a + b)
        fc = f(c)
        nfev += 1
        if fc == 0.0:
            a = b = c
            break
        if (fc > 0) == (flo > 0):
            a, flo = c, fc
        else:
            b, fhi = c, fc
    return 0.5 * (a + b), abs(b - a), nfev, f(0.5 * (a + b))


def solve_root(f, lo=0.5, hi=0.8, xtol=1e-16):
    """Root of f on [lo,hi] with a hard absolute tolerance.

    Uses scipy brentq when available (same tolerance, ~15 evaluations instead of
    ~52), otherwise plain bisection.  Every path uses this same helper so the
    root-finder contribution is identical across paths.
    """
    try:
        from scipy.optimize import brentq
    except Exception:
        return bisect(f, lo, hi, xtol=xtol)
    flo, fhi = f(lo), f(hi)
    if (flo > 0) == (fhi > 0):
        raise ValueError("no sign change")
    r = brentq(f, lo, hi, xtol=xtol, rtol=4.0 * 2.0 ** -52, maxiter=300)
    return r, float(xtol), -1, f(r)
