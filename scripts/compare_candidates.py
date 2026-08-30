#!/usr/bin/env python3
"""Compare common closed forms against the accepted square-site p_c."""

from __future__ import annotations

import math

PC = 0.59274605079210

CANDIDATES = {
    "(sqrt(5)-1)/2": (math.sqrt(5) - 1) / 2,
    "1/sqrt(3)": 1 / math.sqrt(3),
    "2/pi": 2 / math.pi,
    "pi/sqrt(28)": math.pi / math.sqrt(28),
    "sqrt(pi)/3": math.sqrt(math.pi) / 3,
    "sin(pi/5)": math.sin(math.pi / 5),
    "1/sqrt(e)": 1 / math.sqrt(math.e),
    "2*sin(pi/18)": 2 * math.sin(math.pi / 18),
    "1-2*sin(pi/18)": 1 - 2 * math.sin(math.pi / 18),
}


def main() -> None:
    print(f"p_c = {PC:.14f}")
    print(f"{'name':22} {'value':18} {'pc-value':14} {'rel_ppm':10}")
    for name, value in CANDIDATES.items():
        delta = PC - value
        print(f"{name:22} {value:18.12f} {delta:14.6e} {delta / PC * 1e6:10.3f}")


if __name__ == "__main__":
    main()
