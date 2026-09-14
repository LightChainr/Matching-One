#!/usr/bin/env python3
"""Check the magnetic 5/48 gap in oblique cylinder coordinates.

For primitive u=(a,b), a Bezout transverse step has physical normal height
1/|u|.  A transfer excitation I per t-step therefore corresponds to physical
energy |u|*I, while the circumference is ell=n|u|.  The dimensionless magnetic
gap is

    n*|u|^2*I.

At criticality / the nearby charge root this should tend to 2*pi*(5/48).
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from scipy.optimize import brentq

from oblique_charge_transfer import ObliqueSafeTransfer


def run(a: int, b: int, n: int, state_cap: int) -> dict[str, object]:
    u = (a, b)
    g4 = ObliqueSafeTransfer(n, u, matching=False, state_cap=state_cap)
    g8 = ObliqueSafeTransfer(n, u, matching=True, state_cap=state_cap)

    def equation(p: float) -> float:
        return math.log(g4.lambda0(p)) - math.log(g8.lambda0(1.0 - p))

    root = brentq(equation, 0.5, 0.8, xtol=3e-11)
    lam = g4.lambda0(root)
    I0 = -math.log(lam)
    scaled = n * (a * a + b * b) * I0
    target = 2 * math.pi * 5 / 48
    return {
        "direction": [a, b],
        "n": n,
        "charge_root_p": root,
        "I0_per_transverse_integer_step": I0,
        "scaled_physical_gap_n_normu2_I0": scaled,
        "target_2pi_5_over_48": target,
        "ratio_to_target": scaled / target,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--a", type=int, required=True)
    parser.add_argument("--b", type=int, required=True)
    parser.add_argument("--n", type=int, required=True)
    parser.add_argument("--state-cap", type=int, default=300000)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = run(args.a, args.b, args.n, args.state_cap)
    text = json.dumps(result, indent=2, sort_keys=True)
    if args.output:
        args.output.write_text(text + "\n", encoding="utf-8")
    print(text)


if __name__ == "__main__":
    main()
