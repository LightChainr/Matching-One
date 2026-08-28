#!/usr/bin/env python3
"""Generate Pell approximants to the hexagonal torus and E4 spin-filter diagnostics.

The period matrix uses columns

    v1 = (2m, 0),  v2 = (m, x)

with x^2 - 3 m^2 = 1, hence

    tau = 1/2 + i x/(2m) -> exp(i*pi/3).

For a modular-covariant spin-4 level-4 amplitude proportional to E4(tau),
the simple zero E4(exp(i*pi/3))=0 supplies an additional O(m^-2)
suppression along this sequence.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import mpmath as mp


def pell_plus(count: int, skip_fundamental: bool = True) -> list[tuple[int, int]]:
    """Return solutions (x,m) of x^2-3m^2=1 in increasing m."""
    if count <= 0:
        return []
    x, m = 2, 1
    out: list[tuple[int, int]] = []
    while len(out) < count + int(skip_fundamental):
        out.append((x, m))
        x, m = 2 * x + 3 * m, x + 2 * m
    return out[1:] if skip_fundamental else out[:count]


def sigma3(n: int) -> int:
    total = 0
    d = 1
    while d * d <= n:
        if n % d == 0:
            q = n // d
            total += d**3
            if q != d:
                total += q**3
        d += 1
    return total


def eisenstein_e4(tau: mp.mpc, terms: int = 80) -> mp.mpc:
    q = mp.e ** (2 * mp.pi * 1j * tau)
    return 1 + 240 * mp.fsum(sigma3(n) * q**n for n in range(1, terms + 1))


def record(x: int, m: int, e4_i: mp.mpc, terms: int) -> dict[str, object]:
    tau = mp.mpf(1) / 2 + 1j * mp.mpf(x) / (2 * m)
    rho = mp.mpf(1) / 2 + 1j * mp.sqrt(3) / 2
    e4 = eisenstein_e4(tau, terms)
    ratio = e4 / e4_i
    return {
        "x": x,
        "m": m,
        "pell_residual": x * x - 3 * m * m,
        "period_matrix_rows": [[2 * m, m], [0, x]],
        "site_count_det": 2 * m * x,
        "tau_real": mp.nstr(mp.re(tau), 40),
        "tau_imag": mp.nstr(mp.im(tau), 40),
        "tau_distance_to_hex": mp.nstr(abs(tau - rho), 30),
        "E4": mp.nstr(e4, 40),
        # Re(tau)=1/2 makes E4(tau)/E4(i) real.  The truncated q-series can
        # retain an O(10^-dps) imaginary roundoff, which must not leak into a
        # field whose machine-readable contract is a real scalar.
        "E4_over_E4_i": mp.nstr(mp.re(ratio), 40),
        "m2_times_E4_ratio": mp.nstr(mp.re(m * m * ratio), 40),
    }


def build(count: int, dps: int, terms: int) -> dict[str, object]:
    mp.mp.dps = dps
    rho = mp.mpf(1) / 2 + 1j * mp.sqrt(3) / 2
    e4_i = eisenstein_e4(1j, terms)
    e4_rho = eisenstein_e4(rho, terms)
    rows = [record(x, m, e4_i, terms) for x, m in pell_plus(count)]
    return {
        "schema": "hexagonal Pell modular spin filter v1",
        "hexagonal_modulus": {
            "real": "0.5",
            "imag": mp.nstr(mp.sqrt(3) / 2, 40),
            "E4_numeric": mp.nstr(e4_rho, 40),
        },
        "square_modulus_E4": mp.nstr(e4_i, 40),
        "pell_equation": "x^2 - 3 m^2 = 1",
        "periods": "v1=(2m,0), v2=(m,x)",
        "records": rows,
        "frozen_scaling_hypotheses": {
            "H4_level4": {
                "generic_matching_residual": "L^(-13/4)",
                "pell_hex_matching_residual": "L^(-21/4)",
                "generic_root_bias": "L^(-4)",
                "pell_hex_root_bias": "L^(-6)",
                "reason": "spin-4 one-point is forbidden at tau=exp(i*pi/3); E4 has a simple zero and Pell shape error is O(L^-2)",
            },
            "same_exponent_H12_alias": {
                "pell_hex_matching_residual": "L^(-13/4)",
                "pell_hex_root_bias": "L^(-4)",
                "reason": "spin 12 is invariant under the hexagonal torus 60-degree automorphism",
            },
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--count", type=int, default=4)
    parser.add_argument("--dps", type=int, default=70)
    parser.add_argument("--terms", type=int, default=100)
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()
    payload = build(args.count, args.dps, args.terms)
    text = json.dumps(payload, indent=2) + "\n"
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(text, encoding="utf-8")
    else:
        print(text, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
