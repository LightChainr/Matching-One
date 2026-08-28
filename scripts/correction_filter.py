#!/usr/bin/env python3
"""Generate finite-size correction-annihilation weights.

Two modes are supported.

zero:
    Build a homogeneous filter for a root condition.  For sizes L_i,
    choose weights w_i such that

        sum_i w_i = 0
        sum_i w_i L_i**(-q) = 0    for each q in --powers

    The last weight is fixed before a harmless overall normalization.
    With no powers and two sizes this is the ordinary first difference.
    Applied to S_L(p)=L^(13/4) M_L(p), this recovers the structure of the
    Mertens-Ziff two-size condition.  Adding powers creates higher
    annihilators.

limit:
    Build Richardson-style weights for an asymptotic limit E_inf:

        sum_i w_i = 1
        sum_i w_i L_i**(-q) = 0    for each q in --powers.

All arithmetic uses mpmath.  The program reports cancellation residuals
and simple weight-amplification diagnostics; it does not claim that the
chosen correction powers are physically correct.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Sequence

import mpmath as mp


def parse_mpf(text: str) -> mp.mpf:
    try:
        return mp.mpf(text)
    except (ValueError, TypeError) as exc:
        raise argparse.ArgumentTypeError(f"invalid real value: {text!r}") from exc


def fmt(x: mp.mpf, digits: int = 30) -> str:
    return mp.nstr(x, digits, strip_zeros=False)


def _basis(size: mp.mpf, power: mp.mpf) -> mp.mpf:
    return mp.power(size, -power)


def zero_weights(sizes: Sequence[mp.mpf], powers: Sequence[mp.mpf]) -> list[mp.mpf]:
    expected = len(powers) + 2
    if len(sizes) != expected:
        raise ValueError(
            f"zero mode with {len(powers)} correction powers needs exactly "
            f"{expected} sizes, got {len(sizes)}"
        )

    # Fix the final weight to one, then solve the remaining square system.
    rows: list[list[mp.mpf]] = [[mp.mpf(1) for _ in sizes[:-1]]]
    rhs: list[mp.mpf] = [-mp.mpf(1)]
    for power in powers:
        rows.append([_basis(size, power) for size in sizes[:-1]])
        rhs.append(-_basis(sizes[-1], power))

    solved = mp.lu_solve(mp.matrix(rows), mp.matrix(rhs))
    weights = [solved[i] for i in range(len(sizes) - 1)] + [mp.mpf(1)]

    # Root equations are scale invariant. Normalize to max |w| = 1.
    scale = max(abs(weight) for weight in weights)
    if scale == 0:
        raise ValueError("degenerate all-zero filter")
    return [weight / scale for weight in weights]


def limit_weights(sizes: Sequence[mp.mpf], powers: Sequence[mp.mpf]) -> list[mp.mpf]:
    expected = len(powers) + 1
    if len(sizes) != expected:
        raise ValueError(
            f"limit mode with {len(powers)} correction powers needs exactly "
            f"{expected} sizes, got {len(sizes)}"
        )

    rows: list[list[mp.mpf]] = [[mp.mpf(1) for _ in sizes]]
    rhs: list[mp.mpf] = [mp.mpf(1)]
    for power in powers:
        rows.append([_basis(size, power) for size in sizes])
        rhs.append(mp.mpf(0))

    solved = mp.lu_solve(mp.matrix(rows), mp.matrix(rhs))
    return [solved[i] for i in range(len(sizes))]


def residuals(
    sizes: Sequence[mp.mpf],
    weights: Sequence[mp.mpf],
    powers: Sequence[mp.mpf],
    mode: str,
) -> dict[str, mp.mpf]:
    out: dict[str, mp.mpf] = {}
    target = mp.mpf(0) if mode == "zero" else mp.mpf(1)
    out["constant"] = mp.fsum(weights) - target
    for power in powers:
        out[f"L^-{fmt(power, 12)}"] = mp.fsum(
            weight * _basis(size, power)
            for size, weight in zip(sizes, weights)
        )
    return out


def diagnostics(weights: Sequence[mp.mpf]) -> dict[str, mp.mpf]:
    return {
        "l1": mp.fsum(abs(weight) for weight in weights),
        "l2": mp.sqrt(mp.fsum(weight * weight for weight in weights)),
        "max_abs": max(abs(weight) for weight in weights),
        "min_abs_nonzero": min(abs(weight) for weight in weights if weight != 0),
    }


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--mode", choices=("zero", "limit"), required=True)
    p.add_argument("--sizes", nargs="+", type=parse_mpf, required=True)
    p.add_argument(
        "--powers",
        nargs="*",
        type=parse_mpf,
        default=[],
        help="correction powers q to annihilate as L^-q",
    )
    p.add_argument("--dps", type=int, default=80)
    p.add_argument("--json", type=Path, default=None)
    p.add_argument(
        "--next-power",
        type=parse_mpf,
        default=None,
        help="optional first surviving correction power for a heuristic root-bias exponent",
    )
    p.add_argument(
        "--derivative-power",
        type=parse_mpf,
        default=None,
        help="optional k if dS/dp scales as L^k; use with --next-power",
    )
    return p


def main() -> int:
    args = parser().parse_args()
    if args.dps < 40:
        raise SystemExit("--dps must be at least 40")
    mp.mp.dps = args.dps

    sizes = list(args.sizes)
    powers = list(args.powers)
    if any(size <= 0 for size in sizes):
        raise SystemExit("all sizes must be positive")
    if len(set(map(str, sizes))) != len(sizes):
        raise SystemExit("sizes must be distinct")
    if any(power <= 0 for power in powers):
        raise SystemExit("all correction powers must be positive")

    try:
        if args.mode == "zero":
            weights = zero_weights(sizes, powers)
        else:
            weights = limit_weights(sizes, powers)
    except (ValueError, ZeroDivisionError) as exc:
        raise SystemExit(str(exc)) from exc

    cancels = residuals(sizes, weights, powers, args.mode)
    diag = diagnostics(weights)

    print(f"mode: {args.mode}")
    print("size                 weight")
    for size, weight in zip(sizes, weights):
        print(f"{fmt(size, 16):>16}  {fmt(weight, 32)}")

    print("\ncancellation residuals")
    for name, value in cancels.items():
        print(f"{name:>18}: {fmt(value, 12)}")

    print("\nweight diagnostics")
    for name, value in diag.items():
        print(f"{name:>18}: {fmt(value, 16)}")

    heuristic = None
    if args.next_power is not None or args.derivative_power is not None:
        if args.mode != "zero":
            raise SystemExit(
                "--next-power/--derivative-power heuristic is defined only for zero mode"
            )
        if args.next_power is None or args.derivative_power is None:
            raise SystemExit("provide both --next-power and --derivative-power")
        heuristic = args.next_power + args.derivative_power
        print(
            "\nheuristic: if the first surviving critical residual is "
            f"L^-{fmt(args.next_power, 12)} and dS/dp ~ L^{fmt(args.derivative_power, 12)}, "
            f"expect root bias roughly L^-{fmt(heuristic, 12)}."
        )
        print("This is an asymptotic scaling heuristic, not a fitted result.")

    if args.json is not None:
        payload = {
            "mode": args.mode,
            "dps": args.dps,
            "sizes": [fmt(value) for value in sizes],
            "powers": [fmt(value) for value in powers],
            "weights": [fmt(value) for value in weights],
            "residuals": {key: fmt(value) for key, value in cancels.items()},
            "diagnostics": {key: fmt(value) for key, value in diag.items()},
            "heuristic_root_bias_power": None if heuristic is None else fmt(heuristic),
        }
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
