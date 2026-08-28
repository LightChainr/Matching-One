#!/usr/bin/env python3
"""Power-budget calculations for orientation-difference discovery runs.

This is deliberately a planning tool, not a statistical inference package.
It extrapolates an observed standard error with the iid 1/sqrt(n) law and a
user-specified signal power.  Use it to reject obviously underpowered lattice
sizes before spending compute.

Examples
--------
Observed root-gap SE 5.98e-5 with 1e6 replicas at physical L~17.  If the
orientation root gap is expected to be C L^-4 with C~1:

    python scripts/orientation_power_budget.py \
        --baseline-L 17 --baseline-samples 1000000 --baseline-se 5.98e-5 \
        --signal-amplitude 1 --signal-power 4 --target-z 3 \
        --sizes 7 8 12 17 29 41 99

For fixed-p matching-function differences use signal-power 13/4 instead.
"""

from __future__ import annotations

import argparse
import math


def projected_se(baseline_se: float, baseline_samples: float, samples: float) -> float:
    return baseline_se * math.sqrt(baseline_samples / samples)


def expected_signal(amplitude: float, length: float, power: float) -> float:
    return abs(amplitude) * length ** (-power)


def required_samples(
    *,
    baseline_se: float,
    baseline_samples: float,
    signal: float,
    target_z: float,
) -> float:
    if signal <= 0:
        return math.inf
    target_se = signal / target_z
    return baseline_samples * (baseline_se / target_se) ** 2


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-L", type=float, required=True)
    parser.add_argument("--baseline-samples", type=float, required=True)
    parser.add_argument("--baseline-se", type=float, required=True)
    parser.add_argument("--signal-amplitude", type=float, required=True)
    parser.add_argument("--signal-power", type=float, required=True)
    parser.add_argument("--target-z", type=float, default=3.0)
    parser.add_argument("--sizes", type=float, nargs="+", required=True)
    parser.add_argument(
        "--variance-power",
        type=float,
        default=0.0,
        help=(
            "optional empirical scaling Var(single-replica estimator)~L^q. "
            "Default q=0; positive q makes large sizes even more expensive"
        ),
    )
    args = parser.parse_args()

    if min(
        args.baseline_L,
        args.baseline_samples,
        args.baseline_se,
        args.signal_amplitude,
        args.signal_power,
        args.target_z,
        *args.sizes,
    ) <= 0:
        raise SystemExit("all lengths/counts/scales/powers must be positive")

    print(
        "L          signal        SE@baseline_n    z@baseline_n    "
        f"samples_for_z{args.target_z:g}"
    )
    for length in args.sizes:
        signal = expected_signal(args.signal_amplitude, length, args.signal_power)
        # Baseline SE transported to the new size using the optional
        # single-replica variance scaling law.
        se_at_baseline_n = args.baseline_se * (length / args.baseline_L) ** (
            args.variance_power / 2
        )
        z_at_baseline = signal / se_at_baseline_n
        n_req = required_samples(
            baseline_se=se_at_baseline_n,
            baseline_samples=args.baseline_samples,
            signal=signal,
            target_z=args.target_z,
        )
        print(
            f"{length:8.3f}  {signal:12.5e}  {se_at_baseline_n:14.5e}  "
            f"{z_at_baseline:12.5g}  {n_req:16.6g}"
        )

    exponent = 2 * args.signal_power + args.variance_power
    print()
    print(
        "Asymptotic sample-count scaling under these assumptions: "
        f"n ~ L^{exponent:g}."
    )
    print(
        "If per-replica connectivity work is O(L^2), total work scales near "
        f"L^{exponent + 2:g}."
    )
    print(
        "Caveat: improved common-random-number coupling can make the variance "
        "of a *difference* decrease with L; measure that empirically and pass "
        "a negative --variance-power if justified."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
