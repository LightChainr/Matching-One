#!/usr/bin/env python3
"""Finite-duration Markov-pulse control for the exact P398 double-pulse result.

Requires mpmath. No sampling. Matrix exponentials use physical G+epsilon*H,
never the signed operator H as a standalone Markov generator. Numerical values
are diagnostics, not rigorous interval certificates. Historical results are
never overwritten by this command.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from mpmath import mp
from p398_double_pulse_exact import build_model


def report() -> dict:
    with mp.workdps(45):
        model = build_model(4)
        n = len(model['states'])
        g, h = mp.matrix(model['G']), mp.matrix(model['H'])
        source = mp.matrix([[mp.mpf(x) / n for x in model['sources_scaled'][2]]])
        readout = mp.matrix([x[2] for x in model['F']])
        tau, duration = mp.mpf(1) / 4, mp.mpf(1) / 32
        waiting = mp.expm(tau * g)
        # An independent Frechet derivative series, not a finite difference:
        # D_k = d/d epsilon (G+epsilon H)^k |0.
        power, derivative = mp.eye(n), mp.zeros(n)
        window_derivative, weight = mp.zeros(n), mp.mpf(1)
        for k in range(1, 36):
            derivative = g * derivative + h * power
            power = g * power
            weight *= duration / k
            window_derivative += weight * derivative
        truth = (source * window_derivative * waiting * window_derivative * readout)[0]
        ideal = (source * h * waiting * h * readout)[0]
        rows = []
        for denominator in (8, 16, 32):
            epsilon = mp.mpf(1) / denominator
            plus = mp.expm(duration * (g + epsilon * h))
            minus = mp.expm(duration * (g - epsilon * h))
            pp = (source * plus * waiting * plus * readout)[0]
            pm = (source * plus * waiting * minus * readout)[0]
            mp_value = (source * minus * waiting * plus * readout)[0]
            mm = (source * minus * waiting * minus * readout)[0]
            mixed = (pp - pm - mp_value + mm) / (4 * epsilon**2)
            rows.append(dict(
                epsilon=f'1/{denominator}',
                centered_mixed_difference=mp.nstr(mixed, 35),
                finite_window_derivative_error=mp.nstr(mixed - truth, 25),
                parity_pair_max_error=mp.nstr(max(abs(pp - mm), abs(pm - mp_value)), 15),
            ))
        return dict(
            schema='matching-one.p398-double-pulse-numerical-control.v1',
            standing='Numerical independent matrix-exponential control, not an exact certificate or simulation.',
            dps=45, source='delta_wrapped_pair', readout='wrap', tau='1/4',
            pulse_duration='1/32', frechet_series_terms=35,
            finite_window_mixed_derivative=mp.nstr(truth, 35),
            short_pulse_kernel=mp.nstr(ideal, 35),
            finite_window_derivative_divided_by_delta_squared=mp.nstr(truth / duration**2, 35),
            ladder=rows,
        )


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path)
    args = parser.parse_args()
    text = json.dumps(report(), indent=2, allow_nan=False) + '\n'
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        with args.out.open('x', encoding='utf-8') as handle:
            handle.write(text)
    else:
        print(text, end='')
