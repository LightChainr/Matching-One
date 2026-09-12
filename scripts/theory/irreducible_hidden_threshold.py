#!/usr/bin/env python3
"""An irreducible product-chain no-go, not a percolation simulation.

Row generators act on terminal observables. B=e0 tensor 1; C=e0^T tensor
uniform initial law. Algebra proves the response; small mpmath checks are
numerical diagnostics, never the evidence for an all-L theorem.
"""
from __future__ import annotations
import json
from mpmath import mp


def generator(length, p, center):
    if length < 2:
        raise ValueError("length must be at least two")
    a, b = mp.exp(p-center), mp.exp(center-p)
    g = mp.zeros(2*length)
    for hand in range(2):
        for i in range(length):
            row = hand*length+i
            g[row, (1-hand)*length+i] = 2
            if i+1 < length:
                g[row, row+1] = a
            if i:
                g[row, row-1] = b
            g[row, row] = -mp.fsum(g[row, j] for j in range(2*length) if j != row)
    return g


def check():
    with mp.workdps(40):
        error = mp.mpf(0)
        count = 0
        for length in (2, 3, 5):
            B = mp.matrix([1]*length+[0]*length)
            C = mp.matrix([[mp.mpf(1)/length]*length+[0]*length])
            for p in (mp.mpf(0), mp.mpf('0.5'), mp.mpf(1)):
                for center in (mp.mpf(1)/2, mp.mpf(1)/3):
                    g = generator(length, p, center)
                    assert mp.norm(g*mp.ones(2*length, 1)) < mp.mpf('1e-35')
                    # All rails in both directions and every rung are positive.
                    for hand in range(2):
                        for i in range(length-1):
                            j = hand*length+i
                            assert g[j, j+1] > 0 and g[j+1, j] > 0
                    for i in range(length):
                        assert g[i, i+length] == g[i+length, i] == 2
                    for t in (mp.mpf(0), mp.mpf('0.2'), mp.mpf(1)):
                        got = (C*mp.expm(t*g)*B)[0]
                        want = (1+mp.exp(-4*t))/2
                        error = max(error, abs(got-want))
                        count += 1
        assert error < mp.mpf('1e-35')
        return {
            'scope': 'Abstract irreducible Markov control; not percolation',
            'response_checks': count,
            'precision_decimal_digits': 40,
            'max_response_error_diagnostic': mp.nstr(error, 8),
            'exact_response': '(1+exp(-4*t))/2',
            'finite_gap': '2*cosh(p-c)-2*cos(pi/L), L>=2, c=1/2 or 1/3',
            'limiting_gap': '4*sinh((p-c)/2)^2',
            'exact_task_order_including_constant_mode': 2,
            'all_checked_graphs_connected': True,
        }


if __name__ == '__main__':
    print(json.dumps(check(), indent=2, allow_nan=False))
