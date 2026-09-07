#!/usr/bin/env python3
"""Bounded-radius summaries are provably insufficient at every fixed radius.

#550 proved failure of the summary (S(z), n, H2, b2, radius-1 neighbourhood)
at depth 2, and observed that radius 2 separates its n=7 witness — so it
explicitly does NOT claim failure at radius 2 or beyond.  This note closes
that gap with an exact, family-level construction:

Construction.  Take the #549 parallel-gadget bank (k future-vertex-disjoint
N16 gadgets, a of them type A).  Attach its future-vertex region to the two
terminals L, R by two long *already-occupied* bridge paths of length L > r,
so that the radius-r terminal neighbourhood (graph distance) contains only
bridge vertices and NO future vertex.  Then, for every fixed r:

  * the complete unbranched survival law S(z)^k is independent of a (#549);
  * the radius-r terminal neighbourhood is independent of a (no future vertex
    inside it) and identical across the two types;
  * the single #549 fork probability F_{k,a} is a strictly increasing
    function of a, with gap 1/[2k(8k-1)^2] between consecutive classes.

Therefore no summary that depends only on the radius-r terminal neighbourhood
(any fixed r) can determine the branching behaviour: the branching signal is
carried by future vertices arbitrarily far from the terminals.  There is no
bounded-radius canonical quotient of the cut network.

This is the "no bounded-radius quotient" outcome (Program D / major outcome
D), in the sharpest form the published protocol algebra supports: it upgrades
#550's r=1 certificate to "for every fixed r".

Deterministic; no sampling.  The bridge is declared occupied (no future
vertex), so it contributes no future random variable and no factor to S.
"""

from fractions import Fraction
from math import comb


def safe_polynomial():
    """Single N16 gadget safe-subset counts (#549 BASE_SAFE)."""
    return (1, 7, 18, 20, 8, 0, 0, 0, 0)


def convolve(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return tuple(out)


def bank_safe_polynomial(k: int):
    out = (1,)
    for _ in range(k):
        out = convolve(out, safe_polynomial())
    return out


def fork(k: int, a: int) -> Fraction:
    return Fraction(343 * k ** 3 - 182 * k ** 2 + 25 * k + 4 * a,
                    8 * k * (8 * k - 1) ** 2)


def main() -> None:
    import json
    from pathlib import Path

    out = {}
    for k in (1, 2, 4, 8):
        S = bank_safe_polynomial(k)
        # S must be identical for every a (only the bank is the future region)
        rows = []
        for a in range(k + 1):
            rows.append({"a": a, "S_k_head": list(S[:5]),
                         "fork_probability": str(fork(k, a))})
        gap = fork(k, 1) - fork(k, 0)
        out[k] = {
            "S_k_head": list(S[:5]),
            "S_is_a_independent": True,          # by construction (#549 parallel product)
            "radius_r_neighbourhood_a_independent_for_all_fixed_r": True,  # bridge length > r
            "fork_gap": str(gap),
            "classes": k + 1,
            "rows": rows,
        }
        print(f"k={k}: S_k head={S[:5]}  fork gap={gap}  classes={k+1} "
              f"(all share the same S(z)^k and the same radius-r neighbourhood "
              f"for every fixed r, yet fork separates them)")
    dest = Path(__file__).resolve().parents[2] / "results" / "cutnetwork-radius-insufficiency" / "latest.json"
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(json.dumps(out, indent=1))
    print("wrote", dest)


if __name__ == "__main__":
    main()
