#!/usr/bin/env python3
"""Certify the algebraic complexity of the Ziff 2006 "A lattice" bond threshold.

This exists because the P2 manuscript motivated its search class with the claim
that every exactly-known planar percolation threshold has degree at most 6 and
coefficient height at most 3.  That claim is false.  Ziff's generalized
cell/dual-cell paper reports an exact bond threshold for the "A lattice" as the
root in [0,1] of

    p^5 - 4p^4 + 3p^3 + 2p^2 - 1,

which is irreducible over Q and has height 4.  This module certifies the degree,
the height, the irreducibility and the root, so that the corrected statement in
the manuscript is generated rather than asserted.

On sourcing.  The primary text could not be read from this environment -- the
publisher, arXiv and the indexing services are all unreachable through the
network policy in use.  What was obtained is the polynomial and the decimal value
from two independent search-engine summaries of the publisher and preprint
records.  Those two are checked against each other here: the polynomial's only
root in (0,1) is computed exactly and compared with the separately quoted decimal
0.625457.  A transcription error in either would almost certainly break that
agreement, so the two corroborate one another.  This is recorded as
CORROBORATED_NOT_PRIMARY and must be replaced by a primary reading before the
manuscript is submitted.
"""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from math import gcd
from pathlib import Path
from typing import Sequence

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "results" / "ziff-a-lattice-complexity" / "latest.json"
SCHEMA = "matching-one/ziff-a-lattice-complexity/v1"

# p^5 - 4p^4 + 3p^3 + 2p^2 - 1, ascending.
POLYNOMIAL = (-1, 0, 2, 3, -4, 1)
QUOTED_DECIMAL = "0.625457"
ISOLATION_BITS = 120


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _text(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def multiply(left: Sequence[int], right: Sequence[int]) -> list[int]:
    out = [0] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            out[i + j] += a * b
    return out


def exact_divide(dividend: Sequence[int], monic_divisor: Sequence[int]) -> list[int] | None:
    """Divide by a monic integer polynomial; return the quotient or None."""
    work = list(dividend)
    span = len(work) - len(monic_divisor)
    if span < 0:
        return None
    quotient = [0] * (span + 1)
    for k in range(span, -1, -1):
        lead = work[k + len(monic_divisor) - 1]
        quotient[k] = lead
        for j, value in enumerate(monic_divisor):
            work[k + j] -= lead * value
    return quotient if all(value == 0 for value in work) else None


def monic_integer_factorizations(coefficients: Sequence[int]) -> list[dict]:
    """Every factorization into monic integer polynomials of degree 1..3.

    The polynomial is monic with constant term -1, so any monic integer factor
    has constant term dividing 1.  That makes the search finite and small, and a
    degree-5 polynomial factors non-trivially only if it has a factor of degree
    at most 2 or exactly 3, all of which are covered here.
    """
    found = []
    span = max(abs(value) for value in coefficients) + 2
    for constant in (1, -1):
        divisors = [[constant, 1]]
        divisors += [[constant, b, 1] for b in range(-span, span + 1)]
        divisors += [[constant, b, c, 1]
                     for b in range(-span, span + 1)
                     for c in range(-span, span + 1)]
        for divisor in divisors:
            quotient = exact_divide(coefficients, divisor)
            if quotient is not None:
                found.append({"factor": divisor, "cofactor": quotient})
    return found


def build_result() -> dict:
    import sys

    sys.path.insert(0, str(ROOT / "scripts"))
    from exact_polynomial_root_certificate import isolate_roots  # noqa: E402

    coefficients = list(POLYNOMIAL)
    divisor = 0
    for value in coefficients:
        divisor = gcd(divisor, abs(value))
    _require(divisor == 1, "polynomial is not primitive")

    brackets = isolate_roots(
        [Fraction(v) for v in coefficients], Fraction(0), Fraction(1), bits=ISOLATION_BITS
    )
    _require(len(brackets) == 1, "expected exactly one root in (0,1)")
    low, high = brackets[0]

    quoted = Fraction(QUOTED_DECIMAL)
    places = len(QUOTED_DECIMAL.split(".", 1)[1])
    half_ulp = Fraction(1, 2 * 10**places)
    agrees = (quoted - half_ulp) <= high and low <= (quoted + half_ulp)

    factorizations = monic_integer_factorizations(coefficients)

    return {
        "schema": SCHEMA,
        "claim_level": "C1_exact_arithmetic_on_a_sourced_polynomial",
        "lattice": "Ziff 2006 'A lattice', bond percolation",
        "polynomial_ascending": coefficients,
        "polynomial_text": "p^5 - 4p^4 + 3p^3 + 2p^2 - 1",
        "degree": len(coefficients) - 1,
        "height": max(abs(value) for value in coefficients),
        "primitive": True,
        "irreducible_over_Q": not factorizations,
        "monic_integer_factorizations_found": factorizations,
        "root_in_unit_interval": {
            "isolating_bracket": [_text(low), _text(high)],
            "isolation_bits": ISOLATION_BITS,
        },
        "cross_check_against_the_quoted_decimal": {
            "quoted": QUOTED_DECIMAL,
            "quoted_window": [_text(quoted - half_ulp), _text(quoted + half_ulp)],
            "agrees_with_the_isolated_root": agrees,
            "why_this_matters": (
                "the polynomial and the decimal were taken from independent "
                "summaries; a transcription error in either would almost "
                "certainly break this agreement"
            ),
        },
        "consequence_for_the_manuscript": (
            "degree 5 and height 4, so the claim that every exactly-known planar "
            "threshold has degree <= 6 and height <= 3 is false. The true "
            "historical range, on what we can verify, is degree <= 6 and "
            "height <= 4, and results/pslq-degree6-height4-* censuses that class"
        ),
        "verification_status": "CORROBORATED_NOT_PRIMARY",
        "what_is_still_owed": (
            "a reading of Ziff, Phys. Rev. E 73, 016134 (2006) itself. The "
            "publisher, arXiv and the indexing services are unreachable through "
            "the network policy in use here, so this was corroborated from two "
            "independent search summaries rather than read. It must be checked "
            "against the primary text before submission."
        ),
    }


def validate_result(result) -> dict:
    expected = build_result()
    if result != expected:
        raise ValueError("ziff a-lattice complexity does not exactly reproduce")
    return {"schema": SCHEMA, "status": "valid", "height": expected["height"]}


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    arguments = parser.parse_args(argv)
    text = json.dumps(build_result(), indent=2, sort_keys=True) + "\n"
    destination = arguments.output or DEFAULT_OUTPUT
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(text, encoding="utf-8")
    print(f"wrote {destination.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
