#!/usr/bin/env python3
"""Configuration-level Euler–Poincaré / Betti lift of the matching identity.

On the implemented square tori the occupied cell-complex Euler characteristic
is configuration-identical to a Betti/wrapping combination

    chi = V - E + F0 = beta0_black - beta0_white - q

which is the same statement as the P34 identity

    C_black - C_white = q + V - E + F0.

Locked consequences, all exhaustive on ``named_tiny_geometries``:

- wrapping-difference channels remain configuration-identical;
- ``q`` is a wrapping-event difference in ``{-1,0,+1}``, not the homology-rank
  difference ``r_black - r_white`` (empty mask: ``q=-1``, ``r_black=0``,
  ``r_white=2``);
- cyclomatic numbers satisfy ``kappa_black = E_primal - V + beta0_black >=
  r_black`` and the matching-white analogue;
- every proposed Betti control reduces to the existing ``(V,E,F0,q,C)``
  statistics, so the variance-reduction branch of issue #111 closes.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

from euler_motif_controls import (
    cluster_difference,
    configuration_identity,
    count_motifs,
    named_tiny_geometries,
)
from integer_period_torus import IntegerTorusGeometry, classify_configuration


# Exhaustive counts of configurations on which q equals r_black - r_white.
# The complementary majority is the reason q is not a homology-rank difference.
Q_EQUALS_RANK_DIFFERENCE = {
    ("axis", 2): (4, 16),
    ("axis", 3): (162, 512),
    ("gaussian-2-1", 0): (10, 32),
    ("diamond", 2): (68, 256),
}


@dataclass(frozen=True)
class BettiRecord:
    mask: int
    chi: int
    q: int
    beta0_black: int
    beta0_white: int
    r_black: int
    r_white: int
    kappa_black: int
    kappa_white: int
    wrapping: dict[str, int]
    residual: int

    @property
    def q_equals_rank_difference(self) -> bool:
        return self.q == self.r_black - self.r_white


def occupied_edges(
    geometry: IntegerTorusGeometry, active: Sequence[bool], *, matching: bool
) -> int:
    edges = geometry.matching_edges if matching else geometry.primal_edges
    return sum(1 for edge in edges if active[edge.i] and active[edge.j])


def betti_record(
    geometry: IntegerTorusGeometry, active: Sequence[bool], mask: int = -1
) -> BettiRecord:
    identity = configuration_identity(geometry, active, mask)
    white = [not value for value in active]
    black_wrap, black_components = classify_configuration(geometry, active)
    white_wrap, white_components = classify_configuration(
        geometry, white, matching=True
    )
    beta0_black = len(black_components)
    beta0_white = len(white_components)
    cluster_black, cluster_white = cluster_difference(geometry, active)
    if (cluster_black, cluster_white) != (beta0_black, beta0_white):
        raise AssertionError("homology component counts disagree with cluster_stats")
    motifs = count_motifs(geometry, active)
    chi = motifs["V"] - motifs["E"] + motifs["F0"]
    residual = chi - (beta0_black - beta0_white - identity.q)
    kappa_black = motifs["E"] - motifs["V"] + beta0_black
    kappa_white = (
        occupied_edges(geometry, white, matching=True)
        - (geometry.n - motifs["V"])
        + beta0_white
    )
    return BettiRecord(
        mask=mask,
        chi=chi,
        q=identity.q,
        beta0_black=beta0_black,
        beta0_white=beta0_white,
        r_black=black_wrap.max_rank,
        r_white=white_wrap.max_rank,
        kappa_black=kappa_black,
        kappa_white=kappa_white,
        wrapping=identity.wrapping,
        residual=residual,
    )


def exhaustive_betti(geometry: IntegerTorusGeometry) -> dict[str, object]:
    if geometry.n > 16:
        raise ValueError("exhaustive Betti identity is limited to N<=16")
    samples = 1 << geometry.n
    identity_failures = 0
    wrapping_not_identical = 0
    cyclomatic_failures = 0
    q_equals_rank = 0
    q_values = set()
    empty = None
    full = None
    for mask in range(samples):
        active = [bool((mask >> vertex) & 1) for vertex in range(geometry.n)]
        record = betti_record(geometry, active, mask)
        q_values.add(record.q)
        if len(set(record.wrapping.values())) != 1:
            wrapping_not_identical += 1
        if record.residual != 0:
            identity_failures += 1
        if record.kappa_black < record.r_black or record.kappa_white < record.r_white:
            cyclomatic_failures += 1
        if record.q_equals_rank_difference:
            q_equals_rank += 1
        if mask == 0:
            empty = record
        if mask == samples - 1:
            full = record
    if empty is None or full is None:
        raise AssertionError("empty/full configurations were not recorded")
    expected = Q_EQUALS_RANK_DIFFERENCE[(geometry.name, geometry.L)]
    return {
        "name": geometry.name,
        "N": geometry.n,
        "L": geometry.L,
        "physical_period": geometry.physical_period,
        "configurations": samples,
        "q_values": sorted(q_values),
        "identity_failures": identity_failures,
        "wrapping_not_identical": wrapping_not_identical,
        "cyclomatic_failures": cyclomatic_failures,
        "q_equals_rank_difference": q_equals_rank,
        "q_equals_rank_difference_locked": list(expected),
        "empty": _boundary_payload(empty),
        "full": _boundary_payload(full),
        "passed": (
            identity_failures == 0
            and wrapping_not_identical == 0
            and cyclomatic_failures == 0
            and q_equals_rank == expected[0]
            and samples == expected[1]
        ),
    }


def _boundary_payload(record: BettiRecord) -> dict[str, int]:
    return {
        "mask": record.mask,
        "chi": record.chi,
        "q": record.q,
        "beta0_black": record.beta0_black,
        "beta0_white": record.beta0_white,
        "r_black": record.r_black,
        "r_white": record.r_white,
        "kappa_black": record.kappa_black,
        "kappa_white": record.kappa_white,
        "residual": record.residual,
    }


def run_betti_suite() -> dict[str, object]:
    rows = [exhaustive_betti(geometry) for geometry in named_tiny_geometries()]
    return {
        "schema": "issue-111 configuration Euler/Betti identity v1",
        "identity": "chi = V - E + F0 = beta0_black - beta0_white - q",
        "equivalent_p34": "C_black - C_white = q + V - E + F0",
        "q_is": "wrapping-event difference, configuration-identical across channels",
        "q_is_not": "homology-rank difference r_black - r_white",
        "variance_reduction_branch": "closed",
        "variance_reduction_reason": (
            "Betti statistics reduce to the existing P34 controls V, E, F0, q, C"
        ),
        "claim_level": "C5",
        "passed": all(row["passed"] for row in rows),
        "exhaustive": rows,
    }


def render_report(payload: dict[str, object]) -> str:
    lines = [
        "# Configuration Euler–Poincaré / Betti identity",
        "",
        "Source: `scripts/configuration_betti_identity.py`.",
        "Claim level: C5 finite identity. Equivalent to P34; not a new control variate.",
        "",
        "On every enumerated square torus configuration",
        "",
        "```text",
        "chi = V - E + F0 = beta0_black - beta0_white - q",
        "```",
        "",
        "The wrapping-difference variable `q` is the common event difference in",
        "`{-1,0,+1}` already used by P34. It is not the homology-rank difference",
        "`r_black - r_white`. Cyclomatic numbers bound wrapping rank from above:",
        "",
        "```text",
        "kappa_black = E_primal - V_black + beta0_black >= r_black",
        "kappa_white = E_matching_white - V_white + beta0_white >= r_white",
        "```",
        "",
        "The issue-#111 variance-reduction branch therefore closes: every proposed",
        "Betti control is an algebraic rewrite of `(V, E, F0, q, C_black, C_white)`.",
        "",
        "## Exhaustive tiny quotients",
        "",
        "| geometry | N | configs | identity fail | cyclo fail | `q = r_b-r_w` | empty `(q,r_b,r_w)` |",
        "|---|---:|---:|---:|---:|---:|---|",
    ]
    for row in payload["exhaustive"]:
        empty = row["empty"]
        lines.append(
            "| {name} | {N} | {configurations} | {identity_failures} | "
            "{cyclomatic_failures} | {q_equals_rank_difference}/{configurations} | "
            "`({q},{r_black},{r_white})` |".format(
                q=empty["q"],
                r_black=empty["r_black"],
                r_white=empty["r_white"],
                **row,
            )
        )
    lines.extend(
        [
            "",
            "Empty-mask counterexample (every listed quotient): `q=-1`, `r_black=0`,",
            "`r_white=2`, `beta0_black=0`, `beta0_white=1`. The Euler identity still",
            "holds because `chi = 0 = 0 - 1 - (-1)`.",
            "",
            "## Boundary",
            "",
            "This does not introduce a production Newman–Ziff Betti statistic and does",
            "not claim a variance reduction relative to P34 motif controls. Averaging",
            "the configuration identity recovers the expected Mertens–Ziff relation",
            "already used in P34.",
            "",
        ]
    )
    return "\n".join(lines)


def json_ready(value: object) -> object:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {str(key): json_ready(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_ready(item) for item in value]
    return value


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exact-json", type=Path)
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()
    payload = run_betti_suite()
    status = "PASS" if payload["passed"] else "FAIL"
    print("Betti identity suite " + status)
    for row in payload["exhaustive"]:
        print(
            "{name} N={N}: identity_failures={identity_failures} "
            "cyclomatic_failures={cyclomatic_failures} "
            "q_equals_rank={q_equals_rank_difference}/{configurations}".format(**row)
        )
    if args.exact_json is not None:
        args.exact_json.parent.mkdir(parents=True, exist_ok=True)
        args.exact_json.write_text(
            json.dumps(json_ready(payload), indent=2) + "\n", encoding="utf-8"
        )
        print("wrote " + str(args.exact_json))
    if args.report is not None:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(render_report(payload), encoding="utf-8")
        print("wrote " + str(args.report))
    return 0 if payload["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
