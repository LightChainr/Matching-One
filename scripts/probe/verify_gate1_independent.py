#!/usr/bin/env python3
"""Verify the P398 ground facts independently (Gate-1 style).

Checks, on widths 4..8:

1. state counts equal Catalan(w): 14, 42, 132, 429, 1430;
2. join/detach maps close inside the state space and are idempotent;
3. every reflection fixes exactly C(w, floor(w/2)) states (Burnside match is
   therefore necessary but never identifying -- this is #598's caveat);
4. exactly one reflection preserves the D0 dictionary plus all four sources:
   the cut reflection i -> (w-1) - i;
5. the coarsest exact strong lumping admissible for D0 and stable for J and D
   separately is *block-for-block equal* to the R-orbit partition;
6. the reported r_positive numbers 10, 26, 76, 232, 750 are reproduced;
7. readout parity classification, including the odd-width odd readout
   ``halves_linked``.
"""

from __future__ import annotations

import json
import math
import sys
import time
from fractions import Fraction
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parent))

from p398core import (  # noqa: E402
    ALL_READOUTS,
    Generator,
    HELD_OUT,
    PRIMARY,
    apply_permutation,
    canonical_rgs,
    detach_point,
    join_adjacent,
    orbit_labels,
    readout_functions,
    reflection_perm,
    source_vectors,
    state_permutation,
    states,
)

sys.setrecursionlimit(20000)

REPORTED_R_POSITIVE = {4: 10, 5: 26, 6: 76, 7: 232, 8: 750}


# --------------------------------------------------------------------------
# exact partition refinement


def observable_labels(g: Generator, names: Sequence[str]) -> List[int]:
    fns = readout_functions(g.width)
    signature: Dict[Tuple, int] = {}
    labels = []
    for s in g.states:
        key = tuple(fns[n](s) for n in names)
        if key not in signature:
            signature[key] = len(signature)
        labels.append(signature[key])
    return labels


def coarsest_lumping(
    g: Generator, rate_sets: Sequence[Dict[Tuple[str, int], Fraction]], labels: List[int]
) -> List[int]:
    """Coarsest strong lumping refining ``labels``, stable for each rate set."""
    operators = [g.rows_exact(rs) for rs in rate_sets]
    labels = list(labels)
    while True:
        signature: Dict[Tuple, int] = {}
        refined = []
        for source in range(g.size):
            key: List = [labels[source]]
            for rows in operators:
                totals: Dict[int, Fraction] = {}
                for dest, value in rows[source]:
                    block = labels[dest]
                    totals[block] = totals.get(block, Fraction(0)) + value
                key.append(tuple(sorted(totals.items())))
            frozen = tuple(key)
            if frozen not in signature:
                signature[frozen] = len(signature)
            refined.append(signature[frozen])
        if refined == labels:
            return labels
        labels = refined


def rate_join(g: Generator) -> Dict[Tuple[str, int], Fraction]:
    return {m: (Fraction(1) if m[0] == "join" else Fraction(0)) for m in g.moves}


def rate_detach(g: Generator) -> Dict[Tuple[str, int], Fraction]:
    return {m: (Fraction(1) if m[0] == "detach" else Fraction(0)) for m in g.moves}


def canonical_labels(labels: Sequence[int]) -> Tuple[int, ...]:
    return canonical_rgs(list(labels))


# --------------------------------------------------------------------------


def main() -> None:
    report: Dict = {}
    for w in range(4, 9):
        t0 = time.time()
        st = states(w)
        g = Generator(w)
        report.setdefault("widths", {})[w] = {
            "state_count": len(st),
            "catalan": math.comb(2 * w, w) // (w + 1),
            "moves": len(g.moves),
        }
        assert len(st) == math.comb(2 * w, w) // (w + 1), f"state count fails at w={w}"

        # moves close + idempotent
        for p in range(w):
            for s in st:
                j = join_adjacent(s, p)
                d = detach_point(s, p)
                assert j in st and d in st
                assert join_adjacent(j, p) == j
                assert detach_point(d, p) == d

        # reflection fixed-point census
        perm = reflection_perm(w, w - 1)
        pi = state_permutation(g, perm)
        fixed = sum(1 for i in range(g.size) if pi[i] == i)
        orbits = (g.size + fixed) // 2
        central = math.comb(w, w // 2)
        assert fixed == central
        report["widths"][w]["cut_reflection"] = {
            "fixed": fixed,
            "central_binomial": central,
            "orbits": orbits,
            "reported_r_positive": REPORTED_R_POSITIVE[w],
            "orbit_matches_reported": orbits == REPORTED_R_POSITIVE[w],
        }
        # every reflection same fixed count
        fix_counts = {
            sum(1 for i in range(g.size) if state_permutation(g, reflection_perm(w, k))[i] == i)
            for k in range(w)
        }
        assert fix_counts == {central}

        # which reflections preserve D0 readouts and sources
        fns = readout_functions(w)
        vals = {n: [f(s) for s in g.states] for n, f in fns.items()}
        srcs = source_vectors(g)
        d0_preserving, all_preserving = [], []
        for k in range(w):
            pk = state_permutation(g, reflection_perm(w, k))
            preserved = [n for n in ALL_READOUTS if all(vals[n][pk[i]] == vals[n][i] for i in range(g.size))]
            src_preserved = [n for n, v in srcs.items() if all(v[pk[i]] == v[i] for i in range(g.size))]
            if set(PRIMARY) <= set(preserved) and len(src_preserved) == len(srcs):
                d0_preserving.append(k)
            if len(preserved) == len(ALL_READOUTS):
                all_preserving.append(k)
        report["widths"][w]["reflections_preserving_D0_plus_sources"] = d0_preserving
        report["widths"][w]["reflections_preserving_all_readouts"] = all_preserving
        assert d0_preserving == [w - 1], f"task does not name the cut reflection at w={w}"

        # Gate-1: coarsest admissible lumping == orbit partition?
        initial = observable_labels(g, PRIMARY)
        coarsest = coarsest_lumping(g, (rate_join(g), rate_detach(g)), initial)
        orbits_labels = orbit_labels(g, pi)
        # is the orbit partition admissible (stable)?
        orbit_again = coarsest_lumping(g, (rate_join(g), rate_detach(g)), list(orbits_labels))
        orbit_admissible = canonical_labels(orbit_again) == canonical_labels(orbits_labels)
        equal = canonical_labels(coarsest) == canonical_labels(orbits_labels)
        report["widths"][w]["gate1"] = {
            "initial_blocks": len(set(initial)),
            "coarsest_blocks": len(set(coarsest)),
            "orbit_blocks": len(set(orbits_labels)),
            "orbit_admissible": orbit_admissible,
            "coarsest_equals_orbit_partition": equal,
        }
        assert equal, f"Gate-1 partition equality fails at w={w}"

        # readout parity defect under R (relative L2 defect from even)
        even_defects = {}
        for n in ALL_READOUTS:
            v = vals[n]
            scale = math.sqrt(sum(x * x for x in v))
            defect = math.sqrt(sum((v[i] - v[pi[i]]) ** 2 for i in range(g.size))) / (2 * scale) if scale else 0.0
            even_defects[n] = defect
        report["widths"][w]["readout_even_defects"] = {
            n: round(d, 12) for n, d in even_defects.items()
        }
        report["widths"][w]["odd_readouts"] = [n for n, d in even_defects.items() if d > 1e-9]
        report["widths"][w]["seconds"] = round(time.time() - t0, 2)

    # widths 9/10 state counts (enumeration only, quick)
    for w in (9, 10):
        t0 = time.time()
        n = len(states(w))
        report.setdefault("widths", {})[w] = {
            "state_count": n,
            "catalan": math.comb(2 * w, w) // (w + 1),
            "seconds": round(time.time() - t0, 2),
        }
        assert n == math.comb(2 * w, w) // (w + 1)

    print(json.dumps(report, indent=2, sort_keys=True))
    verdict = all(
        report["widths"][w]["gate1"]["coarsest_equals_orbit_partition"] for w in range(4, 9)
    ) and all(report["widths"][w]["cut_reflection"]["orbit_matches_reported"] for w in range(4, 9))
    print("GATE1_VERDICT:", "PASS" if verdict else "FAIL")


if __name__ == "__main__":
    main()
