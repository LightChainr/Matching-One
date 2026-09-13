#!/usr/bin/env python3
"""Issue #593 Phase B Deliverables 2 and 3 at widths 9 and 10.

Runs the repository's own ``width_experiment`` transport-scoring pipeline on a
generator shim that is bit-exact against the repository at widths 1..8, so the
frozen configuration (seeds, span ladder, lag grid, thresholds) is the
repository's code, not a reimplementation.

One substitution, declared: ``r_linear`` at widths 9/10 is computed at full
budget by ``p593_phaseB_rlinear_true`` (exact mod-p elimination, numpy
echelon) and injected into the repository's rank cache; D1/D2 stay
budget-truncated at 150 exactly as the repository does at width 8.

Deliverable 3 (out-of-pencil lumping collapse) is computed here for widths
4..10 with the repository's ``exact_lumping`` against the baseline generator
plus the ``single_point_join`` tilt.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import p398_intervention_transport as pit
from p593_noncrossing_fast import WideGenerator, noncrossing_states_fast
from p593_phaseB_rlinear_true import krylov_dimension_numpy, sparse_mod_matrix

OUTPUT = Path(__file__).resolve().parents[1] / "results" / "p593-width9-10-memory-degree" / "raw"

DICTIONARY_SEEDS = {
    "D0_additive_local_counts": ("blocks", "singletons", "wrap"),
    "D1_plus_size_and_extent": (
        "blocks", "singletons", "wrap", "max_block", "linked_pairs", "boundary_span",
    ),
    "D2_plus_nonlocal_topology": (
        "blocks", "singletons", "wrap", "max_block", "linked_pairs",
        "boundary_span", "halves_linked", "covering_depth",
    ),
}


def observable_vector(states, name):
    size = len(states)
    if name == "one":
        return [1] * size
    if name == "blocks":
        return [max(s) + 1 for s in states]
    if name == "wrap":
        return [1 if s[0] == s[-1] else 0 for s in states]
    if name == "singletons":
        return [sum(1 for b in _blocks(s) if len(b) == 1) for s in states]
    if name == "max_block":
        return [max(len(b) for b in _blocks(s)) for s in states]
    if name == "linked_pairs":
        return [sum(len(b) * (len(b) - 1) // 2 for b in _blocks(s)) for s in states]
    if name == "boundary_span":
        return [sum(max(b) - min(b) for b in _blocks(s)) for s in states]
    if name == "halves_linked":
        half = len(states[0]) // 2
        return [
            1 if any(
                any(p < half for p in b) and any(p >= half for p in b) for b in _blocks(s)
            ) else 0
            for s in states
        ]
    if name == "covering_depth":
        out = []
        for s in states:
            blocks = [sorted(b) for b in _blocks(s)]
            best = 0
            for point in range(len(s)):
                depth = sum(1 for b in blocks if b[0] < point < b[-1] and point not in b)
                best = max(best, depth)
            out.append(best)
        return out
    raise KeyError(name)


def _blocks(state):
    out = {}
    for point, label in enumerate(state):
        out.setdefault(label, []).append(point)
    return out.values()


def prepopulate_linear_cache(generator, budget=None):
    """Fill the repository's rank cache with true D0 and capped D1/D2 values."""
    budget = budget if budget is not None else pit.linear_rank_budget(generator.size)
    sparse = sparse_mod_matrix(generator)
    injected = {}
    for name, members in DICTIONARY_SEEDS.items():
        vectors = [observable_vector(generator.states, member) for member in members]
        cap = generator.size if name == "D0_additive_local_counts" else budget
        dim, depth = krylov_dimension_numpy(sparse, generator.size, vectors, budget=cap)
        key = (generator.width, name, budget)
        pit._LINEAR_RANK_CACHE[key] = {
            "dimension": dim,
            "arithmetic": f"exact modulo the prime {pit.RANK_PRIME}",
            "budget": cap,
            "reached_the_whole_state_space": dim >= generator.size,
            "limited_by_budget": dim >= cap and cap < generator.size,
            "is_a_certified_lower_bound_on_the_rational_rank": True,
        }
        injected[name] = dim
    return injected


def lumping_checks(generator):
    baseline = generator.baseline_rates()
    colour = [
        tuple(
            observable_vector(generator.states, member)[state]
            for member in DICTIONARY_SEEDS["D0_additive_local_counts"]
        )
        for state in range(generator.size)
    ]
    joint = pit.exact_lumping(generator, colour, [baseline, generator.coefficients("uniform_join_minus_detach")])
    baseline_only = pit.exact_lumping(generator, colour, [baseline])
    out_of_pencil = pit.exact_lumping(
        generator, colour, [baseline, generator.coefficients("single_point_join")]
    )
    return {
        "r_positive_baseline_lumping": baseline_only["blocks"],
        "r_positive_whole_affine_family": joint["blocks"],
        "out_of_pencil_lumping_blocks": out_of_pencil["blocks"],
        "out_of_pencil_collapses_to_identity": out_of_pencil["blocks"] == generator.size,
        "states": generator.size,
    }


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    mode = sys.argv[1]
    widths = [int(a) for a in sys.argv[2:]]
    if mode == "lumping":
        for width in widths:
            generator = WideGenerator(width, noncrossing_states_fast(width))
            result = lumping_checks(generator)
            result["width"] = width
            (OUTPUT / f"lumping_width{width}.json").write_text(json.dumps(result, indent=2))
            print(json.dumps(result), flush=True)
    elif mode == "transport":
        for width in widths:
            start = time.time()
            generator = WideGenerator(width, noncrossing_states_fast(width))
            injected = prepopulate_linear_cache(generator)
            experiment = pit.width_experiment(
                generator,
                "uniform_join_minus_detach",
                ("krylov",),
                pit.RANKS,
                # r_transport reads only the eta=0 baseline and the PRIMARY_ETAS
                # excess (max over +-1/4); the full repo ladder adds +-1/8, +-1/2,
                # +-1, which do not enter dictionary_scores.
                (0.0,) + pit.PRIMARY_ETAS,
                pit.LAGS,
            )
            rank_rows = experiment["rank_notions"]
            slim = {
                "width": width,
                "states": experiment["states"],
                "injected_r_linear": injected,
                "rank_notions": {
                    name: {
                        "r_linear": row["r_linear"]["dimension"]
                        if isinstance(row["r_linear"], dict) else row["r_linear"],
                        "r_positive_baseline_lumping": row["r_positive_baseline_lumping"]["blocks"],
                        "r_positive_whole_affine_family": row["r_positive_whole_affine_family"]["blocks"],
                        "r_transport": row["r_transport"],
                    }
                    for name, row in rank_rows.items()
                },
                "spans": [
                    {
                        "family": span["family"],
                        "requested_rank": span["requested_rank"],
                        "status": span["status"],
                        "dictionary_scores": span.get("dictionary_scores"),
                    }
                    for span in experiment["spans"]
                ],
                "seconds": round(time.time() - start, 1),
            }
            (OUTPUT / f"transport_width{width}.json").write_text(json.dumps(slim, indent=2))
            print(json.dumps(slim), flush=True)


if __name__ == "__main__":
    main()
