#!/usr/bin/env python3
"""Issue #593 Phase B controls.  Non-negotiable checks from the ticket:

1. transition-table content hash widths 1-8 (repository's own builder);
2. exact row sums zero, rate positivity for |eta| <= 1, strong connectivity;
3. exp(tG) 1 = 1 to machine precision (uniformization evolution);
4. uniformization against an independent dense Taylor exponential at width 5;
5. full-rank projection gives K identically zero and a generator-closed Krylov
   span gives K = 0 to 1e-8 relative (repository's ``exact_controls``);
6. resolvent Schur identity at width 5 (repository's ``exact_controls``).
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from p593_noncrossing_fast import WideGenerator, noncrossing_states_fast, validate_shim
from p398_intervention_transport import (
    Generator,
    evolve_observables,
)
from p398_projected_memory import KERNEL_GRID, exact_controls

OUTPUT = Path(__file__).resolve().parents[1] / "results" / "p593-width9-10-memory-degree" / "raw"


def chain_controls(generator, eta_values=(-0.5, -0.25, 0.0, 0.25, 0.5, 1.0)):
    size = generator.size
    worst_row_sum = 0.0
    min_rate = float("inf")
    for eta in eta_values:
        for intervention in ("uniform_join_minus_detach", "uniform_detach_only", "single_point_join"):
            rates = generator.rates(intervention, eta)
            min_rate = min(min_rate, min(rates.values()))
            rows = generator.rows(rates)
            for row in rows:
                worst_row_sum = max(worst_row_sum, abs(sum(v for _, v in row)))
    # strong connectivity of the transition graph (any rate > 0 edges)
    rates = generator.baseline_rates()
    adjacency = [[] for _ in range(size)]
    for source in range(size):
        for move, rate in rates.items():
            image = generator.target[move][source]
            if image is not None and rate > 0:
                adjacency[source].append(image)
    seen = {0}
    stack = [0]
    while stack:
        node = stack.pop()
        for nxt in adjacency[node]:
            if nxt not in seen:
                seen.add(nxt)
                stack.append(nxt)
    return {
        "worst_abs_row_sum": worst_row_sum,
        "min_rate_over_declared_etas": min_rate,
        "strongly_connected": len(seen) == size,
    }


def exp_constant_control(generator):
    rates = generator.baseline_rates()
    rows = generator.rows(rates)
    rate = generator.exit_rate(rates)
    constant = [[1.0] * generator.size]
    evolved = evolve_observables(rows, generator.size, rate, constant, list(KERNEL_GRID))
    worst = max(abs(value - 1.0) for block in evolved for value in block[0])
    return worst


def uniformization_versus_taylor(width=5):
    """exp(tG) f by uniformization against a dense Taylor/scipy exponential."""
    import scipy.linalg

    generator = Generator(width)
    size = generator.size
    rates = generator.baseline_rates()
    rows = generator.rows(rates)
    rate = generator.exit_rate(rates)
    dense = np.zeros((size, size))
    for source, row in enumerate(rows):
        for destination, value in row:
            dense[source, destination] = value
    observables = [[1.0 if abs(s[0] - s[-1]) == 0 else 0.0 for s in generator.states]]
    worst = 0.0
    for t in (0.25, 1.0, 2.0):
        exact = scipy.linalg.expm(t * dense) @ np.array(observables[0])
        evolved = evolve_observables(rows, size, rate, observables, [t])
        worst = max(worst, float(np.max(np.abs(exact - np.array(evolved[0][0])))))
    return worst


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    result = {}
    result["control_1_transition_table"] = validate_shim(8)
    widths = {}
    for width in (4, 6, 8, 9, 10):
        generator = WideGenerator(width, noncrossing_states_fast(width))
        widths[str(width)] = chain_controls(generator)
    result["control_2_chain_structure"] = widths
    result["control_3_exp_tG_one"] = {
        str(w): exp_constant_control(WideGenerator(w, noncrossing_states_fast(w)))
        for w in (4, 8, 10)
    }
    result["control_4_uniformization_vs_taylor_width5"] = uniformization_versus_taylor(5)
    result["controls_5_6_exact_width4"] = exact_controls(4)
    result["controls_5_6_exact_width5"] = exact_controls(5)
    path = OUTPUT / "controls.json"
    path.write_text(json.dumps(result, indent=2, default=str))
    print(json.dumps(result, indent=2, default=str)[:3000])


if __name__ == "__main__":
    main()
