#!/usr/bin/env python3
"""Three specified geometries and their mechanism responses; not a census.

Reuse only the lifted topology function. All outputs are exact fractions.
No parameter search, permutation enumeration, sampling, or old-data rescore.
"""
from fractions import Fraction as F
from functools import lru_cache
import importlib.util
from itertools import combinations
import json
from math import comb
from pathlib import Path

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent / "birth-completion-geometry-20260929/analyze.py"
SPEC = importlib.util.spec_from_file_location("topology_only", SOURCE)
TOPO = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(TOPO)
TOPO.L, TOPO.N = 4, 16
NEIGHBORS = [[((v % 4 + dx) % 4 + 4*((v//4 + dy) % 4), dx, dy)
              for dx, dy in [(1, 0), (-1, 0), (0, 1), (0, -1)]]
             for v in range(16)]


@lru_cache(None)
def topology(mask):
    return TOPO.topology(mask, NEIGHBORS)


@lru_cache(None)
def survival(mask, h):
    if topology(mask)[0] != 1:
        return F(0)
    if h == 0:
        return F(1)
    vacant = [v for v in range(16) if not mask >> v & 1]
    return sum((survival(mask | (1 << v), h-1) for v in vacant), F(0))/len(vacant)


def geometry(name, vertices):
    mask = sum(1 << v for v in vertices)
    rank, direction = topology(mask)
    vacancies = [v for v in range(16) if not mask >> v & 1]
    completing = [v for v in vacancies if topology(mask | (1 << v))[0] == 2]
    safe = [v for v in vacancies if v not in completing]
    edges = [(v, w) for v, w in combinations(safe, 2)
             if topology(mask | (1 << v) | (1 << w))[0] == 2]
    degrees = [sum(v in pair for pair in edges) for v in safe]
    m, s, c, e = len(vacancies), len(safe), len(completing), len(edges)
    mean_d = F(sum(degrees), s)
    variance_d = F(sum(d*d for d in degrees), s)-mean_d**2
    retention_counts = [sum(topology(sum(1 << v for v in subset))[0] >= 1
                            for subset in combinations(vertices, a))
                        for a in range(len(vertices)+1)]
    # exp(theta)=2: a finite, fully normalized source kick, not a derivative
    # inferred by numerical differencing. Only the first insertion is tilted.
    weights = [2**d for d in degrees]
    tilted_mean_d = F(sum(w*d for w, d in zip(weights, degrees)), sum(weights))
    lag_rows = []
    for h in range(1, m+1):
        tails = [survival(mask | (1 << v), h-1) for v in safe]
        lag_rows.append({
            "lag_insertions": h,
            "uniform_survival": str(survival(mask, h)),
            "first_kick_derivative_at_zero": str(sum(((d-mean_d)*f for d, f in zip(degrees, tails)), F(0))/m),
            "first_kick_exp_theta_2_survival": str(F(s, m)*sum((w*f for w, f in zip(weights, tails)), F(0))/sum(weights)),
        })
    return {
        "name": name, "occupied": vertices, "rank": rank,
        "direction": direction, "c": c, "e": e, "synergy_edges": edges,
        "safe_degree_counts": {str(d): degrees.count(d) for d in sorted(set(degrees))},
        "safe_degree_mean": str(mean_d), "safe_degree_variance": str(variance_d),
        "retaining_subsets_by_size": retention_counts,
        "birth_cdf_at_count_a": [str(F(n, comb(len(vertices), a)))
                                  for a, n in enumerate(retention_counts)],
        "survival_two_uniform": str(F(s*(s-1)-2*e, m*(m-1))),
        "survival_two_first_kick_derivative_at_zero": str(-F(s, m*(m-1))*variance_d),
        "survival_two_first_kick_exp_theta_2": str(F(s, m)*(1-F(c, m-1)-tilted_mean_d/(m-1))),
        "all_lags_single_kick": lag_rows,
        "mean_remaining_steps_to_completion_derivative": str(sum((F(r["first_kick_derivative_at_zero"]) for r in lag_rows), F(0))),
    }


def main():
    rows = [geometry("A", [0, 1, 4, 5, 8, 12]),
            geometry("B", [0, 2, 4, 5, 8, 12]),
            geometry("C", [0, 1, 5, 8, 9, 12])]
    # A and B are the already published L4 geometric witness. C is one
    # specified six-site bent loop, not a best-effect search over subsets.
    assert all((r["rank"], r["direction"], r["c"]) == (1, (0, 1), 0) for r in rows)
    assert rows[0]["birth_cdf_at_count_a"] == rows[1]["birth_cdf_at_count_a"]
    assert rows[0]["e"] == rows[2]["e"] == 2 and rows[1]["e"] == 3
    result = {
        "schema": "matching-one.birth-selection-intervention.v1",
        "date": "2026-09-30", "lattice": "square NN occupied-site L4 torus",
        "clock": "insertion count; k=6; one tilted insertion then original uniform continuation at every later count",
        "source": "P_theta(v|A,safe) proportional to exp(theta*synergy_degree(v)); exit mass c/m unchanged",
        "scope": "exact specified-preparation controls, not uniform-start population inference or L512 evidence",
        "topology_source": str(SOURCE.relative_to(HERE.parent.parent)),
        "topology_source_commit": "8b37c20e5f2b1946cefbed98b09280a0cad25ed2",
        "rows": rows,
        "equal_A_B_entry_ensemble_e_tilt_survival_derivative": str(-F(1, 180)),
    }
    (HERE / "result.json").write_text(json.dumps(result, indent=2)+"\n")
    lines = ["# Three physical mechanism controls", "",
             "Exact fractions on three specified square-L4 preparations; no census or production.", "",
             "| Set | c | e | P(J1<=4 given A6) | Var_safe(degree) | S2(0) | S2'(0) | S2(log 2) |",
             "|---|---:|---:|---:|---:|---:|---:|---:|"]
    for r in rows:
        lines.append(f"| {r['name']} | {r['c']} | {r['e']} | {r['birth_cdf_at_count_a'][4]} | "
                     f"{r['safe_degree_variance']} | {r['survival_two_uniform']} | "
                     f"{r['survival_two_first_kick_derivative_at_zero']} | {r['survival_two_first_kick_exp_theta_2']} |")
    lines += ["", "A/B have the same entire conditional birth law but different futures and source response.",
              "A/C have the same c,e and two-step response but different conditional birth laws.",
              "For iid thinning, their homology-retention polynomials are u^4,u^4,u^6.",
              "The source kick is an explicitly modified source, not the original Bernoulli experiment.", ""]
    lines += ["## Actual longer-window response, not extrapolated curvature", "",
              "| Insertions after kick | A derivative | B derivative | C derivative |",
              "|---:|---:|---:|---:|"]
    for h in range(10):
        lines.append("| " + str(h+1) + " | " + " | ".join(r["all_lags_single_kick"][h]["first_kick_derivative_at_zero"] for r in rows) + " |")
    lines += ["", "Mean remaining steps-to-completion derivatives: " + "; ".join(
        f"{r['name']}={r['mean_remaining_steps_to_completion_derivative']}" for r in rows) + ".", "",
        "Only subsets reachable from the three specified starts are visited; no uniform-start full census.", ""]
    (HERE / "RESULT.md").write_text("\n".join(lines))
    print("\n".join(lines))


if __name__ == "__main__":
    main()
