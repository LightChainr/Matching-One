#!/usr/bin/env python3
"""Fourth-order physical increment and exact <=4 clocks, two real prefixes."""
from fractions import Fraction
import json
from math import comb
from pathlib import Path

from p334_pair_only_survival import contiguous, frac, multiply
from p334_pair_triple_survival import polynomial, strip_one_plus_z, binomial

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "results/p334-quartic-clock"


def main():
    quartics = json.loads((OUTPUT / "full_quartics.json").read_text())["checkpoints"]
    triples = {r["replica_counter"]: r for r in json.loads((ROOT / "results/p334-pair-triple-clock/full_triples.json").read_text())["checkpoints"]}
    prior = {r["counter"]: r for r in json.loads((ROOT / "results/p334-pair-triple-clock/pair_triple_survival.json").read_text())["records"]}
    records = []
    for q in quartics:
        counter = q["replica_counter"]
        old = prior[counter]
        graph = json.loads((ROOT / old["source_graph"]).read_text())
        counts, structure = polynomial(graph, triples[counter]["all_minimal_nonfaces"], q["all_minimal_quartics"])
        if counts[:4] != old["independent_counts"][:4] or counts[4] != q["actual_safe_four_sets"]:
            raise ValueError("polynomial does not reproduce the new physical fourth-order coefficient")
        d = q["d"]
        survival = [Fraction(counts[k], comb(d, k)) for k in range(d + 1)]
        hazard = [None] + [1 - survival[k] / survival[k - 1] if survival[k - 1] else None for k in range(1, d + 1)]
        mean = sum(survival[:-1])
        row = {"counter": counter, "seed": q["seed"], "N": q["N"], "k0": q["k0"], "d": d,
               "source_graph": old["source_graph"], "minimal_quartics": q["minimal_nonfaces_size4"],
               "physical_safe_four_sets": q["actual_safe_four_sets"], **structure,
               "independent_counts": counts, "survival": [frac(s) for s in survival],
               "hazard": [frac(h) if h is not None else None for h in hazard],
               "maximum_safe_k": max(k for k, count in enumerate(counts) if count),
               "mean_first_trigger_step": frac(mean),
               "mean_shortening_from_quartics": frac(Fraction(old["mean_first_pair_triple_trigger_step"]["exact"]) - mean),
               "first_trigger_quantiles": {str(q): next(k for k in range(d + 1) if survival[k] <= 1 - q)
                                           for q in (Fraction(1, 10), Fraction(1, 2), Fraction(9, 10))},
               "removed_triple_safe_counts": [a - b for a, b in zip(old["independent_counts"], counts)]}
        records.append(row)
        print(counter, "quartics", row["minimal_quartics"], "states", row["memoized_states"],
              "groups", len(row["false_twin_groups"]), "mean", row["mean_first_trigger_step"],
              "shortening", row["mean_shortening_from_quartics"], "max_k", row["maximum_safe_k"], flush=True)
    A, B = records
    difference = [b - a for a, b in zip(A["independent_counts"], B["independent_counts"])]
    loss_difference = [b - a for a, b in zip(A["removed_triple_safe_counts"], B["removed_triple_safe_counts"])]
    crosses = {k: A["independent_counts"][k] * B["independent_counts"][k - 1]
                  - B["independent_counts"][k] * A["independent_counts"][k - 1]
               for k in range(1, d + 1) if A["independent_counts"][k - 1] and B["independent_counts"][k - 1]}
    factors = [strip_one_plus_z(row["independent_counts"]) for row in records]
    common = min(power for power, _ in factors)
    residuals = [multiply(p, binomial(power - common)) for power, p in factors]
    width = max(map(len, residuals))
    for p in residuals:
        p += [0] * (width - len(p))
    residual_difference = [b - a for a, b in zip(*residuals)]
    z_power = next(k for k, value in enumerate(residual_difference) if value)
    new_gap = Fraction(B["mean_first_trigger_step"]["exact"]) - Fraction(A["mean_first_trigger_step"]["exact"])
    old_gap = Fraction(prior[B["counter"]]["mean_first_pair_triple_trigger_step"]["exact"]) - Fraction(prior[A["counter"]]["mean_first_pair_triple_trigger_step"]["exact"])
    comparison = {"count_B_minus_A": difference,
                  "survival_B_above": contiguous([k for k, v in enumerate(difference) if v > 0]),
                  "survival_B_below": contiguous([k for k, v in enumerate(difference) if v < 0]),
                  "survival_equal": contiguous([k for k, v in enumerate(difference) if v == 0]),
                  "hazard_crossproducts_B_minus_A": crosses,
                  "hazard_B_above": contiguous([k for k, v in crosses.items() if v > 0]),
                  "hazard_B_below": contiguous([k for k, v in crosses.items() if v < 0]),
                  "hazard_equal": contiguous([k for k, v in crosses.items() if v == 0]),
                  "mean_first_trigger_B_minus_A": frac(new_gap),
                  "fraction_of_previous_mean_gap_retained": frac(new_gap / old_gap),
                  "quartic_shortening_B_minus_A": frac(old_gap - new_gap),
                  "quartic_loss_count_B_minus_A": loss_difference,
                  "quartic_survival_loss_B_above": contiguous([k for k, v in enumerate(loss_difference) if v > 0]),
                  "quartic_survival_loss_B_below": contiguous([k for k, v in enumerate(loss_difference) if v < 0]),
                  "dominance_factor": {"common_1_plus_z_power": common, "z_power": z_power,
                                       "remaining_coefficients": residual_difference[z_power:],
                                       "all_remaining_coefficients_strictly_positive": all(v > 0 for v in residual_difference[z_power:])},
                  "physical_four_step_survival_B_minus_A": frac(Fraction(difference[4], comb(d, 4))),
                  "physical_four_step_hazard_B_minus_A": frac(Fraction(B["hazard"][4]["exact"]) - Fraction(A["hazard"][4]["exact"]))}
    output = {"parent_commit": "d5d2cc89e77ebb2ec6252df75dc858e9c240e6ce", "new_samples": 0,
              "records": records, "comparison": comparison,
              "boundary": "Complete physical rank-one survival through k=4 on two fixed real checkpoints. Full <=4 hypergraph clocks thereafter are upper bounds on true survival, not full physical clocks or hazard bounds."}
    (OUTPUT / "quartic_survival.json").write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print({k: v for k, v in comparison.items() if k not in ("count_B_minus_A", "hazard_crossproducts_B_minus_A", "quartic_loss_count_B_minus_A", "dominance_factor")})


if __name__ == "__main__":
    main()
