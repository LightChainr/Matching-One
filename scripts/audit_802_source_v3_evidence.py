#!/usr/bin/env python3
"""verify802src -- assemble evidence.json on the cloud (no local arithmetic).

Merges v1_exact.json + v2_c4.json and adds the derived quantities that the
erratum needs (corrected phi, corrected Nhat/S_g, raw-vs-corrected comparison).

Outputs /workspace/v802src/out/evidence.json
"""
from __future__ import annotations
import json

IN = "/workspace/v802src/out"
PC = 0.5927460507921


def main():
    v1 = json.load(open(IN + "/v1_exact.json"))
    v2 = json.load(open(IN + "/v2_c4.json"))
    resp = json.load(open("/workspace/v802src/in/root-response.json"))
    at = resp["sources"]

    derived = []
    for w in [4, 5, 6, 7, 8]:
        w = str(w)
        fB = int(w) * PC * PC
        fW = int(w) * (1 - PC) ** 2
        row = {}
        for src, fg in (("black", fB), ("white", fW)):
            e = at[src]["at_pc"][w]
            Dg = e["Delta_g"]
            Sg_raw = e["S_g"]
            Sg_phys = Sg_raw + 2 * fg
            Sp = resp["delta_table"][w]["0.5927460507921"]["S_p"]
            Dp = resp["delta_table"][w]["0.5927460507921"]["Delta_p"]
            row[src] = {
                "f_g": fg,
                "S_g_unrenormalised": Sg_raw,
                "S_g_physical": Sg_phys,
                "Nhat_raw": e["Nhat_per_row"],
                "Nhat_physical": e["Nhat_per_row"] + fg,
                "phi_raw": (Sp * Dg) / (Dp * Sg_raw),
                "phi_physical": (Sp * Dg) / (Dp * Sg_phys),
                "one_over_abs_phi_raw": abs(1.0 / ((Sp * Dg) / (Dp * Sg_raw))),
                "one_over_abs_phi_physical": abs(1.0 / ((Sp * Dg) / (Dp * Sg_phys))),
            }
        row["w"] = int(w)
        row["T_black"] = at["black"]["at_pc"][w]["T"]
        row["T_white"] = at["white"]["at_pc"][w]["T"]
        row["T_white_minus_T_black"] = row["T_white"] - row["T_black"]
        row["two_p_one_minus_p"] = 2 * PC * (1 - PC)
        derived.append(row)
        print("w=%s phi_black_raw=%.6e phi_black_phys=%.6e  "
              "1/|phi| raw=%.1f -> phys=%.1f"
              % (w, row["black"]["phi_raw"], row["black"]["phi_physical"],
                 row["black"]["one_over_abs_phi_raw"],
                 row["black"]["one_over_abs_phi_physical"]), flush=True)

    out = {
        "meta": {
            "p_c": PC,
            "team_json": "sector-root-response-20260914.json (commit 723f629)",
            "analysis_note": "docs/research-bridges-post-compass-20260914.md / "
                             "/tmp/concurrent-out/two-observable-response-and-angular-alias-20260914.md "
                             "(§7)",
            "engine_weight_reading": "sector802_lib.build_matrices: R_ij = sum p^k (1-p)^(w-k) "
                                     "exp(gamma h)  -> UNNORMALISED row weight (no /Z_row)",
            "local_computation": "none; every number produced on DevEnvC_NePnUn",
        },
        "C1a_identity_full_enumeration": v1["C1a_identity"],
        "C1b_measure_equivalence_exact": v1["C1b_measure_equivalence"],
        "C1c_T_identity_on_team_data": v1["C1c_T_identity"],
        "C1d_finite_torus_independent": v1["C1d_finite_torus"],
        "C2f_row_partition_exact": v1["C2f_row_partition_exact"],
        "C2_recompute": v1["C2_recompute"],
        "C3_counterexample": v1["C3_counterexample"],
        "derived_phi_and_Sg": derived,
        "C4_torus_cases": v2["cases"],
    }
    with open(IN + "/evidence.json", "w") as f:
        json.dump(out, f, indent=1, default=str)
    print("-> %s/evidence.json" % IN)


if __name__ == "__main__":
    main()
