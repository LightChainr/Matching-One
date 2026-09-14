#!/usr/bin/env python3
"""STEP 0 -- convention audit of the shipped pair.

Reads the two shipped artefacts and decides, with numbers, what the reported
8.6e-9 difference in Omega(axis, l=8) is made of:

  file A (oblique)  : oblique-spin4-controls.json   (rev769 oblique producer)
  file B (closure)  : closure-amplitude-raw.json    (sector802b producer)

Also re-checks the *definition* of every printed field
(`A_estimate`, `root_minus_pc`, `shift_times_ell4`) against the two candidate
definitions, so that a definition mismatch is never reported as a numeric one.
"""
from __future__ import annotations
import sys, os, json, math

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import fl_common as C

IN = "/workspace/dpfloor/in"


def main():
    obl_f = sys.argv[1] if len(sys.argv) > 1 else os.path.join(IN, "oblique-spin4-controls.json")
    clo_f = sys.argv[2] if len(sys.argv) > 2 else os.path.join(IN, "closure-amplitude-raw.json")

    O = json.load(open(obl_f))
    K = json.load(open(clo_f))

    pc_obl = O["reference_pc"]
    pc_clo = K["p_c"]
    rec8 = [r for r in O["records"] if r["direction"] == [1, 0] and r["n"] == 8][0]
    root_obl = rec8["charge_root_p"]
    ell = rec8["physical_circumference"]
    ell4 = ell ** 4

    w8 = K["root"]["8"]
    root_clo = w8["p_ch"]
    Delta_clo = w8["Delta"]
    Delta_p_clo = w8["Delta_p"]

    # ---- rebuild every printed field from the two candidate definitions ----
    def_checks = {}
    # candidate 1:  Omega = -(p - pc) * ell^4                (the task's Omega)
    # candidate 2:  A_est = -(p - pc) * ell^4 / cos(4theta)   (n1105mix claim)
    c4 = rec8["cos_4theta"]
    om_c1 = -(root_obl - pc_obl) * ell4
    a_c2 = -(root_obl - pc_obl) * ell4 / c4
    def_checks["oblique"] = {
        "shift_times_ell4_printed": rec8["shift_times_ell4"],
        "Omega_from_def1": om_c1,
        "Omega_minus_printed": om_c1 - rec8["shift_times_ell4"],
        "A_est_from_def2": a_c2,
        "A_est_printed": rec8["A_estimate"],
        "A_est_minus_printed": a_c2 - rec8["A_estimate"],
        "root_minus_pc_recomputed": root_obl - pc_obl,
        "root_minus_pc_printed": rec8["root_minus_pc"],
    }
    om_clo = -(root_clo - pc_clo) * ell4
    def_checks["closure"] = {
        "w4_times_offset_printed": w8["w4_times_offset"],
        "Omega_from_def1": om_clo,
        "Omega_minus_printed": om_clo - w8["w4_times_offset"],
        "p_c_minus_p_ch_recomputed": pc_clo - root_clo,
        "p_c_minus_p_ch_printed": w8["p_c_minus_p_ch"],
    }

    # ---- the decomposition of the shipped Omega difference ----
    d_omega = om_clo - om_c1
    d_root = root_clo - root_obl
    d_pc = pc_clo - pc_obl
    part_pc = d_pc * ell4
    part_root = -d_root * ell4
    resid = d_omega - (part_pc + part_root)

    conv = {
        "ell": ell, "ell4": ell4,
        "p_c_oblique_file": pc_obl,
        "p_c_closure_file": pc_clo,
        "p_c_difference_closure_minus_oblique": d_pc,
        "root_oblique_file": root_obl,
        "root_closure_file": root_clo,
        "root_difference_closure_minus_oblique": d_root,
        "Omega_oblique_file": om_c1,
        "Omega_closure_file": om_clo,
        "Omega_difference_closure_minus_oblique": d_omega,
        "Omega_difference_relative": d_omega / om_clo,
        "decomposition": {
            "from_different_p_c": part_pc,
            "from_different_root": part_root,
            "residual": resid,
            "share_from_p_c": part_pc / d_omega,
            "share_from_root": part_root / d_omega,
        },
        "delta_p_of_the_Omega_difference": d_omega / ell4,
        "p_c_difference_in_ulp_of_pc": d_pc / (2.0 ** -52 * pc_clo),
    }

    # ---- does anything amplitude-like enter Omega? ----
    amp = {
        "closure_A_open_G4_at_pc_w8": K["amplitude"]["8"]["G4_at_pc"]["A_open"],
        "closure_A_open_G8_at_1mpc_w8": K["amplitude"]["8"]["G8_at_1mpc"]["A_open"],
        "closure_Delta_w8": Delta_clo,
        "closure_Delta_p_w8": Delta_p_clo,
        "closure_p_ch_from_root_solve": root_clo,
        "closure_lin_offset": w8.get("lin_offset"),
        "closure_p_c_minus_p_ch": w8["p_c_minus_p_ch"],
        "statement": ("Delta_w is built from Perron ROOTS only (sector802b's own "
                      "note, C3(a)); A_open / trace coefficients are diagnostics "
                      "and cannot enter Omega.  The p_c the closure file used is "
                      "stored as its own 'p_c' field."),
    }

    print(json.dumps({"definition_checks": def_checks, "convention": conv,
                      "amplitude_irrelevant": amp}, indent=2))
    out = {"definition_checks": def_checks, "convention": conv,
           "amplitude_irrelevant": amp,
           "files": {"oblique": obl_f, "closure": clo_f}}
    op = sys.argv[3] if len(sys.argv) > 3 else "/workspace/dpfloor/out/step0_convention.json"
    with open(op, "w") as fh:
        json.dump(out, fh, indent=2)
    print("->", op)


if __name__ == "__main__":
    main()
