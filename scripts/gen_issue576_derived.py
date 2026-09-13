#!/usr/bin/env python3
"""Generate derived JSON for issue #576 Part 1b (exact enumeration) and Part 2
(Pinson/Arguin continuum wrapping at tau=i,2i,4i)."""
import json
import mpmath as mp
import sys

sys.path.insert(0, "scripts")
import pinson_arguin_primitive as P

mp.mp.dps = 90

res = {
    "formula": "Pinson 1994 / Arguin 2002 (hep-th/0111193) continuum critical Q=1 FK wrapping probability pi_tau({1,0})",
    "coupling_note": "percolation specialization e0=2/3 (Coulomb-gas g such that A=2*pi|u|^2/(3*Im tau)); tau=i*r below",
    "values": [],
}
for r in (1, 2, 4):
    tau = mp.mpc(0, r)
    d = P.primitive_probability_direct(1, 0, tau, dps=80)
    t = P.primitive_probability_theta(1, 0, tau, dps=80)
    res["values"].append({
        "r": r,
        "tau": "%di" % r,
        "pi_direct": mp.nstr(d, 55),
        "pi_theta": mp.nstr(t, 55),
        "abs_diff_direct_theta": mp.nstr(abs(d - t), 3),
        "dps": 80,
    })
with open("results/issue576/wrapping-grounding/derived/pinson_arguin_r1r2r4.json", "w") as fh:
    json.dump(res, fh, indent=2)
print("pinson values:")
for v in res["values"]:
    print("  r=%d  pi({1,0})(%si) = %s  (direct~theta diff %s)" % (v["r"], v["r"], v["pi_direct"], v["abs_diff_direct_theta"]))

enum = {
    "definition": "NN square SITE percolation on L x L torus; horizontal (period-1) winding. Exact brute force 2^(L^2) via scripts/exact_wrapping_enum.cpp (method: double the grid horizontally, detect component containing (r,0) and (r,L)).",
    "values": [
        {"L": 2, "count": 7, "fraction": "7/16", "p_half": 0.4375},
        {"L": 3, "count": 175, "fraction": "175/512", "p_half": 0.341796875},
        {"L": 4, "count": 19571, "fraction": "19571/65536", "p_half": 0.2986297607421875},
        {"L": 5, "count": 8853291, "fraction": "8853291/33554432", "p_half": 0.263848632574},
    ],
    "reachable": "L<=5 exact (2^25 = 33,554,432 configs, ~20s in C++). L=6 (2^36) not attempted.",
}
with open("results/issue576/wrapping-grounding/derived/exact_enum_l2l5.json", "w") as fh:
    json.dump(enum, fh, indent=2)
print("wrote exact_enum json")
