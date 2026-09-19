"""Issue 814: exact bookkeeping only; no lattice enumeration or fitting."""
from fractions import Fraction as F
from pathlib import Path
import json

def main():
    cases = {
        "identity_spin4_root": (F(4)-F(5,4), F(11,4)),
        "thermal_level4_root": (F(21,4)-F(5,4), F(4)),
        "row_gap_exponent": (F(21,4)-1, F(17,4)),
        "row_slope_exponent": (F(5,4)-1, F(1,4)),
        "row_ratio_exponent": (F(17,4)-F(1,4), F(4)),
        "ward_real_pair_coefficient_over_pi4": (2*3*F(5,8)/45, F(1,12)),
    }
    for name, (actual, expected) in cases.items():
        assert actual == expected, (name, actual, expected)
    result = {"issue":814,"decision":"c","E4_downstream":"frozen",
              "T6_executed":False,"checks":{k:str(v[0]) for k,v in cases.items()},
              "passed":len(cases),"quoted_not_recomputed":{"Jacobsen_Delta1":"4.0001(2)"},
              "scope":"Exact rational bookkeeping only; not a lattice or Ward-interface proof."}
    Path(__file__).with_name("result.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
    print(f"{len(cases)} exact rational checks passed")

if __name__ == "__main__":
    main()
