from __future__ import annotations

import json
import math
from pathlib import Path

lam = math.log(2.0)
p0 = math.exp(-lam)
p1 = lam * p0
pge2 = 1.0 - p0 - p1

out = {
    "lambda": lam,
    "P_K_0": p0,
    "P_K_1": p1,
    "P_K_ge_2": pge2,
    "P_K_ge_2_given_K_ge_1": pge2 / (1.0 - p0),
    "P_K_1_given_K_ge_1": p1 / (1.0 - p0),
    "E_K_given_K_ge_1": lam / (1.0 - p0),
    "same_parameter_lower_window_joint_law": {
        "(W4,W8)=(0,1)": p0,
        "(W4,W8)=(1,1)": p1,
        "(W4,W8) has k>=2 equal counts": pge2,
    },
    "E_W8_lower_window_if_W4_Poisson": lam + p0,
}

if __name__ == "__main__":
    text = json.dumps(out, indent=2)
    print(text)
    Path("results/geometric-consistency/poisson-topology-constants-20260914.json").write_text(text + "\n")
