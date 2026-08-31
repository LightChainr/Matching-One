#!/usr/bin/env python3
"""Reduce the regular pair counterterm Gram exactly, without a parameter scan."""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import sympy as sp


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", default="analysis/regular_pair_counterterm_contract.json")
    parser.add_argument("--output-dir", default="results/regular-pair-counterterm")
    args = parser.parse_args()
    started = time.perf_counter()
    contract_path = Path(args.contract)
    contract = json.loads(contract_path.read_text())
    q = sp.Symbol("Q")
    delta = q-1
    # These three closed contractions are inputs from the pinned source note.
    g22 = q*(q-3)*(3*q*q-9*q+8)/(8*(q-1)*(q-2))
    g20 = (q-3)/(4*(q-1))
    g00 = (2*q*q-4*q+3)/(2*q*(q-1))
    regular = sp.factor(g22+2*g20+g00)
    gram_q = sp.Matrix([
        [regular, sp.factor(delta*(g20+g00))],
        [sp.factor(delta*(g20+g00)), sp.factor(delta**2*g00)],
    ])
    endpoint = gram_q.applyfunc(lambda x:sp.limit(x,q,1))
    gram = gram_q.applyfunc(lambda x:sp.limit(sp.diff(x,q),q,1))
    g,a = sp.symbols("g a", real=True)
    quadratic = sp.expand((sp.Matrix([g,a]).T*gram*sp.Matrix([g,a]))[0])
    schur = sp.factor(gram[0,0]-gram[0,1]**2/gram[1,1])
    completed = schur*g*g+gram[1,1]*(a+gram[0,1]*g/gram[1,1])**2
    alpha_x,alpha_y = sp.symbols("alpha_x alpha_y", real=True)
    cross = sp.expand((sp.Matrix([1,alpha_x]).T*gram*sp.Matrix([1,alpha_y]))[0])
    vx,vy = sp.symbols("v_x v_y", positive=True)
    beta = sp.Symbol("beta", real=True)
    # beta is an arbitrary next Taylor coefficient, not a selected extra model.
    c = 1+a*delta+beta*delta**2
    general_closed = sp.cancel(g22+2*c*g20+c*c*g00)
    arbitrary_second_jet = sp.factor(sp.limit(general_closed/delta,q,1))
    result = {
        "schema":"matching-one.regular-pair-counterterm.v1",
        "source":contract["source"],
        "basis":contract["basis"],
        "rational_Gram_Q":[[str(x) for x in row] for row in gram_q.tolist()],
        "endpoint_Gram":[[str(x) for x in row] for row in endpoint.tolist()],
        "first_Q_Gram":[[str(x) for x in row] for row in gram.tolist()],
        "determinant":str(gram.det()),
        "quadratic_form":str(quadratic),
        "completed_square":str(completed),
        "positive_definite_reason":"leading principal minors 13/8 and 3/4 are positive",
        "unit_pair_lower_bound":str(schur),
        "same_counterterm_first_Q_two_copy":str(sp.factor(quadratic.subs(g,1))),
        "arbitrary_second_Taylor_coefficient_result":str(arbitrary_second_jet),
        "different_counterterm_cross":str(cross),
        "conditional_connected_first_Q_coefficient":str(quadratic/((1+vx)*(1+vy))),
        "conditional_unit_pair_lower_bound":str(schur/((1+vx)*(1+vy))),
        "single_mark_four_free_Q_activation":"g; aL has zero one-Q derivative in the four-singleton exterior",
        "single_U_response_family":"W(g,a)=g*W_canonical-a*V_old; previously scored values are inputs, no new U score",
        "decisions":[
            "no_nonzero_first_jet_direction_in_this_regular_two_dimensional_family_has_zero_identical_two_copy_Q_activation",
            "all_real_uniform_counterterms_with_unit_pair_coefficient_have_shared_four_line_Q_activation_at_least_3_over_2",
            "a_single_N25_U_mixed_coefficient_is_not_completion_invariant",
        ],
        "scope":[
            "same real counterterm prescription at both insertions for the positive lower bound",
            "specified physical exterior and Q derivative, not ordinary Q1 probability or continuum norm",
            "no bound for arbitrary distinct alpha_x,alpha_y cross contractions",
            "no assertion about occupation-summed homogeneous global U or asymptotic scaling",
            "completed-square minimizer is an algebraic bound, not a recommendation to select that counterterm",
        ],
    }
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    output = out/"latest.json"
    output.write_text(json.dumps(result,indent=2)+"\n")
    report = out/"REPORT.md"
    report.write_text("""# The regular pair interaction survives every uniform finite singlet counterterm

Use the already specified basis `Kreg=K2+K0`, `L=(Q−1)K0`.
The complete first-Q four-line Gram matrix is

```text
G1 = [[13/8, -1/4], [-1/4, 1/2]],     det(G1)=3/4.
```

It is strictly positive definite. For a real regular first-jet direction
`K=g Kreg+a L`, the identical two-copy activation is exactly

`(g,a) G1 (g,a)^T = (3/2)g²+(1/2)(a−g/2)²`.

For `K2+c(Q)K0`, with `c(1)=1` and arbitrary real `alpha=c'(1)`, this gives

`H_alpha=13/8−alpha/2+alpha²/2 ≥ 3/2`.

Higher Taylor coefficients of c do not enter. This is a bound over the
whole stated family, not a coefficient fit or scan. It excludes an additive
first-Q interaction in the fixed shared-four-line exterior for every such
uniform counterterm, even though one N25 mixed U response can be shifted
by `W_alpha=W_canonical−alpha V_old`.

On the existing contractible L17 two-hole exterior, summing both holes'
vacant/occupied states divides the connected mixed Q coefficient by
`(1+v_x)(1+v_y)`. Its unit-pair lower bound is therefore
`3/[2(1+v_x)(1+v_y)]`, strictly positive for positive activities.

Different counterterms at the two insertions instead give
`13/8−(alpha_x+alpha_y)/4+alpha_x alpha_y/2`; no positive lower bound is
claimed for that cross response. Nor is this a universal field norm or a
result for the occupation-summed homogeneous global observable.

The old canonical mixed U and pure-pair U scores were not rerun. The
calculation reduces three pinned rational contractions symbolically,
without an occupation enumeration, Monte Carlo, new roots or tests.
See latest.json and run.json for exact formulas, definitions and provenance.
""")
    source = contract["source"]
    source_bytes = subprocess.check_output(["git","show",source["commit"]+":"+source["path"]])
    sha = lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
    receipt = {
        "schema":"matching-one.regular-pair-counterterm.run.v1",
        "definition_commit":subprocess.check_output(["git","rev-parse","HEAD"],text=True).strip(),
        "command":[sys.executable,*sys.argv],
        "created_utc":datetime.now(timezone.utc).isoformat(),
        "python":sys.version,"sympy":sp.__version__,"machine":platform.machine(),
        "elapsed_seconds":time.perf_counter()-started,
        "source_sha256":hashlib.sha256(source_bytes).hexdigest(),
        "script_sha256":sha(__file__),"contract_sha256":sha(contract_path),
        "output_sha256":{p.name:sha(p) for p in [output,report]},
        "new_occupation_enumerations":0,"new_random_samples":0,"new_U_scores":0,
        "parameter_scans":0,"root_searches":0,"tests_run":0,"cloud_jobs":0,
    }
    (out/"run.json").write_text(json.dumps(receipt,indent=2)+"\n")
    print(json.dumps({"first_Q_Gram":result["first_Q_Gram"],"determinant":result["determinant"],"completed_square":result["completed_square"],"unit_pair_lower_bound":result["unit_pair_lower_bound"],"higher_jet_result":result["arbitrary_second_Taylor_coefficient_result"],"elapsed_seconds":receipt["elapsed_seconds"]},indent=2))


if __name__ == "__main__":
    main()
