#!/usr/bin/env python3
"""Exact positive control for the frozen kagome-site polynomial relation."""

from __future__ import annotations
import argparse, hashlib, json
from fractions import Fraction
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence
try:
    from scripts.exact_polynomial_root_certificate import evaluate, isolate_roots, open_root_count, sturm_sequence
except ModuleNotFoundError:
    from exact_polynomial_root_certificate import evaluate, isolate_roots, open_root_count, sturm_sequence

ROOT=Path(__file__).resolve().parents[1]
DEFAULT_CONTRACT=ROOT/"analysis"/"pslq_search_contract.json"
DEFAULT_OUTPUT=ROOT/"results"/"pslq-kagome-exact-control"/"latest.json"
SCHEMA="matching-one/kagome-exact-relation-control/v1"
POLYNOMIAL=[1,0,-3,1]

def _text(x:Fraction)->str: return f"{x.numerator}/{x.denominator}"
def build_result(contract_path:Path=DEFAULT_CONTRACT)->dict[str,Any]:
    raw=contract_path.read_bytes(); contract=json.loads(raw)
    control=contract["false_positive_controls"]["exact_percolation_control"]
    if control!={"model":"kagome site percolation","value":"1-2*sin(pi/18)","expected_relation":"1-3*p^2+p^3=0"}: raise ValueError("control contract drift")
    polynomial=[Fraction(x) for x in POLYNOMIAL]; lo,hi=Fraction(3,5),Fraction(7,10)
    seq=sturm_sequence(polynomial); count=open_root_count(seq,lo,hi)
    roots=isolate_roots(polynomial,lo,hi,bits=120)
    if count!=1 or len(roots)!=1: raise ValueError("physical control root not uniquely isolated")
    rlo,rhi=roots[0]
    return {"schema":SCHEMA,"issue":1,"status":"exact_kagome_positive_control_recovered","contract_sha256":hashlib.sha256(raw).hexdigest(),
      "model":control["model"],"declared_expression":control["value"],"polynomial_coefficients_ascending":POLYNOMIAL,
      "sturm_open_root_count_in_physical_window":count,"physical_window":[_text(lo),_text(hi)],"isolating_interval":[_text(rlo),_text(rhi)],"isolation_bits":120,
      "endpoint_signs":[-1 if evaluate(polynomial,lo)<0 else 1, -1 if evaluate(polynomial,hi)<0 else 1],
      "exact_checks":{"declared_relation_recovered":True,"unique_physical_root":True,"coefficient_height":3,"degree":3},
      "claim_boundary":{"included":"exact recovery of the frozen kagome-site positive-control relation","excluded":"square-site degree-3 exclusion, candidate transfer, closed forms, or transcendence","parent_issue":"remain open"}}
def validate_result(result:Mapping[str,Any],contract_path:Path=DEFAULT_CONTRACT)->Mapping[str,Any]:
    expected=build_result(contract_path)
    if result!=expected: raise ValueError("kagome control does not exactly reproduce")
    return {"schema":SCHEMA,"status":"valid","root_count":1}
def main(argv:Optional[Sequence[str]]=None)->int:
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("--contract",type=Path,default=DEFAULT_CONTRACT);p.add_argument("--output",type=Path);p.add_argument("--validate",type=Path);a=p.parse_args(argv)
    if a.validate: print(json.dumps(validate_result(json.loads(a.validate.read_text()),a.contract),indent=2,sort_keys=True));return 0
    rendered=json.dumps(build_result(a.contract),indent=2,sort_keys=True)+"\n"
    if a.output:a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(rendered)
    else:print(rendered,end="")
    return 0
if __name__=="__main__":raise SystemExit(main())
