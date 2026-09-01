#!/usr/bin/env python3
"""Exact hypothesis-count ledger for every frozen first-pass search family."""

from __future__ import annotations
import argparse, hashlib, json
from functools import lru_cache
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence

ROOT=Path(__file__).resolve().parents[1]
DEFAULT_CONTRACT=ROOT/"analysis"/"pslq_search_contract.json"
DEFAULT_OUTPUT=ROOT/"results"/"pslq-look-elsewhere-ledger"/"latest.json"
SCHEMA="matching-one/pslq-look-elsewhere-ledger/v1"

def mobius(n:int)->int:
    value=n;prime_count=0;p=2
    while p*p<=value:
        if value%p==0:
            value//=p;prime_count+=1
            if value%p==0:return 0
            while value%p==0:value//=p
        p+=1
    if value>1:prime_count+=1
    return -1 if prime_count%2 else 1
def primitive_polynomial_count(degree:int,height:int)->int:
    if degree<1 or height<1:raise ValueError("degree and height must be positive")
    return sum(mobius(g)*(height//g)*(2*(height//g)+1)**degree for g in range(1,height+1))
def primitive_pairwise_constant_count(height:int)->int:
    # a+b*p+c*C with b>0, c!=0, |a|,|c|<=H, gcd(a,b,c)=1.
    total=primitive_polynomial_count(2,height)
    no_constant=primitive_polynomial_count(1,height)
    return total-no_constant
@lru_cache(maxsize=2)
def build_result(contract_path:Path=DEFAULT_CONTRACT)->dict[str,Any]:
    raw=contract_path.read_bytes();contract=json.loads(raw);poly=contract["search_stages"]["algebraic_polynomial"];height=poly["coefficient_height_max"]
    degree_counts={str(d):primitive_polynomial_count(d,height) for d in range(poly["degree_min"],poly["degree_max"]+1)}
    constants=len(contract["search_stages"]["standard_constant_pairwise"]["library"]);per_constant=primitive_pairwise_constant_count(height)
    lattice=len(contract["search_stages"]["lattice_native_candidates"]["library"]);synthetic=contract["false_positive_controls"]["synthetic_random_constants"]["count"]
    method_intervals=len(contract["intervals"])
    families={"algebraic_polynomials_per_interval":sum(degree_counts.values()),"standard_constant_relations_per_interval":constants*per_constant,"lattice_native_candidates_per_interval":lattice,"synthetic_controls":synthetic}
    return {"schema":SCHEMA,"issue":1,"status":"exact_frozen_hypothesis_count_ledger","contract_sha256":hashlib.sha256(raw).hexdigest(),"coefficient_height_max":height,
      "primitive_polynomial_counts_by_degree":degree_counts,"standard_constant":{"constant_count":constants,"primitive_relations_per_constant":per_constant,"relations_per_method_interval":constants*per_constant},
      "lattice_native_candidate_count":lattice,"synthetic_control_count":synthetic,"method_interval_count":method_intervals,"family_counts":families,
      "total_declared_interval_comparisons":method_intervals*(families["algebraic_polynomials_per_interval"]+families["standard_constant_relations_per_interval"]+lattice),
      "exact_checks":{"degree1_matches_committed":degree_counts["1"]==12175,"degree2_matches_committed":degree_counts["2"]==3355121,"contract_requires_look_elsewhere":contract["result_policy"]["look_elsewhere_count_required"] is True},
      "claim_boundary":{"included":"exact cardinalities of the frozen finite search families","excluded":"a multiplicity correction, null distribution, p-value, near-hit promotion, closed form, or transcendence","parent_issue":"remain open"}}
def validate_result(result:Mapping[str,Any],contract_path:Path=DEFAULT_CONTRACT)->Mapping[str,Any]:
    expected=build_result(contract_path)
    if result!=expected:raise ValueError("look-elsewhere ledger does not exactly reproduce")
    return {"schema":SCHEMA,"status":"valid","interval_comparisons":expected["total_declared_interval_comparisons"]}
def main(argv:Optional[Sequence[str]]=None)->int:
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("--contract",type=Path,default=DEFAULT_CONTRACT);p.add_argument("--output",type=Path);p.add_argument("--validate",type=Path);a=p.parse_args(argv)
    if a.validate:print(json.dumps(validate_result(json.loads(a.validate.read_text()),a.contract),indent=2,sort_keys=True));return 0
    rendered=json.dumps(build_result(a.contract),indent=2,sort_keys=True)+"\n"
    if a.output:a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(rendered)
    else:print(rendered,end="")
    return 0
if __name__=="__main__":raise SystemExit(main())
