#!/usr/bin/env python3
"""Seeded exact-decimal false-positive calibration for the degree-one search."""

from __future__ import annotations
import argparse, hashlib, json, random
from fractions import Fraction
from functools import lru_cache
from math import gcd
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence

ROOT=Path(__file__).resolve().parents[1]
DEFAULT_CONTRACT=ROOT/"analysis"/"pslq_search_contract.json"
DEFAULT_OUTPUT=ROOT/"results"/"pslq-synthetic-false-positive-calibration"/"latest.json"
SCHEMA="matching-one/pslq-synthetic-false-positive-calibration/v1"

def _text(x:Fraction)->str:return f"{x.numerator}/{x.denominator}"
def synthetic_values(count:int,seed:int)->list[Fraction]:
    rng=random.Random(seed);den=10**15;lower=55*10**13;width=10**14
    return [Fraction(lower+rng.randrange(width+1),den) for _ in range(count)]
def closest_degree_one(value:Fraction,height:int)->tuple[Fraction,tuple[int,int]]:
    best=None
    for a1 in range(1,height+1):
        for a0 in range(-height,height+1):
            if gcd(abs(a0),a1)!=1:continue
            root_distance=abs(Fraction(-a0,a1)-value)
            candidate=(root_distance,(a0,a1))
            if best is None or candidate<best:best=candidate
    assert best is not None;return best
def _decimal(value:Fraction)->str:
    scale=10**15;scaled=value.numerator*scale//value.denominator
    return f"{scaled//scale}.{scaled%scale:015d}"
@lru_cache(maxsize=2)
def build_result(contract_path:Path=DEFAULT_CONTRACT)->dict[str,Any]:
    raw=contract_path.read_bytes();contract=json.loads(raw);spec=contract["false_positive_controls"]["synthetic_random_constants"]
    if spec!={"count":100,"seed":20260830,"domain":"uniform decimal strings on [0.55,0.65]"}:raise ValueError("synthetic-control contract drift")
    height=contract["search_stages"]["algebraic_polynomial"]["coefficient_height_max"];floor=Fraction(contract["input_policy"]["comparison_resolution_floor"])
    rows=[]
    for index,value in enumerate(synthetic_values(spec["count"],spec["seed"])):
        distance,coefficients=closest_degree_one(value,height)
        rows.append({"index":index,"value":_decimal(value),"closest_coefficients_ascending":list(coefficients),"root_distance":_text(distance),"within_resolution_floor":distance<=floor})
    distances=sorted(Fraction(row["root_distance"]) for row in rows)
    return {"schema":SCHEMA,"issue":1,"status":"seeded_degree1_false_positive_calibration_complete","contract_sha256":hashlib.sha256(raw).hexdigest(),
      "generator":{"seed":spec["seed"],"count":spec["count"],"domain":spec["domain"],"decimal_places":15,"algorithm":"random.Random.randrange over the closed 10^14 grid"},
      "search":{"degree":1,"coefficient_height_max":height,"primitive_only":True,"resolution_floor":contract["input_policy"]["comparison_resolution_floor"]},
      "rows":rows,"summary":{"resolution_floor_hits":sum(row["within_resolution_floor"] for row in rows),"minimum_root_distance":_text(distances[0]),"median_root_distance":_text(distances[len(distances)//2]),"maximum_root_distance":_text(distances[-1])},
      "claim_boundary":{"included":"seeded 100-value degree-one finite-height false-positive calibration","excluded":"probability model for p_c, degree-2/3/4 calibration, p-values, near-hit promotion, closed forms, or transcendence","parent_issue":"remain open"}}
def validate_result(result:Mapping[str,Any],contract_path:Path=DEFAULT_CONTRACT)->Mapping[str,Any]:
    expected=build_result(contract_path)
    if result!=expected:raise ValueError("synthetic calibration does not exactly reproduce")
    return {"schema":SCHEMA,"status":"valid","controls":len(expected["rows"]),"hits":expected["summary"]["resolution_floor_hits"]}
def main(argv:Optional[Sequence[str]]=None)->int:
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("--contract",type=Path,default=DEFAULT_CONTRACT);p.add_argument("--output",type=Path);p.add_argument("--validate",type=Path);a=p.parse_args(argv)
    if a.validate:print(json.dumps(validate_result(json.loads(a.validate.read_text()),a.contract),indent=2,sort_keys=True));return 0
    rendered=json.dumps(build_result(a.contract),indent=2,sort_keys=True)+"\n"
    if a.output:a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(rendered)
    else:print(rendered,end="")
    return 0
if __name__=="__main__":raise SystemExit(main())
