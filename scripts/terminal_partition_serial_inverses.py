#!/usr/bin/env python3
"""Exact generalized-inverse census for the typed serial monoid."""
from __future__ import annotations
import argparse, json
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence
try:
    from scripts.terminal_partition_canonical import enumerate_rgs
    from scripts.terminal_partition_serial_category import serial_compose
except ModuleNotFoundError:
    from terminal_partition_canonical import enumerate_rgs
    from terminal_partition_serial_category import serial_compose
ROOT=Path(__file__).resolve().parents[1]; SCHEMA="matching-one/terminal-partition-serial-inverses/v1"

def table():
    states=enumerate_rgs(4); index={s:i for i,s in enumerate(states)}
    return tuple(tuple(index[serial_compose(a,b)] for b in states) for a in states)

def inverse_sets(product):
    return [tuple(b for b in range(15) if product[product[a][b]][a]==a and product[product[b][a]][b]==b) for a in range(15)]

def build_artifact()->dict[str,Any]:
    product=table(); inverses=inverse_sets(product); regular=[i for i,v in enumerate(inverses) if v]
    idempotents=[i for i in range(15) if product[i][i]==i]
    units=[i for i in range(15) if any(product[i][j]==6==product[j][i] for j in range(15))]
    return {"schema":SCHEMA,"issue":13,"status":"exact_generalized_inverse_census",
      "inverse_sets":[list(v) for v in inverses],"inverse_count_profile":[len(v) for v in inverses],
      "regular_element_indices":regular,"idempotent_indices":idempotents,"unit_indices":units,
      "exact_checks":{
        "all_elements_are_regular":regular==list(range(15)),
        "inverse_relation_is_symmetric":all((b in inverses[a])==(a in inverses[b]) for a in range(15) for b in range(15)),
        "inverse_count_profile_is_exact":[len(v) for v in inverses]==[9,6,6,4,4,6,1,4,1,6,4,4,4,4,4],
        "units_form_the_order_two_maximal_subgroup":units==[6,8] and product[8][8]==6,
        "exactly_twelve_elements_are_idempotent":len(idempotents)==12},
      "claim_boundary":{"included":"generalized inverses, regularity, idempotents, and units of the finite typed serial monoid","excluded":"endomorphism or congruence lattices, planar realization, reliability, or thresholds","parent_issue":"remain open"}}

def validate_artifact(value:Mapping[str,Any]):
    expected=build_artifact()
    if value!=expected: raise ValueError("inverse artifact does not exactly reproduce")
    if set(expected["exact_checks"].values())!={True}: raise ValueError("all exact checks must pass")
    return {"schema":SCHEMA,"status":"valid","regular_elements":15,"units":2}

def main(argv:Optional[Sequence[str]]=None)->int:
    p=argparse.ArgumentParser(description=__doc__); p.add_argument("--output",type=Path); p.add_argument("--validate",type=Path); a=p.parse_args(argv)
    if a.validate: print(json.dumps(validate_artifact(json.loads(a.validate.read_text())),indent=2,sort_keys=True)); return 0
    rendered=json.dumps(build_artifact(),indent=2,sort_keys=True)+"\n"
    if a.output: a.output.parent.mkdir(parents=True,exist_ok=True); a.output.write_text(rendered,encoding="utf-8")
    else: print(rendered,end="")
    return 0
if __name__=="__main__": raise SystemExit(main())
