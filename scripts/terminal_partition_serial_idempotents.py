#!/usr/bin/env python3
"""Exact idempotent-generated sector of the typed serial monoid."""
from __future__ import annotations
import argparse,json
from itertools import combinations
from pathlib import Path
from typing import Any,Mapping,Optional,Sequence
try:
 from scripts.terminal_partition_canonical import enumerate_rgs
 from scripts.terminal_partition_serial_category import serial_compose
except ModuleNotFoundError:
 from terminal_partition_canonical import enumerate_rgs
 from terminal_partition_serial_category import serial_compose
ROOT=Path(__file__).resolve().parents[1];SCHEMA="matching-one/terminal-partition-serial-idempotents/v1";IDENTITY=6
def table():
 states=enumerate_rgs(4);index={s:i for i,s in enumerate(states)}
 return tuple(tuple(index[serial_compose(a,b)] for b in states) for a in states)
def closure(seed,product):
 if any(type(x) is not int or not 0<=x<15 for x in seed): raise ValueError("generator index outside table")
 out={IDENTITY,*seed};changed=True
 while changed:
  changed=False
  for a in tuple(out):
   for b in tuple(out):
    if product[a][b] not in out: out.add(product[a][b]);changed=True
 return frozenset(out)
def build_artifact()->dict[str,Any]:
 product=table(); idem=[i for i in range(15) if product[i][i]==i]; sector=closure(idem,product); minimum=[];rank=None
 for k in range(len(idem)+1):
  minimum=[s for s in combinations(idem,k) if closure(s,product)==sector]
  if minimum: rank=k;break
 return {"schema":SCHEMA,"issue":13,"status":"exact_idempotent_generated_sector","idempotent_indices":idem,
  "idempotent_generated_sector":sorted(sector),"excluded_indices":sorted(set(range(15))-sector),"idempotent_rank":rank,
  "minimum_idempotent_generating_sets":[list(s) for s in minimum],
  "exact_checks":{"exactly_twelve_idempotents":len(idem)==12,"idempotents_generate_exactly_fourteen_states":len(sector)==14,
   "nontrivial_unit_is_not_idempotent_generated":set(range(15))-sector=={8},"idempotent_rank_is_three":rank==3,
   "exactly_three_minimum_idempotent_sets":len(minimum)==3},
  "claim_boundary":{"included":"idempotent set, its generated submonoid, idempotent rank, and all minimum idempotent generating sets","excluded":"all-element rank beyond the committed certificate, endomorphisms, congruences, planarity, reliability, or thresholds","parent_issue":"remain open"}}
def validate_artifact(value:Mapping[str,Any]):
 expected=build_artifact()
 if value!=expected: raise ValueError("idempotent artifact does not exactly reproduce")
 if set(expected["exact_checks"].values())!={True}: raise ValueError("all exact checks must pass")
 return {"schema":SCHEMA,"status":"valid","idempotents":12,"generated_states":14,"rank":3}
def main(argv:Optional[Sequence[str]]=None)->int:
 p=argparse.ArgumentParser(description=__doc__);p.add_argument("--output",type=Path);p.add_argument("--validate",type=Path);a=p.parse_args(argv)
 if a.validate: print(json.dumps(validate_artifact(json.loads(a.validate.read_text())),indent=2,sort_keys=True));return 0
 rendered=json.dumps(build_artifact(),indent=2,sort_keys=True)+"\n"
 if a.output:a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(rendered,encoding="utf-8")
 else:print(rendered,end="")
 return 0
if __name__=="__main__":raise SystemExit(main())
