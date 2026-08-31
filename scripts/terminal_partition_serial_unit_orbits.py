#!/usr/bin/env python3
"""Exact unit-group orbit decompositions of typed serial partitions."""
from __future__ import annotations
import argparse,json
from pathlib import Path
from typing import Any,Mapping,Optional,Sequence
try:
 from scripts.terminal_partition_canonical import enumerate_rgs
 from scripts.terminal_partition_serial_category import serial_compose
except ModuleNotFoundError:
 from terminal_partition_canonical import enumerate_rgs
 from terminal_partition_serial_category import serial_compose
ROOT=Path(__file__).resolve().parents[1];SCHEMA="matching-one/terminal-partition-serial-unit-orbits/v1"
def build_artifact()->dict[str,Any]:
 st=enumerate_rgs(4);idx={s:i for i,s in enumerate(st)};p=[[idx[serial_compose(a,b)] for b in st] for a in st];u=(6,8)
 def orbits(fn):
  seen=set();out=[]
  for x in range(15):
   if x in seen:continue
   o={fn(a,x,b) for a in u for b in u};seen|=o;out.append(sorted(o))
  return sorted(out,key=lambda x:(len(x),x))
 left=orbits(lambda a,x,b:p[a][x]);right=orbits(lambda a,x,b:p[x][b]);double=orbits(lambda a,x,b:p[p[a][x]][b]);conj=orbits(lambda a,x,b:p[p[a][x]][a])
 return {"schema":SCHEMA,"issue":13,"status":"exact_unit_action_orbits","unit_indices":list(u),"left_orbits":left,"right_orbits":right,"double_orbits":double,"conjugation_orbits":conj,
 "exact_checks":{"left_and_right_have_eleven_orbits":len(left)==len(right)==11,"double_action_has_nine_orbits":len(double)==9,"conjugation_has_eleven_orbits":len(conj)==11,"four_decompositions_are_not_all_equal":len({json.dumps(x) for x in (left,right,double,conj)})==4},
 "claim_boundary":{"included":"left, right, double, and inner-conjugation actions of the C2 unit group","excluded":"full automorphism action already certified, semigroup conjugacy, planarity, reliability, or thresholds","parent_issue":"remain open"}}
def validate_artifact(v:Mapping[str,Any]):
 e=build_artifact()
 if v!=e:raise ValueError("unit orbit artifact mismatch")
 if set(e["exact_checks"].values())!={True}:raise ValueError("checks failed")
 return {"schema":SCHEMA,"status":"valid","double_orbits":9}
def main(argv:Optional[Sequence[str]]=None)->int:
 q=argparse.ArgumentParser(description=__doc__);q.add_argument("--output",type=Path);q.add_argument("--validate",type=Path);a=q.parse_args(argv)
 if a.validate:print(json.dumps(validate_artifact(json.loads(a.validate.read_text())),indent=2,sort_keys=True));return 0
 s=json.dumps(build_artifact(),indent=2,sort_keys=True)+"\n";a.output.write_text(s) if a.output else print(s,end="");return 0
if __name__=="__main__":raise SystemExit(main())
