#!/usr/bin/env python3
"""Exact principal left, right, and two-sided divisibility posets."""
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
SCHEMA="matching-one/terminal-partition-serial-divisibility-posets/v1"
def data():
 st=enumerate_rgs(4);idx={s:i for i,s in enumerate(st)};p=tuple(tuple(idx[serial_compose(a,b)] for b in st) for a in st)
 left=[frozenset(p[s][a] for s in range(15)) for a in range(15)];right=[frozenset(p[a][s] for s in range(15)) for a in range(15)];two=[]
 for a in range(15):two.append(frozenset(p[p[s][a]][t] for s in range(15) for t in range(15)))
 return p,left,right,two
def poset(name,ideals):
 carriers=sorted(set(ideals),key=lambda s:(len(s),sorted(s)));classes=[[i for i,x in enumerate(ideals) if x==s] for s in carriers];covers=[]
 for i,a in enumerate(carriers):
  for j,b in enumerate(carriers):
   if a<b and not any(a<c<b for c in carriers):covers.append([i,j])
 return {"name":name,"class_count":len(carriers),"classes":classes,"principal_ideals":[sorted(s) for s in carriers],"cover_relations":covers}
def build_artifact()->dict[str,Any]:
 _,l,r,j=data();lp,rp,jp=poset("left",l),poset("right",r),poset("two_sided",j)
 return {"schema":SCHEMA,"issue":13,"status":"complete_principal_divisibility_posets","left":lp,"right":rp,"two_sided":jp,
 "exact_checks":{"left_and_right_have_six_classes":lp["class_count"]==rp["class_count"]==6,"two_sided_has_three_classes":jp["class_count"]==3,"left_and_right_cover_profiles_match":lp["cover_relations"]==rp["cover_relations"]==[[0,2],[1,3],[1,4],[2,5],[3,5],[4,5]],"two_sided_is_three_element_chain":jp["cover_relations"]==[[0,1],[1,2]]},
 "claim_boundary":{"included":"principal left, right, and two-sided ideal classes and their inclusion covers","excluded":"full one-sided ideal lattices already certified, congruences, planarity, reliability, or thresholds","parent_issue":"remain open"}}
def validate_artifact(v:Mapping[str,Any]):
 e=build_artifact()
 if v!=e:raise ValueError("divisibility artifact mismatch")
 if set(e["exact_checks"].values())!={True}:raise ValueError("checks failed")
 return {"schema":SCHEMA,"status":"valid","class_counts":[6,6,3]}
def main(argv:Optional[Sequence[str]]=None)->int:
 q=argparse.ArgumentParser(description=__doc__);q.add_argument("--output",type=Path);q.add_argument("--validate",type=Path);a=q.parse_args(argv)
 if a.validate:print(json.dumps(validate_artifact(json.loads(a.validate.read_text())),indent=2,sort_keys=True));return 0
 s=json.dumps(build_artifact(),indent=2,sort_keys=True)+"\n"
 if a.output:a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(s,encoding="utf-8")
 else:print(s,end="")
 return 0
if __name__=="__main__":raise SystemExit(main())
