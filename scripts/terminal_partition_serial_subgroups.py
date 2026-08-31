#!/usr/bin/env python3
"""Complete subgroup census of the typed serial monoid."""
from __future__ import annotations
import argparse,json
from collections import Counter
from pathlib import Path
from typing import Any,Mapping,Optional,Sequence
try:
 from scripts.terminal_partition_canonical import enumerate_rgs
 from scripts.terminal_partition_serial_category import serial_compose
except ModuleNotFoundError:
 from terminal_partition_canonical import enumerate_rgs
 from terminal_partition_serial_category import serial_compose
ROOT=Path(__file__).resolve().parents[1];SCHEMA="matching-one/terminal-partition-serial-subgroups/v1"
def table():
 states=enumerate_rgs(4);idx={s:i for i,s in enumerate(states)}
 return tuple(tuple(idx[serial_compose(a,b)] for b in states) for a in states)
def build_artifact()->dict[str,Any]:
 p=table();groups=[]
 for mask in range(1,1<<15):
  s={i for i in range(15) if mask>>i&1}
  if not all(p[a][b] in s for a in s for b in s):continue
  identities=[e for e in s if all(p[e][x]==x==p[x][e] for x in s)]
  if len(identities)!=1:continue
  e=identities[0]
  if all(any(p[x][y]==e==p[y][x] for y in s) for x in s):groups.append({"carrier":sorted(s),"identity":e})
 groups.sort(key=lambda x:(len(x["carrier"]),x["carrier"]))
 hist=Counter(len(x["carrier"]) for x in groups)
 return {"schema":SCHEMA,"issue":13,"status":"complete_subgroup_census","subgroups":groups,"subgroup_count":len(groups),
  "size_histogram":{str(k):v for k,v in sorted(hist.items())},
  "exact_checks":{"exactly_twelve_trivial_subgroups":hist[1]==12,"unique_nontrivial_subgroup":[x["carrier"] for x in groups if len(x["carrier"])>1]==[[6,8]],
   "nontrivial_subgroup_is_c2":p[8][8]==6,"every_group_identity_is_idempotent":all(p[x["identity"]][x["identity"]]==x["identity"] for x in groups)},
  "claim_boundary":{"included":"all internal subgroups with their own identities","excluded":"Schutzenberger groups, external representations, congruences, planarity, reliability, or thresholds","parent_issue":"remain open"}}
def validate_artifact(v:Mapping[str,Any]):
 e=build_artifact()
 if v!=e:raise ValueError("subgroup artifact does not reproduce")
 if set(e["exact_checks"].values())!={True}:raise ValueError("checks failed")
 return {"schema":SCHEMA,"status":"valid","subgroups":13}
def main(argv:Optional[Sequence[str]]=None)->int:
 q=argparse.ArgumentParser(description=__doc__);q.add_argument("--output",type=Path);q.add_argument("--validate",type=Path);a=q.parse_args(argv)
 if a.validate:print(json.dumps(validate_artifact(json.loads(a.validate.read_text())),indent=2,sort_keys=True));return 0
 s=json.dumps(build_artifact(),indent=2,sort_keys=True)+"\n"
 if a.output:a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(s,encoding="utf-8")
 else:print(s,end="")
 return 0
if __name__=="__main__":raise SystemExit(main())
