#!/usr/bin/env python3
"""Complete subsemigroup census of typed serial composition."""
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
ROOT=Path(__file__).resolve().parents[1];SCHEMA="matching-one/terminal-partition-serial-subsemigroups/v1"
def table():
 states=enumerate_rgs(4);idx={s:i for i,s in enumerate(states)}
 return tuple(tuple(idx[serial_compose(a,b)] for b in states) for a in states)
def build_artifact()->dict[str,Any]:
 p=table();values=[]
 for mask in range(1<<15):
  s=frozenset(i for i in range(15) if mask>>i&1)
  if all(p[a][b] in s for a in s for b in s):values.append(s)
 values=sorted(values,key=lambda x:(len(x),sorted(x)));proper=[s for s in values if len(s)<15]
 maximal=sorted((sorted(s) for s in proper if not any(s<t and len(t)<15 for t in values)),key=lambda x:(len(x),x));idem={i for i in range(15) if p[i][i]==i}
 hist=Counter(map(len,values));without_identity=sum(6 not in s for s in values)
 return {"schema":SCHEMA,"issue":13,"status":"complete_subsemigroup_census","subsemigroup_count":len(values),"subsemigroups":[sorted(s) for s in values],
  "size_histogram":{str(k):v for k,v in sorted(hist.items())},"without_wire_identity_count":without_identity,"with_wire_identity_count":len(values)-without_identity,
  "maximal_proper_subsemigroups":maximal,
  "exact_checks":{"all_32768_subsets_yield_416_subsemigroups":len(values)==416,"exactly_228_contain_wire_identity":len(values)-without_identity==228,
   "exactly_188_exclude_wire_identity":without_identity==188,"exactly_five_maximal_proper_subsemigroups":len(maximal)==5,
   "every_nonempty_subsemigroup_contains_an_idempotent":all(not s or bool(s&idem) for s in values)},
  "claim_boundary":{"included":"all multiplication-closed subsets, including empty and identity-free subsemigroups, plus maximal proper sectors","excluded":"submonoid symmetry orbits already certified, congruences, representations, planarity, reliability, or thresholds","parent_issue":"remain open"}}
def validate_artifact(v:Mapping[str,Any]):
 e=build_artifact()
 if v!=e:raise ValueError("subsemigroup artifact does not reproduce")
 if set(e["exact_checks"].values())!={True}:raise ValueError("checks failed")
 return {"schema":SCHEMA,"status":"valid","subsemigroups":416}
def main(argv:Optional[Sequence[str]]=None)->int:
 q=argparse.ArgumentParser(description=__doc__);q.add_argument("--output",type=Path);q.add_argument("--validate",type=Path);a=q.parse_args(argv)
 if a.validate:print(json.dumps(validate_artifact(json.loads(a.validate.read_text())),indent=2,sort_keys=True));return 0
 s=json.dumps(build_artifact(),indent=2,sort_keys=True)+"\n"
 if a.output:a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(s,encoding="utf-8")
 else:print(s,end="")
 return 0
if __name__=="__main__":raise SystemExit(main())
