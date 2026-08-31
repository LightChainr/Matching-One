#!/usr/bin/env python3
"""Exact symmetry orbits of all typed serial subsemigroups."""
from __future__ import annotations
import argparse,json
from collections import Counter
from pathlib import Path
from typing import Any,Mapping,Optional,Sequence
try:
 from scripts.terminal_partition_canonical import enumerate_rgs,apply_permutation
 from scripts.terminal_partition_serial_category import serial_compose
 from scripts.terminal_partition_serial_reversal import reverse_ports
except ModuleNotFoundError:
 from terminal_partition_canonical import enumerate_rgs,apply_permutation
 from terminal_partition_serial_category import serial_compose
 from terminal_partition_serial_reversal import reverse_ports
ROOT=Path(__file__).resolve().parents[1];SCHEMA="matching-one/terminal-partition-serial-subsemigroup-symmetry/v1"
def data():
 states=enumerate_rgs(4);idx={s:i for i,s in enumerate(states)};p=tuple(tuple(idx[serial_compose(a,b)] for b in states) for a in states)
 sems=[]
 for mask in range(1<<15):
  s=frozenset(i for i in range(15) if mask>>i&1)
  if all(p[a][b] in s for a in s for b in s):sems.append(s)
 rev=tuple(idx[reverse_ports(s)] for s in states);lane=tuple(idx[apply_permutation(s,(1,0,3,2))] for s in states);both=tuple(rev[lane[i]] for i in range(15))
 return sems,(tuple(range(15)),rev,lane,both),rev,lane
def image(s,m):return frozenset(m[x] for x in s)
def build_artifact()->dict[str,Any]:
 sems,maps,rev,lane=data();universe=set(sems);seen=set();orbits=[]
 for s in sorted(sems,key=lambda x:(len(x),sorted(x))):
  if s in seen:continue
  o={image(s,m) for m in maps};seen|=o;orbits.append(sorted((sorted(x) for x in o),key=lambda x:(len(x),x)))
 hist=Counter(map(len,orbits));rs=sum(image(s,rev)==s for s in sems);ls=sum(image(s,lane)==s for s in sems)
 return {"schema":SCHEMA,"issue":13,"status":"complete_subsemigroup_symmetry_orbits","orbit_count":len(orbits),"orbit_size_histogram":{str(k):v for k,v in sorted(hist.items())},"reversal_stable_count":rs,"lane_swap_stable_count":ls,"orbits":orbits,
 "exact_checks":{"orbits_partition_all_416_subsemigroups":seen==universe,"orbit_profile_is_exact":hist=={1:34,2:71,4:60},"stable_counts_are_exact":rs==54 and ls==128},
 "claim_boundary":{"included":"lane-swap and port-reversal action on all subsemigroups","excluded":"submonoid orbit certificate already committed, congruences, planarity, reliability, or thresholds","parent_issue":"remain open"}}
def validate_artifact(v:Mapping[str,Any]):
 e=build_artifact()
 if v!=e:raise ValueError("symmetry artifact mismatch")
 if set(e["exact_checks"].values())!={True}:raise ValueError("checks failed")
 return {"schema":SCHEMA,"status":"valid","orbits":165}
def main(argv:Optional[Sequence[str]]=None)->int:
 q=argparse.ArgumentParser(description=__doc__);q.add_argument("--output",type=Path);q.add_argument("--validate",type=Path);a=q.parse_args(argv)
 if a.validate:print(json.dumps(validate_artifact(json.loads(a.validate.read_text())),indent=2,sort_keys=True));return 0
 s=json.dumps(build_artifact(),indent=2,sort_keys=True)+"\n";a.output.write_text(s) if a.output else print(s,end="");return 0
if __name__=="__main__":raise SystemExit(main())
