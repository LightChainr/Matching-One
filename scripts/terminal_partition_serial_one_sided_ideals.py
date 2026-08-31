#!/usr/bin/env python3
"""Complete left and right ideal lattices of the typed serial monoid."""
from __future__ import annotations
import argparse,json
from collections import Counter
from pathlib import Path
from typing import Any,Mapping,Optional,Sequence
try:
 from scripts.terminal_partition_canonical import enumerate_rgs
 from scripts.terminal_partition_serial_category import serial_compose
 from scripts.terminal_partition_serial_reversal import reverse_ports
except ModuleNotFoundError:
 from terminal_partition_canonical import enumerate_rgs
 from terminal_partition_serial_category import serial_compose
 from terminal_partition_serial_reversal import reverse_ports
ROOT=Path(__file__).resolve().parents[1];SCHEMA="matching-one/terminal-partition-serial-one-sided-ideals/v1"
def table():
 states=enumerate_rgs(4);idx={s:i for i,s in enumerate(states)}
 return tuple(tuple(idx[serial_compose(a,b)] for b in states) for a in states)
def covers(items):return [(a,b) for a in items for b in items if a<b and not any(a<c<b for c in items)]
def height(items):
 ordered=sorted(items,key=lambda x:(len(x),sorted(x)));dp={x:1 for x in ordered}
 for x in ordered:
  for y in ordered:
   if y<x:dp[x]=max(dp[x],dp[y]+1)
 return max(dp.values())
def build_artifact()->dict[str,Any]:
 states=enumerate_rgs(4);idx={s:i for i,s in enumerate(states)};p=table();rev=tuple(idx[reverse_ports(s)] for s in states)
 left=[];right=[]
 for mask in range(1<<15):
  s=frozenset(i for i in range(15) if mask>>i&1)
  if all(p[a][x] in s for a in range(15) for x in s):left.append(s)
  if all(p[x][a] in s for a in range(15) for x in s):right.append(s)
 left=sorted(left,key=lambda x:(len(x),sorted(x)));right=sorted(right,key=lambda x:(len(x),sorted(x)));lc=covers(left);rc=covers(right)
 return {"schema":SCHEMA,"issue":13,"status":"complete_one_sided_ideal_lattices","left_ideals":[sorted(x) for x in left],"right_ideals":[sorted(x) for x in right],
  "left_size_histogram":{str(k):v for k,v in sorted(Counter(map(len,left)).items())},"right_size_histogram":{str(k):v for k,v in sorted(Counter(map(len,right)).items())},
  "left_cover_relations":[[sorted(a),sorted(b)] for a,b in lc],"right_cover_relations":[[sorted(a),sorted(b)] for a,b in rc],"left_height":height(left),"right_height":height(right),
  "exact_checks":{"exactly_sixteen_left_and_right_ideals":len(left)==len(right)==16,"each_cover_graph_has_26_edges":len(lc)==len(rc)==26,
   "both_lattices_have_height_seven":height(left)==height(right)==7,"reversal_bijects_left_to_right":{frozenset(rev[x] for x in s) for s in left}==set(right),
   "size_histograms_agree":Counter(map(len,left))==Counter(map(len,right))},
  "claim_boundary":{"included":"all left/right ideals, inclusion covers, heights, and reversal anti-isomorphism","excluded":"congruence lattice, principal-factor representations, planarity, reliability, or thresholds","parent_issue":"remain open"}}
def validate_artifact(v:Mapping[str,Any]):
 e=build_artifact()
 if v!=e:raise ValueError("one-sided ideal artifact does not reproduce")
 if set(e["exact_checks"].values())!={True}:raise ValueError("checks failed")
 return {"schema":SCHEMA,"status":"valid","left_ideals":16,"right_ideals":16}
def main(argv:Optional[Sequence[str]]=None)->int:
 q=argparse.ArgumentParser(description=__doc__);q.add_argument("--output",type=Path);q.add_argument("--validate",type=Path);a=q.parse_args(argv)
 if a.validate:print(json.dumps(validate_artifact(json.loads(a.validate.read_text())),indent=2,sort_keys=True));return 0
 s=json.dumps(build_artifact(),indent=2,sort_keys=True)+"\n"
 if a.output:a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(s,encoding="utf-8")
 else:print(s,end="")
 return 0
if __name__=="__main__":raise SystemExit(main())
