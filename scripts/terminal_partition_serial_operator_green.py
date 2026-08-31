#!/usr/bin/env python3
"""Exact Green relations of the two-sided operator semigroup."""
from __future__ import annotations
import argparse,json
from collections import Counter
from pathlib import Path
from typing import Any,Mapping,Optional,Sequence
try:
 from scripts.terminal_partition_serial_two_sided_actions import table,operator,compose
except ModuleNotFoundError:
 from terminal_partition_serial_two_sided_actions import table,operator,compose
SCHEMA="matching-one/terminal-partition-serial-operator-green/v1"
def data():
 p=table();ops=sorted({operator(a,b,p) for a in range(15) for b in range(15)});ix={f:i for i,f in enumerate(ops)};m=tuple(tuple(ix[compose(f,g)] for g in ops) for f in ops);return ops,m
def classes(keys):
 d={}
 for i,k in enumerate(keys):d.setdefault(k,[]).append(i)
 return sorted((sorted(x) for x in d.values()),key=lambda x:(len(x),x))
def profile(c):return {str(k):v for k,v in sorted(Counter(map(len,c)).items())}
def build_artifact()->dict[str,Any]:
 _,m=data();n=len(m);L=[frozenset(m[s][a] for s in range(n)) for a in range(n)];R=[frozenset(m[a][s] for s in range(n)) for a in range(n)];J=[]
 for a in range(n):J.append(frozenset(m[m[s][a]][t] for s in range(n) for t in range(n)))
 lc,rc,hc,jc=classes(L),classes(R),classes([(L[i],R[i]) for i in range(n)]),classes(J);parent=list(range(n))
 def find(x):
  while parent[x]!=x:parent[x]=parent[parent[x]];x=parent[x]
  return x
 def union(a,b):
  a,b=find(a),find(b)
  if a!=b:parent[b]=a
 for c in lc+rc:
  for x in c[1:]:union(c[0],x)
 dc=classes([find(i) for i in range(n)])
 return {"schema":SCHEMA,"issue":13,"status":"complete_two_sided_operator_green_relations","L_classes":lc,"R_classes":rc,"H_classes":hc,"J_classes":jc,"D_classes":dc,"class_counts":{"L":len(lc),"R":len(rc),"H":len(hc),"J":len(jc),"D":len(dc)},"class_size_histograms":{"L":profile(lc),"R":profile(rc),"H":profile(hc),"J":profile(jc),"D":profile(dc)},
 "exact_checks":{"class_counts_are_exact":(len(lc),len(rc),len(hc),len(jc),len(dc))==(21,24,112,7,7),"L_profile_is_exact":profile(lc)=={"2":4,"4":2,"6":6,"9":9},"R_profile_is_exact":profile(rc)=={"1":4,"2":4,"4":1,"6":6,"9":9},"H_profile_is_exact":profile(hc)=={"1":93,"2":18,"4":1},"D_equals_J":{frozenset(x) for x in dc}=={frozenset(x) for x in jc} and profile(jc)=={"4":4,"18":2,"81":1}},
 "claim_boundary":{"included":"complete L, R, H, J, and D classes of the 133 two-sided operators","excluded":"base-monoid Green relations already certified, congruences, physical transfer operators, reliability, or thresholds","parent_issue":"remain open"}}
def validate_artifact(v:Mapping[str,Any]):
 e=build_artifact()
 if v!=e:raise ValueError("operator Green artifact mismatch")
 if set(e["exact_checks"].values())!={True}:raise ValueError("checks failed")
 return {"schema":SCHEMA,"status":"valid","J_classes":7}
def main(argv:Optional[Sequence[str]]=None)->int:
 q=argparse.ArgumentParser(description=__doc__);q.add_argument("--output",type=Path);q.add_argument("--validate",type=Path);a=q.parse_args(argv)
 if a.validate:print(json.dumps(validate_artifact(json.loads(a.validate.read_text())),indent=2,sort_keys=True));return 0
 s=json.dumps(build_artifact(),indent=2,sort_keys=True)+"\n"
 if a.output:a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(s,encoding="utf-8")
 else:print(s,end="")
 return 0
if __name__=="__main__":raise SystemExit(main())
