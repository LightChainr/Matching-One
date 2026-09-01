#!/usr/bin/env python3
"""Exact inclusion-lattice invariants of all typed serial subsemigroups."""
from __future__ import annotations
import argparse,json
from collections import deque
from pathlib import Path
from typing import Any,Mapping,Optional,Sequence
try:
 from scripts.terminal_partition_canonical import enumerate_rgs
 from scripts.terminal_partition_serial_category import serial_compose
except ModuleNotFoundError:
 from terminal_partition_canonical import enumerate_rgs
 from terminal_partition_serial_category import serial_compose
ROOT=Path(__file__).resolve().parents[1];SCHEMA="matching-one/terminal-partition-serial-subsemigroup-lattice/v1"
def subsemigroups():
 st=enumerate_rgs(4);idx={s:i for i,s in enumerate(st)};p=[[idx[serial_compose(a,b)] for b in st] for a in st];out=[]
 for mask in range(1<<15):
  s=frozenset(i for i in range(15) if mask>>i&1)
  if all(p[a][b] in s for a in s for b in s):out.append(s)
 return out
def width(items):
 adj={u:[v for v in items if u<v] for u in items};pu={u:None for u in items};pv={v:None for v in items};match=0;INF=10**9
 while True:
  d={};q=deque()
  for u in items:d[u]=0 if pu[u] is None else INF;q.append(u) if pu[u] is None else None
  found=False
  while q:
   u=q.popleft()
   for v in adj[u]:
    w=pv[v]
    if w is None:found=True
    elif d[w]==INF:d[w]=d[u]+1;q.append(w)
  if not found:break
  def dfs(u):
   for v in adj[u]:
    w=pv[v]
    if w is None or (d[w]==d[u]+1 and dfs(w)):pu[u]=v;pv[v]=u;return True
   d[u]=INF;return False
  for u in items:
   if pu[u] is None and dfs(u):match+=1
 return len(items)-match,match
def build_artifact()->dict[str,Any]:
 items=subsemigroups();covers=sum(a<b and not any(a<c<b for c in items) for a in items for b in items);ordered=sorted(items,key=len);dp={x:1 for x in ordered}
 for x in ordered:
  for y in ordered:
   if y<x:dp[x]=max(dp[x],dp[y]+1)
 w,m=width(items)
 return {"schema":SCHEMA,"issue":13,"status":"exact_subsemigroup_inclusion_lattice","elements":len(items),"cover_count":covers,"height":max(dp.values()),"width":w,"dilworth_matching_size":m,
 "exact_checks":{"lattice_has_416_elements":len(items)==416,"cover_count_is_1400":covers==1400,"height_is_twelve":max(dp.values())==12,"width_is_eighty_two":w==82 and m==334},
 "claim_boundary":{"included":"inclusion covers, longest chains, and exact Dilworth width","excluded":"congruence lattice, symmetry orbits, planarity, reliability, or thresholds","parent_issue":"remain open"}}
def validate_artifact(v:Mapping[str,Any]):
 e=build_artifact()
 if v!=e:raise ValueError("lattice artifact mismatch")
 if set(e["exact_checks"].values())!={True}:raise ValueError("checks failed")
 return {"schema":SCHEMA,"status":"valid","height":12,"width":82}
def main(argv:Optional[Sequence[str]]=None)->int:
 q=argparse.ArgumentParser(description=__doc__);q.add_argument("--output",type=Path);q.add_argument("--validate",type=Path);a=q.parse_args(argv)
 if a.validate:print(json.dumps(validate_artifact(json.loads(a.validate.read_text())),indent=2,sort_keys=True));return 0
 s=json.dumps(build_artifact(),indent=2,sort_keys=True)+"\n";a.output.write_text(s) if a.output else print(s,end="");return 0
if __name__=="__main__":raise SystemExit(main())
