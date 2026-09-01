#!/usr/bin/env python3
"""Faithful left and right regular actions of typed serial composition."""
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
ROOT=Path(__file__).resolve().parents[1];SCHEMA="matching-one/terminal-partition-serial-regular-actions/v1"
def table():
 states=enumerate_rgs(4);idx={s:i for i,s in enumerate(states)}
 return tuple(tuple(idx[serial_compose(a,b)] for b in states) for a in states)
def compose(f,g):
 if len(f)!=len(g):raise ValueError("transformation widths differ")
 return tuple(f[g[i]] for i in range(len(f)))
def build_artifact()->dict[str,Any]:
 p=table();left=[tuple(p[a][x] for x in range(15)) for a in range(15)];right=[tuple(p[x][a] for x in range(15)) for a in range(15)]
 lh=Counter(len(set(x)) for x in left);rh=Counter(len(set(x)) for x in right);intersection=[(i,j) for i,a in enumerate(left) for j,b in enumerate(right) if a==b]
 return {"schema":SCHEMA,"issue":13,"status":"faithful_left_right_regular_actions","left_transformations":[list(x) for x in left],"right_transformations":[list(x) for x in right],
  "left_rank_histogram":{str(k):v for k,v in sorted(lh.items())},"right_rank_histogram":{str(k):v for k,v in sorted(rh.items())},"left_right_intersection":[list(x) for x in intersection],
  "exact_checks":{"left_action_is_faithful":len(set(left))==15,"right_action_is_faithful":len(set(right))==15,
   "left_action_is_homomorphic":all(left[p[a][b]]==compose(left[a],left[b]) for a in range(15) for b in range(15)),
   "right_action_is_anti_homomorphic":all(right[p[a][b]]==compose(right[b],right[a]) for a in range(15) for b in range(15)),
   "rank_profiles_agree":lh==rh=={2:4,5:9,15:2},"only_identity_translation_is_shared":intersection==[(6,6)]},
  "claim_boundary":{"included":"complete left/right regular transformation actions, faithfulness, ranks, and composition laws","excluded":"endomorphism census, Schutzenberger actions, probabilistic dynamics, reliability, or thresholds","parent_issue":"remain open"}}
def validate_artifact(v:Mapping[str,Any]):
 e=build_artifact()
 if v!=e:raise ValueError("regular-action artifact does not reproduce")
 if set(e["exact_checks"].values())!={True}:raise ValueError("checks failed")
 return {"schema":SCHEMA,"status":"valid","left":15,"right":15}
def main(argv:Optional[Sequence[str]]=None)->int:
 q=argparse.ArgumentParser(description=__doc__);q.add_argument("--output",type=Path);q.add_argument("--validate",type=Path);a=q.parse_args(argv)
 if a.validate:print(json.dumps(validate_artifact(json.loads(a.validate.read_text())),indent=2,sort_keys=True));return 0
 s=json.dumps(build_artifact(),indent=2,sort_keys=True)+"\n"
 if a.output:a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(s,encoding="utf-8")
 else:print(s,end="")
 return 0
if __name__=="__main__":raise SystemExit(main())
