#!/usr/bin/env python3
"""Exact stabilizers, translation ranks, and cancellation in the serial monoid."""
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
ROOT=Path(__file__).resolve().parents[1];SCHEMA="matching-one/terminal-partition-serial-stabilizers/v1"
def table():
 states=enumerate_rgs(4);index={s:i for i,s in enumerate(states)}
 return tuple(tuple(index[serial_compose(a,b)] for b in states) for a in states)
def build_artifact()->dict[str,Any]:
 states=enumerate_rgs(4);index={s:i for i,s in enumerate(states)};product=table();rev=tuple(index[reverse_ports(s)] for s in states)
 left=[tuple(s for s in range(15) if product[s][a]==a) for a in range(15)]
 right=[tuple(s for s in range(15) if product[a][s]==a) for a in range(15)]
 left_ranks=[len(set(product[a])) for a in range(15)];right_ranks=[len({product[a][b] for a in range(15)}) for b in range(15)]
 left_cancel=[a for a,r in enumerate(left_ranks) if r==15];right_cancel=[a for a,r in enumerate(right_ranks) if r==15]
 return {"schema":SCHEMA,"issue":13,"status":"exact_stabilizer_and_cancellation_census",
  "left_stabilizers":[list(v) for v in left],"right_stabilizers":[list(v) for v in right],
  "left_stabilizer_size_histogram":{str(k):v for k,v in sorted(Counter(map(len,left)).items())},
  "right_stabilizer_size_histogram":{str(k):v for k,v in sorted(Counter(map(len,right)).items())},
  "left_translation_ranks":left_ranks,"right_translation_ranks":right_ranks,"left_cancellative_indices":left_cancel,"right_cancellative_indices":right_cancel,
  "exact_checks":{"left_and_right_cancellative_elements_are_exactly_units":left_cancel==right_cancel==[6,8],
   "translation_rank_profiles_agree":left_ranks==right_ranks==[5,5,5,2,2,5,15,5,15,5,5,5,5,2,2],
   "reversal_exchanges_left_and_right_stabilizers":all({rev[x] for x in left[a]}==set(right[rev[a]]) for a in range(15)),
   "identity_has_trivial_stabilizers":left[6]==right[6]==(6,),"nontrivial_unit_has_trivial_stabilizers":left[8]==right[8]==(6,)},
  "claim_boundary":{"included":"element stabilizers, regular-translation image ranks, cancellation, and reversal covariance","excluded":"shortlex dynamics, endomorphisms, congruences, physical cancellation, reliability, or thresholds","parent_issue":"remain open"}}
def validate_artifact(value:Mapping[str,Any]):
 expected=build_artifact()
 if value!=expected:raise ValueError("stabilizer artifact does not exactly reproduce")
 if set(expected["exact_checks"].values())!={True}:raise ValueError("all exact checks must pass")
 return {"schema":SCHEMA,"status":"valid","left_cancellative":2,"right_cancellative":2}
def main(argv:Optional[Sequence[str]]=None)->int:
 p=argparse.ArgumentParser(description=__doc__);p.add_argument("--output",type=Path);p.add_argument("--validate",type=Path);a=p.parse_args(argv)
 if a.validate:print(json.dumps(validate_artifact(json.loads(a.validate.read_text())),indent=2,sort_keys=True));return 0
 rendered=json.dumps(build_artifact(),indent=2,sort_keys=True)+"\n"
 if a.output:a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(rendered,encoding="utf-8")
 else:print(rendered,end="")
 return 0
if __name__=="__main__":raise SystemExit(main())
