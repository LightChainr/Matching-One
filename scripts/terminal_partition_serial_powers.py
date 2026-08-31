#!/usr/bin/env python3
"""Exact element power dynamics of the typed serial monoid."""
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
ROOT=Path(__file__).resolve().parents[1];SCHEMA="matching-one/terminal-partition-serial-powers/v1"
def table():
 states=enumerate_rgs(4);index={s:i for i,s in enumerate(states)}
 return tuple(tuple(index[serial_compose(a,b)] for b in states) for a in states)
def power_profile(element:int,product):
 if type(element) is not int or not 0<=element<len(product):raise ValueError("element index outside table")
 seen={};sequence=[];value=element
 while value not in seen:
  seen[value]=len(sequence);sequence.append(value);value=product[value][element]
 return {"element":element,"index":seen[value]+1,"period":len(sequence)-seen[value],"distinct_powers":sequence,"repeat_target":value}
def build_artifact()->dict[str,Any]:
 product=table();profiles=[power_profile(i,product) for i in range(15)];hist=Counter((x["index"],x["period"]) for x in profiles)
 return {"schema":SCHEMA,"issue":13,"status":"exact_element_power_dynamics","profiles":profiles,
  "index_period_histogram":{f"{k[0]},{k[1]}":v for k,v in sorted(hist.items())},"global_eventual_index":max(x["index"] for x in profiles),"global_eventual_period":2,
  "exact_checks":{"twelve_idempotent_fixed_points":hist[(1,1)]==12,"two_index_two_collapses":hist[(2,1)]==2,
   "unique_period_two_element_is_nontrivial_unit":[x["element"] for x in profiles if x["period"]==2]==[8],
   "elements_ten_and_eleven_square_to_fourteen":product[10][10]==product[11][11]==14,"period_two_cycle_is_8_6":profiles[8]["distinct_powers"]==[8,6]},
  "claim_boundary":{"included":"all positive power sequences, eventual indices, periods, and repeat targets","excluded":"shortlex representatives, random walks, endomorphisms, congruences, planar realization, or thresholds","parent_issue":"remain open"}}
def validate_artifact(value:Mapping[str,Any]):
 expected=build_artifact()
 if value!=expected:raise ValueError("power artifact does not exactly reproduce")
 if set(expected["exact_checks"].values())!={True}:raise ValueError("all exact checks must pass")
 return {"schema":SCHEMA,"status":"valid","profiles":15,"global_period":2}
def main(argv:Optional[Sequence[str]]=None)->int:
 p=argparse.ArgumentParser(description=__doc__);p.add_argument("--output",type=Path);p.add_argument("--validate",type=Path);a=p.parse_args(argv)
 if a.validate:print(json.dumps(validate_artifact(json.loads(a.validate.read_text())),indent=2,sort_keys=True));return 0
 rendered=json.dumps(build_artifact(),indent=2,sort_keys=True)+"\n"
 if a.output:a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(rendered,encoding="utf-8")
 else:print(rendered,end="")
 return 0
if __name__=="__main__":raise SystemExit(main())
