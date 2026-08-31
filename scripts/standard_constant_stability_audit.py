#!/usr/bin/env python3
"""Endpoint and 80/160-digit stability audit for closest pairwise witnesses."""

from __future__ import annotations
import argparse, hashlib, json
from bisect import bisect_left
from fractions import Fraction
from functools import lru_cache
from math import gcd
import mpmath as mp
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence

ROOT=Path(__file__).resolve().parents[1]
DEFAULT_CONTRACT=ROOT/"analysis"/"pslq_search_contract.json"
DEFAULT_OUTPUT=ROOT/"results"/"pslq-standard-constant-stability"/"latest.json"
SCHEMA="matching-one/standard-constant-stability-audit/v1"
INITIAL={
"pi":("3.14159265358979323846264338327950288419716939937510","3.14159265358979323846264338327950288419716939937511"),"e":("2.71828182845904523536028747135266249775724709369995","2.71828182845904523536028747135266249775724709369996"),"log2":("0.69314718055994530941723212145817656807550013436025","0.69314718055994530941723212145817656807550013436026"),"sqrt2":("1.41421356237309504880168872420969807856967187537694","1.41421356237309504880168872420969807856967187537695"),"sqrt3":("1.73205080756887729352744634150587236694280525381038","1.73205080756887729352744634150587236694280525381039"),"sqrt5":("2.23606797749978969640917366873127623544061835961152","2.23606797749978969640917366873127623544061835961153")}
CONFIRM={
"pi":("3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679","3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170680"),"e":("2.7182818284590452353602874713526624977572470936999595749669676277240766303535475945713821785251664274","2.7182818284590452353602874713526624977572470936999595749669676277240766303535475945713821785251664275"),"log2":("0.6931471805599453094172321214581765680755001343602552541206800094933936219696947156058633269964186875","0.6931471805599453094172321214581765680755001343602552541206800094933936219696947156058633269964186876"),"sqrt2":("1.4142135623730950488016887242096980785696718753769480731766797379907324784621070388503875343276415727","1.4142135623730950488016887242096980785696718753769480731766797379907324784621070388503875343276415728"),"sqrt3":("1.7320508075688772935274463415058723669428052538103806280558069794519330169088000370811461867572485756","1.7320508075688772935274463415058723669428052538103806280558069794519330169088000370811461867572485757"),"sqrt5":("2.2360679774997896964091736687312762354406183596115257242708972454105209256378048994144144083787822749","2.2360679774997896964091736687312762354406183596115257242708972454105209256378048994144144083787822750")}
IV_CONSTANTS={"pi":lambda:mp.iv.pi,"e":lambda:mp.iv.e,"log2":lambda:mp.iv.log(2),"sqrt2":lambda:mp.iv.sqrt(2),"sqrt3":lambda:mp.iv.sqrt(3),"sqrt5":lambda:mp.iv.sqrt(5)}
def verify_bounds(bounds:Mapping[str,tuple[str,str]],digits:int)->None:
 mp.iv.dps=digits
 for cid,(lo,hi) in bounds.items():
  if IV_CONSTANTS[cid]() not in mp.iv.mpf([lo,hi]):raise ValueError(f"{cid} enclosure is not outward at {digits} digits")
def _ceil(x:Fraction)->int:return -((-x.numerator)//x.denominator)
def _distance(lo:Fraction,hi:Fraction)->Fraction:return Fraction(0) if lo<=0<=hi else min(abs(lo),abs(hi))
def _text(x:Fraction)->str:return f"{x.numerator}/{x.denominator}"
def _scale(c:int,lo:Fraction,hi:Fraction):return (c*lo,c*hi) if c>=0 else (c*hi,c*lo)
def closest(p_lo:Fraction,p_hi:Fraction,c_lo:Fraction,c_hi:Fraction,height:int=100):
 best=None
 for b in range(1,height+1):
  for c in range(-height,height+1):
   if c==0:continue
   allowed=[a for a in range(-height,height+1) if gcd(gcd(abs(a),b),abs(c))==1]
   cl,ch=_scale(c,c_lo,c_hi);bl,bh=b*p_lo+cl,b*p_hi+ch;idx=bisect_left(allowed,_ceil(-bh))
   for j in (idx-1,idx,idx+1):
    if 0<=j<len(allowed):
     a=allowed[j];lo,hi=a+bl,a+bh;candidate=(_distance(lo,hi),(a,b,c))
     if best is None or candidate<best:best=candidate
 assert best is not None;return best
@lru_cache(maxsize=2)
def build_result(contract_path:Path=DEFAULT_CONTRACT)->dict[str,Any]:
 raw=contract_path.read_bytes();contract=json.loads(raw);verify_bounds(INITIAL,80);verify_bounds(CONFIRM,160);rows=[]
 for cid in INITIAL:
  ilo,ihi=map(Fraction,INITIAL[cid]);clo,chi=map(Fraction,CONFIRM[cid])
  for interval in contract["intervals"]:
   plo,phi=Fraction(interval["lower"]),Fraction(interval["upper"]);initial=closest(plo,phi,ilo,ihi);confirmation=closest(plo,phi,clo,chi)
   coefficients=initial[1];checks=[]
   for precision,bounds in ((80,(ilo,ihi)),(160,(clo,chi))):
    for point_id,p in (("lower",plo),("midpoint",(plo+phi)/2),("upper",phi)):
     a,b,c=coefficients;xlo,xhi=_scale(c,*bounds);rlo,rhi=a+b*p+xlo,a+b*p+xhi
     checks.append({"precision_digits":precision,"point":point_id,"residual_interval":[_text(rlo),_text(rhi)],"excludes_zero":not(rlo<=0<=rhi)})
   rows.append({"constant_id":cid,"interval_id":interval["id"],"closest_coefficients_for_1_p_constant":list(coefficients),"closest_witness_same_at_confirmation_precision":initial[1]==confirmation[1],"initial_minimum_residual":_text(initial[0]),"confirmation_minimum_residual":_text(confirmation[0]),"endpoint_checks":checks,"all_endpoint_checks_exclude_zero":all(x["excludes_zero"] for x in checks)})
 return {"schema":SCHEMA,"issue":1,"status":"standard_constant_stability_audit_complete","contract_sha256":hashlib.sha256(raw).hexdigest(),"initial_enclosure_digits":50,"confirmation_enclosure_digits":100,"declared_protocol_labels":[80,160],"constant_enclosures_verified_with_mpmath_iv":True,"rows":rows,
 "conclusion":{"all_closest_witnesses_stable":all(r["closest_witness_same_at_confirmation_precision"] for r in rows),"all_144_endpoint_checks_exclude_zero":all(r["all_endpoint_checks_exclude_zero"] for r in rows)},
 "claim_boundary":{"included":"lower/midpoint/upper and two-enclosure confirmation of the frozen closest pairwise witnesses","excluded":"formal proof about arbitrary precision, expanded bases, near-hit promotion, p-values, closed forms, or transcendence","parent_issue":"remain open"}}
def validate_result(result:Mapping[str,Any],contract_path:Path=DEFAULT_CONTRACT)->Mapping[str,Any]:
 expected=build_result(contract_path)
 if result!=expected:raise ValueError("stability audit does not exactly reproduce")
 return {"schema":SCHEMA,"status":"valid","rows":len(expected["rows"]),"endpoint_checks":sum(len(r["endpoint_checks"]) for r in expected["rows"])}
def main(argv:Optional[Sequence[str]]=None)->int:
 p=argparse.ArgumentParser(description=__doc__);p.add_argument("--contract",type=Path,default=DEFAULT_CONTRACT);p.add_argument("--output",type=Path);p.add_argument("--validate",type=Path);a=p.parse_args(argv)
 if a.validate:print(json.dumps(validate_result(json.loads(a.validate.read_text()),a.contract),indent=2,sort_keys=True));return 0
 rendered=json.dumps(build_result(a.contract),indent=2,sort_keys=True)+"\n"
 if a.output:a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(rendered)
 else:print(rendered,end="")
 return 0
if __name__=="__main__":raise SystemExit(main())
