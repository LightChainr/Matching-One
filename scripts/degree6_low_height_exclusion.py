#!/usr/bin/env python3
"""Certified exclusion at the historical complexity of exact planar thresholds.

Every exactly-known planar percolation threshold is algebraic of degree at
most 6 and coefficient height at most 3, and exactly one of them -- the
(3,12^2) site value ``sqrt(1-2 sin(pi/18))``, a root of ``x^6-3x^4+1`` -- lies
outside the degree-4 height-100 class already censused.  This module closes
that gap by exhausting ``C(d, 3)`` for ``d = 1..6`` against each frozen method
interval.

The class is small enough that every member is evaluated exactly at both
interval endpoints, in integer arithmetic on the common denominator, so no
statistical or floating-point step enters at any stage.  Two certified
consequences of ``|P'| <= D`` on ``[0,1]``, with ``D = 3*d*(d+1)/2`` bounding
``sum k|a_k|`` over the class, do the work:

* **screen** -- if ``|P(l)| > D*(u-l)`` then ``P`` cannot vanish on ``[l, u]``,
  since ``|P(x)| >= |P(l)| - D*(x-l)``.  Only polynomials failing this test
  reach an exact Sturm decision.
* **approach floor** -- for a root ``xi`` of ``P`` outside ``[l, u]``, either
  ``l - xi >= |P(l)|/D`` or ``xi - u >= |P(u)|/D``, so
  ``min(|P(l)|, |P(u)|)/D`` is a certified lower bound on the distance from the
  interval to the nearest root.  This extends the boundary-degree table of the
  P2 draft to degrees 5 and 6 without assuming monotonicity.
"""

from __future__ import annotations
import argparse, hashlib, json
from fractions import Fraction
from functools import lru_cache
from math import gcd
from pathlib import Path
from typing import Any, Iterator, Mapping, Optional, Sequence

try:
    from scripts.exact_polynomial_root_certificate import isolate_roots
    from scripts.pslq_look_elsewhere_ledger import primitive_polynomial_count
except ModuleNotFoundError:
    from exact_polynomial_root_certificate import isolate_roots
    from pslq_look_elsewhere_ledger import primitive_polynomial_count

ROOT=Path(__file__).resolve().parents[1]
DEFAULT_CONTRACT=ROOT/"analysis"/"pslq_search_contract.json"
SCHEMA="matching-one/degree6-low-height-exclusion/v1"
ISSUE=559
HEIGHT=3
DEGREE_MAX=6
ISOLATION_BITS=120


def _require(condition:bool,message:str)->None:
    if not condition:raise ValueError(message)


def _text(value:Fraction)->str:
    return f"{value.numerator}/{value.denominator}"


def _digest(path:Path)->str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def output_path(interval_id:str)->Path:
    return ROOT/"results"/f"pslq-degree6-low-height-{interval_id}"/"latest.json"


def derivative_bound(degree:int,height:int=HEIGHT)->int:
    """Bound on sum k|a_k|, hence on |P'| over [0,1], for all of C(degree, height)."""
    return height*degree*(degree+1)//2


@lru_cache(maxsize=DEGREE_MAX+1)
def _class_members(degree:int,height:int)->tuple[tuple[int,...],...]:
    return tuple(_generate_class(degree,height))


def enumerate_class(degree:int,height:int=HEIGHT)->Iterator[tuple[int,...]]:
    """Primitive, sign-normalized integer tuples of exactly this degree.

    The membership is fixed by (degree, height), so it is built once and
    replayed; every interval and every control trial scans the same class.
    """
    return iter(_class_members(degree,height))


def _generate_class(degree:int,height:int)->Iterator[tuple[int,...]]:
    span=list(range(-height,height+1))
    def rec(prefix:tuple[int,...])->Iterator[tuple[int,...]]:
        if len(prefix)==degree:
            for leading in range(1,height+1):
                candidate=prefix+(leading,)
                common=0
                for value in candidate:
                    common=gcd(common,value)
                    if common==1:break
                if common==1:yield candidate
            return
        for value in span:
            yield from rec(prefix+(value,))
    yield from rec(())


def _scan_degree(degree:int,lower:Fraction,upper:Fraction)->dict[str,Any]:
    denominator=lower.denominator*upper.denominator//gcd(lower.denominator,upper.denominator)
    low_numerator=lower.numerator*(denominator//lower.denominator)
    high_numerator=upper.numerator*(denominator//upper.denominator)
    scale=denominator**degree
    low_weights=[low_numerator**k*denominator**(degree-k) for k in range(degree+1)]
    high_weights=[high_numerator**k*denominator**(degree-k) for k in range(degree+1)]
    bound=derivative_bound(degree)
    # |P(l)| > D*(u-l)  =>  no root in [l, u];  scaled by denominator**degree.
    screen_limit=bound*(high_numerator-low_numerator)*denominator**(degree-1)

    counted=0
    screened=0
    root_containing=0
    distinct_roots=0
    witnesses:list[dict[str,Any]]=[]
    best:tuple[int,tuple[int,...],int,int]|None=None

    for coefficients in enumerate_class(degree):
        counted+=1
        at_low=0
        at_high=0
        for index,coefficient in enumerate(coefficients):
            if coefficient:
                at_low+=low_weights[index]*coefficient
                at_high+=high_weights[index]*coefficient
        nearest=min(abs(at_low),abs(at_high))
        if abs(at_low)<=screen_limit:
            screened+=1
            polynomial=[Fraction(value) for value in coefficients]
            roots=isolate_roots(polynomial,lower,upper,bits=ISOLATION_BITS)
            if roots:
                root_containing+=1
                distinct_roots+=len(roots)
                nearest=0
                witnesses.append({"coefficients_ascending":list(coefficients),
                                  "root_brackets":[[_text(lo),_text(hi)] for lo,hi in roots],
                                  "isolation_bits":ISOLATION_BITS})
        candidate=(nearest,coefficients,at_low,at_high)
        if best is None or candidate[:2]<best[:2]:
            best=candidate

    _require(counted==primitive_polynomial_count(degree,HEIGHT),
             f"degree {degree} enumeration size disagrees with the committed counter")
    assert best is not None
    minimum=Fraction(best[0],scale)
    floor=minimum/bound
    return {
        "degree":degree,
        "polynomials_in_class":counted,
        "screen_candidates_exactly_decided":screened,
        "root_containing_polynomials":root_containing,
        "distinct_roots_in_interval":distinct_roots,
        "excluded":root_containing==0,
        "root_witnesses":witnesses,
        "derivative_bound_on_unit_interval":bound,
        "closest_polynomial":{
            "coefficients_ascending":list(best[1]),
            "height":max(abs(value) for value in best[1]),
            "minimum_absolute_endpoint_residual":_text(minimum),
            "polynomial_endpoint_values":[_text(Fraction(best[2],scale)),_text(Fraction(best[3],scale))],
        },
        "root_distance_lower_bound_text":_text(floor),
        "floor_to_interval_width_ratio_text":_text(floor/(upper-lower)),
    }


def run_search(interval:Mapping[str,Any])->dict[str,Any]:
    lower,upper=Fraction(interval["lower"]),Fraction(interval["upper"])
    _require(0<lower<upper<1,"interval must be a nonempty subinterval of (0,1)")
    degrees=[_scan_degree(degree,lower,upper) for degree in range(1,DEGREE_MAX+1)]
    return {
        "interval_id":interval["id"],"source_id":interval["source_id"],
        "lower":interval["lower"],"upper":interval["upper"],
        "width_text":_text(upper-lower),
        "polynomials_per_interval":sum(row["polynomials_in_class"] for row in degrees),
        "by_degree":degrees,
        "excluded":all(row["excluded"] for row in degrees),
        "degrees_excluded":[row["degree"] for row in degrees if row["excluded"]],
        "degrees_with_roots":[row["degree"] for row in degrees if not row["excluded"]],
    }


def build_result(interval_id:str,contract_path:Path=DEFAULT_CONTRACT)->dict[str,Any]:
    return _build_result(interval_id,Path(contract_path).resolve())


@lru_cache(maxsize=8)
def _build_result(interval_id:str,contract_path:Path)->dict[str,Any]:
    raw=contract_path.read_bytes()
    contract=json.loads(raw)
    provenance=contract["provenance"]
    provenance_digest=_digest(ROOT/provenance["path"])
    _require(provenance_digest==provenance["sha256"],"provenance digest drift")
    rows=[row for row in contract["intervals"] if row["id"]==interval_id]
    _require(len(rows)==1,"interval id is not uniquely frozen")
    return {
        "schema":SCHEMA,"issue":ISSUE,"status":"degree6_low_height_exclusion_complete",
        "contract_sha256":hashlib.sha256(raw).hexdigest(),"provenance_sha256":provenance_digest,
        "search":{
            "degree_min":1,"degree_max":DEGREE_MAX,"coefficient_height_max":HEIGHT,
            "primitive_coefficients_only":True,"nonzero_leading_coefficient":True,
            "sign_normalization":"leading_positive",
            "decision":"exact integer endpoint evaluation, certified derivative screen, exact Sturm isolation",
            "motivation":"every exactly-known planar percolation threshold has degree <= 6 and height <= 3; "
                         "the (3,12^2) site value x^6-3x^4+1 is the one such form outside C(<=4, <=100)",
        },
        "interval_result":run_search(rows[0]),
        "claim_boundary":{
            "included":f"exhaustive degree-1..{DEGREE_MAX} height-{HEIGHT} exclusion on {interval_id} only",
            "excluded":"other method intervals, higher degree or height, near-hit promotion, p-values, "
                       "closed forms, or transcendence",
            "parent_issue":"remain open",
        },
    }


def validate_result(result:Mapping[str,Any],contract_path:Path=DEFAULT_CONTRACT)->Mapping[str,Any]:
    interval_id=result["interval_result"]["interval_id"]
    expected=build_result(interval_id,contract_path)
    if result!=expected:raise ValueError("degree-6 low-height exclusion does not exactly reproduce")
    return {"schema":SCHEMA,"status":"valid","interval_id":interval_id,
            "excluded":expected["interval_result"]["excluded"]}


def main(argv:Optional[Sequence[str]]=None)->int:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("interval",nargs="?",help="frozen interval id; omit with --all")
    parser.add_argument("--all",action="store_true",help="run every frozen interval")
    parser.add_argument("--contract",type=Path,default=DEFAULT_CONTRACT)
    parser.add_argument("--output",type=Path)
    parser.add_argument("--validate",type=Path)
    args=parser.parse_args(argv)
    if args.validate:
        print(json.dumps(validate_result(json.loads(args.validate.read_text(encoding="utf-8")),args.contract),indent=2,sort_keys=True))
        return 0
    contract=json.loads(Path(args.contract).read_bytes())
    targets=[row["id"] for row in contract["intervals"]] if args.all else [args.interval]
    _require(all(targets),"an interval id or --all is required")
    for interval_id in targets:
        rendered=json.dumps(build_result(interval_id,args.contract),indent=2,sort_keys=True)+"\n"
        destination=args.output if (args.output and not args.all) else output_path(interval_id)
        destination.parent.mkdir(parents=True,exist_ok=True)
        destination.write_text(rendered,encoding="utf-8")
        print(f"wrote {destination.relative_to(ROOT)}")
    return 0


if __name__=="__main__":raise SystemExit(main())
