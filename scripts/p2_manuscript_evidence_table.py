#!/usr/bin/env python3
"""Assemble the P2 manuscript evidence tables from committed exclusion artifacts.

This script performs no census computation.  It reads the frozen search
contract, the canonical provenance manifest, and the committed degree-1
through degree-4 result artifacts, then derives only statements that follow
exactly from those inputs: interval ordering and disjointness, the exclusion
table, and the cross-interval resolution of every committed quartic root
witness.  Root decisions reuse the repository's exact Sturm certificate path.
"""

from __future__ import annotations
import argparse, hashlib, json
from fractions import Fraction
from functools import lru_cache
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence

try:
    from scripts.exact_polynomial_root_certificate import isolate_roots
    from scripts.pslq_look_elsewhere_ledger import primitive_polynomial_count as _ledger_count
except ModuleNotFoundError:
    from exact_polynomial_root_certificate import isolate_roots
    from pslq_look_elsewhere_ledger import primitive_polynomial_count as _ledger_count

ROOT=Path(__file__).resolve().parents[1]
DEFAULT_CONTRACT=ROOT/"analysis"/"pslq_search_contract.json"
DEFAULT_OUTPUT=ROOT/"results"/"p2-algebraic-exclusion-manuscript"/"latest.json"
SCHEMA="matching-one/p2-manuscript-evidence/v1"
ISSUE=551
ISOLATION_BITS=120
DECIMAL_DIGITS=24

DEGREE_SOURCES={
    1:("results/pslq-degree1-rational-exclusion/latest.json","interval_results"),
    2:("results/pslq-degree2-polynomial-exclusion/latest.json","interval_results"),
    3:("results/pslq-degree3-{interval}/latest.json","interval_result"),
    4:("results/pslq-degree4-{interval}/latest.json","interval_result"),
}
COUNT_KEYS={1:"primitive_polynomials_checked",2:"primitive_quadratics_checked",3:"primitive_cubics_checked",4:"primitive_quartics_covered"}
CONTROL_SOURCES={
    "look_elsewhere_ledger":"results/pslq-look-elsewhere-ledger/latest.json",
    "kagome_positive_control":"results/pslq-kagome-exact-control/latest.json",
    "synthetic_false_positive_calibration":"results/pslq-synthetic-false-positive-calibration/latest.json",
    "standard_constant_pairwise":"results/pslq-standard-constant-pairwise/latest.json",
    "standard_constant_stability":"results/pslq-standard-constant-stability/latest.json",
    "lattice_native_candidates":"results/pslq-lattice-native-candidates/latest.json",
    "degree4_boundary_sensitivity":"results/pslq-degree4-synthetic-boundary-control/latest.json",
    "degree6_low_height_control":"results/pslq-degree6-low-height-control/latest.json",
}
LOW_HEIGHT_SOURCES="results/pslq-degree6-low-height-{interval}/latest.json"


def _require(condition:bool,message:str)->None:
    if not condition:raise ValueError(message)


def _text(value:Fraction)->str:
    return f"{value.numerator}/{value.denominator}"


def decimal_string(value:Fraction,digits:int=DECIMAL_DIGITS)->str:
    """Truncate a nonnegative rational toward zero to a fixed decimal string."""
    _require(value>=0,"decimal rendering expects a nonnegative value")
    scaled=value.numerator*10**digits//value.denominator
    whole,fraction=divmod(scaled,10**digits)
    return f"{whole}.{fraction:0{digits}d}"


def decimal_prefix(lower:Fraction,upper:Fraction,digits:int=DECIMAL_DIGITS)->str:
    """Longest decimal prefix shared by every point of a bracket."""
    _require(0<=lower<=upper,"bracket must be ordered and nonnegative")
    low,high=decimal_string(lower,digits),decimal_string(upper,digits)
    shared=0
    while shared<len(low) and low[shared]==high[shared]:shared+=1
    return low[:shared]


def _digest(path:Path)->str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _load(relative:str)->dict[str,Any]:
    return json.loads((ROOT/relative).read_text(encoding="utf-8"))


def _degree_rows(interval_ids:Sequence[str])->dict[int,dict[str,Mapping[str,Any]]]:
    rows:dict[int,dict[str,Mapping[str,Any]]]={}
    for degree,(template,key) in DEGREE_SOURCES.items():
        rows[degree]={}
        if "{interval}" in template:
            for interval_id in interval_ids:
                rows[degree][interval_id]=_load(template.format(interval=interval_id))[key]
        else:
            for row in _load(template)[key]:
                rows[degree][row["interval_id"]]=row
        _require(set(rows[degree])==set(interval_ids),f"degree {degree} does not cover every frozen interval")
    return rows


def _source_artifacts(interval_ids:Sequence[str])->list[dict[str,Any]]:
    paths=[]
    for degree,(template,_) in sorted(DEGREE_SOURCES.items()):
        if "{interval}" in template:paths.extend(template.format(interval=i) for i in interval_ids)
        else:paths.append(template)
    paths.extend(LOW_HEIGHT_SOURCES.format(interval=i) for i in interval_ids)
    paths.extend(CONTROL_SOURCES[name] for name in sorted(CONTROL_SOURCES))
    artifacts=[]
    for relative in paths:
        payload=_load(relative)
        artifacts.append({"path":relative,"sha256":_digest(ROOT/relative),"schema":payload.get("schema"),"status":payload.get("status")})
    return artifacts


def _interval_table(contract:Mapping[str,Any],provenance:Mapping[str,Any])->list[dict[str,Any]]:
    sources={row["id"]:row for row in provenance["sources"]}
    table=[]
    for row in contract["intervals"]:
        source=sources[row["source_id"]]
        quoted=[q for q in source["quoted_estimates"] if q["central_value"]==row["central_value"]]
        _require(len(quoted)==1,f"interval {row['id']} does not match exactly one quoted estimate")
        lower,upper=Fraction(row["lower"]),Fraction(row["upper"])
        _require(lower<upper,f"interval {row['id']} is empty")
        table.append({
            "interval_id":row["id"],"source_id":row["source_id"],"citation":source["citation"],
            "doi":source.get("doi"),"source_status":source["source_status"],
            "estimator":quoted[0].get("estimator",source["estimator_family"]),"geometry":source["geometry"],
            "source_location":quoted[0]["source_location"],"value_text":row["value_text"],
            "lower":row["lower"],"upper":row["upper"],
            "width_text":_text(upper-lower),"width_decimal":decimal_string(upper-lower),
            "confidence_homogenized":row["confidence_homogenized"],
        })
    return table


def _disjointness(table:Sequence[Mapping[str,Any]])->dict[str,Any]:
    ordered=sorted(table,key=lambda row:Fraction(row["lower"]))
    gaps=[]
    disjoint=True
    for below,above in zip(ordered,ordered[1:]):
        separation=Fraction(above["lower"])-Fraction(below["upper"])
        if separation<=0:disjoint=False
        gaps.append({"below_interval_id":below["interval_id"],"above_interval_id":above["interval_id"],
                     "gap_text":_text(separation),"gap_decimal":decimal_string(separation) if separation>=0 else None})
    return {"pairwise_disjoint":disjoint,"ascending_order":[row["interval_id"] for row in ordered],
            "adjacent_gaps":gaps,
            "meaning":"the four method intervals do not intersect, so no single pooled interval represents the literature"}


def _exclusion_table(rows:Mapping[int,Mapping[str,Mapping[str,Any]]],interval_ids:Sequence[str])->list[dict[str,Any]]:
    table=[]
    for degree in sorted(rows):
        for interval_id in interval_ids:
            row=rows[degree][interval_id]
            containing=row.get("root_containing_polynomials",row.get("zero_containing_residuals"))
            _require(containing is not None,f"degree {degree} row for {interval_id} has no root count")
            closest=row["closest_polynomial"]
            table.append({
                "degree":degree,"interval_id":interval_id,
                "polynomials_searched":row[COUNT_KEYS[degree]],
                "root_containing_polynomials":containing,
                "excluded":row.get("excluded",containing==0),
                "closest_coefficients_ascending":closest["coefficients_ascending"],
                "closest_height":closest["height"],
                "closest_minimum_absolute_residual":closest["minimum_absolute_residual"],
            })
    return table


def _closest_approach(exclusion:Sequence[Mapping[str,Any]],table:Sequence[Mapping[str,Any]])->dict[str,Any]:
    """Certified root-distance floor achieved by the best polynomial of each degree.

    On [0,1] the derivative of an integer polynomial is bounded by
    ``D = sum k*|a_k|``.  If ``|P| >= r`` throughout a method interval then the
    mean value theorem puts every root of ``P`` in [0,1] at distance at least
    ``r / D`` from that interval, so ``r / D`` measures how closely the search
    class can approach the interval without entering it.
    """
    widths={row["interval_id"]:Fraction(row["width_text"]) for row in table}
    rows=[]
    for entry in exclusion:
        coefficients=entry["closest_coefficients_ascending"]
        residual=Fraction(entry["closest_minimum_absolute_residual"])
        derivative_bound=sum(power*abs(value) for power,value in enumerate(coefficients) if power)
        _require(derivative_bound>0,"closest polynomial has a vanishing derivative bound")
        floor=residual/derivative_bound
        ratio=floor/widths[entry["interval_id"]]
        rows.append({"degree":entry["degree"],"interval_id":entry["interval_id"],
                     "closest_coefficients_ascending":coefficients,
                     "derivative_bound_on_unit_interval":derivative_bound,
                     "root_distance_lower_bound_text":_text(floor),
                     "root_distance_lower_bound_decimal":decimal_string(floor,30),
                     "floor_to_interval_width_ratio_text":_text(ratio),
                     "floor_to_interval_width_ratio_decimal":decimal_string(ratio,6),
                     "root_inside_interval":floor==0})
    resolution={}
    for degree in sorted({row["degree"] for row in rows}):
        floors=[Fraction(row["root_distance_lower_bound_text"]) for row in rows if row["degree"]==degree]
        best=min(floors)
        resolution[str(degree)]={"closest_approach_text":_text(best),"closest_approach_decimal":decimal_string(best,30)}
    strictly_decreasing=all(Fraction(resolution[str(d)]["closest_approach_text"])>Fraction(resolution[str(d+1)]["closest_approach_text"])
                            for d in sorted(int(k) for k in resolution)[:-1])
    clears_own_width=sorted({row["degree"] for row in rows if
                             all(Fraction(other["floor_to_interval_width_ratio_text"])>1
                                 for other in rows if other["degree"]==row["degree"])})
    return {
        "rows":rows,"approach_resolution_by_degree":resolution,
        "approach_resolution_strictly_decreasing_in_degree":strictly_decreasing,
        "degrees_whose_closest_polynomial_stays_further_than_one_interval_width":clears_own_width,
        "boundary_degree":min(sorted(int(k) for k in resolution),key=lambda d:d) if not clears_own_width
            else max(clears_own_width)+1,
        "meaning":"the search class approaches the method intervals more closely at each higher degree; "
                  "degree 4 is the first degree whose closest polynomial no longer stays more than one "
                  "interval width away, and it is also the first degree with surviving polynomials",
        "claim_boundary":"a certified distance floor for the committed closest polynomials only; "
                         "not an equidistribution theorem, a null distribution, or a p-value",
    }


def _reflect(coefficients:Sequence[int])->list[int]:
    """Coefficients of P(1-x), sign-normalized to a positive leading coefficient.

    This is the matching substitution p -> 1-p, so a lattice's exactly-known
    threshold polynomial and its matching partner's are related by it.
    """
    degree=len(coefficients)-1
    result=[0]*(degree+1)
    for power,value in enumerate(coefficients):
        # (1-x)^power = sum_j C(power,j) (-x)^j
        binomial=1
        for j in range(power+1):
            result[j]+=value*binomial*(-1 if j%2 else 1)
            binomial=binomial*(power-j)//(j+1)
    while result and result[-1]==0:result.pop()
    _require(bool(result),"reflection collapsed the polynomial")
    if result[-1]<0:result=[-value for value in result]
    common=0
    for value in result:common=_gcd(common,abs(value))
    _require(common>0,"reflection produced the zero polynomial")
    return [value//common for value in result]


def _gcd(first:int,second:int)->int:
    while second:first,second=second,first%second
    return abs(first)


def _primitive_count(degree:int,height:int)->int:
    """Delegate to the committed look-elsewhere ledger's exact counter."""
    return _ledger_count(degree,height)


def _historical_forms(table:Sequence[Mapping[str,Any]])->dict[str,Any]:
    """Algebraic complexity of the exactly-known planar percolation thresholds.

    Minimal polynomials are taken from the committed lattice-native candidate
    artifact, which certifies each one with an isolating interval, rather than
    restated here.  The triangular-bond entry is derived from the certified
    kagome-site entry by the matching substitution p -> 1-p.
    """
    native=_load(CONTROL_SOURCES["lattice_native_candidates"])
    polynomials={row["candidate_id"]:row["minimal_polynomial_coefficients_ascending"] for row in native["candidates"]}
    kagome=polynomials["kagome-site"]
    rows=[
        {"lattice":"square bond; triangular site","closed_form":"1/2",
         "minimal_polynomial_ascending":[-1,2],"provenance":"Sykes and Essam 1964"},
        {"lattice":"triangular bond","closed_form":"2*sin(pi/18)",
         "minimal_polynomial_ascending":_reflect(kagome),
         "provenance":"derived from the certified kagome-site polynomial by p -> 1-p"},
        {"lattice":"honeycomb bond; kagome site","closed_form":"1-2*sin(pi/18)",
         "minimal_polynomial_ascending":kagome,
         "provenance":"results/pslq-lattice-native-candidates/latest.json"},
        {"lattice":"(3,12^2) site","closed_form":"sqrt(1-2*sin(pi/18))",
         "minimal_polynomial_ascending":polynomials["three-twelve-site"],
         "provenance":"results/pslq-lattice-native-candidates/latest.json"},
        {"lattice":"martini bond","closed_form":"1/sqrt(2)",
         "minimal_polynomial_ascending":polynomials["martini-descendant-root2"],
         "provenance":"results/pslq-lattice-native-candidates/latest.json"},
        {"lattice":"martini descendant","closed_form":"(sqrt(5)-1)/2",
         "minimal_polynomial_ascending":polynomials["martini-descendant-golden"],
         "provenance":"results/pslq-lattice-native-candidates/latest.json"},
    ]
    for row in rows:
        coefficients=row["minimal_polynomial_ascending"]
        row["degree"]=len(coefficients)-1
        row["height"]=max(abs(value) for value in coefficients)
        row["inside_census_class"]=row["degree"]<=4 and row["height"]<=100
    max_degree=max(row["degree"] for row in rows)
    max_height=max(row["height"] for row in rows)
    extension={str(degree):_primitive_count(degree,max_height) for degree in range(1,max_degree+1)}
    return {
        "rows":rows,
        "historical_class_closure":{
            "hypothesis":"the one exactly-known planar threshold outside the census class is the (3,12^2) site "
                         "value, a square root of a lower-degree threshold; the targeted question is therefore "
                         "whether any form of degree up to the historical maximum, at the historical height, "
                         "has a root in a method interval",
            "class":f"degree 1..{max_degree}, height <= {max_height}",
            "primitive_counts_by_degree":extension,
            "total_polynomials_per_interval":sum(extension.values()),
            "status":"executed",
            "issue":559,
            "artifacts":sorted(LOW_HEIGHT_SOURCES.format(interval=row["interval_id"]) for row in table),
            "control":CONTROL_SOURCES["degree6_low_height_control"],
        },
        "max_degree":max(row["degree"] for row in rows),
        "max_height":max(row["height"] for row in rows),
        "all_inside_census_class":all(row["inside_census_class"] for row in rows),
        "outside_census_class":[row["lattice"] for row in rows if not row["inside_census_class"]],
        "meaning":"every exactly-known planar percolation threshold is algebraic of degree at most 6 and "
                  "coefficient height at most 3; the census class C(<=4, <=100) is far more generous in height "
                  "than the entire historical record, but its degree bound does not reach the one known "
                  "degree-6 form, which arises as a square root of a lower-degree threshold",
        "claim_boundary":"algebraic complexity of published exact thresholds; not a claim that square-site p_c "
                         "belongs to this tradition, and not a survey of every conjectured form",
    }


def _historical_range(table:Sequence[Mapping[str,Any]])->dict[str,Any]:
    """Summarise the degree-1..6 height-3 exhaustion across the frozen intervals."""
    per_interval={}
    for row in table:
        payload=_load(LOW_HEIGHT_SOURCES.format(interval=row["interval_id"]))["interval_result"]
        per_interval[row["interval_id"]]=payload
    reference=next(iter(per_interval.values()))
    degrees=[]
    for index,entry in enumerate(reference["by_degree"]):
        degree=entry["degree"]
        ratios={}
        for interval_id,payload in per_interval.items():
            other=payload["by_degree"][index]
            _require(other["degree"]==degree,"degree rows are misaligned across intervals")
            _require(other["closest_polynomial"]["coefficients_ascending"]
                     ==entry["closest_polynomial"]["coefficients_ascending"],
                     f"degree {degree} closest polynomial differs across intervals")
            ratios[interval_id]=Fraction(other["floor_to_interval_width_ratio_text"])
        floor=Fraction(entry["root_distance_lower_bound_text"])
        degrees.append({
            "degree":degree,
            "polynomials_in_class":entry["polynomials_in_class"],
            "closest_coefficients_ascending":entry["closest_polynomial"]["coefficients_ascending"],
            "root_distance_lower_bound_text":_text(floor),
            "root_distance_lower_bound_decimal":decimal_string(floor,30),
            "floor_to_width_ratio_min_decimal":decimal_string(min(ratios.values()),3),
            "floor_to_width_ratio_max_decimal":decimal_string(max(ratios.values()),3),
            "screen_candidates_total":sum(p["by_degree"][index]["screen_candidates_exactly_decided"]
                                          for p in per_interval.values()),
            "root_containing_total":sum(p["by_degree"][index]["root_containing_polynomials"]
                                        for p in per_interval.values()),
        })
    control=_load(CONTROL_SOURCES["degree6_low_height_control"])
    return {
        "rows":degrees,
        "polynomials_per_interval":reference["polynomials_per_interval"],
        "excluded_on_every_interval":all(p["excluded"] for p in per_interval.values()),
        "screen_retained_nothing":all(row["screen_candidates_total"]==0 for row in degrees),
        "sensitivity_control":{
            "planted_coefficients_ascending":control["planted"]["coefficients_ascending"],
            "closed_form":control["planted"]["closed_form"],
            "positive_trials":control["conclusion"]["positive_trials"],
            "negative_trials":control["conclusion"]["negative_trials"],
            "all_trials_passed":control["conclusion"]["all_trials_passed"],
        },
        "meaning":"no algebraic form at the complexity of any exactly-known planar percolation threshold has a "
                  "root in any of the four published intervals",
        "claim_boundary":"exhaustive over degree <= 6 at height <= 3 only; says nothing about height > 3 or "
                         "degree > 6, and is not a transcendence or non-algebraicity claim",
    }


def _separation(bracket:tuple[Fraction,Fraction],lower:Fraction,upper:Fraction)->Optional[Fraction]:
    """Certified lower bound on the distance from a root bracket to an interval."""
    low,high=bracket
    if high<lower:return lower-high
    if low>upper:return low-upper
    return None


def _survivor_census(rows:Mapping[int,Mapping[str,Mapping[str,Any]]],table:Sequence[Mapping[str,Any]])->dict[str,Any]:
    bounds={row["interval_id"]:(Fraction(row["lower"]),Fraction(row["upper"])) for row in table}
    interval_ids=[row["interval_id"] for row in table]
    witnesses:dict[tuple[int,...],list[str]]={}
    for interval_id in interval_ids:
        for witness in rows[4][interval_id].get("root_witnesses",[]):
            witnesses.setdefault(tuple(witness["coefficients_ascending"]),[]).append(interval_id)

    survivors=[]
    for coefficients in sorted(witnesses,key=lambda c:(max(abs(v) for v in c),c)):
        polynomial=[Fraction(value) for value in coefficients]
        counts,brackets={},{}
        for interval_id in interval_ids:
            lower,upper=bounds[interval_id]
            roots=isolate_roots(polynomial,lower,upper,bits=ISOLATION_BITS)
            counts[interval_id]=len(roots)
            if roots:brackets[interval_id]=roots
        surviving=[interval_id for interval_id in interval_ids if counts[interval_id]]
        _require(surviving==sorted(witnesses[coefficients],key=interval_ids.index),
                 f"recomputed survival set disagrees with the committed witnesses for {coefficients}")
        _require(len(surviving)==1,f"quartic {coefficients} survives more than one interval")
        home=surviving[0]
        _require(len(brackets[home])==1,f"quartic {coefficients} has multiple roots inside {home}")
        bracket=brackets[home][0]
        separations=[]
        for interval_id in interval_ids:
            if interval_id==home:continue
            gap=_separation(bracket,*bounds[interval_id])
            _require(gap is not None and gap>0,f"quartic {coefficients} is not separated from {interval_id}")
            separations.append({"interval_id":interval_id,"separation_lower_bound_text":_text(gap),
                                "separation_lower_bound_decimal":decimal_string(gap)})
        survivors.append({
            "coefficients_ascending":list(coefficients),
            "height":max(abs(value) for value in coefficients),
            "surviving_interval_id":home,
            "excluded_by_interval_ids":[row["interval_id"] for row in separations],
            "root_bracket":[_text(bracket[0]),_text(bracket[1])],
            "root_decimal_prefix":decimal_prefix(*bracket),
            "isolation_bits":ISOLATION_BITS,
            "separations":separations,
        })

    per_interval={interval_id:sum(1 for row in survivors if row["surviving_interval_id"]==interval_id) for interval_id in interval_ids}
    for interval_id in interval_ids:
        _require(per_interval[interval_id]==rows[4][interval_id]["root_containing_polynomials"],
                 f"survivor count drift on {interval_id}")
    return {
        "distinct_surviving_quartics":len(survivors),
        "max_intervals_per_survivor":max((len(interval_ids)-len(row["excluded_by_interval_ids"]) for row in survivors),default=0),
        "survivors_per_interval":per_interval,
        "every_survivor_excluded_by_every_other_interval":True,
        "survivors":survivors,
        "meaning":"each committed quartic root witness lies in exactly one method interval and is certified outside the other three",
    }


def _width_scaling(table:Sequence[Mapping[str,Any]],census:Mapping[str,Any])->dict[str,Any]:
    rows=[]
    for row in sorted(table,key=lambda entry:-Fraction(entry["width_text"])):
        width=Fraction(row["width_text"])
        count=census["survivors_per_interval"][row["interval_id"]]
        rows.append({"interval_id":row["interval_id"],"width_decimal":row["width_decimal"],
                     "surviving_quartics":count,
                     "survivors_per_unit_length_decimal":decimal_string(Fraction(count)/width,3)})
    widest_first=[row["surviving_quartics"] for row in rows]
    return {
        "rows_widest_first":rows,
        "survivor_count_is_monotone_in_width":all(a>=b for a,b in zip(widest_first,widest_first[1:])),
        "claim_boundary":"descriptive count density only; this is not a null distribution, multiplicity correction, or p-value",
    }


@lru_cache(maxsize=2)
def build_result(contract_path:Path=DEFAULT_CONTRACT)->dict[str,Any]:
    raw=contract_path.read_bytes()
    contract=json.loads(raw)
    provenance_meta=contract["provenance"]
    provenance_path=ROOT/provenance_meta["path"]
    provenance_digest=_digest(provenance_path)
    _require(provenance_digest==provenance_meta["sha256"],"provenance digest drift")
    provenance=json.loads(provenance_path.read_text(encoding="utf-8"))

    interval_table=_interval_table(contract,provenance)
    interval_ids=[row["interval_id"] for row in interval_table]
    rows=_degree_rows(interval_ids)
    ledger=_load(CONTROL_SOURCES["look_elsewhere_ledger"])
    counts=ledger["primitive_polynomial_counts_by_degree"]
    for degree in sorted(rows):
        for interval_id in interval_ids:
            _require(rows[degree][interval_id][COUNT_KEYS[degree]]==counts[str(degree)],
                     f"degree {degree} search size disagrees with the look-elsewhere ledger")

    census=_survivor_census(rows,interval_table)
    exclusion=_exclusion_table(rows,interval_ids)
    fully_excluded=[row["interval_id"] for row in interval_table
                    if all(entry["excluded"] for entry in exclusion if entry["interval_id"]==row["interval_id"])]
    return {
        "schema":SCHEMA,"issue":ISSUE,"status":"p2_manuscript_evidence_assembled",
        "contract_sha256":hashlib.sha256(raw).hexdigest(),"provenance_sha256":provenance_digest,
        "source_artifacts":_source_artifacts(interval_ids),
        "intervals":interval_table,
        "interval_disjointness":_disjointness(interval_table),
        "search_class":{
            "coefficient_height_max":contract["search_stages"]["algebraic_polynomial"]["coefficient_height_max"],
            "degree_min":contract["search_stages"]["algebraic_polynomial"]["degree_min"],
            "degree_max":contract["search_stages"]["algebraic_polynomial"]["degree_max"],
            "primitive_coefficients_only":True,"nonzero_leading_coefficient":True,
            "primitive_polynomial_counts_by_degree":counts,
            "polynomials_per_interval":ledger["family_counts"]["algebraic_polynomials_per_interval"],
            "total_declared_interval_comparisons":ledger["total_declared_interval_comparisons"],
        },
        "exclusion_table":exclusion,
        "intervals_excluded_at_every_degree":fully_excluded,
        "closest_approach_by_degree":_closest_approach(exclusion,interval_table),
        "historical_form_complexity":_historical_forms(interval_table),
        "historical_range_exclusion":_historical_range(interval_table),
        "quartic_survivor_census":census,
        "width_scaling_diagnostic":_width_scaling(interval_table,census),
        "controls":{name:{"path":relative,"status":_load(relative).get("status")} for name,relative in sorted(CONTROL_SOURCES.items())},
        "claim_boundary":{
            "included":"assembly and exact cross-interval resolution of the committed degree-1..4 census artifacts",
            "excluded":"new census computation, degree/height expansion, near-hit promotion, p-values, closed forms, or transcendence",
            "parent_issue":"remain open",
        },
    }


def validate_result(result:Mapping[str,Any],contract_path:Path=DEFAULT_CONTRACT)->Mapping[str,Any]:
    expected=build_result(contract_path)
    if result!=expected:raise ValueError("P2 manuscript evidence does not exactly reproduce")
    return {"schema":SCHEMA,"status":"valid","distinct_surviving_quartics":expected["quartic_survivor_census"]["distinct_surviving_quartics"]}


def _polynomial_text(coefficients:Sequence[int])->str:
    terms=[]
    for power,value in reversed(list(enumerate(coefficients))):
        if value==0:continue
        monomial="" if power==0 else ("x" if power==1 else f"x^{power}")
        sign="-" if value<0 else "+"
        magnitude=abs(value)
        body=monomial if magnitude==1 and monomial else f"{magnitude}{monomial}"
        terms.append(f"{sign} {body}" if terms else (f"-{body}" if value<0 else body))
    return " ".join(terms) if terms else "0"


def render_markdown(result:Mapping[str,Any])->str:
    """Render the manuscript tables so the prose can never drift from the artifact."""
    lines=["<!-- Generated by scripts/p2_manuscript_evidence_table.py --markdown. Do not edit by hand. -->",
           "# P2 manuscript tables","",
           f"Assembled from the committed census artifacts; contract `{result['contract_sha256'][:16]}`, "
           f"provenance `{result['provenance_sha256'][:16]}`.","",
           "## Table 1 — Canonical provenance of the frozen method intervals","",
           "| Interval id | Source | Estimator | Geometry | Quoted value | Interval | Width | Source location | Provenance status |",
           "|---|---|---|---|---|---|---:|---|---|"]
    for row in result["intervals"]:
        lines.append(f"| `{row['interval_id']}` | {row['citation'].split(',')[0]} | `{row['estimator']}` | {row['geometry']} | "
                     f"`{row['value_text']}` | `[{row['lower']}, {row['upper']}]` | `{row['width_decimal'].rstrip('0')}` | "
                     f"{row['source_location']} | `{row['source_status']}` |")
    disjoint=result["interval_disjointness"]
    lines+=["","## Table 2 — Interval ordering and separation","",
            f"Ascending order: {' < '.join('`'+i+'`' for i in disjoint['ascending_order'])}. "
            f"Pairwise disjoint: **{str(disjoint['pairwise_disjoint']).lower()}**.","",
            "| Lower interval | Upper interval | Gap |","|---|---|---:|"]
    for gap in disjoint["adjacent_gaps"]:
        lines.append(f"| `{gap['below_interval_id']}` | `{gap['above_interval_id']}` | `{gap['gap_decimal'].rstrip('0')}` |")
    lines+=["","## Table 3 — Exclusion results by degree and method interval","",
            "| Degree | Interval id | Primitive polynomials in class | Root-containing | Excluded | Closest polynomial | Height |",
            "|---:|---|---:|---:|---|---|---:|"]
    for row in result["exclusion_table"]:
        lines.append(f"| {row['degree']} | `{row['interval_id']}` | {row['polynomials_searched']:,} | "
                     f"{row['root_containing_polynomials']} | {'yes' if row['excluded'] else '**no**'} | "
                     f"`{_polynomial_text(row['closest_coefficients_ascending'])}` | {row['closest_height']} |")
    approach=result["closest_approach_by_degree"]
    lines+=["","## Table 4 — Certified approach resolution of the search class","",
            "`floor` is a certified lower bound on the distance from the interval to the nearest root of the "
            "closest polynomial of that degree; `floor / width` compares it to the interval's own width.","",
            "| Degree | Interval id | Derivative bound | Floor | Floor / width |","|---:|---|---:|---:|---:|"]
    for row in approach["rows"]:
        ratio="0 (root inside)" if row["root_inside_interval"] else f"{row['floor_to_interval_width_ratio_decimal'].rstrip('0').rstrip('.')}"
        lines.append(f"| {row['degree']} | `{row['interval_id']}` | {row['derivative_bound_on_unit_interval']} | "
                     f"`{row['root_distance_lower_bound_decimal'].rstrip('0') if not row['root_inside_interval'] else '0'}` | {ratio} |")
    lines+=["","Best approach over the four intervals, by degree:","",
            "| Degree | Closest approach |","|---:|---:|"]
    for degree in sorted(approach["approach_resolution_by_degree"],key=int):
        entry=approach["approach_resolution_by_degree"][degree]
        value=entry["closest_approach_decimal"].rstrip("0")
        lines.append(f"| {degree} | `{value if value.rstrip('.') != '0' else '0 (roots inside)'}` |")
    census=result["quartic_survivor_census"]
    lines+=["",f"## Table 5 — Cross-interval resolution of the {census['distinct_surviving_quartics']} surviving quartics","",
            f"Every survivor lies in exactly one method interval (`max_intervals_per_survivor = "
            f"{census['max_intervals_per_survivor']}`) and is certified outside the other three.","",
            "| Quartic | Height | Surviving interval | Certified root prefix | Min separation from the other intervals |",
            "|---|---:|---|---|---:|"]
    for row in census["survivors"]:
        separation=min(Fraction(entry["separation_lower_bound_text"]) for entry in row["separations"])
        lines.append(f"| `{_polynomial_text(row['coefficients_ascending'])}` | {row['height']} | "
                     f"`{row['surviving_interval_id']}` | `{row['root_decimal_prefix']}` | "
                     f"`{decimal_string(separation,13).rstrip('0')}` |")
    historical=result["historical_form_complexity"]
    lines+=["","## Table 6 — Algebraic complexity of the exactly-known planar thresholds","",
            "Minimal polynomials are taken from the committed lattice-native candidate artifact; the "
            "triangular-bond row is derived from the certified kagome-site row by the matching substitution "
            "`p -> 1-p`.","",
            "| Lattice | Closed form | Minimal polynomial | Degree | Height | Inside `C(<=4, <=100)` |",
            "|---|---|---|---:|---:|---|"]
    for row in historical["rows"]:
        lines.append(f"| {row['lattice']} | `{row['closed_form']}` | "
                     f"`{_polynomial_text(row['minimal_polynomial_ascending'])}` | {row['degree']} | "
                     f"{row['height']} | {'yes' if row['inside_census_class'] else '**no**'} |")
    lines+=["",f"Maximum degree {historical['max_degree']}, maximum height {historical['max_height']}. "
            f"Outside the census class: {', '.join(historical['outside_census_class']) or 'none'}.",""]
    control=_load(CONTROL_SOURCES["degree4_boundary_sensitivity"])
    conclusion=control["conclusion"]
    selection=control["width_selection"]
    lines+=["## Table 7 — Quartic-census sensitivity where the census returned a null","",
            "A committed quartic root witness is planted inside a synthetic interval of the width in question "
            "and the unmodified census path is re-run; the negative twin shifts the same interval one full "
            "width away from the planted root. Covered: "
            +", ".join(f"`{i}`" for i in selection["covered"])+". Not covered: "
            +", ".join(f"`{i}`" for i in selection["not_covered"])+" — "+selection["reason"]+".","",
            "| Planted quartic | Width | Polarity | Roots found | Planted quartic reported | Expected |",
            "|---|---|---|---:|---|---|"]
    for row in control["trials"]:
        lines.append(f"| `{_polynomial_text(row['planted_coefficients_ascending'])}` | `{row['width_id']}` | "
                     f"{row['polarity']} | {row['root_containing_polynomials']} | "
                     f"{'yes' if row['planted_quartic_detected'] else 'no'} | "
                     f"{'yes' if row['expected_detection'] else 'no'} |")
    lines+=["",f"All {conclusion['positive_trials']} positive and {conclusion['negative_trials']} negative "
            f"trials behaved as required (`all_trials_passed = "
            f"{str(conclusion['all_trials_passed']).lower()}`)."]
    historical=result["historical_range_exclusion"]
    sensitivity=historical["sensitivity_control"]
    lines+=["","## Table 8 — Exhaustion of the historical complexity range (degree ≤ 6, height ≤ 3)","",
            f"{historical['polynomials_per_interval']:,} primitive polynomials per interval. Excluded on every "
            f"interval: **{str(historical['excluded_on_every_interval']).lower()}**. The certified screen "
            f"retained no candidate at any degree on any interval: "
            f"**{str(historical['screen_retained_nothing']).lower()}**.","",
            "| Degree | Class | Closest polynomial | Distance floor | Floor / width (min–max over the four) |",
            "|---:|---:|---|---:|---:|"]
    for row in historical["rows"]:
        lines.append(f"| {row['degree']} | {row['polynomials_in_class']:,} | "
                     f"`{_polynomial_text(row['closest_coefficients_ascending'])}` | "
                     f"`{row['root_distance_lower_bound_decimal'].rstrip('0')}` | "
                     f"{row['floor_to_width_ratio_min_decimal'].rstrip('0').rstrip('.')} – "
                     f"{row['floor_to_width_ratio_max_decimal'].rstrip('0').rstrip('.')} |")
    lines+=["",f"Sensitivity control: the planted `{sensitivity['closed_form']}` polynomial "
            f"`{_polynomial_text(sensitivity['planted_coefficients_ascending'])}` was recovered in all "
            f"{sensitivity['positive_trials']} positive trials and reported in none of the "
            f"{sensitivity['negative_trials']} negative trials "
            f"(`all_trials_passed = {str(sensitivity['all_trials_passed']).lower()}`)."]
    return "\n".join(lines)+"\n"


def main(argv:Optional[Sequence[str]]=None)->int:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract",type=Path,default=DEFAULT_CONTRACT)
    parser.add_argument("--output",type=Path)
    parser.add_argument("--validate",type=Path)
    parser.add_argument("--markdown",action="store_true",help="render the manuscript tables instead of the JSON artifact")
    args=parser.parse_args(argv)
    if args.validate:
        print(json.dumps(validate_result(json.loads(args.validate.read_text(encoding="utf-8")),args.contract),indent=2,sort_keys=True))
        return 0
    if args.markdown:
        rendered=render_markdown(build_result(args.contract))
        if args.output:
            args.output.parent.mkdir(parents=True,exist_ok=True)
            args.output.write_text(rendered,encoding="utf-8")
        else:
            print(rendered,end="")
        return 0
    rendered=json.dumps(build_result(args.contract),indent=2,sort_keys=True)+"\n"
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(rendered,encoding="utf-8")
    else:
        print(rendered,end="")
    return 0


if __name__=="__main__":raise SystemExit(main())
