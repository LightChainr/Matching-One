#!/usr/bin/env python3
"""Fixed all-147 archived-prefix expansion; existing twelve outputs are reused."""
import csv
from fractions import Fraction
import hashlib
import json
from math import comb
from pathlib import Path
import signal
import time

from p334_checkpoint_scalar_collision import archived_permutation
from p334_contracted_birth_network import build
from p334_full_birth_reliability import safety_polynomial
from p334_pair_only_survival import frac, multiply

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "results/p334-all147-prefix-clocks"
MANIFEST = ROOT / "analysis/p334_all147_clock_manifest.json"


def timeout(_signum, _frame):
    raise TimeoutError("frozen 20-second per-prefix wall limit")


def main():
    manifest = json.loads(MANIFEST.read_text())
    original = json.loads((ROOT / manifest["source_selection_manifest"]).read_text())
    csv_path = ROOT / manifest["source_archive"]
    metadata = json.loads((ROOT / manifest["source_metadata"]).read_text())
    with csv_path.open(newline="") as stream:
        rows = [r for r in csv.DictReader(stream)
                if int(r["replica"]) not in original["exclude_counters"]
                and all(r[k] == str(v) for k,v in original["eligibility"].items() if k != "checkpoint_rank")
                and int(r["k1"]) <= int(r["k0"]) < int(r["k2"])]
    rows.sort(key=lambda r:int(r["replica"]))
    if len(rows) != manifest["expected_eligible_rows"]:
        raise ValueError("fixed source count differs; no outcome-based repair")
    old_clocks = {r["counter"]:r for r in json.loads((ROOT / "results/p334-twelve-prefix-clocks/full_clocks.json").read_text())["records"]}
    OUTPUT.mkdir(parents=True,exist_ok=True)
    (OUTPUT / "prefixes").mkdir(exist_ok=True)
    (OUTPUT / "selection.json").write_text(json.dumps({"manifest_commit":"8d7ac0e9", "manifest_sha256":hashlib.sha256(MANIFEST.read_bytes()).hexdigest(),
                                                        "source_sha256":hashlib.sha256(csv_path.read_bytes()).hexdigest(),
                                                        "selected_original_rows":rows, "reused_counters":sorted(old_clocks)},indent=2,sort_keys=True)+"\n")
    signal.signal(signal.SIGALRM,timeout)
    results = []
    batch_started = time.monotonic()
    for index,row in enumerate(rows):
        counter = int(row["replica"])
        if counter in old_clocks:
            r = dict(old_clocks[counter])
            r.update(reused_original_twelve=True, new_evaluation_seconds=0,
                     factor_polynomial_artifact=f"results/p334-component-birth-race/{counter}.json",
                     mapping_artifact=f"results/p334-twelve-prefix-clocks/maps/{counter}.json")
            results.append(r)
            continue
        artifact = OUTPUT / "prefixes" / f"{counter}.json"
        if artifact.exists():
            saved = json.loads(artifact.read_text())
            results.append(saved["clock"])
            continue
        r = {"counter":counter,"original_row":row,"status":"not_started","reused_original_twelve":False,
             "prefix_artifact":str(artifact.relative_to(ROOT))}
        saved = {"mapping":None,"factors":[],"clock":r}
        started = time.monotonic()
        signal.setitimer(signal.ITIMER_REAL,manifest["per_prefix_limits"]["wall_seconds"])
        try:
            replay = {"replica_counter":counter,"seed":metadata["seed"],"N":425,"k0":252,
                      "occupied_prefix_labels":archived_permutation(425,metadata["seed"],counter)[:252],
                      "period_matrix":[[425,268],[0,1]],"ell":[12,-19]}
            mapping = build(replay)
            saved["mapping"] = mapping
            components = mapping["port_components"]
            networks = [c for c in components if "two_terminal_network" in c]
            r["structure"] = {"address_class_counts":[len(c["addresses"]) for c in components],
                              "occupied_components":len(mapping["occupied_components"]),
                              "essential_carriers":len(mapping["essential_component_roots"]),
                              "parallel_two_port_factors":len(networks)}
            if not networks or any(len(c["addresses"])>2 for c in components):
                r["status"] = "not_solved_non_two_port"
            else:
                coefficients,core_sites,peak_states,width,joins = [1],set(),0,0,0
                for ci,c in enumerate(networks):
                    network = c["two_terminal_network"]
                    sites = set(network["vacant_sites"])
                    if core_sites & sites:
                        raise ValueError("factor site overlap")
                    core_sites |= sites
                    fc,stats = safety_polynomial(network,sites)
                    coefficients = multiply(coefficients,fc)
                    saved["factors"].append({"component":ci,"site_labels":sorted(sites),"n_sites":len(sites),
                                             "port_addresses":c["addresses"],"safe_coefficients":fc,
                                             "dp":{k:v for k,v in stats.items() if k!="bag_state_counts"}})
                    peak_states=max(peak_states,stats["maximum_states"])
                    width=max(width,stats["treewidth_upper_bound"])
                    joins+=stats["join_pairs"]
                free=173-len(core_sites)
                coefficients=multiply(coefficients,[comb(free,k) for k in range(free+1)])
                if coefficients[:3] != [1,173-int(row["H2"]),int(row["checkpoint_b2_safe_pairs"])]:
                    raise ValueError("full-clock coefficients mismatch original exact singleton/pair fields")
                survival=[Fraction(coefficients[k],comb(173,k)) for k in range(174)]
                mean=sum(survival[:-1])
                r.update(status="solved_full_physical",true_safe_counts=coefficients,
                         true_survival=[frac(s) for s in survival],mean_true_birth_step=frac(mean),
                         variance_true_birth_step=frac(sum((2*k+1)*survival[k] for k in range(173))-mean*mean),
                         tail_after_20=frac(survival[20]),tail_after_40=frac(survival[40]),tail_after_65=frac(survival[65]),
                         maximum_true_safe_k=max(k for k,c in enumerate(coefficients) if c),
                         birth_quantiles={str(q):next(k for k in range(174) if survival[k]<=1-q) for q in (Fraction(1,10),Fraction(1,2),Fraction(9,10))},
                         dp={"maximum_states":peak_states,"join_pairs":joins})
                r["structure"].update(core_random_sites=len(core_sites),off_core_random_sites=free,treewidth_upper_bound=width)
        except TimeoutError as error:
            r.update(status="not_solved_time_limit",reason=str(error))
        except (ValueError,RuntimeError) as error:
            r.update(status="not_solved_structure_or_state_limit",reason=str(error))
        finally:
            signal.setitimer(signal.ITIMER_REAL,0)
        r["new_evaluation_seconds"]=time.monotonic()-started
        artifact.write_text(json.dumps(saved,separators=(",",":"),sort_keys=True)+"\n")
        results.append(r)
        if index%10==0 or r["status"]!="solved_full_physical":
            print(index+1,"/147",counter,r["status"],"mean",r.get("mean_true_birth_step",{}).get("value"),
                  "width",r.get("structure",{}).get("treewidth_upper_bound"),"seconds",r["new_evaluation_seconds"],flush=True)
    output={"manifest_commit":"8d7ac0e9","new_samples":0,"selected_rows":len(rows),"reused_rows":len(old_clocks),
            "new_rows":len(rows)-len(old_clocks),"batch_wall_seconds":time.monotonic()-batch_started,"records":results}
    (OUTPUT/"full_clocks.json").write_text(json.dumps(output,indent=2,sort_keys=True)+"\n")
    print("DONE",len(results),"solved",sum(r["status"]=="solved_full_physical" for r in results),"seconds",output["batch_wall_seconds"],flush=True)


if __name__=="__main__":
    main()
