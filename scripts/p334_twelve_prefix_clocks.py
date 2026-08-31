#!/usr/bin/env python3
"""Frozen deterministic twelve-prefix full physical clock extraction; no MC."""
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
OUTPUT = ROOT / "results/p334-twelve-prefix-clocks"
MANIFEST = ROOT / "analysis/p334_twelve_prefix_selection_manifest.json"


def timeout(_signum, _frame):
    raise TimeoutError("frozen 20-second per-prefix wall limit")


def main():
    manifest = json.loads(MANIFEST.read_text())
    source_path = ROOT / manifest["source_archive"]
    metadata = json.loads((ROOT / manifest["source_metadata"]).read_text())
    with source_path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    filters = manifest["eligibility"]
    eligible = [r for r in rows if int(r["replica"]) not in manifest["exclude_counters"]
                and all(r[k] == str(v) for k,v in filters.items() if k != "checkpoint_rank")
                and int(r["k1"]) <= int(r["k0"]) < int(r["k2"])]
    eligible.sort(key=lambda r: int(r["replica"]))
    if len(eligible) < manifest["requested_count"]:
        raise RuntimeError(f"only {len(eligible)} eligible prefixes; selection not widened")
    chosen = eligible[:manifest["requested_count"]]
    OUTPUT.mkdir(parents=True, exist_ok=True)
    (OUTPUT / "maps").mkdir(exist_ok=True)
    selection = {"manifest_commit": "b9cbe13e", "manifest_sha256": hashlib.sha256(MANIFEST.read_bytes()).hexdigest(),
                 "source_sha256": hashlib.sha256(source_path.read_bytes()).hexdigest(),
                 "eligible_prefixes": len(eligible), "selected_original_rows": chosen,
                 "checkpoint_rank_definition": "k1 <= k0 < k2; source metadata restricts geometry-pilot rows to rank one"}
    (OUTPUT / "selected_prefixes.json").write_text(json.dumps(selection, indent=2, sort_keys=True) + "\n")
    signal.signal(signal.SIGALRM, timeout)
    previous_path = OUTPUT / "full_clocks.json"
    previous = json.loads(previous_path.read_text()) if previous_path.exists() else {"records": []}
    previous_by_counter = {r["counter"]: r for r in previous["records"]}
    if previous["records"] and not (OUTPUT / "single_factor_initial_pass.json").exists():
        (OUTPUT / "single_factor_initial_pass.json").write_text(json.dumps(previous, indent=2, sort_keys=True)+"\n")
    results = []
    for row in chosen:
        counter, k0 = int(row["replica"]), int(row["k0"])
        if previous_by_counter.get(counter, {}).get("status") == "solved_full_physical":
            result = previous_by_counter[counter]
            result["reused_completed_single_factor_result"] = True
            result["structure"]["parallel_two_port_factors"] = 1
            results.append(result)
            continue
        prefix = archived_permutation(425, metadata["seed"], counter)[:k0]
        replay = {"replica_counter": counter, "seed": metadata["seed"], "N": 425,
                  "k0": k0, "occupied_prefix_labels": prefix,
                  "period_matrix": [[425,268],[0,1]], "ell": [12,-19]}
        result = {"counter": counter, "original_row": row, "status": "not_started"}
        started = time.monotonic()
        signal.setitimer(signal.ITIMER_REAL, manifest["per_prefix_limits"]["wall_seconds"])
        try:
            map_path = OUTPUT / "maps" / f"{counter}.json"
            if map_path.exists():
                mapping = json.loads(map_path.read_text())
            else:
                mapping = build(replay)
                map_path.write_text(json.dumps(mapping, indent=2, sort_keys=True) + "\n")
            components = mapping["port_components"]
            networks = [c for c in components if "two_terminal_network" in c]
            result["structure"] = {"occupied_components": len(mapping["occupied_components"]),
                                   "essential_carriers": len(mapping["essential_component_roots"]),
                                   "address_class_counts": [len(c["addresses"]) for c in components],
                                   "contracted_gain_edges": len(mapping["contracted_edges_with_transverse_gain"])}
            if not networks or any(len(c["addresses"]) > 2 for c in components):
                result["status"] = "not_solved_non_two_port"
            else:
                counts, core_sites, shared, factor_stats = [1], set(), set(), []
                for component in networks:
                    network = component["two_terminal_network"]
                    factor_sites = set(network["vacant_sites"])
                    if core_sites & factor_sites:
                        raise ValueError("parallel factors are not site-disjoint")
                    core_sites |= factor_sites
                    port_sets = [set(v for v,a in component["ports"] if a == address) for address in component["addresses"]]
                    shared |= set.intersection(*port_sets)
                    factor_counts, stats = safety_polynomial(network, factor_sites)
                    counts = multiply(counts, factor_counts)
                    factor_stats.append({"core_sites": len(factor_sites), "port_addresses": component["addresses"],
                                         "port_site_counts": [len(s) for s in port_sets],
                                         "dp": {k:v for k,v in stats.items() if k != "bag_state_counts"}})
                free = 173-len(core_sites)
                counts = multiply(counts, [comb(free,k) for k in range(free+1)])
                result["structure"].update({"parallel_two_port_factors": len(networks),
                                            "factor_details": factor_stats, "shared_port_sites": sorted(shared),
                                            "core_random_sites": len(core_sites), "off_core_random_sites": free,
                                            "treewidth_upper_bound": max(f["dp"]["treewidth_upper_bound"] for f in factor_stats)})
                stats = {"factor_count": len(networks), "maximum_states": max(f["dp"]["maximum_states"] for f in factor_stats),
                         "join_pairs": sum(f["dp"]["join_pairs"] for f in factor_stats),
                         "elapsed_seconds": sum(f["dp"]["elapsed_seconds"] for f in factor_stats)}
                if counts[:3] != [1, 173-int(row["H2"]), int(row["checkpoint_b2_safe_pairs"])]:
                    raise ValueError("whole-event coefficients do not match original exact singleton/pair fields")
                survival = [Fraction(counts[k], comb(173,k)) for k in range(174)]
                hazard = [None] + [1-survival[k]/survival[k-1] if survival[k-1] else None for k in range(1,174)]
                mean = sum(survival[:-1])
                second = sum((2*k+1)*survival[k] for k in range(173))
                result.update({"status": "solved_full_physical", "true_safe_counts": counts,
                               "true_survival": [frac(s) for s in survival],
                               "true_hazard": [frac(h) if h is not None else None for h in hazard],
                               "mean_true_birth_step": frac(mean), "variance_true_birth_step": frac(second-mean*mean),
                               "tail_after_20": frac(survival[20]), "tail_after_40": frac(survival[40]),
                               "tail_after_65": frac(survival[65]),
                               "birth_quantiles": {str(q): next(k for k in range(174) if survival[k] <= 1-q)
                                                   for q in (Fraction(1,10),Fraction(1,2),Fraction(9,10))},
                               "maximum_true_safe_k": max(k for k,c in enumerate(counts) if c),
                               "dp": {k:v for k,v in stats.items() if k != "bag_state_counts"}})
        except TimeoutError as error:
            result.update(status="not_solved_time_limit", reason=str(error))
        except (ValueError, RuntimeError) as error:
            result.update(status="not_solved_structural_or_state_limit", reason=str(error))
        finally:
            signal.setitimer(signal.ITIMER_REAL, 0)
        result["single_process_seconds"] = time.monotonic()-started
        results.append(result)
        (OUTPUT / "full_clocks.json").write_text(json.dumps({"manifest_commit": "b9cbe13e", "new_samples": 0,
                                                             "method_increment": "product of vertex-disjoint two-port safety factors; original four solved single-factor outputs reused",
                                                             "records": results}, indent=2, sort_keys=True)+"\n")
        print(counter, result["status"], "H2", row["H2"], "mean", result.get("mean_true_birth_step",{}).get("value"),
              "tail40", result.get("tail_after_40",{}).get("value"), "width", result.get("structure",{}).get("treewidth_upper_bound"),
              "seconds", result["single_process_seconds"], flush=True)
    (OUTPUT / "full_clocks.json").write_text(json.dumps({"manifest_commit": "b9cbe13e", "new_samples": 0,
                                                         "method_increment": "product of vertex-disjoint two-port safety factors; original four solved single-factor outputs reused",
                                                         "records": results}, indent=2, sort_keys=True)+"\n")


if __name__ == "__main__":
    main()
