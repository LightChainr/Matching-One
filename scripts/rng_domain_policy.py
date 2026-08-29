#!/usr/bin/env python3
"""Derive and record production RNG domains for Matching One.

The default keeps common random numbers inside one same-N orientation pair and
domain-separates distinct sizes.  Cross-size coupling is an explicit exception
that must name the residual whose variance it is intended to reduce.
"""

from __future__ import annotations

import argparse
import hashlib
import json


SCHEMA = "matching-one-rng-domain-v1"


def derive_size_seed(base_seed: int, experiment_tag: str, size: int) -> int:
    if not experiment_tag or "\x00" in experiment_tag:
        raise ValueError("experiment_tag must be nonempty and contain no NUL")
    if not 0 <= base_seed < 2**64:
        raise ValueError("base_seed must fit uint64")
    if size <= 0:
        raise ValueError("size must be positive")
    payload = f"{SCHEMA}\0{experiment_tag}\0N={size}\0seed={base_seed}".encode()
    return int.from_bytes(hashlib.sha256(payload).digest()[:8], "big")


def domain_record(base_seed: int, experiment_tag: str, size: int, *,
                  mode: str = "domain_separated", coupled_residual: str | None = None):
    if mode not in {"domain_separated", "intentional_cross_size_coupling"}:
        raise ValueError("unknown RNG domain mode")
    if mode == "domain_separated":
        if coupled_residual:
            raise ValueError("coupled_residual is forbidden for domain-separated runs")
        effective_seed = derive_size_seed(base_seed, experiment_tag, size)
        covariance = "independent_across_sizes"
    else:
        if not coupled_residual:
            raise ValueError("intentional coupling must name its prespecified residual")
        effective_seed = base_seed
        covariance = "aligned_batches_full_cross_size_covariance_required"
    return {
        "schema": SCHEMA,
        "mode": mode,
        "experiment_tag": experiment_tag,
        "size": size,
        "base_seed": base_seed,
        "effective_seed": effective_seed,
        "within_same_size_orientation_pair": "shared_common_field",
        "cross_size_covariance_contract": covariance,
        "coupled_residual": coupled_residual,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-seed", type=int, required=True)
    parser.add_argument("--experiment-tag", required=True)
    parser.add_argument("--size", type=int, required=True)
    parser.add_argument("--mode", choices=["domain_separated", "intentional_cross_size_coupling"],
                        default="domain_separated")
    parser.add_argument("--coupled-residual")
    args = parser.parse_args()
    print(json.dumps(domain_record(
        args.base_seed, args.experiment_tag, args.size,
        mode=args.mode, coupled_residual=args.coupled_residual,
    ), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
