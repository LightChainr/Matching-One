#!/usr/bin/env python3
"""Issue #593 Phase B Deliverable 1: projected memory degree, widths 4-10.

Frozen configuration is the repository's own code end to end
(``frozen_span``, ``memory_kernel``, ``kernel_statistics``, ``block_hankel``
from ``p398_projected_memory``), applied to a width-9/10 generator shim that
is bit-exact against the repository at widths 1..8.  Everything is at
``eta = 0`` on the declared kernel grid ``t = 0, 0.25, ..., 4.0``.

One deliberate departure, declared on the ticket: singular values of the
block Hankel matrix come from a true SVD (numpy) instead of the
eigendecomposition of ``H^T H``, exactly as the ticket asks.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from p593_noncrossing_fast import WideGenerator, noncrossing_states_fast
from p398_projected_memory import (
    KERNEL_GRID,
    block_hankel,
    frozen_span,
    kernel_statistics,
    memory_kernel,
    unresolved_columns,
)

OUTPUT = Path(__file__).resolve().parents[1] / "results" / "p593-width9-10-memory-degree" / "raw"


def analyse(width: int) -> dict:
    start = time.time()
    generator = WideGenerator(width, noncrossing_states_fast(width))
    size = generator.size
    rates = generator.baseline_rates()
    rows = generator.rows(rates)
    rate = generator.exit_rate(rates)

    basis = frozen_span(generator, rows)
    kernels = memory_kernel(basis, rows, size, rate, KERNEL_GRID)
    stats = kernel_statistics(kernels, KERNEL_GRID)

    # rank C from a true SVD of the coupling matrix C = Q G Phi
    columns = unresolved_columns(basis, rows)
    c_matrix = np.array(columns, dtype=float)
    c_singular = np.linalg.svd(c_matrix, compute_uv=False)
    c_rank = int(sum(1 for v in c_singular if v > 1e-10 * c_singular[0]))

    # block Hankel from samples 1..16 (drops t=0), 8x8 blocks of size 6, true SVD
    hankel = np.array(block_hankel(kernels), dtype=float)
    singular = np.linalg.svd(hankel, compute_uv=False)
    largest = float(singular[0])
    numerical_rank = int(sum(1 for v in singular if v > 1e-6 * largest))
    total = float(np.sum(singular**2))
    effective = {}
    for level in (0.99, 0.999):
        running, count = 0.0, 0
        for value in singular:
            running += value * value
            count += 1
            if running >= level * total:
                break
        effective[f"{level:g}"] = count

    return {
        "width": width,
        "states": size,
        "eta": 0.0,
        "frozen_rank": len(basis),
        "rank_C": c_rank,
        "C_singular_values": [float(v) for v in c_singular],
        "hankel_shape": list(hankel.shape),
        "hankel_numerical_rank_tol_1e-6": numerical_rank,
        "hankel_effective_orders": effective,
        "hankel_singular_values_normalized_top16": [
            float(v / largest) for v in singular[:16]
        ],
        "integrated_norm": stats["integrated_weight"],
        "decay_time": stats["decay_time"],
        "tail_mass_fraction": stats["tail_mass_fraction"],
        "norm_profile": stats["norm_profile"],
        "seconds": round(time.time() - start, 1),
    }


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    widths = [int(a) for a in sys.argv[1:]] or [4, 5, 6, 7, 8, 9, 10]
    for width in widths:
        result = analyse(width)
        path = OUTPUT / f"width{width}.json"
        path.write_text(json.dumps(result, indent=2))
        print(
            json.dumps(
                {
                    "width": width,
                    "rank_C": result["rank_C"],
                    "hankel_rank": result["hankel_numerical_rank_tol_1e-6"],
                    "orders": result["hankel_effective_orders"],
                    "integrated": result["integrated_norm"],
                    "decay": result["decay_time"],
                    "tail": result["tail_mass_fraction"],
                    "seconds": result["seconds"],
                }
            ),
            flush=True,
        )


if __name__ == "__main__":
    main()
