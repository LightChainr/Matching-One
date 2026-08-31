#!/usr/bin/env python3
"""Fixed 01/10 gate/clock half-difference products; production reader pending."""
import numpy as np

CELLS = ("01", "10")
FIELDS = (
    "prevalence", "mean_gate_R0_first", "mean_gate_R1_second",
    "same_label_joint_gate", "gate_cov_R0_R0", "gate_cov_R0_R1", "gate_cov_R1_R1",
    "gate_R0_cov_A", "gate_R0_cov_E", "gate_R1_cov_A", "gate_R1_cov_E",
    "paired_mean_A", "paired_mean_E",
)


def quartet_vector(ranks, next_ranks, births, n, delta_cos4):
    """Return (cell, 13-vector), or (None, zeros) outside 01/10.

    ranks: original (first,second) ranks.
    next_ranks: [label U/V, original orientation first/second].
    births: [label U/V, suffix 0/1, original orientation first/second, K1/K2].
    These are the common next-label ranks, not terminal ranks or future marks.
    """
    ranks = np.asarray(ranks, dtype=int)
    cell = "".join(str(x) for x in ranks)
    if cell not in CELLS:
        return None, np.zeros(len(FIELDS))
    next_ranks = np.asarray(next_ranks, dtype=int)
    births = np.asarray(births, dtype=float)
    if next_ranks.shape != (2, 2) or births.shape != (2, 2, 2, 2):
        raise ValueError("A quartet needs two labels by two paired suffixes")
    low, high = (0, 1) if cell == "01" else (1, 0)
    g = np.column_stack((next_ranks[:, low] >= 1, next_ranks[:, high] == 2)).astype(float)
    a = 1-(births[..., 0]+births[..., 1])/(n+1)
    e = 1-(births[..., 1]-births[..., 0])/(n+1)
    observations = np.stack(((a[..., 0]-a[..., 1])/delta_cos4,
                             (e[..., 0]-e[..., 1])/delta_cos4), axis=-1)
    m = observations.mean(axis=1)
    dg, dm = g[0]-g[1], m[0]-m[1]
    gate_cov = .5*np.outer(dg, dg)
    response_cov = .5*np.outer(dg, dm)
    return cell, np.r_[1., g.mean(axis=0), np.prod(g, axis=1).mean(),
                       gate_cov[0, 0], gate_cov[0, 1], gate_cov[1, 1],
                       response_cov.ravel(), m.mean(axis=0)]
