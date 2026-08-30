# Linear evaluator for the frozen two-cell block event

Status: deterministic engineering slice for Issue 112. No probability or threshold bound is claimed.

## Algorithm

The input is an explicit Boolean occupation flag for every site of an open `2s x s` rectangle.
The evaluator constructs two disjoint-set forests in one scan over the microscopic edges:

1. the **full** forest includes every edge and is used only for the final connection test;
2. the **halves** forest suppresses edges crossing `x=s` and is used to select the unique largest
   open component independently inside each `s x s` half.

An empty half fails. A tie at the largest size fails. Otherwise the selected left and right
components succeed exactly when representative sites share a root in the full forest. This preserves
the ordering frozen by the tiny exact oracle.

The rectangle has `2s^2` sites. The square graph scans `4s^2-3s` full edges and `4s(s-1)` within-half
edges. Adding both diagonals gives `8s^2-9s+2` full edges and `8s^2-12s+4` within-half edges. Thus
time and memory are both `O(s^2)`.

## Validation gate

The implementation is compared configuration-by-configuration with the independent set/BFS oracle:

- both graphs at `s=1`: `2 * 4` cases;
- both graphs at `s=2`: `2 * 256` cases;
- total: 520 event decisions, including reasons and the selected clusters.

The required mismatch count is zero. Deterministic all-open controls at `s=64` additionally exercise
8,192 sites and require both selected halves to contain exactly 4,096 sites.

The square nonmonotonicity witness remains explicit: mask 6 succeeds, while adding site 4 to obtain
mask 22 creates a largest-cluster tie and fails. The matching graph connects that diagonal pair.

## Why there is no sampler here

A reproducible PRNG stream is useful engineering, but reproducibility is not a proof that final trials
are independent Bernoulli observations. Mixing event semantics, RNG policy, exploratory reuse, and the
confidence calculation in one step would blur the auditable boundary. This slice therefore accepts
only caller-supplied deterministic configurations.

Still open are the sampling protocol, randomness and domain-separation policy, production probability
measurement, and any certified statement about `p_c`.
