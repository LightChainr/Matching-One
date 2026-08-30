# P250 common-counter radius selector

Primary gate: all four two-point denominators plus a nondegenerate cubic covariance at at least two separations.

| radius | d=1 min z / usable | d=2 min z / usable | d=3 min z / usable | pass |
|---:|---:|---:|---:|---|
| R1 | 4.51 / True | 0.48 / False | 0.00387 / False | False |
| R2 | 4.6 / True | 0.642 / False | 0.122 / False | False |
| R3 | 4.52 / True | 0.756 / False | 0.43 / False | False |
| R4 | 6.09 / True | 1.11 / False | 0.214 / False | False |

Decision: `no_production_candidate_in_radius_1_to_4_local_landing_h4_family`.

The R4/d3 cubic zero score is descriptive only: its two-point denominator fails the frozen gate, so it cannot select R4.

Next selector: freeze a different charged insertion, such as a leg-defect or mesoscopic row, before acquiring more cubic replicas
