# P275 post-reveal mean-J_D4 discovery score

This reuses the frozen nine streams but replaces `Cov(A_top,J_D4)/B` with
transported `E[J_D4]`. It does not change the Phase-1 result.

| geometry | p_N | Re N^(13/8)E[J_D4] | Im | Re E[J_D4] | Im |
|---|---:|---:|---:|---:|---:|
| N50/i | 0.592615702394 | 25.728448 | -0.027244878 | 0.044625936 | -4.725618e-05 |
| N50/2i | 0.593162621815 | -23.057641 | -1.4390164 | -0.039993426 | -0.0024959706 |
| N50/5i_over_2 | 0.593010690818 | -8.8545335 | -1.8357623 | -0.015358168 | -0.0031841256 |
| N130/i | 0.59273769904 | 30.025465 | -0.22801429 | 0.011023812 | -8.3715161e-05 |
| N130/2i | 0.592688170172 | 50.12825 | -2.5969661 | 0.018404524 | -0.00095347285 |
| N130/5i_over_2 | 0.592830160517 | -45.8191 | -3.0980715 | -0.016822425 | -0.0011374531 |
| N170/i | 0.592721625884 | 67.617506 | -0.14017138 | 0.016053867 | -3.3279735e-05 |
| N170/2i | 0.592726207363 | 49.446256 | -5.5218263 | 0.011739617 | -0.0013110017 |
| N170/5i_over_2 | 0.592709301529 | 46.079873 | -4.1163181 | 0.010940364 | -0.00097730345 |

## Same Phase-1 model bases

| model | chi2 | df | survival p |
|---|---:|---:|---:|
| Q4_epsilon_ordinary | 25516 | 16 | 0 |
| Q4_energy_Jordan | 17329.9 | 12 | 0 |
| generic_allowed_H4_pure | 2678.75 | 12 | 0 |
| generic_allowed_H4_affine_log | 868.748 | 6 | 2.13982e-184 |
| zero_response | 33410.4 | 18 | 0 |

## Reading

No Phase-1 basis survives for the single-quotient mean source either. The least
bad basis is generic affine-log H4 (`chi2=868.75/6`), still decisively rejected.
This does not contradict commit `634040d`: its `0.33085` q2 ratio belongs to a
paired-orientation P4 difference on an exact N65/N130 parent-child chain. These
nine rows are individual quotients at N50/N130/N170 and contain neither that
orientation subtraction nor its exact radial chain. Therefore the q2 mean-source
hint remains a separate candidate, while the present nine-geometry archive cannot
promote transported single-quotient `E[J_D4]` to the field-identity covector.

## Boundary

- This score was chosen after the connected Gamma reveal and is discovery-only.
- It reuses the frozen samples, roots, transport, scaling, covariance and Phase-1 model bases.
- It does not alter the frozen Phase-1 rejection of the Gamma scaling map.
- The nine rows are single quotients, not the paired-orientation P4 q2 parent/child difference scored in 634040d.
- A surviving basis nominates a preregistered local-source follow-up; it does not identify a field by itself.
