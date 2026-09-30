# Directional next-insertion hazard: one retrospective check

Same completed 140k block, unchanged cutoffs, cohorts and exact-D overlap weights.

| Quantity | Estimate | Aligned batch SE |
|---|---:|---:|
| pooled_nu_early_minus_late | -2.184732729 | 0.4794961673 |
| within_D_nu_early_minus_late | 1.629168057 | 0.3909761647 |
| within_D_survival_early_minus_late | -0.004149992812 | 0.004511941724 |
| within_D_next_insertion_exit_difference | 1.526024089e-05 | 3.662231425e-06 |
| within_D_window_scaled_exit_difference | 0.01121627705 | 0.002691740097 |

Positive nu contrast means greater early-cohort next-insertion exit probability.
Window scaling is not a finite-window exit estimate. The last survival contrast is the unchanged original result.

This necessity check was selected after the original report; it is not independent validation.
No new observations, subgroup searches, thresholds, cutoffs or samples are added.
