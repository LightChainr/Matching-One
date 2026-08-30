# P418 paired-anchor 5k pilot

Decision: `one_anchor_noise_visible_but_no_mask_specific_or_current_vs_full_bias`.

| channel | current masked d2 (p) | independent masked d2 (p) | full masked d2 (p) | channel decision |
|---|---:|---:|---:|---|
| plus_r1 | 75.185 (0.0198) | 85.239 (0.009901) | 0 (1) | one_anchor_only_tension |
| plus_r2 | 96.6315 (0.009901) | 40.1255 (0.495) | 0 (1) | one_anchor_only_tension |
| minus_r1 | 49.5738 (0.1881) | 38.4754 (0.5545) | 0 (1) | no_estimator_rejects_at_5k |
| minus_r2 | 46.3893 (0.2178) | 65.6649 (0.0198) | 0 (1) | no_estimator_rejects_at_5k |

The current estimator replays 25 historical radius-4 batches with maximum summed-coordinate error `1.24e-14`.

Maximum absolute masked-minus-raw distance across all estimators/channels: `0.00180318`. Current-minus-independent has Hotelling p<0.01 only in `['plus_r2']`; current-minus-full has none: `[]`.

The full-anchor rows have 41 resolved covariance modes and masked-design rank 41, so their exact zero distance is structurally saturated rather than an independent acceptance of the CRT mask. The causal observation is instead that switching anchor streams changes ordinary finite-sample cone distance but creates essentially no masked-minus-raw penalty.

For scale only, the committed archive score has masked-minus-raw increments at least `513.144`, more than `2.85e+05` times this pilot maximum. This is descriptive because the sample sizes and block assembly differ, but it rules out reproducing the archive-specific penalty in the paired radius-4 gate.

All 100 batch rows retain the joint current/independent/full × hand × charge coordinates, so the complete paired covariance remains reconstructible. This is the only authorized 5k pilot; no production extension was run.
