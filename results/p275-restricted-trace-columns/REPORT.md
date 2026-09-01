# #275 restricted-trace candidate-column audit

Status: `PARTIAL_COLUMNS_COMPLETE_EXISTING_COVARIANCES_NOT_DIRECTLY_SCOREABLE`.

This analysis used the existing theory and covariance assets only.  It generated
zero random samples and did not start the PR #546 queue.

## Result

After quotienting the two common-normalizer directions, the aligned
restricted-trace jet has four physical coordinates

```text
(eta_B, xi_B, partial_p eta_B, partial_p xi_B).
```

The vacuum/Ward identity available on the critical surface fixes
`eta_B=0`, leaving a rank-3 family.  The stronger neighbourhood version would
also fix `partial_p eta_B=0`, leaving rank 2.  The current thermal-Q4/Jordan
assets fix the downstream `E4hat` modulus shape and Jordan `log N` slope, but
no same-source restricted-sector relation; their honest current envelope has
rank 4.  Therefore the critical Ward image is a rank-3 subspace contained in
the rank-4 thermal envelope.  The present theory is short by one restricted
sector direction at the critical surface, followed by a numeric source scale
and the map through the original-U normalizer.

The covariance matrices are not the bottleneck.  The rho-child covariance has
rank 9/9 and the global K1/K2 covariance has rank 12/12, but neither archive
contains any of the six same-`B` source-jet coordinates.  Their direct
candidate-forward-map rank is therefore zero until the missing theory column
is supplied.

As a conditional diagnostic only, the existing K1/K2 covariance distinguishes
the fixed transfer spaces algebraically: fixed semisimple `kappa=0.5` and
Jordan `kappa=1` each have rank 8 in 12 observations, their ordinary
intersection has dimension 4, and the combined rank is 12.  The GLS replay is

```text
fixed semisimple proxy: chi2 = 11.7924660378 on 4 residual df
Jordan proxy:           chi2 =  6.43252655289 on 4 residual df
```

These are proxy transfer scores, not physical vacuum/Ward versus
thermal-Q4/Jordan field scores.  The literal candidate decision remains
`UNIDENTIFIABLE_WITH_CURRENT_ASSETS`.

## Execution provenance

XPk2PZ was requested and started from Ready, but its saved key failed with
`Permission denied (publickey)`.  No key reset was attempted, no remote file
was written, and no remote analysis command ran.  The owned local tunnel was
stopped and XPk2PZ was sent back toward Ready.  The deterministic analysis then
ran on the local Mac as an explicitly recorded fallback.

```text
python3 -m unittest discover -s tests -p 'test_p275_restricted_trace_columns.py' -v
python3 scripts/analyze_p275_restricted_trace_columns.py \
  --output results/p275-restricted-trace-columns/latest.json
```

Four focused tests passed.  The analysis wall time recorded inside the result
was 0.0217 s; the complete command used 0.12 s wall time and about 36.5 MB
maximum RSS on local Python 3.13.7 / NumPy 2.4.6.

## Next direct calculation

For each named microscopic source, fix `partial_p eta_B`, `xi_B`, and
`partial_p xi_B` with its source units, rank-1 denominator, and pooled-root
original-U map.  Then apply one existing-covariance score; no new samples are
needed for that score.
