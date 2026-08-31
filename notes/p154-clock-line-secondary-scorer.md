# P154 two clock lines: implementation receipt before outcome reading

This thin scorer implements rule commit
`83f3eba88d7f1290704f82610c28669dc5e12f3c`, using the exact flat-value
gain formula from `c2828e3430fe1ac7e02fbe0e5ddc0e6a24c99847`.
It is written and committed before its author reads any fresh P154 value.
This receipt does not run the scorer or modify the official primary result.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
/Users/lc/python-envs/research-py311/bin/python \
  scripts/p154_clock_line_secondary.py \
  --official-package OFFICIAL_PACKAGE_DIRECTORY \
  --official-result OFFICIAL_PACKAGE_DIRECTORY/PROSPECTIVE_RESULT.json \
  --raw-source-commit FULL_RAW_DELIVERY_SHA \
  --official-result-commit FULL_OFFICIAL_RESULT_SHA \
  --output NEW_SECONDARY_JSON
```

Raw shards default to `OFFICIAL_PACKAGE_DIRECTORY/production`; `--raw-dir`
may name their completed delivery directory. Inputs are never overwritten.
An existing output is refused. Commit IDs record caller-declared provenance;
actual input SHA256 hashes and official receipt correspondence are checked
and retained. The metadata also records implementation checkout SHA and
script SHA256, production freeze and secondary-rule commit.

Only the exact nine completed shards for production freeze
`0820b8d203e2dc534bb883d6fdb4d6d1e0acb11f` are accepted. The seven official
contract/code/vendor files are pinned by SHA256 to that commit. Every raw
receipt must match the receipt used by the completed official result,
including its digest, freeze, domain and samples. Missing/altered inputs,
bad brackets, nonfinite jets or zero gain norms exit with `UNSCORABLE`;
there is no replacement bracket, model or partial score.

Only raw `sum_q,sum_e` and row identities/counts are numerically converted.
The other source/event fields are not used. The numerical `moments` function
is imported from the verified official `archive_channel_split.py` without
calling its guarded main. Source centering, event kernels and the official
source scorer are never executed.

For each N, the same block's unmarked pooled root and first/second p-jets
give the M10/M11 entry/completion gains. They are recomputed for every
original batch omission and paired with the **already saved** official
source response omission. The residual is exactly

```
(C_entry*v_completion - C_completion*v_entry)/hypot(C_entry,C_completion).
```

The central and omission values therefore propagate correlated baseline
and source uncertainty. The output retains four residuals in frozen order
`N85.M10,N85.M11,N340.M10,N340.M11`, their4x4 covariance, and the10-coordinate
joint factors/covariance with official entry/completion/net at both N.
Each N has its own200-row factor `sqrt(199/200)*(loo-mean(loo))`.
Baseline central jets and gain/root omissions are retained as the direct
reproduction path. The official primary decision is copied unchanged;
the reconstructed primary covariance is checked against its saved value.

Only the fixed four-residual Bonferroni normal family is scored, using
`Phi^-1(1-.05/8)`. Either size excluding zero rejects that pure flat-jet
restriction. Otherwise the label is `not_excluded`. No equivalence test,
forced winner, alpha/mixture fit, new mode or derivative scan is present.
These are correlated interpretations of the same fresh block, not extra
independent evidence. Implementation checks use synthetic jets and frozen
code/contract bytes only; no production values or old gain table are needed.
