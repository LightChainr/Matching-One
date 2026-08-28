# MERGE_BLOCKER_PRECISION: set mpmath dps before CSV parse

Queue item `MERGE_BLOCKER_PRECISION` from
`experiments/server_compute_queue_overrides_20260828.yaml` on
`research/second-wave-20260828` (`25dd9c9`). Target: PR #21. This branch
starts from `origin/server/huawei-analysis-20260828` at
`107ddc948c8ae8e5c8a30549c9b123ac9be091b7`.

Wall time: 0.115 s (18:55:55 Asia/Shanghai, 2026-08-28). Single-thread CPU,
no GPU. The 54 Huawei grid jobs were **not** recomputed.

## Bug

`scripts/summarize_finite_size_grid.py` called `load_observations()` before
`mp.mp.dps` was set to `chosen_dps`. `mp.mpf(row["value"])` therefore ran at
the default mpmath precision (15 decimal digits). Lost digits cannot be
recovered by raising `dps` later.

The original decimal CSV used by the Huawei 54-grid summary is
`data/jacobsen_2015_square_site_cylinder.csv`
(sha256 `d73191d771e5ca7d95a6a2f6dc36c444f40eacbe7819804c442497222d562210`).
Widths 1--21 are Jacobsen 2015 Table 2 cylinder values; several entries have
~40 decimal digits, which default dps truncates. Example, width 19:

```text
CSV:              0.5927438107312915517933469085441350226515
old JSON actual:  0.59274381073129156494161406953935511
new JSON actual:  0.59274381073129155179334690854413502
```

The old `actual` is exactly `mp.mpf(csv)` at dps=15.

## Fix

1. Load grid payloads and read `dps` from JSON fields.
2. Set `mp.mp.dps = chosen_dps` (160 for this campaign).
3. Only then parse the original decimal CSV strings into `mp.mpf`.
4. Refit the already-selected configuration on those high-dps observations.

`load_decimal_rows()` keeps the original strings; conversion happens in
`observations_from_decimals()` after precision is set.

## Rerun

Same inputs as the Huawei summary:

```text
CSV:       data/jacobsen_2015_square_site_cylinder.csv
raw grid:  results/server-20260828/issue-5-grid/raw   (54 JSON files)
final-tail: 3  (knowledge cutoff n=18; withheld 19, 20, 21)
```

Selected configuration is unchanged: powers `4,6,8,10,12`, `n_min=8`,
holdout 2, dps 160.

### Old vs new intercept

```text
old (default-dps parse):  0.59274605094603179333564156705784225
new (dps=160 parse):      0.59274605094603206266439366806726549
shift (new - old):        +2.6932875210100942324e-16
```

Final-tail RMSE moved from `1.3775763431649e-11` to `1.3775861250986e-11`.
Maximum absolute error moved from `1.9374735274555e-11` to
`1.9374870639209e-11`. Signed errors on 19--21 remain positive and
monotone. The shift is about `2.7e-16` and does **not** change the current
scientific conclusion: this is still an out-of-sample finite-width
prediction with same-sign drift, not a `1e-11` claim on the infinite-lattice
threshold.

The original Huawei summary is preserved as superseded provenance:

- original: `results/server-20260828/issue-5-summary.json`
- copy: `results/server-20260828/MERGE_BLOCKER_PRECISION/issue-5-summary.superseded.json`

The corrected summary is `issue-5-summary.json` in this directory.

## Regression test

`tests/test_finite_size_toolchain.py::HighPrecisionCsvParseTests.test_csv_mpf_matches_value_parsed_after_high_dps`

Uses CSV value `0.31415926535897932384626433832795028841971693993751`
(50 non-binary-friendly decimal digits) with a grid JSON at dps=80. The
summarizer's stored `actual` matches `mp.nstr(mp.mpf(value), 35)` constructed
**after** setting dps=80, and differs from the dps=15 truncation.

```text
Ran 4 tests in 0.006s
OK
test_csv_mpf_matches_value_parsed_after_high_dps ... PASS
```

## Scope

C05 files and the sibling worktree `/workspace/Matching-One-c05` were not
touched. This branch is not merged and no pull request was opened.
