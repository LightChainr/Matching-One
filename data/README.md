# Canonical Literature Data

This directory contains source data used by Matching One analyses. Numerical values are evidence attached to a specific estimator, geometry, method, and source; they are not definitions of the infinite-lattice threshold.

## Threshold provenance

`literature_threshold_sources.json` is the canonical manifest for the square-site threshold literature currently used by the project. It records:

- exact citation and DOI/arXiv identifiers;
- estimator and geometry definitions;
- finite-size coverage;
- whether the primary source has actually been checked;
- machine-readable data paths and SHA-256 hashes where transcription is complete;
- quoted estimates without collapsing incompatible methods into one synthetic confidence interval;
- explicit pending transcription tasks.

Current machine-readable tables:

- `jacobsen_2015_square_site_cylinder.csv` — periodic-cylinder eigenvalue-identity sequence, `n=1..21`;
- `mertens_2022_square_site_estimators.csv` — exact `p_med` (`n=1..24`) and `p_cell` (`n=2..24`) finite-square estimators.

The Yang–Zhou 2024 publisher abstract is sufficient to verify their quoted corrected estimate, but its full `h(n)` and cylindrical `p_c(n)` tables are intentionally left pending until they can be transcribed from a primary full-text source. Likewise, a numerical value attributed by secondary compilations to the 2024 Jacobsen Reply is not canonical until the primary full text is checked.

## Rules for adding data

1. Preserve printed decimal strings. Do not parse through binary floating point before committing them.
2. Record the exact table/equation/source location and estimator definition.
3. Prefer primary sources. If only secondary evidence is available, mark it as pending rather than silently promoting it.
4. Add a SHA-256 to the manifest and to `SHA256SUMS` for every canonical machine-readable source table. `tests/test_literature_provenance.py` requires the two to agree with each other and with the files.
5. Add or update a regression test that locks row count, endpoint values, and any source-specific invariant.
6. Do not delete a superseded estimate. Change its status and preserve its provenance.
7. Do not average quoted uncertainties from different methods unless a separate statistical model justifies that operation.
8. Check the digits against the primary source before committing, and record that the check happened. A SHA-256 pin detects drift after a file is committed; it cannot detect an error made while transcribing. Both are needed.

## Transcription corrections

Two rows of `jacobsen_2015_square_site_cylinder.csv` did not reproduce the strings printed in Table 2 of arXiv:1507.03027v1 and were corrected on 2026-09-04: `n=4` carried `0` where the source prints `9` in the 20th decimal, and `n=1` carried two trailing zeros beyond the 40 decimals printed. The remaining 19 rows reproduce the printed strings exactly, and `mertens_2022_square_site_estimators.csv` was checked the same way and is clean in all 47 cells.

The `n=4` error stood under a matching file digest for the whole history of the file, because every check the repository ran was a check for *drift* and the file had not drifted. The digit was wrong the day it was written. The comparison that found it was a digit-by-digit read of all 21 rows against the paper.

Nothing downstream moved. Every committed fit over this table trains on `n_min >= 5`, so neither row is in any training window, and the frozen pre-registration `predictions/polynomial_widths_22_24.yaml` regenerates byte-identically apart from its `input_sha256` line. `tests/test_preregister_width_predictions.py` pins the digest of that file with the `input_sha256` line removed, so a future input correction can move the input digest but not a prediction. Full records are in `literature_threshold_sources.json` (`transcription_verification`, `transcription_corrections`, `correction_impact`) and in `analysis/rational_stage_b_quarantine_manifest.json` (`transcription_correction`).

Correcting the table also moves the digest of `literature_threshold_sources.json`, which `analysis/pslq_search_contract.json` pins and every PSLQ result artifact re-verifies, so those artifacts were regenerated. That is the seal working as intended: the numbers in them are unchanged and only the recorded digests moved.

Issue #4 tracks provenance completion. Issue #1 (bounded integer-relation/PSLQ work) is downstream of this dataset and must not select a preferred rounded value independently.
