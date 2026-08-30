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
4. Add a SHA-256 to the manifest for every canonical machine-readable source table.
5. Add or update a regression test that locks row count, endpoint values, and any source-specific invariant.
6. Do not delete a superseded estimate. Change its status and preserve its provenance.
7. Do not average quoted uncertainties from different methods unless a separate statistical model justifies that operation.

Issue #4 tracks provenance completion. Issue #1 (bounded integer-relation/PSLQ work) is downstream of this dataset and must not select a preferred rounded value independently.
