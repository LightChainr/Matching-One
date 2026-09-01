# P537 frozen-witness joint-incidence refinement

This directory is the bounded joint-incidence refinement of the committed
N25 radius-one collar result.  It keeps the original root, Schur projection,
P4 geometry pooling and preferred four-cell filter fixed.

- `REPORT.md` gives the scientific decision and its limits.
- `result.json` is the exact interval score and sector-mixing decomposition.
- `sector-scores.csv` gives every sector determinant and four P4 cell
  midpoints.
- `joint-index.json` records the canonical 12-port identities, terminal
  incidence, `z -> source` roles and exact cell support.
- `schur-aggregates.csv` contains the integer sufficient statistics.  Its
  sector sum reproduces the parent coarse aggregate exactly.
- `RUN-RECEIPT.json` records commands, runtimes, hashes and the base commit.

The paired columns come from the same counterfactual fibre.  For each source
site `y`, the producer freezes the background with `x=y=z=0`, records the
ordered `x[N,E,S,W],y[N,E,S,W],z[N,E,S,W]` component identity at `z=0` and
after `z=1`, and only then evaluates the actual `y=0/1` columns.  No absent
sentinel is matched retrospectively to an unrelated present partition.
