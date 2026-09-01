# P537 exact radius-one finite collar

## Decision

The first frozen finite-collar minor is strictly positive.  The computation
stopped after that single minor, as required:

\[
\det L\in
\left[
\frac{2690418846144177725678447286445842349}{10^{50}},
\frac{2690418846144177725678447286445842359}{10^{50}}
\right],
\]

with midpoint `2.6904188461441777e-14`.  The interval excludes zero.

Rows are the encoded rank births `0->1` and `1->2`; columns are
`axial2:absent` and `axial2:present`.  Both outer attachment bits are one,
all sixteen corner words are summed, and `local_source_contact_mask=0`.

## Exact collar state

The collar is the injective nine-site `B_inf(z,1)`.  In the fixed frame
`x=West(z)`, alternation forces `arm_mask=5` and labels the arms cyclically
as `B_N,W_E,B_S,W_W`.  Bits `NE,SE,SW,NW` record the complete assignment of
the diagonal boundary sites before any outer reconnection.

The two global connectivity bits have only the semantics of outer
attachments:

- `J_B`: the two occupied collar arms reconnect outside;
- `J_W`: the two vacant matching separators reconnect in the off-`z` graph.

Every alternating fibre was checked against

```text
rank1-rank0 = J_B+J_W-1.
```

There were zero identity failures and zero `J_B=J_W=0` fibres in either
geometry.  This separates the finite landing labels from the global
reconnection that caused the previous empty-cell obstruction.

## Computation

The exact producer enumerated `2^23` off-`z` backgrounds and `192,937,984`
source-pair fibres per geometry.  Axis took 22.21 seconds and tilted 17.35
seconds when run concurrently.  Raw CSVs are ignored; committed outputs are
the 386-row Schur aggregate, compact collar counts, interval score, preferred
minor, and summary.

Global `R`, means, and all six C4 source-component beta coefficients were
computed with the existing exact scorer.  Axis and tilted laws were scored
separately and combined only in the final P4 Schur matrix.  The absent source
column was retained with its Schur term.

## Interpretation boundary

This is an exact rank-two certificate for the frozen N25 radius-one collar.
It rejects rank one for that finite ordinary block.  It does not assert a
macroscopic arm probability, scaling law, or asymptotic cancellation result.

## Reproduction

```sh
clang++ -std=c++17 -O3 -Wall -Wextra -pedantic \
  experiments/p537-finite-collar-20260901/siteflip_collar_exact.cpp \
  -o /tmp/siteflip_collar_exact

/tmp/siteflip_collar_exact 5 0 \
  experiments/p537-landing-matrix-preflight-20260901/kernel.tsv \
  results/p537-finite-collar/axis.csv
/tmp/siteflip_collar_exact 4 3 \
  experiments/p537-landing-matrix-preflight-20260901/kernel.tsv \
  results/p537-finite-collar/tilted.csv

python3 experiments/p537-finite-collar-20260901/aggregate_collar_for_schur.py \
  --axis results/p537-finite-collar/axis.csv \
  --tilted results/p537-finite-collar/tilted.csv \
  --output results/p537-finite-collar/schur-aggregates.csv

python3 experiments/p537-siteflip-landing-20260901/score_siteflip_schur.py \
  --aggregates results/p537-finite-collar/schur-aggregates.csv \
  --baseline-axis experiments/p537-landing-matrix-preflight-20260901/baseline-axis.csv \
  --baseline-tilted experiments/p537-landing-matrix-preflight-20260901/baseline-tilted.csv \
  --baseline-root experiments/p537-landing-matrix-preflight-20260901/baseline-root.json \
  --a-raw-denominator 16 --z-orbit-multiplicity 4 --first-nonzero-only \
  --output results/p537-finite-collar/result.json
```
