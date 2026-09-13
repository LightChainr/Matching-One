# p593 Phase A — pinning down the definition of `r_linear` (2026-09-13)

Issue: #593. Context: the #599 probe measured block-Hankel / McMillan orders of
9/8, 4/5, 2/2 at widths 4..6 for the P398 projected-memory object, while the
rank tables in #593/#588 report `r_linear` = 10, 26, 72.  Phase A asks which
candidate definition of `r_linear` reproduces 10/26/72, and whether the probe
numbers contradict the reported ones.

## Method

All quantities computed at `eta = 0` from the repository's own exact generator
(`p398_intervention_transport.Generator`, widths 4..8), with integer-valued
readouts and exact elimination over the Mersenne prime 2147483647 (the
repository's `RANK_PRIME`; no pivot tolerance anywhere).  Script:
`scripts/p593_phaseA_rlinear_candidates.py`; raw output:
`results/p593-phaseA-rlinear/raw/candidates.json`.

## Result table

| quantity | w=4 | w=5 | w=6 | w=7 | w=8 |
|---|---|---|---|---|---|
| states | 14 | 42 | 132 | 429 | 1430 |
| **C1: joint Krylov `dim span{G0^k f}`, D0 seeds** | **10** | **26** | **72** | **218** | **689** |
| C1 variant: seeds + constant `1` | 10 | 26 | 72 | 218 | 689 |
| C1 variant: D2 dictionary seeds (8 readouts) | 10 | 42 | 76 | 415 | 750 |
| C1 variant: `k >= 1` only (drop the seeds) | 9 | 25 | 71 | 217 | 688 |
| C1 variant: Krylov under `H` instead of `G0` | 10 | 22 | 76 | 214 | 750 |
| C1 variant: floating-point Krylov | 10 | 27 | 81 | – | – |
| C2: coarsest exact lumping vs `G0` only (block count) | 10 | 26 | 76 | 232 | 750 |
| C2: coarsest exact lumping vs `G0` and `H` jointly | 10 | 26 | 76 | 232 | 750 |
| C2: distinct D0 readout value-tuples | 8 | 12 | 18 | 24 | 32 |
| C3: rank `G0` mod p | 13 | 41 | 131 | 428 | 1429 |
| C3: rank `H` mod p | 13 | 32 | 130 | 392 | 1429 |
| C3: rank `J` / rank `D` mod p | 13 | 41 | 131 | 428 | 1429 |
| C3: distinct rows of `G0` | 14 | 42 | 132 | 429 | 1430 |
| repository `observable_reachable_dimension` (cross-check) | 10 | 26 | 72 | 218 | not run (budget) |

## Conclusion

**Candidate 1 is pinned, exactly and only.**  `r_linear` is the joint Krylov
dimension

```text
r_linear = dim span{ G0^k f : k >= 0, f in {blocks, singletons, wrap} }
```

computed by exact elimination over the prime 2147483647.  It reproduces
10/26/72 at widths 4..6 and 218 at width 7, and agrees digit-for-digit with
`rank_notions.D0_additive_local_counts.r_linear` stored in
`results/p398-intervention-transport/latest.json`.

- **Candidate 2 (lumping degrees of freedom) is falsified**: the coarsest
  exact strong lumping gives 76/232/750 at widths 6..8 — that is `r_positive`,
  not `r_linear`.  (It coincides with `r_linear` only at widths 4 and 5.)
- **Candidate 3 (pencil ranks / combinatorial structures) is falsified**:
  every pencil rank is `size - 1` or near it, and no combinatorial count
  matches.
- The #599 probe numbers 9/8, 4/5, 2/2 are **not** a competing measurement of
  `r_linear`; they are block-Hankel / McMillan orders of the *projected
  memory kernel* on the frozen rank-6 span — a strictly smaller and different
  object (its `k >= 1`-flavoured count at width 4 is 9, which is the probe's
  block-Hankel value, but the sequence diverges immediately at width 5).

## New value produced

- **True `r_linear` at width 8 = 689** (Krylov depth 621).  The repository
  entry was budget-truncated (`>= 150`, `LINEAR_RANK_BUDGET_ABOVE`); the
  exact elimination at full budget is now affordable and the number is a
  certified lower bound on the rational rank.

## Notes

- The floating-point Krylov variant overestimates (27 at width 5, 81 at
  width 6 against true 26/72), reproducing the failure mode the ticket
  documents for float Gram-Schmidt.
- `dim span{G0^k f}` with the constant function added to the seeds is
  unchanged at every width: `1` is already in the Krylov span of the D0
  readouts.
- Full-budget `r_linear` at widths 9 and 10 is computed in Phase B
  (`results/p593-width9-10-memory-degree/`).

Full Matching-One repository CI has not been run for this commit.
