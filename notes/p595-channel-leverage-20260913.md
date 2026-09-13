# #595 — Channel leverage: which threshold-rank readout maximizes angular resolution per sample

**Branch:** `analysis/p595-channel-leverage-20260913`
**Issue:** #595 (Channel leverage)
**Status:** design recommendation delivered; pre-registration controls applied; exact S/D/Sp/Dp obs() channels declared a buy-back.

## Question (pre-registered)

Over the archived threshold-rank histograms in
`results/server-20260829/P205-norm5-conjugate-coalescence/raw/` (zero new
sampling), which readout maximizes `|A4| / se(A4)` at fixed sample count?
`A4` is the angular amplitude (response to `cos 4θ`); the two archived
orientations (first matrix `[15,-10;10,15]`, second `[17,-6;6,17]`) supply the
angular lever. Because the orientation geometry is identical for every readout,
`|A4|/se(A4)` is proportional to the SNR of the orientation difference
`dR = R_first − R_second`; the constant `1/Δcos 4θ` is common, so the *ranking*
and the *amplification relative to M* are exact.

## Readout family scanned (all reconstructable from the archive)

1. `M` — exact matching observable `M(p_ref)`, via the same recurrence as
   `scripts/analyze_threshold_ranks.py`. **Baseline.**
2. `M` on a **grid of `p`** (15 points, 0.50–0.64), not only `p_ref`.  The
   archived integer-period histogram supports `M(p)` as a function of `p`, so
   the "same channel at a grid of p" family is genuinely in the archive.
3. **Higher tail derivatives** `M'(p)`, `M''(p)` of the same histogram.
4. **Histogram-native functionals** (readable surrogates for the S/D/Sp/Dp
   family): mean threshold rank `<K>`, gap `K_minus − K_plus`, variance `var K`.
5. **Linear combinations** solved as a generalized-eigenvalue problem against
   the delete-one jackknife covariance `Σ`: `w* = argmax_w (w' μ μ' w)/(w' Σ w)`.
   Solved, not searched.

**Buy-back (not fabricable here, declared prospectively).** The exact `S`, `D`,
`Sp`, `Dp` `obs()` outputs are emitted by the C++ threshold-rank engine and are
**not** independently reconstructable from the archived threshold-rank
histograms alone. They are therefore a forward buy-back, not invented. This is
exactly the issue's own governance: select a channel prospectively, then buy it
with new samples before it may carry a verdict.

## Governance (GOVERNANCE §2D/§2E) — enforced

- The scan reports **only** `se(A4)`-equivalent (SE of `dR`) and
  **amplification relative to M**. It never reports `delta`, `rho`, a z-score,
  or any Smith-contamination bound for a channel it selected — those are the
  numbers selection already saw.
- The selected combination is declared and then validated on a **held-out half**
  (split-half control, 40 random splits for stability).
- The honest deliverable is a **design recommendation with a stated
  selection-optimism penalty**, not a measurement.

## Results (fixed `n_batches = 100` per archive file)

Best apparent readout across the family: **`varKminus`** (variance of the
threshold rank), with apparent SNR up to 5.35 (325/425 × C_B) and apparent
amplification vs `M` up to **65×** (apparent, 325_C_B).

| cell | best readout | apparent SNR | apparent amp vs M | transport (single) | transport (optimal combo) | optimism penalty |
|---|---|---:|---:|---:|---:|---:|
| 325_C_A | varKminus | 3.04 | 3.10× | 0.83 | 0.57 | 0.43 |
| 325_C_B | varKminus | 5.35 | 65.5× | 0.92 | 0.74 | 0.26 |
| 425_C_A | varKminus | 4.87 | 18.0× | 0.99 | 0.63 | 0.37 |
| 425_C_B | varKminus | 3.49 | 3.98× | 0.89 | 0.37 | 0.63 |

`transport = snr_heldout / snr_selected` averaged over 40 random splits.

## Interpretation (the pre-registered control fires)

- The **single best readable channel transports** (`0.83–1.00` → real, modest
  leverage). So the spread is **not pure noise**: a single channel
  (variance of the threshold rank) carries a transportable angular lever and is
  the honest design recommendation.
- The **data-snooped generalized-eigenvalue combination does NOT transport**
  (`0.37–0.74`). Solving the GE problem on the same finite sample finds leverage
  that is substantially selection optimism. This is exactly the warning the
  issue pre-registered: do not trust the solved combination as a measurement.

## Conclusion

- **Forward recommendation (with selection-optimism penalty stated):** score the
  variance-of-threshold-rank channel (`varKminus`) prospectively with new
  samples; it is the readable channel with the best *transportable* leverage.
  Do **not** deploy the GE-optimized combination without held-out buy-back — its
  apparent amplification is mostly optimism.
- **Limitation / buy-back:** the issue's specific `Sp`-vs-`M` 4.5× claim and the
  exact `S/D/Sp/Dp` obs() channels cannot be confirmed or refuted from this
  archive (they are not in the threshold-rank histograms). That comparison
  remains a prospective buy-back. The issue therefore does **not** close negative
  on "the spread is noise" (a single channel does transport); it yields a design
  recommendation plus a declared data gap.
- No contamination bound, `delta`, `rho`, or z-score is reported for any selected
  channel, per governance.
