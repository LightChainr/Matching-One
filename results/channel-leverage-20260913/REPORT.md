# #595 Channel-leverage scan — REPORT

**Issue:** #595 (Channel leverage: which observable maximizes angular resolution per sample)
**Branch:** `analysis/p595-channel-leverage-20260913`
**Inputs:** `results/server-20260829/P205-norm5-conjugate-coalescence/raw/` (archived threshold-rank histograms, zero new sampling)
**Method script:** `scripts/channel_leverage_scan.py`

## What was scanned

The readable readout family reconstructable from the archived threshold-rank
histograms, all at zero new sampling cost:

1. `M(p_ref)` — exact matching observable (baseline).
2. `M(p)` on a 15-point grid of `p` (0.50–0.64).
3. Tail derivatives `M'(p)`, `M''(p)`.
4. Histogram-native functionals: mean threshold rank, gap, variance `varKminus`.
5. Linear combinations solved as a generalized-eigenvalue problem against the
   delete-one jackknife covariance.

The exact `S`/`D`/`Sp`/`Dp` `obs()` channels are **not** reconstructable from the
archived threshold-rank histograms and are declared a prospective buy-back (not
fabricated).

## Governance compliance

- Only `se(A4)`-equivalent (SE of orientation difference `dR`) and amplification
  relative to `M` are reported. No `delta`, `rho`, z-score, or contamination bound
  for any selected channel.
- Selected combination validated on held-out batches (split-half control, 40
  random splits).
- Deliverable is a design recommendation with a stated selection-optimism penalty.

## Headline numbers

- Best apparent readout: **`varKminus`** (variance of threshold rank), apparent
  SNR up to 5.35, apparent amplification vs `M` up to **65×** (325_C_B).
- **Transport (snr_heldout / snr_selected), 40 splits:**
  - single best readout: **0.83–1.00** → real, transportable leverage.
  - GE-optimal combination: **0.37–0.74** → dominated by selection optimism.

| cell | best readout | SNR | amp vs M | transport (single) | transport (combo) | optimism |
|---|---|---:|---:|---:|---:|---:|
| 325_C_A | varKminus | 3.04 | 3.10× | 0.83 | 0.57 | 0.43 |
| 325_C_B | varKminus | 5.35 | 65.5× | 0.92 | 0.74 | 0.26 |
| 425_C_A | varKminus | 4.87 | 18.0× | 0.99 | 0.63 | 0.37 |
| 425_C_B | varKminus | 3.49 | 3.98× | 0.89 | 0.37 | 0.63 |

## Conclusion

The spread is **not pure noise** — the single readable channel `varKminus` carries
a transportable angular lever (transport ~0.83–1.00). But the **data-snooped
GE-optimal combination does not transport** (0.37–0.74): its apparent amplification
is mostly selection optimism, exactly the pre-registered risk. Recommendation:
score `varKminus` prospectively with new samples; do not ship the optimized
combination without held-out buy-back. The exact `Sp`-vs-`M` 4.5× claim is a
declared buy-back (channels absent from the archive).

Full Matching-One repository CI has not been run for this commit.
