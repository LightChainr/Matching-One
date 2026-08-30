# Issue #154 Phase-E mixed local-singlet pilot

The pilot stores `B`, `B^2`, `I0*B`, and `I2*B` on the same stream,
so `J_top=Cov(I2-I0,B)` and `J_bulk=Cov(I2+I0,B)` are directly scoreable.

| N | P4[B] | P4[J_top] | P4[J_bulk] | P4[Var(B)] |
|---:|---:|---:|---:|---:|
| 65 | -6.601562e-05 +/- 0.00014 (z=-0.48) | -0.0002954513 +/- 0.00014 (z=-2.08) | 7.824877e-05 +/- 0.00011 (z=0.73) | -2.043047e-05 +/- 7.6e-06 (z=-2.68) |
| 130 | -0.0001335829 +/- 0.00011 (z=-1.19) | -0.0001454797 +/- 0.00011 (z=-1.39) | -2.638123e-05 +/- 6.5e-05 (z=-0.41) | 3.307401e-06 +/- 3.8e-06 (z=0.88) |

## Common-plane score

| candidate | scale N130/N65 | chi2 / 2 df | p |
|---|---:|---:|---:|
| A/E/C | 0.32043 | 2.274 | 0.3208 |
| A/E/J_bulk | 0.098113 | 0.29264 | 0.8639 |
| A/E/J_top | 0.2265 | 0.93489 | 0.6266 |
| A/E/B | 0.18857 | 1.1637 | 0.5589 |

## Decision

`mixed_local_plane_not_selected_at_20k`. The A/E/J_bulk improvement over A/E/C is 1.9814; resolved at both sizes: False.

This is a finite-volume mixed-observation decision. It does not name `B` as the
continuum thermal energy or turn a surviving plane into a field identity.

The stronger visual fit of `A/E/J_bulk` is not promoted: its improvement is below
the frozen threshold, `J_bulk` is unresolved at both sizes, and its point estimate
changes sign. `J_top` is negative at both sizes but reaches only `2.08` and `1.39`
standard errors. The declared local-singlet definition therefore stops at 20k.

## Provenance

- Freeze/authorization: `0578105` / `d165b1e`.
- N65: `DevEnvC_ZyTrST`, 20,000 replicas, 100 batches, 0.0510 s.
- N130: `DevEnvC_XPk2PZ`, 20,000 replicas, 100 batches, 0.0748 s.
- Remote and local SHA-256 hashes of all four raw files agree; see `SHA256SUMS`.
