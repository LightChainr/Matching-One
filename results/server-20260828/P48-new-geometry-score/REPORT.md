# P48 prospective new-geometry four-channel score

## Why this score exists

The N=185/265 full curves are genuinely new target geometries. They can therefore score the original intrinsic-center P48 pure-power amplitudes frozen from N=65,85,130 without using the targets to fit a parameter.

This score is separate from the Issue #43 fixed-coordinate wrapping-channel contract. PR #134 already records that the original Issue #43 even prediction used an `either/even` source against a `cross/even` target and supplies the exact channel-map erratum. This report does not reopen or replace that erratum.

## Frozen source

Source amplitudes and uncertainties are read directly from

`results/server-20260828/P48-retrospective/summary.json`.

They were trained on `N=65,85,130` before the N=185/265 target values existed.

## Target and covariance

Target values are the independent 500M-per-size intrinsic-center analyses

- `P43-heldout-fullcurve-500m/analysis/n185.p48.json`;
- `P43-heldout-fullcurve-500m/analysis/n265.p48.json`.

N=185 and N=265 use disjoint counter domains. For each channel the residual covariance is

```text
C = diag(target sampling variance) + Var(frozen source amplitude) * 11^T.
```

No target parameter is fit.

## Result

| channel | frozen pure law | N=185 scaled | N=265 scaled | frozen chi2 / 2 | zero chi2 / 2 |
|:---|:---|---:|---:|---:|---:|
| `P4[S]` | `N^-1` | -0.009549 +/- 0.001270 | -0.010808 +/- 0.001444 | **1.139** | 112.540 |
| `P4[D]` | `N^-13/8` | +0.280620 +/- 0.070214 | +0.319722 +/- 0.087227 | **0.281** | 29.408 |
| `P4[D']` | `N^-5/8` | -0.023707 +/- 0.004402 | -0.022663 +/- 0.004111 | **0.088** | 59.393 |
| `P4[S']` | `N^-5/4` | +2.579708 +/- 0.100202 | +2.858439 +/- 0.115193 | **52.716** | 1278.555 |

The N=185/265 marginal residuals from the frozen `P4[S]` amplitude are only `-0.57` and `-1.04` combined standard errors. `D` and `D'` are even closer to their frozen amplitudes.

`S'` is the unique clear pure-law failure among these four intrinsic-center channels on the prospective new geometries.

## Interpretation

The prospective data strengthen the empirical P48 parity pattern:

```text
P4[S]   pure law survives
P4[D]   pure law survives
P4[D']  pure law survives
P4[S']  pure law fails; finite-size correction required
```

The already-frozen q=2 and log/Jordan correction models for `S'` remain the first alternatives. This score does not identify which correction mechanism is correct.

The Issue #43 channel-map erratum remains a different statement: its original registered even score failed because source and target wrapping channels were mismatched, while the exact `either -> cross` correction gives a compatible lower-status protocol-repair diagnostic. Neither result should be conflated with this intrinsic-center P48 score.

## Reproduction

```bash
python3 scripts/score_p48_new_geometry_channels.py \
  --output results/server-20260828/P48-new-geometry-score/score.json
```
