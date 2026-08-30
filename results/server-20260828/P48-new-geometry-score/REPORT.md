# P48 prospective new-geometry four-channel score

## Why this score exists

The N=185/265 full curves in P43 are genuinely new target geometries. They can therefore score the original intrinsic-center P48 pure-power amplitudes frozen from N=65,85,130 without using the targets to fit a parameter.

This score also fixes an important observable distinction:

- P48 `P4[S]` is the intrinsic-center thermal-even projector;
- P43 `DeltaS` is the separate fixed-coordinate P31 `either/even` observable.

The failure of the P43 positive `DeltaS` forecast must not be described as a failure of the P48 `P4[S] ~ N^-1` law.

## Frozen source

The source amplitudes and uncertainties are read directly from

`results/server-20260828/P48-retrospective/summary.json`.

They were trained on `N=65,85,130` before the N=185/265 target values existed.

## Target

The target values are the independent 500M-per-size intrinsic-center analyses

- `analysis/n185.p48.json`;
- `analysis/n265.p48.json`

from `P43-heldout-fullcurve-500m`.

N=185 and N=265 use disjoint counter domains. The residual covariance is therefore

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

The N=185/265 marginal residuals from the frozen P48 `S` amplitude are only `-0.57` and `-1.04` combined standard errors. `D` and `D'` are even closer to their frozen amplitudes.

`S'` is the unique clear pure-law failure among these four channels on the prospective new geometries.

## Interpretation

The prospective data sharpen rather than destroy the P48 parity picture:

```text
P4[S]   pure law survives
P4[D]   pure law survives
P4[D']  pure law survives
P4[S']  pure law fails; finite-size correction required
```

This is compatible with the already-frozen q=2/log correction program for `S'`. It does not identify which correction mechanism is correct.

Separately, the P43 fixed-coordinate `either/even DeltaS` positive forecast fails strongly. The two observables should be analyzed as distinct finite-size constructions until a derivation relates their corrections.

## Reproduction

```bash
python3 scripts/score_p48_new_geometry_channels.py \
  --output results/server-20260828/P48-new-geometry-score/score.json
```
