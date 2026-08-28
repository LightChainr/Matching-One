# P48 S-prime: minimal relative 1/N correction

## Status

Prospective correction model frozen before N=185/265 target scoring.

The original pure-power package in Issue #48 remains the primary test and must be
scored first.  This note defines the **single predeclared correction alternative**
for the channel that failed retrospectively, `P4[S']`.

## Motivation

The retrospective P33 score supports the matching-parity pattern but finds upward
drift in the scaled quantity

```text
Z_Sp(N) = N^(5/4) P4[S'].
```

Do not immediately replace `x=21/4` by a free/lower exponent.  The repository now
has independent evidence that relative `L^-2 = N^-1` finite-size corrections occur
in center-slope quantities (P49 / PR #83), and the observed operator algebra also
permits a matching-even scalar correction multiplying the matching-odd H4 field.

Therefore freeze the minimal model

```text
P4[S'] = N^(-5/4) * (A + B/N),
```

before trying logs, mixed powers or free exponents.

## Source fit

Use only the already-known retrospective training sizes

```text
N = 65, 85, 130
```

and the synchronized cross-size covariance committed in
`results/server-20260828/P48-retrospective/cross_size_covariance.csv`.

The covariance-aware GLS fit to the scaled variable `Z_Sp=N^(5/4)P4[S']` gives

```text
A =  3.152857099000578
B = -88.03321330713617
SE(A) ~= 0.393648631743133
SE(B) ~= 28.105453649642055
```

These values are derived from data that predate the model and are therefore only
source calibration, not confirmatory evidence.

The model retrospectively predicts the old held-out points 145/170 without the
systematic upward failure of the pure `N^-5/4` law; this observation motivated the
freeze but must not be counted as prospective support.

## Frozen N=185/265 source predictions

Propagating **source-fit covariance only** gives:

```text
N=185:
  scaled Z_Sp = 2.6770018919349767 +/- 0.24641162774047082
  P4[S']      = 0.003923593601402982 +/- 0.00036115741599830305

N=265:
  scaled Z_Sp = 2.8206562940679882 +/- 0.2903523906239357
  P4[S']      = 0.0026381095543789495 +/- 0.0002715614154949255
```

Target sampling uncertainty and its covariance with other channels must be added
by the prospective scorer.

## Scoring order

1. Score the original Issue #48 four pure powers exactly as frozen.
2. Score the no-fit operator-ratio competitors in
   `predictions/derivative_operator_ratio_test_20260828.yaml`.
3. For `P4[S']` only, score this fixed `N^-5/4(A+B/N)` alternative using the
   source covariance above plus target covariance.
4. Only then inspect a logarithmic correction, mixed fields, or a free exponent.

## Interpretation

- Pure `N^-5/4` passes prospectively: the old drift was statistical/preasymptotic.
- Pure fails but this `1/N` correction passes: retain `x=21/4`; the drift is
  consistent with an ordinary finite-size correction rather than a new leading field.
- A lower-x operator competitor wins both central `D` and derivative `S'`: revise
  the operator assignment rather than hiding the failure in a correction term.
- Neither correction nor one-field competitors work: require a genuine multi-field
  or logarithmic/Jordan analysis.

This model does not identify the `1/N` correction operator.  It is deliberately the
smallest correction structure already motivated elsewhere in the repository.
