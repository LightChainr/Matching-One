# Prospective correction test for the P48 `P4[S']` drift

## Status

This note freezes a design response to a **retrospective** anomaly before the prospective `N=185,265` full-curve results are revealed. It is not evidence for a logarithmic correction.

P48 predicted, for the matching-odd thermal-family spin-4 hypothesis,

\[
P_4[S'] \sim N^{-5/4}.
\]

The parity pattern survived retrospectively, but the scaled quantity

\[
y_N = N^{5/4} P_4[S']
\]

rose across the old sizes and the pure constant-amplitude law failed its declared held-out conjunction. The existing five-size design has weak power to distinguish radial correction families, so this note freezes only two physically motivated competitors.

## Training boundary

Use only the original P48 training sizes

```text
N = 65, 85, 130.
```

The already-known `N=145,170` values are deliberately excluded from all coefficient fits. The new target sizes are `N=185,265` from Issue #43.

The full synchronized-jackknife covariance of the three training values is used in GLS.

## Frozen competitors

### Rank-2 logarithmic/Jordan form

\[
y_N=A+B\log(N/100).
\]

Training-only GLS gives

```text
A = 2.2764144236775126
B = 1.092329622789324
```

and prospectively predicts

```text
N=185: 2.9483999207703566 +/- 0.3231005705113109  [source-fit uncertainty]
N=265: 3.3409547876223700 +/- 0.4431673354383646
```

A positive logarithmic slope is compatible with the possibility that a `c=0` thermal logarithmic multiplet contaminates the nominal `N^-5/4` law. This is only a phenomenological rank-2 form; it is not an LCFT derivation.

### Ordinary relative `L^-2` correction

Because `L^2=N`, a relative `L^-2` correction becomes `N^-1` after the leading `N^-5/4` factor is removed:

\[
y_N=A+B/N.
\]

Training-only GLS gives

```text
A = 3.1528570990005758
B = -88.03321330713636
```

and predicts

```text
N=185: 2.6770018919349736 +/- 0.24641162774047082
N=265: 2.8206562940679856 +/- 0.2903523906239357
```

## Scoring order

1. retain the original pure-power P48 law as the baseline;
2. score the frozen log/Jordan competitor;
3. score the frozen ordinary `N^-1` correction;
4. compare zero effect;
5. only then consider free exponents or additional correction terms.

Target covariance must be combined with the source-fit prediction covariance committed in `predictions/p48_sprime_correction_competitors_20260828.yaml`.

## Interpretation

A prospective preference for the logarithmic form would justify a dedicated logarithmic/Jordan study, but would not identify a specific logarithmic partner. Preference for `N^-1` would instead favor a conventional analytic/irrelevant correction. If both fail, the parity model can remain viable while the naive one-field derivative asymptotics are rejected.

This discriminator is deliberately secondary to Issue #43's original frozen central-amplitude predictions and to the norm-5 harmonic test. It must not be used to retune those primary predictions.
