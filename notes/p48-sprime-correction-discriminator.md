# Prospective correction discriminator for the P48 `S'` channel

Status: experiment-design note. The target `N=185,265` full-curve production was already in progress when this note was committed, but its target values had not been revealed. The original Issue #43 `DeltaM/DeltaS` score remains primary and must be reported first.

## Why this channel matters

The matching/thermal parity rule predicts that a matching-odd field with scaling dimension `x=21/4` contributes at the intrinsic center as

\[
P_4[S'] \sim L^{2-x+y_t}
          = L^{-5/2}
          = N^{-5/4},
\qquad y_t=3/4.
\]

The exponent arithmetic is therefore not the problem in the retrospective P48 result. The problem is that

\[
Y_N=N^{5/4}P_4[S']
\]

rises across the old sizes and the pure-constant law fails its held-out conjunction, while the parity signal itself is clearly nonzero.

## Two pre-target explanations

Only `N=65,85,130` are used to fix coefficients. `N=145,170` remain excluded from the fit.

### A. Rank-2 / logarithmic correction

Because the percolation thermal sector is logarithmic at `c=0`, the first correction tested prospectively is

\[
Y_N=A+B\log N.
\]

GLS with the committed P48 covariance gives

\[
A=-2.7539493884634396,\qquad
B=1.0923296227893287.
\]

Frozen targets:

\[
Y_{185}=2.9483999207704414\pm0.3231005705113089,
\]

\[
Y_{265}=3.340954787622457\pm0.4431673354383616.
\]

This is a phenomenological rank-2/log test. Passing it would not by itself identify a particular LCFT Jordan partner.

### B. Ordinary analytic correction

A relative `L^-2=N^-1` correction gives

\[
Y_N=A+B/N.
\]

The same GLS training fit gives

\[
A=3.1528570990005798,\qquad
B=-88.03321330713665.
\]

Frozen targets:

\[
Y_{185}=2.6770018919349763\pm0.2464116277404708,
\]

\[
Y_{265}=2.8206562940679887\pm0.2903523906239357.
\]

## Baseline

The original P48 pure-power law remains the first score:

\[
Y_N=A_0,\qquad
A_0=1.9434247576878727\pm0.0766048795577632.
\]

Thus the new full-curve targets separate the three predictions without any extra simulation.

## Scoring discipline

1. Report Issue #43's original `DeltaM/DeltaS` endpoints first.
2. Score the original P48 pure `N^-5/4` law.
3. Score the frozen log correction.
4. Score the frozen `1/N` analytic correction.
5. Only then inspect free exponents or richer mixtures.

Target covariance must be included. Source-prediction covariance is stored explicitly in `predictions/p48_sprime_correction_20260828.yaml`.

## Interpretation boundary

- pure law passes: the old P48 drift was a low-stat fluctuation;
- log beats pure and analytic on fresh targets: evidence for a logarithmic correction in this channel, not yet a Jordan-block identification;
- analytic beats log: ordinary irrelevant/composite contamination is sufficient;
- all fail: revisit the assumed leading field/parity assignment before adding more flexible fits.
