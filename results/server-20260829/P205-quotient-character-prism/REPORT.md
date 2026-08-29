# P205 quotient-character prism: frozen 12M reveal

## Outcome

The quotient-character prism selected H4 cleanly.  All three measured
equal-area contrasts were positive, exactly the H4 sign code and incompatible
with the coordinate-specific sign flips frozen for H8 and H12.

| frozen model | chi-square / df | p-value | delta chi-square from H4 |
|:---|---:|---:|---:|
| H4 | 1.087585 / 2 | 0.580542 | 0 |
| H8 | 753.259290 / 2 | 2.703e-164 | 752.171705 |
| H12 | 26.417646 / 2 | 1.834e-6 | 25.330061 |

H4 is not merely the least bad row.  Its three standardized residuals are
`+0.238, -0.934, -0.398`, so a single frozen `N^-13/8` H4 amplitude describes
all three cross-Smith quotient contrasts without tension.

The fitted H4 amplitude is

```text
A4 = 0.8039184543 +/- 0.0131464871.
```

The planning-only value `0.7885`, recorded before these targets were acquired,
is only 1.17 standard errors away.  It set the sample scale but was not fixed in
the score.

## The mechanism selector fired by sign

The raw fixed-p contrasts are

| N | quotient pair | Delta M first-minus-second | SE | frozen signs H4/H8/H12 |
|---:|:---|---:|---:|:---:|
| 25 | `(5,0)/(4,3)` | +0.0079591560 | 0.0001334756 | `+ / + / +` |
| 50 | `(7,1)/(5,5)` | +0.0024043631 | 0.0001774468 | `+ / - / +` |
| 125 | `(11,2)/(10,5)` | +0.0002374592 | 0.0002193785 | `+ / + / -` |

N50 supplies the H8 veto: after its best global amplitude fit, the N50 H8
residual is `+26.46` standard errors.  N125 supplies the H12 veto: its H12
residual is `+5.04` standard errors.  The separation therefore does not depend
on selecting a delicate correction exponent or on a marginal likelihood
preference.  It is the exact angular character code doing the work.

The observed signed ratios relative to N25 are approximately `+0.3021` at N50
and `+0.02983` at N125.  The frozen H4 ratios were `+0.3242` and `+0.04096`;
their deviations are already included in the covariance-aware H4 score above.

## What changed scientifically

The completed N325/N425 coalescence block established that the low-dimensional
angular transfer survived a same-parent cyclic-to-noncyclic quotient control,
but it did not choose among H4/H8/H12.  Its H4-H8 and H12-H8 chi-square gaps
were only `2.253` and `0.242`.

The prism changes that conclusion.  Across three different equal-area Smith
changes,

- `(5,5) -> (1,25)`;
- `(1,50) -> (5,10)`;
- `(1,125) -> (5,25)`;

one common H4 amplitude survives while the two higher-harmonic character lines
fail.  Within the frozen three-model mechanism competition, this is decisive
selection of H4 rather than generic evidence for an unspecified low-rank
angular law.

The small sizes were intentional: they maximize exact angular discrimination
per CPU.  Their successful one-amplitude collapse means the predeclared
N338/N400/N500 and N400/N450/N500 bridges are not needed as rescue runs.  They
remain useful only if the project wants a separate finite-size transport study,
not to establish which of H4/H8/H12 won this acquisition.

## Execution

The acquisition ran on Huawei DevEnv `DevEnvC_ZyTrST`, ID
`f415a4bcbd9a438b85f5f29e4a507ea4`, using two 8-thread lanes.  N125 occupied
one lane while N50 then N25 occupied the other.

| N | seed | replica interval | samples / batches | elapsed seconds |
|---:|---:|:---|:---|---:|
| 25 | 2026105525 | `[9400000000,9412000000)` | 12M / 100 | 9.2431 |
| 50 | 2026105550 | `[9400000000,9412000000)` | 12M / 100 | 17.0201 |
| 125 | 2026105625 | `[9400000000,9412000000)` | 12M / 100 | 38.7633 |

The two-lane wall time was 39 seconds, essentially identical to the frozen
38.83-second estimate.  All stderr files are empty.  The aarch64 binary SHA256
is `ee9010f524935099ba22f1820fddc05a79dd309d98784dfd5ba7da28129b6856`,
byte-identical to the binary used for the preceding P205 block.  Metadata for
all three jobs records runner commit
`4d82c754a04635de0561b221771099e4ec5a5f88`, the exact seeds and the exact
counter domain.

The scorer was committed as `8114fb5` before local target reveal.  It fits
exactly one amplitude for each of H4, H8 and H12 in frozen order, holds the
radial exponent at `13/8`, and fits no quotient offset or correction term.
Machine-readable output is `analysis/score.json`; copied raw/log artifacts have
matching remote and local SHA256 values in `remote-local-sha256.txt`, and all
archived artifacts are covered by `checksums.sha256`.
